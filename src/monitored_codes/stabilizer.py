"""Sign-free mixed stabiliser states over GF(2).

Entropies of stabiliser (mixed) states depend only on the binary check matrix,
not on Pauli signs, so signs are dropped.  The state on ``n`` qubits is the
row space of a ``(m, 2n)`` binary matrix ``[X | Z]`` with ``m <= n``
independent generators (``m < n`` means the state is mixed).  For a region
``A``, ``S_A = |A| - dim(G_A)`` where ``G_A`` is the subgroup supported in ``A``.
"""
from __future__ import annotations

from functools import lru_cache
from typing import Iterable

import numpy as np


def gf2_rank(mat: np.ndarray) -> int:
    """Rank of a binary matrix over GF(2)."""
    a = mat.astype(np.uint8).copy()
    rank = 0
    rows, cols = a.shape
    for c in range(cols):
        if rank == rows:
            break
        piv = np.nonzero(a[rank:, c])[0]
        if piv.size == 0:
            continue
        p = piv[0] + rank
        if p != rank:
            a[[rank, p]] = a[[p, rank]]
        mask = a[:, c].astype(bool)
        mask[rank] = False
        a[mask] ^= a[rank]
        rank += 1
    return rank


@lru_cache(maxsize=1)
def symplectic_group() -> np.ndarray:
    """All 720 elements of Sp(4,2) acting on columns ``(x_i, x_j, z_i, z_j)``."""
    eye = np.eye(4, dtype=np.uint8)
    h0 = eye[[2, 1, 0, 3]]                                # Hadamard, qubit 0
    h1 = eye[[0, 3, 2, 1]]                                # Hadamard, qubit 1
    s0 = eye.copy(); s0[0, 2] = 1                         # phase gate, qubit 0
    s1 = eye.copy(); s1[1, 3] = 1
    cx = eye.copy(); cx[0, 1] = 1; cx[3, 2] = 1           # CNOT 0 -> 1
    gens = [h0, h1, s0, s1, cx]
    seen = {eye.tobytes(): eye}
    frontier = [eye]
    while frontier:
        new = []
        for m in frontier:
            for g in gens:
                p = ((m.astype(int) @ g.astype(int)) % 2).astype(np.uint8)
                k = p.tobytes()
                if k not in seen:
                    seen[k] = p
                    new.append(p)
        frontier = new
    return np.stack(list(seen.values()))


class MixedStabilizer:
    """Mixed stabiliser state with Clifford gates, measurements and noise."""

    def __init__(self, n: int, gens: np.ndarray | None = None) -> None:
        self.n = n
        self.g = np.zeros((0, 2 * n), dtype=np.uint8) if gens is None else gens.astype(np.uint8)

    @classmethod
    def zero_with_reference(cls, n_sys: int) -> "MixedStabilizer":
        """|0>^n_sys plus a reference qubit (index n_sys) Bell-paired with qubit 0."""
        n = n_sys + 1
        rows = []
        for i in range(1, n_sys):
            r = np.zeros(2 * n, dtype=np.uint8); r[n + i] = 1; rows.append(r)
        r = np.zeros(2 * n, dtype=np.uint8); r[0] = 1; r[n_sys] = 1; rows.append(r)           # XX
        r = np.zeros(2 * n, dtype=np.uint8); r[n] = 1; r[n + n_sys] = 1; rows.append(r)       # ZZ
        return cls(n, np.array(rows))

    # -- gates -----------------------------------------------------------
    def apply_two_qubit(self, i: int, j: int, m: np.ndarray) -> None:
        """Apply the Clifford with symplectic matrix ``m`` to qubits (i, j)."""
        cols = [i, j, self.n + i, self.n + j]
        self.g[:, cols] = (self.g[:, cols].astype(np.int64) @ m) % 2

    # -- helpers ---------------------------------------------------------
    def _eliminate(self, col: int) -> int | None:
        """Leave at most one row with a 1 in ``col``; return its index (or None)."""
        idx = np.nonzero(self.g[:, col])[0]
        if idx.size == 0:
            return None
        p = idx[0]
        for r in idx[1:]:
            self.g[r] ^= self.g[p]
        return int(p)

    # -- measurements and channels --------------------------------------
    def measure_z(self, i: int) -> None:
        """Recorded projective Z measurement (outcome known; purifies)."""
        p = self._eliminate(i)
        z = np.zeros(2 * self.n, dtype=np.uint8); z[self.n + i] = 1
        if p is not None:
            self.g[p] = z
        elif gf2_rank(np.vstack([self.g, z])) > self.g.shape[0]:
            self.g = np.vstack([self.g, z])

    def dephase_z(self, i: int) -> None:
        """Unrecorded Z measurement: keep the subgroup commuting with Z_i."""
        p = self._eliminate(i)
        if p is not None:
            self.g = np.delete(self.g, p, axis=0)

    def trace_out(self, i: int, reset: bool) -> None:
        """Discard qubit i; ``reset`` re-prepares |0> (heralded erasure), else maximally mixed."""
        for col in (i, self.n + i):
            p = self._eliminate(col)
            if p is not None:
                self.g = np.delete(self.g, p, axis=0)
        if reset:
            z = np.zeros((1, 2 * self.n), dtype=np.uint8); z[0, self.n + i] = 1
            self.g = np.vstack([self.g, z])

    # -- entropies -------------------------------------------------------
    def entropy(self, region: Iterable[int]) -> int:
        """Von Neumann entropy (bits) of the reduced state on ``region``."""
        region = list(region)
        m = self.g.shape[0]
        if m == 0:
            return len(region)
        inside = set(region)
        outside = [q for q in range(self.n) if q not in inside]
        cols = outside + [self.n + q for q in outside]
        rank_out = gf2_rank(self.g[:, cols]) if cols else 0
        return len(region) - (m - rank_out)
