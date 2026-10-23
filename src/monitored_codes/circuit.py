"""Monitored brickwork Clifford circuit with optional noise."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from .stabilizer import MixedStabilizer, symplectic_group

Noise = Literal["none", "dephasing", "erasure_unheralded", "erasure_heralded"]


@dataclass(frozen=True)
class TrajectoryConfig:
    """Parameters of one trajectory.

    p: probability of a recorded Z measurement per qubit per layer.
    q: probability of a noise event per qubit per layer.
    noise: channel applied on noise events.
    depth_factor: number of layers is ``depth_factor * L``.
    """

    L: int
    p: float
    q: float = 0.0
    noise: Noise = "none"
    depth_factor: int = 2


def coherent_information(st: MixedStabilizer, L: int) -> float:
    """I_c = S(system) - S(system + reference) in bits; 1 = reference info fully retained."""
    return float(st.entropy(range(L)) - st.entropy(range(L + 1)))


def run_trajectory(cfg: TrajectoryConfig, rng: np.random.Generator) -> float:
    """Evolve for ``depth_factor * L`` layers and return the final coherent information."""
    group = symplectic_group()
    st = MixedStabilizer.zero_with_reference(cfg.L)
    for t in range(cfg.depth_factor * cfg.L):
        for i in range(t % 2, cfg.L - 1, 2):
            st.apply_two_qubit(i, i + 1, group[rng.integers(len(group))])
        for i in range(cfg.L):
            if rng.random() < cfg.p:
                st.measure_z(i)
        if cfg.noise != "none" and cfg.q > 0:
            for i in range(cfg.L):
                if rng.random() < cfg.q:
                    if cfg.noise == "dephasing":
                        st.dephase_z(i)
                    else:
                        st.trace_out(i, reset=cfg.noise == "erasure_heralded")
    return coherent_information(st, cfg.L)
