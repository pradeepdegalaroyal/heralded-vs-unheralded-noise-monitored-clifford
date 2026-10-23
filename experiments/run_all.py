"""Reproduce all results: p-sweeps of coherent information per noise channel, then figures.

Usage: python experiments/run_all.py [--quick]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from monitored_codes import TrajectoryConfig, run_trajectory  # noqa: E402


def crossing(ps: np.ndarray, a: np.ndarray, b: np.ndarray) -> float:
    """Linear-interpolated p where curve a (small L) and b (large L) cross; NaN if none."""
    d = a - b
    for k in range(len(ps) - 1):
        if d[k] * d[k + 1] < 0:
            return float(ps[k] - d[k] * (ps[k + 1] - ps[k]) / (d[k + 1] - d[k]))
    return float("nan")


def half_point(ps: np.ndarray, y: np.ndarray, level: float = 0.5) -> float:
    """Smallest p (interpolated) where y drops below ``level``; NaN if never / always."""
    for k in range(len(ps) - 1):
        if y[k] >= level > y[k + 1]:
            return float(ps[k] + (y[k] - level) * (ps[k + 1] - ps[k]) / (y[k] - y[k + 1]))
    return float("nan")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    quick = args.quick

    sizes = [6, 8] if quick else [8, 12, 16]
    ps = np.round(np.arange(0.04, 0.31, 0.08 if quick else 0.025), 3)
    seeds, per_seed = (2, 3) if quick else (4, 15)
    qs = [0.02] if quick else [0.002, 0.005, 0.02]
    conds = [("none", 0.0)] + [(n, q) for n in ("dephasing", "erasure_unheralded", "erasure_heralded") for q in qs]

    out = ROOT / ("results_quick" if quick else "results")
    out.mkdir(exist_ok=True)
    rows = []
    t0 = time.time()
    for noise, q in conds:
        for L in sizes:
            for p in ps:
                for s in range(seeds):
                    rng = np.random.default_rng(1000 * s + L)
                    vals = [run_trajectory(TrajectoryConfig(L, float(p), q, noise), rng) for _ in range(per_seed)]
                    rows.append(dict(noise=noise, q=q, L=L, p=float(p), seed=s, ic=float(np.mean(vals))))
        print(f"{noise} q={q} done ({time.time() - t0:.0f}s)", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out / "sweep_raw.csv", index=False)

    agg = df.groupby(["noise", "q", "L", "p"]).ic.agg(["mean", "std"]).reset_index()
    agg.to_csv(out / "sweep_summary.csv", index=False)

    # Threshold estimates per condition: crossing of smallest/largest L, and half-point of largest L.
    summ = []
    for (noise, q), g in agg.groupby(["noise", "q"]):
        small = g[g.L == sizes[0]].sort_values("p")["mean"].to_numpy()
        large = g[g.L == sizes[-1]].sort_values("p")["mean"].to_numpy()
        # per-seed uncertainty of half-point of the largest size
        hp = []
        for s in range(seeds):
            y = df[(df.noise == noise) & (df.q == q) & (df.L == sizes[-1]) & (df.seed == s)].sort_values("p").ic.to_numpy()
            hp.append(half_point(ps, y))
        summ.append(dict(
            noise=noise, q=q,
            p_cross=crossing(ps, small, large),
            p_half_mean=float(np.nanmean(hp)) if not np.all(np.isnan(hp)) else float("nan"),
            p_half_std=float(np.nanstd(hp)) if not np.all(np.isnan(hp)) else float("nan"),
            ic_at_p0_largest_L=float(large[0]),
        ))
    pd.DataFrame(summ).to_csv(out / "thresholds.csv", index=False)
    json.dump(dict(sizes=sizes, p_grid=ps.tolist(), seeds=seeds, trajectories_per_seed=per_seed,
                   runtime_s=time.time() - t0), open(out / "meta.json", "w"), indent=1)

    make_figures(agg, pd.DataFrame(summ), sizes, ROOT / ("figures_quick" if quick else "figures"))
    print(f"total {time.time() - t0:.0f}s")


def make_figures(agg: pd.DataFrame, thr: pd.DataFrame, sizes: list[int], figdir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir.mkdir(exist_ok=True)
    # Fig 1: noiseless scaling
    g0 = agg[agg.noise == "none"]
    fig, ax = plt.subplots(figsize=(5, 3.6))
    for L in sizes:
        d = g0[g0.L == L].sort_values("p")
        ax.errorbar(d.p, d["mean"], d["std"], marker="o", ms=3, capsize=2, label=f"L={L}")
    ax.set_xlabel("measurement rate p"); ax.set_ylabel("coherent information $I_c$ (bits)")
    ax.set_title("Noiseless monitored Clifford circuit"); ax.legend()
    fig.tight_layout(); fig.savefig(figdir / "fig1_noiseless_scaling.png", dpi=200); plt.close(fig)

    # Fig 2: channels at each q, largest L
    qs = sorted(q for q in agg.q.unique() if q > 0)
    fig, axes = plt.subplots(1, len(qs), figsize=(4.2 * len(qs), 3.6), sharey=True, squeeze=False)
    Lmax = sizes[-1]
    for ax, q in zip(axes[0], qs):
        d = agg[(agg.noise == "none") & (agg.L == Lmax)].sort_values("p")
        ax.plot(d.p, d["mean"], "k--", label="no noise")
        for noise in ("dephasing", "erasure_unheralded", "erasure_heralded"):
            d = agg[(agg.noise == noise) & (agg.q == q) & (agg.L == Lmax)].sort_values("p")
            ax.errorbar(d.p, d["mean"], d["std"], marker="o", ms=3, capsize=2, label=noise.replace("_", " "))
        ax.set_title(f"noise rate q={q}, L={Lmax}"); ax.set_xlabel("measurement rate p")
    axes[0][0].set_ylabel("coherent information $I_c$ (bits)"); axes[0][0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(figdir / "fig2_noise_channels.png", dpi=200); plt.close(fig)

    # Fig 3: half-point threshold vs q
    fig, ax = plt.subplots(figsize=(5, 3.6))
    for noise in ("dephasing", "erasure_unheralded", "erasure_heralded"):
        d = thr[thr.noise == noise].sort_values("q")
        ax.errorbar(d.q, d.p_half_mean, d.p_half_std, marker="o", capsize=2, label=noise.replace("_", " "))
    n0 = thr[thr.noise == "none"].p_half_mean.iloc[0]
    ax.axhline(n0, color="k", ls="--", label="no noise")
    ax.set_xlabel("noise rate q"); ax.set_ylabel(f"p where $I_c$ falls below 0.5 (L={Lmax})")
    ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(figdir / "fig3_threshold_vs_noise.png", dpi=200); plt.close(fig)


if __name__ == "__main__":
    main()
