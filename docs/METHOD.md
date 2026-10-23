# Method

- **Simulator:** sign-free mixed stabiliser states over GF(2) (own implementation, `src/monitored_codes/stabilizer.py`); entropies via S_A = |A| − dim(G_A). Two-qubit gates are drawn uniformly from Sp(4,2) (720 elements) on a brickwork chain with open boundaries.
- **Layer:** random Clifford gates → recorded Z measurement per qubit with prob. p → noise event per qubit with prob. q. Depth 2L.
- **Noise channels:** dephasing (unrecorded Z measurement); unheralded erasure (qubit traced out, left maximally mixed); heralded erasure (traced out, reset to |0>).
- **Baselines:** (i) the noiseless monitored circuit (standard measurement-induced-transition setting), (ii) the p = 0 unitary-only circuit (I_c = 1 by construction); the unheralded channel is the classical "unknown error location" baseline against the heralded one.
- **Metric:** coherent information I_c; threshold proxy p_half (L = 16); noiseless crossing of L = 8 and 16 curves.
- **Ablation:** heralded vs unheralded erasure (only difference: the reset); dephasing vs erasure (channel strength).
- **Noise sweep:** q ∈ {0.002, 0.005, 0.02}; **scaling:** L ∈ {8, 12, 16}; p from 0.04 to 0.29 in steps of 0.025.
- **Seeds:** 4 independent seeds × 15 trajectories each; error bars are std over seeds.
