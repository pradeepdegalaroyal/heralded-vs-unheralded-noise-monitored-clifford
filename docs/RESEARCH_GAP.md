# Research gap

## State of the art (from the supplied literature)
- Random and monitored circuits are a central model of entanglement dynamics; measurements compete with scrambling and can drive entanglement phase transitions (Fisher et al., 2023, doi:10.1146/annurev-conmatphys-031720-030658).
- Experiments now run real-time error correction near break-even (Sivak et al., 2023, doi:10.1038/s41586-023-05782-6), report the open question of whether quantum error correction (QEC) can extend coherence in practice, and demonstrate mid-circuit erasure conversion, i.e. turning errors into *located* ones (Ma et al., 2023, doi:10.1038/s41586-023-06438-1), with parallel 99.5% two-qubit gates on neutral atoms (Evered et al., 2023, doi:10.1038/s41586-023-06481-y).
- Quantum-information protection across a network is limited by loss and noise (Azuma et al., 2023, doi:10.1103/RevModPhys.95.045006).

## Gap
The monitored-circuit literature mostly treats ideal measurements, whereas hardware gives a mixture of *unheralded* errors (unknown location, e.g. dephasing) and *heralded* erasures (known location). How much of the dynamically generated "coding phase" survives, and where its measurement-rate threshold sits, for each error type at equal error rate is not quantified in the provided papers. We study this at simulation scale (stabiliser/Clifford circuits, L = 8-16 qubits) — a toy model, not a hardware claim.

## Hypothesis and metric
Metric: coherent information I_c = S(system) − S(system+reference) of a reference qubit initially Bell-paired with the chain, after 2L brickwork layers (1 bit = fully retained). Quantity of interest: p_half, the measurement rate per qubit-layer at which I_c (L = 16) falls below 0.5.

**H1 (falsifiable):** at equal event rate q, a heralded erasure (reset to |0>) is less harmful than an unheralded one (qubit left maximally mixed), i.e. p_half(heralded) > p_half(unheralded) by more than the seed standard deviation. Null: no difference. A secondary question is whether unheralded dephasing (a weaker channel) shifts p_half less than full erasure.
