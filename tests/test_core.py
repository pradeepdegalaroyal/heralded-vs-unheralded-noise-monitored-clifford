import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from monitored_codes import MixedStabilizer, TrajectoryConfig, run_trajectory, symplectic_group  # noqa: E402


def test_group_size():
    assert len(symplectic_group()) == 720


def test_initial_state_entropies():
    st = MixedStabilizer.zero_with_reference(4)
    assert st.entropy([4]) == 1          # reference is maximally mixed
    assert st.entropy(range(5)) == 0     # global state pure
    assert st.entropy([1]) == 0


def test_unitary_preserves_purity():
    rng = np.random.default_rng(0)
    st = MixedStabilizer.zero_with_reference(6)
    grp = symplectic_group()
    for _ in range(20):
        st.apply_two_qubit(2, 3, grp[rng.integers(720)])
    assert st.entropy(range(7)) == 0


def test_dephasing_and_erasure():
    st = MixedStabilizer.zero_with_reference(3)
    st.trace_out(0, reset=False)         # destroys the Bell pair
    assert st.entropy([3]) == 1 and st.entropy(range(4)) == 2   # qubit 0 and reference both mixed
    st2 = MixedStabilizer.zero_with_reference(3)
    st2.measure_z(0)                     # collapses reference too
    assert st2.entropy([3]) == 0


def test_noiseless_p0_keeps_information():
    assert run_trajectory(TrajectoryConfig(8, 0.0), np.random.default_rng(1)) == 1.0


def test_heavy_measurement_destroys_information():
    vals = [run_trajectory(TrajectoryConfig(8, 0.9), np.random.default_rng(s)) for s in range(5)]
    assert max(vals) <= 0.0
