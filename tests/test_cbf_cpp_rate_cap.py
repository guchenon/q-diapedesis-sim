"""
Regression test for H2 (ELN #002): the ZCBF-QP must never return a recanalization
rate above the staged CPP cap, including after the distal pH has been buffered
(pH >= 7.15), where the acidic penalty no longer holds the rate down.
"""

import pytest
from src.neurocritical_cbf import NeurocriticalCBFSafetyFilter


@pytest.mark.parametrize("distal_pH", [6.65, 6.90, 7.10, 7.15, 7.25, 7.35])
@pytest.mark.parametrize("recanalized_fraction", [0.0, 0.5, 0.9])
def test_recanalization_rate_never_exceeds_cpp_cap(distal_pH, recanalized_fraction):
    cbf = NeurocriticalCBFSafetyFilter()
    _, _, cap = cbf.compute_staged_cpp_target(recanalized_fraction, distal_pH)

    res = cbf.filter_control(
        nominal_b_mT=18.0,
        nominal_freq_Hz=80.0,
        nominal_recanal_rate_inv_s=0.35,
        hydrodynamic_wall_shear_Pa=0.70,
        current_delta_T_C=0.0,
        current_glycocalyx_ratio=1.0,
        recanalized_fraction=recanalized_fraction,
        distal_pH=distal_pH,
    )

    surge = cbf.Surge_mmHg_per_recanal_rate * res.safe_recanalization_rate_inv_s
    assert res.safe_recanalization_rate_inv_s <= cap + 1e-9
    assert surge <= cbf.max_cpp_surge_mmHg_s + 1e-9
    assert res.h_cpp_mmHg_s >= -1e-9
