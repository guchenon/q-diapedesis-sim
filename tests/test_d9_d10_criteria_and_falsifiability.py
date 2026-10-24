"""
Tests for D9 and D10 (ELN #003).

D9. The barren-plateau criterion was `gradient_variance >= 1e-3`, a fixed floor. Freedom from
barren plateaus means the variance decays only POLYNOMIALLY in N, so a fixed floor tests the
wrong property - any polynomially decaying quantity falls below any constant at large enough N
(here from about N = 24). The criterion is now the measured decay exponent plus the margin over
the unstructured global-cost baseline. OSF Hypothesis 4 registers "without barren-plateau
gradient vanishing" and not the constant, so this is a criterion correction, not a deviation.

D10. `dual_task_fidelity` had no term penalizing an inadequate field, so it was monotonically
non-increasing in field amplitude - a field at 10% of nominal scored HIGHER than full nominal -
and the effective step-out frequency was clipped by `np.clip(scale, 0.85, 1.25)`, making the
score blind to any reduction below 0.85x. Per Mahoney et al. 2014 the step-out frequency is
proportional to the field amplitude, so the clip is removed and the drilling term now carries
the actual rotation rate against the target times the available magnetic torque fraction.

NOTE: D10 changes the operational definition of a PREREGISTERED endpoint
(`target_dual_task_fidelity_min: 0.92` is named in the ELN #001 lock list and asserted by
tests/test_config_sanity.py). It is recorded as DEVIATION-01 in ELN #003 and must be filed on
OSF. The threshold value itself is untouched.
"""

import numpy as np
import pytest
from scipy.linalg import expm, logm
from src.quantum_controller import BimodalFloquetController, LindbladEquivariantTTN


# --------------------------------------------------------------------------- D9
def test_barren_plateau_verdict_uses_the_decay_exponent_not_a_fixed_floor():
    ttn = LindbladEquivariantTTN()
    report = ttn.benchmark_barren_plateau_scaling()
    p = report["equivariant_polynomial_decay_exponent_p"]
    assert 0.0 < p < ttn.max_polynomial_decay_exponent
    assert report["variance_improvement_orders_of_magnitude"] >= ttn.min_improvement_orders
    assert report["equivariant_log_log_fit_max_residual_log2"] <= ttn.max_log_log_fit_residual
    # The verdict must now be true even though the measured variance is far below the old floor.
    assert ttn.barren_plateau_free is True
    assert ttn.measured_gradient_variance < 1.0e-3


def test_a_fixed_floor_would_have_rejected_a_plateau_free_ansatz():
    # Guard against reintroducing the mistake: the polynomial decay crosses 1e-3 well before
    # the configured N = 200, so a constant floor and the real criterion disagree here.
    ttn = LindbladEquivariantTTN()
    assert ttn.measured_gradient_variance < ttn.legacy_gradient_variance_floor
    assert ttn.barren_plateau_free is True


# -------------------------------------------------------------------------- D10
def test_stepout_frequency_is_proportional_to_field_amplitude():
    f = BimodalFloquetController()
    lo = f.evaluate_control_step(0.0, b_scale_factor=0.5, f_alpha_override_Hz=80.0, f_beta_override_Hz=5.0)
    hi = f.evaluate_control_step(0.0, b_scale_factor=1.0, f_alpha_override_Hz=80.0, f_beta_override_Hz=5.0)
    # At half field the alpha step-out drops to 47.5 Hz, below the 80 Hz drive, so the rotor
    # must slip and its rotation rate must fall well short of the drive.
    assert lo.alpha_rotation_rate_Hz < 0.5 * 80.0
    assert hi.alpha_rotation_rate_Hz == pytest.approx(80.0, rel=1e-9)


def test_dual_task_score_now_penalizes_an_inadequate_field():
    f = BimodalFloquetController()
    scores = [
        f.evaluate_control_step(0.0, b_scale_factor=s, f_alpha_override_Hz=80.0, f_beta_override_Hz=5.0).dual_task_fidelity
        for s in (0.1, 0.3, 0.5, 0.8, 0.9, 1.0)
    ]
    assert all(np.diff(scores) > 0.0), f"score must rise with field amplitude, got {scores}"
    assert scores[0] < 0.2, "a field at 10% of nominal must not score well"


def test_nominal_design_point_still_meets_the_preregistered_threshold():
    # The registered 0.92 is untouched; the reformulated score must still clear it at the
    # nominal field, otherwise the redefinition would have moved the design point.
    f = BimodalFloquetController()
    nominal = f.evaluate_control_step(0.0, b_scale_factor=1.0, f_alpha_override_Hz=80.0, f_beta_override_Hz=5.0)
    assert nominal.dual_task_fidelity >= 0.92
    assert nominal.passes_fidelity_target is True


def test_second_order_magnus_term_improves_on_first_order():
    """
    The old second-order term carried an invented coefficient 0.15 and had no measurable effect
    (orders 1 and 2 agreed to eight decimal places). The replacement is the van Vleck term
        H^(2) = ( [[H_+1, H_0], H_-1] + [[H_-1, H_0], H_+1] ) / (2 omega^2)
    and this test checks it against the numerically exact propagator.

    The comparison is on the QUASI-ENERGY SPECTRUM, not the full matrix. The van Vleck effective
    Hamiltonian and the Floquet-Magnus Hamiltonian read off i log U(T) / T at t0 = 0 differ by a
    kick (gauge) operator at O(1/omega), so the matrices never converge to each other however
    many orders are kept - comparing them directly hides the improvement, which is what happened
    on the first attempt here. The quasi-energies are gauge-invariant and do converge.
    """
    f = BimodalFloquetController()

    def quasi_energy_error(f_alpha_Hz: float, order: int) -> float:
        period = 1.0 / f_alpha_Hz
        u_exact = f.exact_one_period_propagator(8.5, 6.0, f_alpha_Hz, n_substeps=8192)
        h_exact = 1.0j * logm(u_exact) / period
        e_exact = np.sort(np.real(np.linalg.eigvalsh(0.5 * (h_exact + h_exact.conj().T))))
        before = f.truncation_order
        f.truncation_order = order
        try:
            h = f.compute_floquet_magnus_hamiltonian(8.5, 6.0, f_alpha_Hz, 5.0)
        finally:
            f.truncation_order = before
        e = np.sort(np.real(np.linalg.eigvalsh(0.5 * (h + h.conj().T))))
        return float(np.linalg.norm(e - e_exact))

    for f_alpha_Hz in (12.0, 20.0, 40.0, 80.0):
        e1 = quasi_energy_error(f_alpha_Hz, 1)
        e2 = quasi_energy_error(f_alpha_Hz, 2)
        assert e2 < 0.5 * e1, f"at {f_alpha_Hz} Hz: order1={e1:.3e}, order2={e2:.3e}"


def test_magnus_truncation_errors_follow_the_expected_order_in_frequency():
    # A correct expansion has order-1 error ~ omega^-2 and order-2 error ~ omega^-3, so the
    # order-2 error must fall off faster with drive frequency than the order-1 error does.
    f = BimodalFloquetController()

    def err(f_alpha_Hz: float, order: int) -> float:
        period = 1.0 / f_alpha_Hz
        h_exact = 1.0j * logm(f.exact_one_period_propagator(8.5, 6.0, f_alpha_Hz, n_substeps=8192)) / period
        e_exact = np.sort(np.real(np.linalg.eigvalsh(0.5 * (h_exact + h_exact.conj().T))))
        before = f.truncation_order
        f.truncation_order = order
        try:
            h = f.compute_floquet_magnus_hamiltonian(8.5, 6.0, f_alpha_Hz, 5.0)
        finally:
            f.truncation_order = before
        e = np.sort(np.real(np.linalg.eigvalsh(0.5 * (h + h.conj().T))))
        return float(np.linalg.norm(e - e_exact))

    freqs = np.array([12.0, 20.0, 40.0, 80.0])
    s1 = np.polyfit(np.log(freqs), np.log([err(x, 1) for x in freqs]), 1)[0]
    s2 = np.polyfit(np.log(freqs), np.log([err(x, 2) for x in freqs]), 1)[0]
    assert s1 < -1.5, f"order-1 error should fall at least as fast as omega^-1.5, got {s1:.2f}"
    assert s2 < s1 - 0.5, f"order-2 should converge faster: slopes {s1:.2f} vs {s2:.2f}"
