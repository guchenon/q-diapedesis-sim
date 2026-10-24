"""
Tests for D3 Stage 3 / 6e (ELN #003):

1. `compute_floquet_magnus_hamiltonian` returned an H_eff that `evaluate_control_step`
   computed and then discarded - `dual_task_fidelity` was a weighted score of two step-out
   ratios and never touched the Hamiltonian. There is now a real, computed fidelity obtained
   by evolving a state under the true time-dependent drive and comparing against the
   propagator H_eff predicts. It is reported alongside, not instead of, the preregistered
   design score, because the two measure different things.

2. The closed-loop pipeline recorded fidelity from the NOMINAL Floquet evaluation, so the CBF
   filter's reduction of the field could not show up in it at all. The recorded fidelity must
   come from the field actually applied.
"""

import numpy as np
import pytest
from src.quantum_controller import BimodalFloquetController, QDiapedesisClosedLoopPipeline


def test_drive_hamiltonian_fourier_components_match_the_magnus_expansion():
    f = BimodalFloquetController()
    b_a, b_b, fa = 8.5, 6.0, 80.0
    period = 1.0 / fa
    ts = np.linspace(0.0, period, 4001)[:-1]
    h_t = np.array([f.bimodal_drive_hamiltonian(t, b_a, b_b, fa) for t in ts])

    # H_0: the time average must equal the zeroth-order Magnus term.
    h_0_numeric = h_t.mean(axis=0)
    h_0_expected = 0.5 * f.b_static_mT * f.sigma_z + 0.5 * b_b * f.sigma_x
    assert np.allclose(h_0_numeric, h_0_expected, atol=1e-6)

    # H_{+1}: the first Fourier coefficient must equal (B_alpha/4)(sigma_x - i sigma_y).
    weights = np.exp(-1.0j * 2.0 * np.pi * fa * ts)
    h_plus_numeric = (h_t * weights[:, None, None]).mean(axis=0)
    h_plus_expected = 0.25 * b_a * (f.sigma_x - 1.0j * f.sigma_y)
    assert np.allclose(h_plus_numeric, h_plus_expected, atol=1e-6)


def test_exact_propagator_is_unitary_and_converged():
    f = BimodalFloquetController()
    u = f.exact_one_period_propagator(8.5, 6.0, n_substeps=2048)
    assert np.allclose(u.conj().T @ u, np.eye(2), atol=1e-10)
    coarse = f.exact_one_period_propagator(8.5, 6.0, n_substeps=512)
    assert np.linalg.norm(u - coarse) < 1e-6, "propagator not converged in n_substeps"


def test_magnus_fidelity_is_a_real_fidelity_bounded_by_one():
    f = BimodalFloquetController()
    fid = f.magnus_propagator_fidelity(8.5, 6.0)
    assert 0.0 <= fid <= 1.0 + 1e-12
    # Perfect agreement when H_eff IS the generator: a zero drive amplitude leaves only H_0.
    assert f.magnus_propagator_fidelity(0.0, 6.0) == pytest.approx(1.0, abs=1e-9)


def test_magnus_fidelity_degrades_as_the_high_frequency_expansion_is_pushed():
    # The expansion is controlled by ||H|| / omega, so lowering the drive frequency at fixed
    # amplitude must make the truncated H_eff a worse description of the exact propagator.
    f = BimodalFloquetController()
    fast = f.magnus_propagator_fidelity(8.5, 6.0, f_alpha_Hz=80.0)
    slow = f.magnus_propagator_fidelity(8.5, 6.0, f_alpha_Hz=2.0)
    assert slow < fast


def test_pipeline_fidelity_responds_to_the_cbf_filtered_field():
    # With the CBF filter on, the applied field is at or below nominal, so the recorded
    # fidelity must differ from the unfiltered run. Before this change both runs recorded the
    # same nominal number and the CBF's effect on fidelity was invisible by construction.
    # Use the 15 um vessel, where the shear ceiling actually binds and the filter has to cut
    # the field. At 30 um there is shear headroom and the filter is almost inactive.
    on = QDiapedesisClosedLoopPipeline(15.0).run_episode(total_time_s=4.0, dt_s=0.05)
    off = QDiapedesisClosedLoopPipeline(15.0).run_episode(
        total_time_s=4.0, dt_s=0.05, enable_cbf_filter=False
    )
    assert on.mean_applied_b_mT < on.mean_nominal_b_mT, "the filter should be binding at 15 um"
    assert on.dual_task_floquet_fidelity != off.dual_task_floquet_fidelity, (
        "the recorded fidelity must depend on the field actually applied"
    )
    # D10 (ELN #003): the direction has reversed, and is now the physical one. Cutting the field
    # LOWERS the score (0.9064 filtered vs 0.9485 unfiltered at 15 um), because a starved field
    # drives the 80 Hz sub-swarm toward step-out. Before D10 the filtered arm scored higher,
    # which made the CBF's effect on actuation fidelity impossible to falsify.
    assert on.dual_task_floquet_fidelity < off.dual_task_floquet_fidelity


# D10 (ELN #003): this xfail has been REMOVED - the defect it recorded is fixed. The score used
# to be monotonically non-increasing in field amplitude (a field at 10% of nominal scored
# 0.994371 against 0.992558 at full nominal) and was blind to any reduction below 0.85x because
# the effective step-out frequency was clipped. The clip is gone and the drilling term now
# carries the available magnetic torque, so the score falls when the field is starved.
def test_dual_task_score_should_penalize_an_inadequate_field():
    f = BimodalFloquetController()
    full = f.evaluate_control_step(0.0, b_scale_factor=1.0, f_alpha_override_Hz=80.0, f_beta_override_Hz=5.0)
    starved = f.evaluate_control_step(0.0, b_scale_factor=0.1, f_alpha_override_Hz=80.0, f_beta_override_Hz=5.0)
    assert starved.dual_task_fidelity < full.dual_task_fidelity
