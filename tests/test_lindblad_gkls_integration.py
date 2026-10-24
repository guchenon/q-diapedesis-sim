"""
Tests for D3 Stage 3 (ELN #003): the Lindblad step must be an integration of the GKLS
master equation and nothing else.

Before this change the step added two terms that are not in the master equation:
a blend `rho <- (1 - a) rho + a rho_attractor` that dragged the state onto the target
attractor directly, and a factor 0.05 on the Hamiltonian commutator. With the blend in
place `dissipative_attractor_overlap > 0.80` was a property of the blend, not of the
dissipator, and the Hamiltonian was effectively scaled away.

These tests require the integrator to be CPTP on its own: trace and Hermiticity preserved
without renormalization, positive semidefinite, purity bounded, the correct analytic fixed
point when H = 0, a Hamiltonian that actually acts, and fourth-order convergence in dt.
"""

import numpy as np
from src.quantum_controller import BimodalFloquetController, LindbladEquivariantTTN

KET_PLUS = np.array([[1.0], [1.0]], dtype=complex) / np.sqrt(2.0)
RHO_PLUS = KET_PLUS @ KET_PLUS.conj().T


def _evolve(ttn, h_eff, n_steps: int, dt_s: float, mu: float = 3.8):
    for _ in range(n_steps):
        state = ttn.step_lindblad_density_matrix(dt_s=dt_s, h_floquet_2x2=h_eff, local_viscosity_mPa_s=mu)
    return state


def test_trace_and_hermiticity_are_preserved_by_the_integrator_itself():
    ttn = LindbladEquivariantTTN()
    h = BimodalFloquetController().compute_floquet_magnus_hamiltonian(b_alpha_mT=8.5, b_beta_mT=6.0)
    for _ in range(25):
        ttn.step_lindblad_density_matrix(dt_s=0.05, h_floquet_2x2=h, local_viscosity_mPa_s=3.8)
        assert np.isclose(np.real(np.trace(ttn.rho)), 1.0, atol=1e-9)
        assert np.allclose(ttn.rho, ttn.rho.conj().T, atol=1e-12)


def test_density_matrix_stays_positive_semidefinite_and_physical():
    ttn = LindbladEquivariantTTN()
    h = BimodalFloquetController().compute_floquet_magnus_hamiltonian(b_alpha_mT=8.5, b_beta_mT=6.0)
    for _ in range(50):
        state = ttn.step_lindblad_density_matrix(dt_s=0.05, h_floquet_2x2=h, local_viscosity_mPa_s=3.8)
        assert np.min(np.real(np.linalg.eigvalsh(ttn.rho))) > -1e-9
    assert 0.5 - 1e-9 <= state.purity_tr_rho2 <= 1.0 + 1e-9


def test_zero_hamiltonian_relaxes_to_the_analytic_gkls_fixed_point():
    # With H = 0 and the single jump operator L = |+><-|, the unique steady state is |+><+|.
    ttn = LindbladEquivariantTTN()
    state = _evolve(ttn, np.zeros((2, 2), dtype=complex), n_steps=4000, dt_s=0.05)
    assert state.dissipative_attractor_overlap > 0.999
    assert np.allclose(ttn.rho, RHO_PLUS, atol=1e-3)


def test_hamiltonian_is_not_scaled_away_and_changes_the_steady_state():
    h = BimodalFloquetController().compute_floquet_magnus_hamiltonian(b_alpha_mT=8.5, b_beta_mT=6.0)
    with_h = _evolve(LindbladEquivariantTTN(), h, n_steps=2000, dt_s=0.05)
    without_h = _evolve(LindbladEquivariantTTN(), np.zeros((2, 2), dtype=complex), 2000, 0.05)
    assert abs(with_h.dissipative_attractor_overlap - without_h.dissipative_attractor_overlap) > 1e-3


def test_zero_dissipation_is_unitary_so_purity_is_conserved():
    # With gamma = 0 the GKLS equation reduces to the von Neumann equation, whose flow is
    # unitary: Tr(rho^2) is an exact invariant and the state must NOT drift to the attractor.
    # Any extra relaxation term in the step shows up here.
    h = BimodalFloquetController().compute_floquet_magnus_hamiltonian(b_alpha_mT=8.5, b_beta_mT=6.0)
    ttn = LindbladEquivariantTTN()
    ttn.gamma_viscous = 0.0
    purity_0 = float(np.real(np.trace(ttn.rho @ ttn.rho)))
    overlap_0 = float(np.real(np.trace(ttn.rho @ RHO_PLUS)))
    # dt small enough that RK4 truncation is far below the tolerance: at dt = 0.01 the residual
    # drift over 2 s is ~1.6e-6, and it falls as dt^4, so 0.002 leaves ~1e-8.
    state = _evolve(ttn, h, n_steps=1000, dt_s=0.002)
    assert abs(state.purity_tr_rho2 - purity_0) < 1e-6, "unitary flow must conserve purity"
    assert abs(state.dissipative_attractor_overlap - overlap_0) < 0.5, (
        "with no dissipation the state must not be pulled onto the attractor"
    )


def test_integrator_is_fourth_order_in_dt():
    h = BimodalFloquetController().compute_floquet_magnus_hamiltonian(b_alpha_mT=8.5, b_beta_mT=6.0)
    # Same physical time, three resolutions. RK4 error ~ dt^4, so halving dt should cut the
    # step-to-step difference by far more than the factor 2 a first-order scheme would give.
    rho = {}
    for dt, n in ((0.02, 50), (0.01, 100), (0.005, 200)):
        ttn = LindbladEquivariantTTN()
        _evolve(ttn, h, n_steps=n, dt_s=dt)
        rho[dt] = ttn.rho.copy()
    coarse = np.linalg.norm(rho[0.02] - rho[0.005])
    fine = np.linalg.norm(rho[0.01] - rho[0.005])
    assert fine < coarse / 4.0, f"not high order: coarse={coarse:.3e}, fine={fine:.3e}"
