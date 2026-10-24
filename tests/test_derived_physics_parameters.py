"""
Tests for D7 / D8 (ELN #003): the Kuramoto coupling and the Lindblad dissipation rate must be
DERIVED from the configured particle physics, not free numbers in the YAML.

D7. `kuramoto_coupling_K: 4.50` and D8 `viscous_dissipation_gamma: 0.15` were free numbers. Both
are present in the v0.1.0-alpha snapshot but neither is named in the ELN #001 lock list nor
asserted by the preregistered tests/test_config_sanity.py, so both are revisable with a
recorded rationale rather than an OSF hypothesis deviation.

The derivation chain, with no new free parameter anywhere:
  zeta_r = 8 pi eta a^3                        rotational drag of the configured particle
  m      = 2 pi f_so zeta_r / B                moment implied by the configured step-out frequency
  K      = (mu_0 m^2 / 4 pi zeta_r) * sum_j (1.5 sin^2 psi_ij - 1) / r_ij^3
  gamma  = m B / zeta_r                        overdamped magnetic-viscous relaxation rate

The geometric factor in K is the part that matters: +0.5 for a neighbour in the rotation plane,
-1 for one along the drive axis, and its isotropic average is exactly zero. These tests pin all
three limits, because an overlapping (non-hard-sphere) sample makes 1/r^3 diverge and silently
fabricates a large positive coupling - a mistake made and caught during this work.
"""

import numpy as np
import pytest
from src.physics import (
    dipolar_kuramoto_coupling,
    magnetic_moment_from_stepout,
    magnetic_relaxation_rate_inv_s,
    rotational_diffusion_rad2_s,
    rotational_drag_coefficient,
    sample_near_wall_shell,
)
from src.biomimicry_mechanics import HydroMagneticSolitonSwarm
from src.quantum_controller import BimodalFloquetController, LindbladEquivariantTTN

A_M = 0.75e-6
ETA = 1.2e-3


def test_rotational_drag_and_diffusion_match_closed_form():
    assert rotational_drag_coefficient(A_M, ETA) == pytest.approx(8 * np.pi * ETA * A_M**3, rel=1e-12)
    d_r = rotational_diffusion_rad2_s(A_M, ETA, 310.0)
    assert d_r == pytest.approx(1.380649e-23 * 310.0 / (8 * np.pi * ETA * A_M**3), rel=1e-12)


def test_moment_from_stepout_inverts_the_mahoney_balance():
    # m is defined so that m B / zeta_r equals 2 pi f_so exactly.
    m = magnetic_moment_from_stepout(95.0, 8.5e-3, A_M, ETA)
    assert magnetic_relaxation_rate_inv_s(m, 8.5e-3, A_M, ETA) == pytest.approx(2 * np.pi * 95.0, rel=1e-12)
    # And the implied effective magnetization must be physically plausible for a composite bead.
    m_per_volume = m / ((4 / 3) * np.pi * A_M**3)
    assert 50.0 < m_per_volume < 5.0e5


@pytest.mark.parametrize(
    "label, builder, expect_sign",
    [
        ("axial chain along the spin axis", lambda: np.column_stack(
            [np.zeros(60), np.zeros(60), np.arange(60) * 2.0 * A_M]), -1),
        ("single layer normal to the spin axis", lambda: (lambda g: np.column_stack(
            [g[0].ravel(), g[1].ravel(), np.zeros(g[0].size)]))(
            np.meshgrid(np.arange(40) * 2.0 * A_M, np.arange(40) * 2.0 * A_M)), +1),
    ],
)
def test_geometric_factor_sign_matches_the_analytic_limits(label, builder, expect_sign):
    m = magnetic_moment_from_stepout(95.0, 8.5e-3, A_M, ETA)
    res = dipolar_kuramoto_coupling(m, builder(), A_M, ETA, n_reference=30, seed=1)
    assert np.sign(res.coupling_K_rad_s) == expect_sign, f"{label}: K={res.coupling_K_rad_s}"


def test_isotropic_hard_sphere_cloud_has_essentially_no_phase_coupling():
    # <1.5 sin^2 psi - 1> = 0 for an isotropic neighbour distribution, so a cubic lattice must
    # give a coupling far below the in-plane value. A non-hard-sphere sample would not.
    m = magnetic_moment_from_stepout(95.0, 8.5e-3, A_M, ETA)
    n = 16
    grid = np.meshgrid(*[np.arange(n) * 2.0 * A_M] * 3, indexing="ij")
    cubic = np.column_stack([g.ravel() for g in grid])
    k_iso = dipolar_kuramoto_coupling(m, cubic, A_M, ETA, n_reference=40, seed=1).coupling_K_rad_s
    sheet = np.meshgrid(np.arange(40) * 2.0 * A_M, np.arange(40) * 2.0 * A_M)
    plane = np.column_stack([sheet[0].ravel(), sheet[1].ravel(), np.zeros(sheet[0].size)])
    k_plane = dipolar_kuramoto_coupling(m, plane, A_M, ETA, n_reference=40, seed=1).coupling_K_rad_s
    assert abs(k_iso) < 0.2 * abs(k_plane)


def test_shell_sampler_refuses_a_dose_that_cannot_physically_fit():
    # 10 000 particles of radius 0.75 um in the near-wall band over only 150 um is above random
    # close packing; the sampler must say so rather than quietly producing overlaps.
    with pytest.raises(ValueError, match="does not fit"):
        sample_near_wall_shell(10000, 15e-6, (0.75, 0.94), 150e-6, A_M, seed=1)


def test_swarm_coupling_is_derived_and_its_sign_is_honoured():
    swarm = HydroMagneticSolitonSwarm()
    assert swarm.coupling_is_derived is True
    assert swarm.derived_coupling_K_rad_s == pytest.approx(swarm.coupling_K, rel=1e-12)
    # The modeled near-wall geometry gives an ANTI-aligning coupling; the code must not clip it
    # to zero or flip it, because that would hide the result.
    assert swarm.coupling_K < 0.0


def test_negative_coupling_does_not_produce_synchronization():
    swarm = HydroMagneticSolitonSwarm(seed=4)
    swarm.reset()
    for _ in range(40):
        state = swarm.step_soliton_swarm(0.05, 12.0, 0.88, True, 0.08, True)
    assert state.kuramoto_order_r < 0.40


def test_dissipation_rate_is_derived_and_commensurate_with_the_hamiltonian():
    ttn = LindbladEquivariantTTN()
    assert ttn.dissipation_is_derived is True
    h = BimodalFloquetController().compute_floquet_magnus_hamiltonian(8.5, 6.0)
    # Both scales carry the same m/zeta_r factor, so the physical ratio is what the dynamics see.
    ratio = (ttn.gamma_viscous * 3.8 / 3.0) / float(np.linalg.norm(h))
    assert ratio == pytest.approx(0.747, abs=0.05), f"physical gamma/||H|| ratio is {ratio}"
    # Far above the free-parameter value it replaces, which was 45x too small.
    assert ttn.gamma_viscous > 20.0 * 0.15
