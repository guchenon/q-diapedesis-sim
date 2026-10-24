"""
Tests for D3 Stage 2 (ELN #003): the Sakaguchi-Kuramoto swarm must be a real
synchronization model, not a pre-synchronized one.

Before this change the model had no natural-frequency dispersion and no rotational
noise, and phases were initialized inside a +/-0.6 rad window, so the order parameter
started near 0.95 and rose trivially. A Kuramoto model with no detuning and no noise
cannot desynchronize, so `kuramoto_order_r >= 0.85` was true by construction rather
than a property of the coupling.

These tests require:
  1. phases spread over the whole circle at reset (incoherent start, r ~ O(1/sqrt(N)));
  2. a non-degenerate distribution of natural frequencies within each sub-swarm;
  3. rotational noise, so two different seeds give different trajectories;
  4. a genuine coupling-dependent transition: strong coupling relative to the frequency
     spread synchronizes, weak coupling does not.
"""

import numpy as np
import pytest
from src.biomimicry_mechanics import HydroMagneticSolitonSwarm


def _evolve(swarm, n_steps: int = 40, dt_s: float = 0.05, coupled: bool = True):
    for _ in range(n_steps):
        state = swarm.step_soliton_swarm(
            dt_s=dt_s,
            axial_rolling_velocity_um_s=12.0,
            wall_proximal_fraction=0.88,
            adhesion_locked=True,
            collateral_skimming_fraction=0.08,
            enable_soliton_coupling=coupled,
        )
    return state


def test_phases_start_incoherent_over_the_full_circle():
    init = HydroMagneticSolitonSwarm(seed=1).reset()
    phases = init.phases_rad
    assert phases.max() - phases.min() > 5.0, "phases must span most of [-pi, pi)"
    # An incoherent start has r of order 1/sqrt(N); N = 200 -> ~0.07.
    assert init.kuramoto_order_r < 0.30


def test_natural_frequencies_are_dispersed_within_each_sub_swarm():
    swarm = HydroMagneticSolitonSwarm(seed=1)
    swarm.reset()
    omegas = np.asarray(swarm.natural_omegas_rad_s, dtype=float)
    assert omegas.shape == (swarm.n_super,)
    for group_id in (0, 1):
        group = omegas[swarm.labels == group_id]
        assert np.std(group) > 0.0, "each sub-swarm needs a non-degenerate frequency spread"


def test_rotational_noise_makes_trajectories_seed_dependent():
    a = _evolve(HydroMagneticSolitonSwarm(seed=11))
    b = _evolve(HydroMagneticSolitonSwarm(seed=12))
    assert not np.allclose(a.phases_rad, b.phases_rad)


def test_synchronization_depends_on_coupling_versus_frequency_spread():
    # Same frequency dispersion, coupling far above vs. far below the critical value
    # K_c ~ sigma_omega * sqrt(8/pi). A small dt keeps the explicit update stable at
    # the large coupling (K * dt << 1).
    def r_at_coupling(multiple_of_sigma: float) -> float:
        swarm = HydroMagneticSolitonSwarm(seed=7)
        swarm.reset()
        sigma = float(np.std(swarm.natural_omegas_rad_s[swarm.labels == 0]))
        swarm.coupling_K = multiple_of_sigma * sigma
        return _evolve(swarm, n_steps=400, dt_s=0.005).alpha_cohesion_r

    r_strong = r_at_coupling(5.0)
    r_weak = r_at_coupling(0.05)
    assert r_strong > r_weak + 0.25, f"no sync transition: r_strong={r_strong}, r_weak={r_weak}"
    assert r_strong > 0.80


def test_uncoupled_swarm_does_not_synchronize():
    free = _evolve(HydroMagneticSolitonSwarm(seed=3), coupled=False)
    assert free.kuramoto_order_r < 0.40
