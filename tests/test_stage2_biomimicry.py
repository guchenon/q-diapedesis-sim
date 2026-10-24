"""
End-to-End Integration & Sensitivity Tests for Stage 2 (Steps 11-15):
Leukocyte Catch-Bond Margination, Sakaguchi-Kuramoto Soliton Synchronization,
and Dual Sub-Swarm Partitioning (< 2.0% Transient Bifurcation Washout Target).
"""

import numpy as np
from src.biomimicry_mechanics import (
    LeukocyteCatchBondMargination,
    HydroMagneticSolitonSwarm,
)


def test_step11_and_step12_margination_and_catch_bond_sweep():
    margination = LeukocyteCatchBondMargination()
    rng = np.random.default_rng(42)

    # Initialize particles near the vessel core (radius 2-8 um inside a 15 um radius lumen)
    r_init = rng.uniform(2.0, 8.0, size=200)
    r_curr = r_init.copy()

    # Step radial margination for 1.5 seconds under 3.5 Pa wall shear stress
    for _ in range(30):
        r_curr = margination.step_radial_margination(
            radial_positions_um=r_curr,
            lumen_radius_um=15.0,
            wall_shear_stress_Pa=3.5,
            dt_s=0.05,
        )

    # Verify outward migration toward the endothelial wall
    assert float(np.mean(r_curr)) > float(np.mean(r_init)) + 3.0

    roll_systole = margination.evaluate_rolling_step(
        radial_positions_um=r_curr,
        lumen_radius_um=15.0,
        wall_shear_stress_Pa=3.5,
        is_systole=True,
        magnetic_drive_freq_Hz=80.0,
    )
    assert roll_systole.wall_proximal_fraction >= 0.80
    assert roll_systole.adhesion_locked is True

    # Verify Audit Table Row 3b sensitivity sweep across F* in {10, 15, 20} pN vs pure slip
    sweep = margination.run_catch_force_sensitivity_sweep(shear_stress_Pa=4.2)
    assert sweep["lifetime_s_at_Fstar_15pN"] > sweep["lifetime_s_pure_slip_baseline"]
    assert sweep["lifetime_s_at_Fstar_10pN"] > 0.0
    assert sweep["lifetime_s_at_Fstar_20pN"] > 0.0


def test_step13_and_step14_soliton_synchronization_and_washout():
    margination = LeukocyteCatchBondMargination()
    swarm = HydroMagneticSolitonSwarm()

    # Verify 60/40 sub-swarm partitioning (120 alpha drills, 80 beta mesh agents)
    init_state = swarm.reset()
    assert init_state.n_super_particles == 200
    assert init_state.n_alpha_drill == 120
    assert init_state.n_beta_mesh == 80

    # Evolve coupled margination + Sakaguchi-Kuramoto soliton synchronization
    r_coords = init_state.radial_positions_um.copy()
    for step in range(40):
        is_systole = (step % 10) < 5
        r_coords = margination.step_radial_margination(
            radial_positions_um=r_coords,
            lumen_radius_um=15.0,
            wall_shear_stress_Pa=3.8,
            dt_s=0.05,
        )
        roll_state = margination.evaluate_rolling_step(
            radial_positions_um=r_coords,
            lumen_radius_um=15.0,
            wall_shear_stress_Pa=3.8,
            is_systole=is_systole,
            magnetic_drive_freq_Hz=80.0,
        )
        swarm_state = swarm.step_soliton_swarm(
            dt_s=0.05,
            axial_rolling_velocity_um_s=roll_state.axial_rolling_velocity_um_s,
            wall_proximal_fraction=roll_state.wall_proximal_fraction,
            adhesion_locked=roll_state.adhesion_locked,
            collateral_skimming_fraction=0.08,
            enable_soliton_coupling=True,
        )

    # Verify Kuramoto phase synchronization increases and transient washout stays < 2.0%
    assert swarm_state.kuramoto_order_r > init_state.kuramoto_order_r
    assert swarm_state.kuramoto_order_r >= 0.85
    assert swarm_state.beta_mesh_coverage_ratio >= 0.90
    assert swarm_state.transient_washout_percent < 2.0
    assert swarm_state.passes_washout_ceiling is True
