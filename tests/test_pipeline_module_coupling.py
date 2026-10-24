"""
Tests for M1 (ELN #002-#003): modules inside the closed loop must actually be coupled to the
simulated state, not run on their library defaults.

Two gaps were open. The PINO observer was invoked as
`self.observer.evaluate_blind_window_tracking()` with no arguments, so it always assumed the
default 280 um/s mean flow no matter what the hemodynamics produced or what vessel diameter
the episode was running at. And `HydroMagneticSolitonSwarm()` was constructed without the
diameter override, so a 15 um or 50 um episode still placed the swarm in a 30 um lumen.
"""

import numpy as np
from src.biomimicry_mechanics import HydroMagneticSolitonSwarm
from src.quantum_controller import QDiapedesisClosedLoopPipeline


def test_swarm_lumen_follows_the_episode_diameter():
    for d in (15.0, 30.0, 50.0):
        pipeline = QDiapedesisClosedLoopPipeline(d)
        assert pipeline.swarm.lumen_radius_um == 0.5 * d, (
            f"swarm lumen radius {pipeline.swarm.lumen_radius_um} does not match diameter {d}"
        )


def test_swarm_accepts_an_explicit_lumen_diameter():
    narrow = HydroMagneticSolitonSwarm(lumen_diameter_um=15.0)
    wide = HydroMagneticSolitonSwarm(lumen_diameter_um=50.0)
    assert narrow.lumen_radius_um == 7.5
    assert wide.lumen_radius_um == 25.0
    # Radial seeding sits inside the margination band of the vessel it was told about.
    assert narrow.reset().radial_positions_um.max() < 7.5
    assert wide.reset().radial_positions_um.max() < 25.0


def test_observer_receives_the_simulated_flow_not_a_library_default():
    # Different diameters give genuinely different simulated velocities, so the observer's
    # tracking error must differ between them. It could not before, since the call took no
    # arguments and always assumed 280 um/s.
    errs = {
        d: QDiapedesisClosedLoopPipeline(d)
        .run_episode(total_time_s=1.0, dt_s=0.05)
        .pino_p95_tracking_error_um
        for d in (15.0, 50.0)
    }
    assert errs[15.0] != errs[50.0], f"observer error identical across diameters: {errs}"
    assert all(np.isfinite(v) and v > 0.0 for v in errs.values())
