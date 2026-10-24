"""
Tests for D4 / C4 (ELN #002-#003): ablations must share one nominal command, toggle
exactly one module each, and support equal-time and iso-recanalization comparisons.
"""

import pytest
from src.quantum_controller import QDiapedesisClosedLoopPipeline


def test_cbf_off_baseline_applies_the_same_nominal_command():
    full = QDiapedesisClosedLoopPipeline(30.0).run_episode(total_time_s=1.0, dt_s=0.05)
    off = QDiapedesisClosedLoopPipeline(30.0).run_episode(total_time_s=1.0, dt_s=0.05, enable_cbf_filter=False)
    # No artificial multipliers: the unfiltered field equals that run's nominal command, and
    # the filtered run can only reduce its own nominal command.
    assert off.mean_applied_b_mT == pytest.approx(off.mean_nominal_b_mT, rel=1e-12)
    assert full.mean_applied_b_mT <= full.mean_nominal_b_mT + 1e-9


def test_ph_buffering_switch_is_independent_of_cbf_switch():
    a = QDiapedesisClosedLoopPipeline(30.0).run_episode(
        total_time_s=2.0, dt_s=0.05, enable_cbf_filter=False, enable_pulsed_pH_buffering=True
    )
    b = QDiapedesisClosedLoopPipeline(30.0).run_episode(
        total_time_s=2.0, dt_s=0.05, enable_cbf_filter=False, enable_pulsed_pH_buffering=False
    )
    assert a.final_distal_pH > b.final_distal_pH


def test_iso_recanalization_checkpoint_stops_at_target_fraction():
    s = QDiapedesisClosedLoopPipeline(30.0).run_episode(
        total_time_s=10.0, dt_s=0.05, enable_cbf_filter=False, stop_at_recanalized_fraction=0.30
    )
    assert 0.30 <= s.final_recanalized_fraction < 0.30 + 0.35 * 0.05 + 1e-9
    assert s.elapsed_time_s < 10.0


def test_default_episode_horizon_comes_from_config():
    s = QDiapedesisClosedLoopPipeline(30.0).run_episode(dt_s=0.05)
    assert s.elapsed_time_s == pytest.approx(10.0)
