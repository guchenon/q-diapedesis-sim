"""
Sanity Check & Unit Tests for Q-DIAPEDESIS Preregistered Configurations.
Verifies physical consistency and physiological Control Barrier Function (CBF) bounds
against the immutable OSF protocol (https://osf.io/htswk).
"""

from pathlib import Path
import yaml
import pytest

ROOT_DIR = Path(__file__).resolve().parent.parent
CONFIG_DIR = ROOT_DIR / "configs"


def load_yaml(filename: str) -> dict:
    filepath = CONFIG_DIR / filename
    assert filepath.exists(), f"Missing configuration file: {filepath}"
    with open(filepath, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def test_vessel_mra_config_sanity():
    cfg = load_yaml("vessel_mra_config.yaml")

    # 1. Verify Penetrating Arteriole Geometry (15-50 um)
    geom = cfg["vessel_geometry"]
    d_min, d_max = geom["diameter_um_range"]
    assert d_min == 15.0 and d_max == 50.0
    assert d_min <= geom["default_diameter_um"] <= d_max
    assert geom["glycocalyx_baseline_thickness_um"] > 0.0

    # 2. Verify Cerebral Perfusion Pressure identity: CPP = MAP - ICP
    hemo = cfg["hemodynamics_pries_secomb"]
    expected_cpp = hemo["map_mmHg"] - hemo["icp_mmHg"]
    assert pytest.approx(hemo["cpp_baseline_mmHg"], rel=1e-5) == expected_cpp
    assert 0.30 <= hemo["discharge_hematocrit"] <= 0.60

    # 3. Verify Post-Occlusive Acidosis & Washout Bounds
    thrombus = cfg["thrombus_and_microenvironment"]
    assert thrombus["post_occlusive_pH_initial"] < 6.80
    assert thrombus["physiological_pH_target"] >= 7.30

    bio = cfg["biomimicry_layer"]
    assert bio["hydro_magnetic_soliton"]["washout_loss_ceiling_percent"] <= 2.0


def test_quantum_floquet_config_sanity():
    cfg = load_yaml("quantum_floquet_config.yaml")

    # 1. Verify Sub-Swarm Partition Conservation (alpha + beta = 1.0)
    scaling = cfg["swarm_scaling"]
    total_ratio = scaling["sub_swarm_alpha_ratio"] + scaling["sub_swarm_beta_ratio"]
    assert pytest.approx(total_ratio, rel=1e-6) == 1.0
    assert scaling["n_super_particles"] >= 50

    # 2. Verify Bimodal Frequency Separation (80 Hz Drill vs 5 Hz Trap)
    drive = cfg["bimodal_floquet_drive"]
    assert drive["alpha_drill_frequency_Hz"] == 80.0
    assert drive["beta_trap_frequency_Hz"] == 5.0
    assert drive["alpha_drill_frequency_Hz"] / drive["beta_trap_frequency_Hz"] >= 10.0
    assert drive["target_dual_task_fidelity_min"] >= 0.92


def test_neurocritical_cbf_config_sanity():
    cfg = load_yaml("neurocritical_cbf_config.yaml")

    # 1. Verify Preregistered Physiological Safety Ceilings
    cbf = cfg["cbf_safety_thresholds"]
    assert cbf["max_wall_shear_stress_Pa"] <= 5.50
    assert cbf["min_glycocalyx_integrity_ratio"] >= 0.92
    assert cbf["max_brain_temp_elevation_C"] <= 0.30

    # 2. Verify Staged CPP Recanalization & Post-Task Retrieval
    staged = cfg["staged_cpp_recanalization"]
    assert staged["enable_staged_stepping"] is True
    assert staged["num_recanalization_stages"] >= 4

    retrieval = cfg["post_task_retrieval"]
    assert retrieval["enable_degaussing_pulse"] is True
    assert retrieval["target_retrieval_efficiency_min"] >= 0.99
