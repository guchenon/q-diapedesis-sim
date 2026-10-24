"""
Stage 1 - Step 08: Endothelial Glycocalyx (Endothelial Surface Layer) Structural Integrity
and Hydraulic Retardation Model.

Grounded in:
- Weinbaum, Tarbell & Damiano, Annu. Rev. Biomed. Eng. 2007 (doi:10.1146/annurev.bioeng.9.060906.151959)
- Pries & Secomb, Am. J. Physiol. Heart Circ. Physiol. 2005 (doi:10.1152/ajpheart.00297.2005)
- Rehm et al., Circulation 2007 (doi:10.1161/CIRCULATIONAHA.106.684852)
- Wang G et al., ATVB 2020 (doi:10.1161/ATVBAHA.119.313399)

Part of the Q-DIAPEDESIS Framework (OSF: https://osf.io/htswk | DOI: 10.5281/zenodo.23599020).
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict
import numpy as np
import yaml


@dataclass
class GlycocalyxLayerState:
    """Spatial and scalar metrics of the modeled endothelial surface layer (ESL)."""
    time_s: float
    baseline_thickness_um: float
    axial_thickness_um: np.ndarray
    integrity_ratio_profile: np.ndarray
    mean_integrity_ratio: float
    min_integrity_ratio: float
    shedding_biomarker_proxy: float
    passes_cbf_threshold: bool


class EndothelialGlycocalyxModel:
    """
    Models the spatial structural thickness L_g(z, t) of the microvascular endothelial
    glycocalyx brush along a cerebral penetrating arteriole (15-50 um).

    Epistemic note (Status B design target):
    `mean_integrity_ratio` represents the modeled retained structural thickness fraction
    L_g(t) / L_{g,0} (target >= 0.92), kept explicitly distinct from `shedding_biomarker_proxy`
    (which tracks cumulative syndecan-1 / heparan sulfate release into plasma).
    Supports baseline thickness sensitivity sweeps across [0.40, 0.50, 1.00] um.
    """

    def __init__(
        self,
        n_axial_nodes: int = 250,
        baseline_thickness_override_um: float | None = None,
        cbf_config_path: Path | str | None = None,
        vessel_config_path: Path | str | None = None,
    ):
        root_dir = Path(__file__).resolve().parent.parent.parent
        if cbf_config_path is None:
            cbf_config_path = root_dir / "configs" / "neurocritical_cbf_config.yaml"
        if vessel_config_path is None:
            vessel_config_path = root_dir / "configs" / "vessel_mra_config.yaml"

        with open(cbf_config_path, "r", encoding="utf-8") as f:
            cbf_cfg = yaml.safe_load(f)
        with open(vessel_config_path, "r", encoding="utf-8") as f:
            vessel_cfg = yaml.safe_load(f)

        default_l0 = float(vessel_cfg["vessel_geometry"]["glycocalyx_baseline_thickness_um"])
        self.baseline_um = (
            float(baseline_thickness_override_um)
            if baseline_thickness_override_um is not None
            else default_l0
        )

        thresh = cbf_cfg["cbf_safety_thresholds"]
        self.tau_crit_Pa = float(thresh["max_wall_shear_stress_Pa"])
        self.min_integrity_target = float(thresh["min_glycocalyx_integrity_ratio"])

        self.n_nodes = int(n_axial_nodes)
        self.k_mech_strip_inv_Pa_s = 0.045
        self.k_oscillatory_fatigue_inv_s = 0.0015
        self.k_enzymatic_acid_inv_s = 0.008

        self.reset()

    def reset(self) -> GlycocalyxLayerState:
        """Resets the spatial glycocalyx thickness profile to intact baseline L_{g,0}."""
        self.current_time_s = 0.0
        self.integrity_profile = np.ones(self.n_nodes, dtype=float)
        self.cumulative_shedding_proxy = 0.0
        return self.get_state()

    def step_spatial_glycocalyx(
        self,
        dt_s: float,
        wall_shear_profile_Pa: np.ndarray,
        pH_profile: np.ndarray | None = None,
    ) -> GlycocalyxLayerState:
        """
        Updates the spatial endothelial glycocalyx profile over timestep `dt_s`
        in response to local total wall shear stress and post-occlusive acidosis.
        """
        tau_Pa = np.asarray(wall_shear_profile_Pa, dtype=float)
        if tau_Pa.shape[0] != self.n_nodes:
            tau_Pa = np.resize(tau_Pa, self.n_nodes)

        # 1. Supra-threshold mechanical stripping when local tau_wall > tau_crit (5.50 Pa)
        excess_shear = np.maximum(0.0, tau_Pa - self.tau_crit_Pa)
        mech_rate = self.k_mech_strip_inv_Pa_s * (excess_shear ** 1.5)

        # 2. Sub-threshold oscillatory shear wear (Wang G et al., ATVB 2020)
        osc_rate = self.k_oscillatory_fatigue_inv_s * ((tau_Pa / self.tau_crit_Pa) ** 2)

        # 3. Acidosis-potentiated heparanase / MMP enzymatic shedding (pH < 6.80)
        if pH_profile is not None:
            pH_arr = np.asarray(pH_profile, dtype=float)
            acid_deficit = np.maximum(0.0, 6.80 - pH_arr)
            enz_rate = self.k_enzymatic_acid_inv_s * acid_deficit
        else:
            enz_rate = np.zeros_like(tau_Pa)

        total_degradation_rate = mech_rate + osc_rate + enz_rate
        delta_integrity = -total_degradation_rate * self.integrity_profile * float(dt_s)

        shed_mass_fraction = float(np.mean(-delta_integrity))
        self.cumulative_shedding_proxy += max(0.0, shed_mass_fraction)

        self.integrity_profile = np.clip(self.integrity_profile + delta_integrity, 0.0, 1.0)
        self.current_time_s += float(dt_s)
        return self.get_state()

    def get_state(self) -> GlycocalyxLayerState:
        """Returns the current spatial and summary glycocalyx integrity state."""
        thickness_um = self. baseline_um * self.integrity_profile
        mean_ratio = float(np.mean(self.integrity_profile))
        min_ratio = float(np.min(self.integrity_profile))
        return GlycocalyxLayerState(
            time_s=self.current_time_s,
            baseline_thickness_um=self.baseline_um,
            axial_thickness_um=thickness_um,
            integrity_ratio_profile=self.integrity_profile.copy(),
            mean_integrity_ratio=mean_ratio,
            min_integrity_ratio=min_ratio,
            shedding_biomarker_proxy=self.cumulative_shedding_proxy,
            passes_cbf_threshold=bool(mean_ratio >= self.min_integrity_target),
        )

    def run_baseline_sensitivity_sweep(
        self,
        wall_shear_profile_Pa: np.ndarray,
        duration_s: float = 2.0,
        dt_s: float = 0.01,
    ) -> Dict[str, float]:
        """
        Evaluates retained integrity across the preregistered baseline thickness
        sensitivity sweep [0.40, 0.50, 1.00] um (Audit Table Row 2b).
        """
        results: Dict[str, float] = {}
        steps = int(duration_s / dt_s)
        for l0_um in (0.40, 0.50, 1.00):
            model = EndothelialGlycocalyxModel(
                n_axial_nodes=self.n_nodes,
                baseline_thickness_override_um=l0_um,
            )
            for _ in range(steps):
                state = model.step_spatial_glycocalyx(dt_s, wall_shear_profile_Pa)
            results[f"retained_um_at_{l0_um:.2f}um_baseline"] = float(
                np.mean(state.axial_thickness_um)
            )
            results[f"integrity_ratio_at_{l0_um:.2f}um_baseline"] = state.mean_integrity_ratio
        return results
