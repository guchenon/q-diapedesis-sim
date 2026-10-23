"""
Stage 1 - Step 06: 3D Cerebral Penetrating Arteriole Geometry, Porous Thrombus Plug,
and Post-Occlusive Acidotic Microenvironment Model.

Part of the Q-DIAPEDESIS Framework (OSF: https://osf.io/htswk | DOI: 10.5281/zenodo.23599020).
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple
import numpy as np
import yaml


@dataclass
class VesselGridState:
    """Container for 1D/3D axial-radial microvascular geometry and chemical state."""
    z_um: np.ndarray                  # Axial coordinate grid [um]
    r_lumen_um: np.ndarray            # Effective open lumen radius [um]
    r_outer_um: np.ndarray            # Baseline vessel radius [um]
    glycocalyx_thickness_um: np.ndarray  # Endothelial glycocalyx thickness [um]
    occlusion_fraction: np.ndarray    # Local thrombus occlusion ratio [0.0 - 1.0]
    darcy_permeability_m2: np.ndarray # Local porous permeability [m^2]
    pH_profile: np.ndarray            # Axial intravascular pH distribution


class PenetratingArterioleGeometry:
    """
    Models a 3D cerebral penetrating arteriole segment (15-50 um diameter) with:
    1. A proximal collateral bifurcation branch (causing particle skimming/washout risk).
    2. A viscoelastic fibrin-NETs thrombus plug governed by Darcy porous permeability.
    3. A distal stagnant acidotic pocket (pH < 6.8) behind the occlusion.
    4. An endothelial glycocalyx protective lining along the vessel wall.
    """

    def __init__(self, config_path: Path | str | None = None, n_axial_nodes: int = 250):
        if config_path is None:
            root_dir = Path(__file__).resolve().parent.parent.parent
            config_path = root_dir / "configs" / "vessel_mra_config.yaml"

        with open(config_path, "r", encoding="utf-8") as f:
            self.cfg = yaml.safe_load(f)

        geom_cfg = self.cfg["vessel_geometry"]
        thrombus_cfg = self.cfg["thrombus_and_microenvironment"]

        self.diameter_um = float(geom_cfg["default_diameter_um"])
        self.radius_um = 0.5 * self.diameter_um
        self.length_um = float(geom_cfg["vessel_length_um"])
        self.bifurcation_angle_deg = float(geom_cfg["bifurcation_angle_deg"])
        self.glyco_baseline_um = float(geom_cfg["glycocalyx_baseline_thickness_um"])

        self.thrombus_center_um = float(thrombus_cfg["thrombus_axial_position_um"])
        self.thrombus_length_um = float(thrombus_cfg["thrombus_length_um"])
        self.initial_occlusion = float(thrombus_cfg["initial_occlusion_fraction"])
        self.darcy_k0_m2 = float(thrombus_cfg["thrombus_porosity_Darcy_m2"])
        self.pH_acidotic = float(thrombus_cfg["post_occlusive_pH_initial"])
        self.pH_physio = float(thrombus_cfg["physiological_pH_target"])

        self.n_nodes = n_axial_nodes
        self.state = self._initialize_grid()

    def _initialize_grid(self) -> VesselGridState:
        """Initializes the axial coordinate mesh, thrombus profile, and acidotic zone."""
        z_um = np.linspace(0.0, self.length_um, self.n_nodes)
        r_outer_um = np.full_like(z_um, self.radius_um)
        glyco_um = np.full_like(z_um, self.glyco_baseline_um)

        # Smooth super-Gaussian spatial profile for the thrombus plug
        half_len = 0.5 * self.thrombus_length_um
        rel_z = (z_um - self.thrombus_center_um) / half_len
        occlusion = self.initial_occlusion * np.exp(-0.5 * (rel_z ** 4))
        occlusion = np.clip(occlusion, 0.0, 0.99)

        # Effective open lumen radius accounting for thrombus and glycocalyx layer
        r_lumen_um = (r_outer_um - glyco_um) * np.sqrt(1.0 - occlusion)

        # Darcy permeability: low inside dense thrombus, high in open lumen
        darcy_m2 = self.darcy_k0_m2 / (occlusion + 1.0e-4)

        # Sigmoidal transition from physiological pH (proximal) to acidotic pH (distal stagnant zone)
        transition_width = 15.0  # [um]
        sigmoid_distal = 1.0 / (1.0 + np.exp(-(z_um - self.thrombus_center_um) / transition_width))
        pH_profile = self.pH_physio - (self.pH_physio - self.pH_acidotic) * sigmoid_distal

        return VesselGridState(
            z_um=z_um,
            r_lumen_um=r_lumen_um,
            r_outer_um=r_outer_um,
            glycocalyx_thickness_um=glyco_um,
            occlusion_fraction=occlusion,
            darcy_permeability_m2=darcy_m2,
            pH_profile=pH_profile,
        )

    def apply_staged_recanalization(self, recanalized_fraction: float, buffer_mix_ratio: float = 0.0) -> None:
        """
        Updates thrombus occlusion and distal pH buffering during micro-channeled recanalization.

        Parameters
        ----------
        recanalized_fraction : float
            Fraction of the thrombus plug cleared by the alpha sub-swarm [0.0 to 1.0].
        buffer_mix_ratio : float
            Degree of fresh arterial blood mixing into the post-occlusive acidotic zone [0.0 to 1.0].
        """
        recanalized_fraction = float(np.clip(recanalized_fraction, 0.0, 1.0))
        buffer_mix_ratio = float(np.clip(buffer_mix_ratio, 0.0, 1.0))

        half_len = 0.5 * self.thrombus_length_um
        rel_z = (self.state.z_um - self.thrombus_center_um) / half_len
        base_profile = np.exp(-0.5 * (rel_z ** 4))

        current_peak_occlusion = self.initial_occlusion * (1.0 - recanalized_fraction)
        self.state.occlusion_fraction = np.clip(current_peak_occlusion * base_profile, 0.0, 0.99)

        self.state.r_lumen_um = (
            (self.state.r_outer_um - self.state.glycocalyx_thickness_um)
            * np.sqrt(1.0 - self.state.occlusion_fraction)
        )
        self.state.darcy_permeability_m2 = self.darcy_k0_m2 / (self.state.occlusion_fraction + 1.0e-4)

        # Update distal pH as pulsed micro-dialysis buffering neutralizes lactic acidosis
        transition_width = 15.0
        sigmoid_distal = 1.0 / (1.0 + np.exp(-(self.state.z_um - self.thrombus_center_um) / transition_width))
        effective_acid_drop = (self.pH_physio - self.pH_acidotic) * (1.0 - buffer_mix_ratio)
        self.state.pH_profile = self.pH_physio - effective_acid_drop * sigmoid_distal

    def sample_3d_cylindrical_mesh(self, n_theta: int = 36) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Generates 3D Cartesian surface coordinates (X, Y, Z) in micrometers for visualization."""
        theta = np.linspace(0.0, 2.0 * np.pi, n_theta)
        Z, Theta = np.meshgrid(self.state.z_um, theta, indexing="ij")
        R = np.tile(self.state.r_lumen_um[:, np.newaxis], (1, n_theta))
        X = R * np.cos(Theta)
        Y = R * np.sin(Theta)
        return X, Y, Z

    def summary_metrics(self) -> Dict[str, float]:
        """Returns key geometric and microenvironmental summary metrics."""
        return {
            "vessel_diameter_um": self.diameter_um,
            "min_open_lumen_radius_um": float(np.min(self.state.r_lumen_um)),
            "peak_occlusion_fraction": float(np.max(self.state.occlusion_fraction)),
            "distal_min_pH": float(np.min(self.state.pH_profile)),
            "proximal_max_pH": float(np.max(self.state.pH_profile)),
            "mean_glycocalyx_um": float(np.mean(self.state.glycocalyx_thickness_um)),
        }
