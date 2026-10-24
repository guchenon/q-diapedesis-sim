"""
Continuous shear dose-response used by the glycocalyx and NETosis models (D3, ELN #003).

Replaces the earlier hard step at the 5.50 Pa CBF ceiling, under which damage was zero by
construction. Reduced-order phenomenological model (Status B): a Hill function of local
wall shear, h(tau) = tau^n / (tau50^n + tau^n). Parameters are fixed a priori in
configs/neurocritical_cbf_config.yaml (implementation_constants.shear_dose_response):
tau50 = 8.0 Pa is the highest dose tested by Baratchi et al., Nat. Commun. 2024
(doi:10.1038/s41467-024-51211-1), assumed half-maximal; n = 6 satisfies the stated criterion
that the Hill fraction at the physiological 3.9 Pa (Rovainen et al., JCBFM 1993) is <= 1.5%.
These parameters are not constrained by data and outcomes are sensitive to them (see ELN #003).
"""

from pathlib import Path
import numpy as np
import yaml


def load_dose_response_constants(cbf_config_path: Path | str | None = None) -> dict:
    if cbf_config_path is None:
        cbf_config_path = Path(__file__).resolve().parent.parent.parent / "configs" / "neurocritical_cbf_config.yaml"
    with open(cbf_config_path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    return dict(cfg["implementation_constants"]["shear_dose_response"])


def hill_fraction(tau_Pa, tau50_Pa: float, n: float):
    """Continuous, monotone activation fraction in [0, 1) for wall shear `tau_Pa`."""
    tau = np.maximum(np.asarray(tau_Pa, dtype=float), 0.0)
    x = (tau / float(tau50_Pa)) ** float(n)
    return x / (1.0 + x)
