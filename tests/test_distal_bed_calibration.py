"""
Test for D1 / C1 (ELN #002-#003): the penetrating arteriole drains into a downstream
capillary bed whose series resistance limits flow once the thrombus is cleared.
R_distal_bed is calibrated once, at the 30 um default diameter, so that the open-vessel
peak hydrodynamic wall shear equals the cited physiological value
(Rovainen et al., J. Cereb. Blood Flow Metab. 1993, 13:359-371: 39 +/- 14 dyn/cm^2 = 3.9 Pa).
"""

import numpy as np
import pytest
from src.env_hemodynamics import PenetratingArterioleGeometry, PriesSecombPulsatileFlow


def _peak_hydro_shear_Pa(recanalized_fraction: float, diameter_um: float | None = None) -> float:
    geom = PenetratingArterioleGeometry(n_axial_nodes=120)
    if diameter_um is not None:
        geom.diameter_um = float(diameter_um)
        geom.radius_um = 0.5 * float(diameter_um)
        geom.state = geom._initialize_grid()
    geom.apply_staged_recanalization(recanalized_fraction, buffer_mix_ratio=1.0)
    flow = PriesSecombPulsatileFlow(geometry=geom)
    times = np.linspace(0.0, 60.0 / flow.heart_rate_bpm, 400)
    return max(float(np.max(flow.solve_flow_step(t).wall_shear_stress_Pa)) for t in times)


def test_open_vessel_peak_shear_matches_cited_physiological_value():
    assert _peak_hydro_shear_Pa(1.0) == pytest.approx(3.9, abs=0.1)


def test_calibration_is_fixed_once_and_not_per_diameter():
    # A fixed downstream bed means smaller arterioles carry different shear than larger ones.
    taus = [_peak_hydro_shear_Pa(1.0, d) for d in (15.0, 30.0, 50.0)]
    assert len({round(t, 2) for t in taus}) == 3


def test_occluded_vessel_shear_stays_low():
    assert _peak_hydro_shear_Pa(0.0) < 1.5
