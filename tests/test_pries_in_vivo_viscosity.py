"""
Regression test for H1 (ELN #002): the viscosity law must match the published
Pries et al. (1994) in-vivo formulation:

    eta*_0.45 = 6 exp(-0.085 D) + 3.2 - 2.44 exp(-0.06 D^0.645)
    eta_rel   = [1 + (eta*_0.45 - 1) ((1-H)^C - 1)/((1-0.45)^C - 1) z] z,  z = (D/(D-1.1))^2

Reference values below were computed independently for H_D = 0.45 and
plasma viscosity 1.20 mPa.s (configs/vessel_mra_config.yaml).
"""

import numpy as np
import pytest
from src.env_hemodynamics import PenetratingArterioleGeometry, PriesSecombPulsatileFlow

REFERENCE_mPa_s = {10.0: 7.047, 20.0: 3.877, 30.0: 3.026, 50.0: 2.738}


@pytest.mark.parametrize("diameter_um, expected", sorted(REFERENCE_mPa_s.items()))
def test_in_vivo_viscosity_matches_pries_1994(diameter_um, expected):
    flow = PriesSecombPulsatileFlow(geometry=PenetratingArterioleGeometry(n_axial_nodes=10))
    mu = float(flow.compute_in_vivo_viscosity_mPa_s(np.array([diameter_um]))[0])
    assert mu == pytest.approx(expected, rel=1e-3)


def test_in_vivo_viscosity_decreases_toward_large_vessels():
    flow = PriesSecombPulsatileFlow(geometry=PenetratingArterioleGeometry(n_axial_nodes=10))
    mu = flow.compute_in_vivo_viscosity_mPa_s(np.array([10.0, 20.0, 30.0, 50.0]))
    assert np.all(np.diff(mu) < 0.0)
