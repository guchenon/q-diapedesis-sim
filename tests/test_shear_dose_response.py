"""
Tests for D3 Stage 1 (ELN #003): glycocalyx shedding and NETosis respond continuously to
wall shear (no hard zero below the 5.50 Pa CBF ceiling), monotonically, and negligibly at
the cited physiological arteriolar shear (Rovainen et al. 1993: 3.9 Pa).
"""

import numpy as np
from src.neurocritical_cbf import EndothelialGlycocalyxModel, GlycocalyxNetosisDynamics


def _glyco_loss_rate(tau_Pa: float) -> float:
    m = EndothelialGlycocalyxModel(n_axial_nodes=10)
    s = m.step_spatial_glycocalyx(dt_s=0.01, wall_shear_profile_Pa=np.full(10, tau_Pa))
    return (1.0 - s.mean_integrity_ratio) / 0.01


def _netosis_after(tau_Pa: float, seconds: float = 10.0) -> float:
    d = GlycocalyxNetosisDynamics()
    for _ in range(int(seconds / 0.05)):
        st = d.step_physiology(0.05, np.full(10, tau_Pa), 0.0, 80.0, 500.0, True)
    return st.netosis_activation_index


def test_damage_is_continuous_below_the_cbf_ceiling():
    assert _glyco_loss_rate(5.4) > _glyco_loss_rate(4.0) > 0.0
    assert _netosis_after(5.4) > _netosis_after(4.0) > 0.0


def test_damage_is_monotone_in_shear():
    taus = [1.0, 3.0, 5.0, 7.0, 9.0, 12.0]
    g = [_glyco_loss_rate(t) for t in taus]
    n = [_netosis_after(t, 2.0) for t in taus]
    assert np.all(np.diff(g) > 0.0)
    assert np.all(np.diff(n) > 0.0)


def test_no_step_at_the_cbf_ceiling():
    below, above = _glyco_loss_rate(5.49), _glyco_loss_rate(5.51)
    assert abs(above - below) < 0.02 * max(above, below)
