import numpy as np
import pytest


@pytest.fixture
def rng():
    return np.random.default_rng(12345)


@pytest.fixture
def sj_like_case():
    """Sensor/condition dictionaries of the sweeping-jet experiment (Robledo et al., 2025)."""
    hfs = dict(s=1e-5, rho=7900.0, cp=460.0, k=14.7 * 1e-5, H=0.276, W=0.1, epsilon=0.95, sides=2,
               s_paint=42e-6, cp_paint=3061.5, rho_paint=1261.175, lambda_paint=1.38)
    cond = dict(V=2.227, I=7.5, Uinf=1.233, Tamb=(294.75, 294.75), L=0.01, dt=1 / 253,
                dx=2.2702e-4, dy=2.2702e-4, rhoinf=1.2, cpinf=1005.0)
    return hfs, cond
