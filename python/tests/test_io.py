"""Tests of the I/O helpers (MATLAB .mat files, CONF_DATA.STR)."""
import numpy as np
import pytest
import scipy.io as sio

from pirt import io

CONF = """-- Resolution of the target
10                   	TargetResolution - Distance between target circles [mm]
223			            TargetXorigin	- X origin target center [pix]
257			            TargetYorigin	- Y origin target center [pix]
--Characteristic length for Nusselt calculation
10e-3					L_char  - [m]
-- Heated Thin Foil parameters
0.01e-3	                s 		- [m] thickness of steel plate
7900	                rho 	- [kg/m^3] density of the material
460 					cp 		- [J/kgK] specific heat
14.7					lambda 	- [W/(mK)] conductive coefficient
0.276					H 		- [m] dimension of the foil
0.1 					W 		- [m] dimension of the foil
0.95					epsilon - [] emisivity of the foil
-- Constants
9.80665 				g 		- [m/s^2] gravity constant
5.67e-8  				sigma   - [W/(m^2 K^4)] Stefan-Boltzmann constant
"""


def test_read_conf_data(tmp_path):
    p = tmp_path / "CONF_DATA.STR"
    p.write_text(CONF)
    conf = io.read_conf_data(p)
    assert conf == {"L": 0.01, "s": 1e-5, "rho": 7900.0, "cp": 460.0, "lambda": 14.7, "H": 0.276,
                    "W": 0.1, "epsilon": 0.95}


def test_load_mat_v5_and_v73(tmp_path, rng):
    A = rng.standard_normal((4, 5, 3))
    sio.savemat(tmp_path / "v5.mat", {"A": A, "s": 2.5})
    d = io.load_mat(tmp_path / "v5.mat")
    np.testing.assert_array_equal(d["A"], A)
    assert d["s"] == 2.5
    h5py = pytest.importorskip("h5py")
    # MATLAB v7.3 stores arrays transposed (column-major); emulate that layout
    with h5py.File(tmp_path / "v73.mat", "w") as f:
        f.create_dataset("A", data=A.T)
        f.create_dataset("s", data=np.array([[2.5]]))
    d73 = io.load_mat(tmp_path / "v73.mat")
    np.testing.assert_array_equal(d73["A"], A)
    assert d73["s"] == 2.5
    only_a = io.load_mat(tmp_path / "v73.mat", variables=["A"], squeeze=False)
    assert set(only_a) == {"A"}


def test_load_sj_case(tmp_path, rng):
    Thot = (300 + rng.standard_normal((6, 7, 5))).astype(np.float32)
    Tcold = 294.3 + rng.standard_normal((6, 7))
    sio.savemat(tmp_path / "Thot.mat", {"Timage_hot": Thot})
    sio.savemat(tmp_path / "Tcold.mat", {"Timage_cold": Tcold})
    sio.savemat(tmp_path / "Resolution.mat", {"dx": 4.40498081, "dy": 4.40498081})
    sio.savemat(tmp_path / "TestConditions.mat", {"V": 2.227, "I": 7.5, "Uinf": 1.233, "Tamb": np.array([[21.6, 21.6]]),
                                                  "dt": 1 / 253})
    (tmp_path / "CONF_DATA.STR").write_text(CONF)
    case = io.load_sj_case(tmp_path, n_snapshots=4, dtype=np.float64)
    assert case["Thot"].shape == (6, 7, 4) and case["Thot"].dtype == np.float64
    np.testing.assert_allclose(case["Conditions"]["Tamb"], [294.75, 294.75])
    assert case["HFS"]["k"] == pytest.approx(14.7e-5)
    assert case["Conditions"]["L"] == 0.01
    assert case["dx"] == pytest.approx(1 / 4.40498081 / 1000)
    assert case["f_acq"] == pytest.approx(253.0)
