"""Tests of the energy balance, derivatives, uncertainty and the PIRT class."""
import numpy as np
import pytest

from pirt import PIRT
from pirt.derivatives import derivative_fd
from pirt.heat_transfer import (STEFAN_BOLTZMANN, air_thermal_conductivity, calculate_heat_transfer,
                                prepare_conditions, prepare_hfs)
from pirt.uncertainty import montecarlo_uncertainty


def _expected_h(hfs, cond, Thot, Tcold):
    qj = cond["V"] * cond["I"] / (hfs["H"] * hfs["W"])
    qr = hfs["sides"] * STEFAN_BOLTZMANN * hfs["epsilon"] * (Thot**4 - cond["Tamb"][1] ** 4)
    return (qj - qr) / (Thot - Tcold * cond["Tamb"][1] / cond["Tamb"][0])


def test_derivative_fd_polynomial_fields():
    dx, dy, dt = 0.2, 0.3, 0.05
    y = np.arange(12)[:, None, None] * dy
    x = np.arange(15)[None, :, None] * dx
    t = np.arange(8)[None, None, :] * dt
    T = 1 + 4 * x**2 - 2 * y**2 + 0.5 * t + 3 * t**2 + x * t
    d2x, d2y, dTdt = derivative_fd(T, spatial=True, temporal=True, dx=dx, dy=dy, dt=dt)
    np.testing.assert_allclose(d2x, 8.0, atol=1e-9)   # x along columns
    np.testing.assert_allclose(d2y, -4.0, atol=1e-9)  # y along rows
    np.testing.assert_allclose(dTdt, np.broadcast_to(0.5 + 6 * t + x, T.shape), atol=1e-9)
    d2x2, d2y2, none = derivative_fd(T[:, :, 0], spatial=True, dx=dx, dy=dy)
    assert none is None and d2x2.shape == T.shape[:2]
    with pytest.raises(ValueError):
        derivative_fd(T, temporal=True)  # dt missing
    with pytest.raises(ValueError):
        derivative_fd(T[:, :, 0], temporal=True, dt=dt)


def test_calculate_heat_transfer_uniform_case(sj_like_case):
    hfs, cond = sj_like_case
    Thot = np.full((5, 6, 4), 300.0)
    Tcold = np.full((5, 6), 294.35)
    res = calculate_heat_transfer(Thot, Tcold, hfs, cond, compute=("h", "Nu", "St"), time_der=True,
                                  spatial_der=True, verbose=False, return_terms=True)
    h_exp = _expected_h(hfs, cond, 300.0, 294.35)
    np.testing.assert_allclose(res["h"], h_exp, rtol=1e-12)
    k_air = air_thermal_conductivity((300.0 + cond["Tamb"][1]) / 2)
    np.testing.assert_allclose(res["Nu"], h_exp * cond["L"] / k_air, rtol=1e-12)
    np.testing.assert_allclose(res["St"], h_exp / (cond["rhoinf"] * cond["cpinf"] * cond["Uinf"]), rtol=1e-12)
    # uniform field: derivative terms vanish
    np.testing.assert_allclose(res["terms"]["q_unsteady"], 0.0, atol=1e-9)
    np.testing.assert_allclose(res["terms"]["q_tangential"], 0.0, atol=1e-9)
    assert res["dTdt_hot"].shape == Thot.shape


def test_calculate_heat_transfer_celsius_custom_q_and_3d_tcold(sj_like_case):
    hfs, cond = sj_like_case
    Thot = np.full((5, 6), 300.0) - 273.15
    Tcold3 = np.full((5, 6, 3), 294.35) - 273.15
    Tcold3[:, :, 0] += 0.3
    Tcold3[:, :, 1] -= 0.3
    res = calculate_heat_transfer(Thot, Tcold3, hfs, cond, compute="h", verbose=False,
                                  custom_q=[100.0, np.full((5, 6), 50.0)])
    qj = cond["V"] * cond["I"] / (hfs["H"] * hfs["W"])
    qr = 2 * STEFAN_BOLTZMANN * 0.95 * (300.0**4 - 294.75**4)
    np.testing.assert_allclose(res["h"], (qj - qr + 150.0) / (300.0 - 294.35), rtol=1e-12)
    with pytest.raises(ValueError):
        calculate_heat_transfer(Thot, Tcold3, hfs, cond, compute="h", verbose=False, custom_q=[np.ones((3, 3))])


def test_unsteady_and_tangential_terms_with_given_derivatives(sj_like_case):
    hfs, cond = sj_like_case
    Thot = np.full((4, 5, 3), 300.0)
    Tcold = np.full((4, 5), 294.35)
    dTdt = np.full(Thot.shape, 2.0)
    d2x = np.full(Thot.shape, 1e4)
    d2y = np.full(Thot.shape, -3e3)
    res = calculate_heat_transfer(Thot, Tcold, hfs, cond, compute="h", time_der=True, spatial_der=True,
                                  dTdt=dTdt, d2Tdx2=d2x, d2Tdy2=d2y, verbose=False)
    qj = cond["V"] * cond["I"] / (hfs["H"] * hfs["W"])
    qr = 2 * STEFAN_BOLTZMANN * 0.95 * (300.0**4 - 294.75**4)
    qu = (hfs["rho"] * hfs["s"] * hfs["cp"] + hfs["rho_paint"] * hfs["cp_paint"] * hfs["s_paint"]) * 2.0
    qk = (hfs["k"] + hfs["s_paint"] * hfs["lambda_paint"]) * (1e4 - 3e3)
    np.testing.assert_allclose(res["h"], (qj - qr - qu - qk) / (300.0 - 294.35), rtol=1e-12)
    # PCB anisotropic model
    hfs_pcb = dict(hfs, Type="PCB", lambdax=1.6, lambday=3.1)
    hfs_pcb.pop("k")
    res_pcb = calculate_heat_transfer(Thot, Tcold, hfs_pcb, cond, compute="h", spatial_der=True,
                                      d2Tdx2=d2x, d2Tdy2=d2y, verbose=False)
    kp = hfs["s_paint"] * hfs["lambda_paint"]
    qk_pcb = (hfs["s"] * 1.6 + kp) * 1e4 + (hfs["s"] * 3.1 + kp) * (-3e3)
    np.testing.assert_allclose(res_pcb["h"], (qj - qr - qk_pcb) / (300.0 - 294.35), rtol=1e-12)


def test_input_validation(sj_like_case):
    hfs, cond = sj_like_case
    with pytest.raises(ValueError):
        prepare_hfs({k: v for k, v in hfs.items() if k != "epsilon"})
    with pytest.raises(ValueError):
        prepare_hfs(dict(hfs, Type="Other"))
    assert prepare_hfs(dict(hfs, Area=0.05))["A"] == 0.05
    assert prepare_hfs({k: v for k, v in hfs.items() if k not in ("H", "W")} | {"Area": 0.02})["A"] == 0.02
    with pytest.raises(ValueError):
        prepare_conditions({k: v for k, v in cond.items() if k != "Tamb"})
    with pytest.raises(ValueError):
        prepare_conditions(dict(cond, Tamb=(300.0,)))
    with pytest.raises(ValueError):
        prepare_conditions({k: v for k, v in cond.items() if k != "cpinf"}, compute_st=True)
    with pytest.warns(UserWarning):
        prepare_hfs({k: v for k, v in hfs.items() if k != "sides"})


def test_montecarlo_zero_uncertainty_recovers_deterministic(sj_like_case, rng):
    hfs, cond = sj_like_case
    Thot = 300.0 + 0.1 * rng.standard_normal((6, 7, 5))
    Tcold = np.full((6, 7), 294.35)
    ref = calculate_heat_transfer(Thot, Tcold, hfs, cond, compute=("h", "Nu"), verbose=False)
    error = {k: 0.0 for k in ("errorT", "errorTamb", "errorV", "errorI", "errorEpsilon", "errorrho", "errorcp",
                              "errors", "errorA", "errorkplate", "errorLchar", "errork", "errors_paint",
                              "errorcp_paint", "errorlambda_paint", "errorrho_paint")}
    error["samples"] = 5
    res = montecarlo_uncertainty(Thot, Tcold, hfs, cond, error, compute=("h", "Nu"), reference=ref,
                                 seed=0, verbose=False)
    # mean mode: h computed from the time-averaged Thot, which is not exactly the time-average of h,
    # but with zero uncertainty all samples are identical
    assert np.ptp(res["errorh"]) == 0.0 and np.ptp(res["errorNu"]) == 0.0
    assert res["errorNu"].shape == (5,)
    assert np.all(res["errorNu_p"] < 0.5)  # h(<T>) vs <h(T)>: second order in the fluctuations
    res_s = montecarlo_uncertainty(Thot, Tcold, hfs, cond, error, compute=("Nu",), reference=ref,
                                   mode="snapshots", seed=0, verbose=False)
    np.testing.assert_allclose(res_s["errorNu"], ref["Nu"].mean(), rtol=1e-12)
    assert np.allclose(res_s["errorNu_p"], 0.0)
    assert "errorNuf" in res_s and np.allclose(res_s["errorNuf_p"], 0.0)


def test_montecarlo_spread_scales_with_uncertainty(sj_like_case):
    hfs, cond = sj_like_case
    Thot = np.full((6, 7, 5), 300.0)
    Tcold = np.full((6, 7), 294.35)
    ref = calculate_heat_transfer(Thot, Tcold, hfs, cond, compute="h", verbose=False)
    base = {k: 0.0 for k in ("errorT", "errorTamb", "errorV", "errorI", "errorEpsilon", "errorrho", "errorcp",
                             "errors", "errorA", "errorkplate", "errorLchar", "errork", "errors_paint",
                             "errorcp_paint", "errorlambda_paint", "errorrho_paint")}
    e1 = dict(base, errorV=0.01, samples=400)
    e2 = dict(base, errorV=0.02, samples=400)
    r1 = montecarlo_uncertainty(Thot, Tcold, hfs, cond, e1, compute="h", reference=ref, seed=1, verbose=False)
    r2 = montecarlo_uncertainty(Thot, Tcold, hfs, cond, e2, compute="h", reference=ref, seed=1, verbose=False)
    # h is linear in V here, so the relative spread must double
    ratio = np.std(r2["errorh"]) / np.std(r1["errorh"])
    assert ratio == pytest.approx(2.0, rel=0.05)
    assert np.std(r1["errorh"]) / np.mean(r1["errorh"]) == pytest.approx(0.01 * ref["h"].mean() / ref["h"].mean() * (1 + 0), rel=0.2)
    with pytest.raises(ValueError):
        montecarlo_uncertainty(Thot, Tcold, hfs, cond, {"errorT": 0.1}, compute="h", verbose=False)


# --------------------------------------------------------------------------- #
# PIRT class
# --------------------------------------------------------------------------- #
def _synthetic_case(rng, nt=40):
    ny, nx = 24, 30
    t = np.arange(nt) / 253.0
    yy, xx = np.mgrid[0:ny, 0:nx]
    base = 300 + 3 * np.exp(-((xx - nx / 2) ** 2 + (yy - ny / 2) ** 2) / 40.0)
    Thot = base[:, :, None] + 0.5 * np.sin(2 * np.pi * 20 * t)[None, None, :] + 0.05 * rng.standard_normal((ny, nx, nt))
    Tcold = 294.35 + 0.02 * rng.standard_normal((ny, nx, 10))
    return Thot, Tcold


def test_pirt_class_full_pipeline(sj_like_case, rng, tmp_path):
    hfs, cond = sj_like_case
    Thot, Tcold = _synthetic_case(rng)
    dx, dy, dt = cond["dx"], cond["dy"], cond["dt"]
    filters = [
        {"Type": "POD", "Parameters": {"Nmod": 5}},                         # MATLAB-like spec
        {"type": "gaussian", "FilterSize": (5, 3, 1), "Sigma": (1.5, 1.5, 0.1)},  # mixed spelling
        {"type": "sgolay32", "kernel_size": (5, 5, 3), "h": (dx, dy, dt)},
    ]
    obj = PIRT(Thot=Thot, Tcold=Tcold, filters=filters, heat_transfer=("Nu", "h"), time_der=True,
               spatial_der=True, HFS=hfs, conditions=cond, verbose=False)
    obj.go()
    r = obj.result
    assert r["Nmod_hot"] == [5, None, None]
    assert r["Thot_new"].shape == (24, 30, 38)  # sgolay crops one snapshot at each end
    assert r["Nu"].shape == (24, 30, 38) and r["h"].shape == (24, 30, 38)
    assert np.isfinite(r["Nu"]).all() and r["Nu"].mean() > 0
    for key in ("dTdt_hot", "d2Tdx2_hot", "d2Tdy2_hot"):
        assert r[key].shape == (24, 30, 38)

    # same computation through the functional API
    from pirt.filters import gaussian_filter3, pod_filter, sgolay32_filter
    T1 = pod_filter(Thot, nmod=5, verbose=False)
    T2 = gaussian_filter3(T1, (1.5, 1.5, 0.1), (5, 3, 1))
    T3, dTdt, d2x, d2y = sgolay32_filter(T2, (5, 5, 3), (dx, dy, dt))
    ref = calculate_heat_transfer(T3, Tcold.mean(axis=2), hfs, cond, compute=("Nu", "h"), time_der=True,
                                  spatial_der=True, dTdt=dTdt, d2Tdx2=d2x, d2Tdy2=d2y, verbose=False)
    np.testing.assert_allclose(r["Nu"], ref["Nu"], rtol=1e-12)

    # output offloaded to disk
    obj2 = PIRT(Thot=Thot, Tcold=Tcold, filters=filters, heat_transfer="Nu", time_der=True,
                HFS=hfs, conditions=cond, verbose=False, output_dir=tmp_path)
    obj2.go()
    assert (tmp_path / "Nu.npy").exists()
    np.testing.assert_allclose(np.asarray(obj2.result["Nu"]), np.load(tmp_path / "Nu.npy"))


def test_pirt_class_crop_and_fd_derivatives(sj_like_case, rng):
    hfs, cond = sj_like_case
    Thot, Tcold = _synthetic_case(rng, nt=12)
    obj = PIRT(Thot=Thot, Tcold=Tcold, crop=((2, 20), (3, 21)), heat_transfer="Nu", time_der=True,
               spatial_der=True, HFS=hfs, conditions=cond, verbose=False)
    obj.go()
    assert obj.result["Nu"].shape == (18, 18, 12)
    ref = calculate_heat_transfer(Thot[3:21, 2:20], Tcold.mean(axis=2)[3:21, 2:20], hfs, cond, compute="Nu",
                                  time_der=True, spatial_der=True, verbose=False)
    np.testing.assert_allclose(obj.result["Nu"], ref["Nu"], rtol=1e-12)
    assert "dTdt_hot" in obj.result and "d2Tdx2_hot" in obj.result


def test_pirt_class_errors_and_aliases(sj_like_case, rng):
    hfs, cond = sj_like_case
    Thot, Tcold = _synthetic_case(rng, nt=8)
    with pytest.raises(ValueError):
        PIRT(verbose=False)
    with pytest.raises(ValueError):
        PIRT(Thot=Thot, heat_transfer="Nu", HFS=hfs, conditions=cond, verbose=False)  # Tcold missing
    with pytest.raises(ValueError):
        PIRT(Thot=Thot, Tcold=Tcold, heat_transfer="Bad", HFS=hfs, conditions=cond, verbose=False)
    with pytest.raises(ValueError):
        PIRT(Thot=Thot, filters=[{"type": "POD", "nonsense": 1}], verbose=False)
    with pytest.raises(ValueError):
        PIRT(Thot=Thot, filters=[{"type": "unknown"}], verbose=False)
    # MATLAB-style keyword aliases
    obj = PIRT(Thot=Thot, Tcold=Tcold, CalculateHeatTransfer=("h",), TimeDer=True, HFS=hfs,
               Conditions=cond, Filter={"Type": "wiener3", "Parameters": {"kernel": (3, 3, 1)}}, verbose=False)
    obj.go()
    assert obj.result["h"].shape == Thot.shape
    assert "noise_hot" in obj.result


def test_pirt_class_montecarlo(sj_like_case, rng):
    hfs, cond = sj_like_case
    Thot, Tcold = _synthetic_case(rng, nt=6)
    error = dict(errorT=0.1, errorTamb=0.1, errorV=0.005, errorI=0.005, errorEpsilon=0.01, errorrho=0.01,
                 errorcp=0.02, errors=0.01, errorA=0.001, errorkplate=0.01, errorLchar=0.01, errork=0.01,
                 errors_paint=0.05, errorcp_paint=0.1, errorlambda_paint=0.0612, errorrho_paint=0.0325,
                 errordTdt=0.1, samples=20)
    obj = PIRT(Thot=Thot, Tcold=Tcold, heat_transfer=("Nu", "h"), time_der=True, HFS=hfs, conditions=cond,
               error=error, error_method="montecarlo", seed=3, verbose=False)
    obj.go()
    assert obj.result["errorNu"].shape == (20,) and obj.result["errorNu_p"].shape == (20,)
    assert 0 < obj.result["errorNu_p"].mean() < 30
    with pytest.raises(ValueError):
        PIRT(Thot=Thot, Tcold=Tcold, heat_transfer="Nu", HFS=hfs, conditions=cond, error_method="montecarlo",
             verbose=False)  # Error data missing


def test_film_temperature_option(sj_like_case):
    hfs, cond = sj_like_case
    Thot = np.full((4, 5, 3), 300.0)
    Tcold = np.full((4, 5), 294.35)
    amb = calculate_heat_transfer(Thot, Tcold, hfs, cond, compute="Nu", verbose=False)
    adi = calculate_heat_transfer(Thot, Tcold, hfs, cond, compute="Nu", film_temperature="adiabatic", verbose=False)
    k_amb = air_thermal_conductivity((300.0 + 294.75) / 2)
    k_adi = air_thermal_conductivity((300.0 + 294.35) / 2)
    np.testing.assert_allclose(adi["Nu"] / amb["Nu"], k_amb / k_adi, rtol=1e-12)
    with pytest.raises(ValueError):
        calculate_heat_transfer(Thot, Tcold, hfs, cond, compute="Nu", film_temperature="bad", verbose=False)
