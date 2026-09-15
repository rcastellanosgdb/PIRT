#!/usr/bin/env python
"""Monte Carlo uncertainty check against ``error_metrics_1000_f.mat`` (Robledo et al., 2025).

The paper's Table 2 values (5.7 % for the mean Nu, 7.8 % for the fluctuating
component) were obtained with ``SJ_main_12_error.m`` and the *paper version* of
``Calculate_HeatTransfer_Error.m`` (``PIRT_paper_SJ_2025/src``), which perturbs
the full snapshot sequence (1000 samples). That MATLAB implementation differs
from the deterministic balance in three points (all corrected in v1.1):

* the radiative term omits the number of exposed sides (``HFS.sides = 2``),
* the y-contribution of the tangential term has the wrong sign,
* the air conductivity is evaluated at ``(T_w + T_cold)/2``.

This script runs (a) the corrected Python Monte Carlo (``mode='snapshots'``) and
(b) an exact replica of the paper's MATLAB sampling formulas, on the same
filtered data (POD + Gaussian + Savitzky-Golay, 'Nu' and 'h' with both
derivative terms), and compares the statistics of the samples with the file.

    python montecarlo_check.py --data ... --ref ... [--samples 60] [--out output/montecarlo]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np

warnings.simplefilter("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from pirt import io  # noqa: E402
from pirt.filters import gaussian_filter3, pod_filter, sgolay32_filter  # noqa: E402
from pirt.heat_transfer import (STEFAN_BOLTZMANN, air_thermal_conductivity,  # noqa: E402
                                calculate_heat_transfer)
from pirt.uncertainty import montecarlo_uncertainty  # noqa: E402

ERROR = dict(errorT=0.1, errorTamb=0.1, errorV=0.005, errorI=0.005, errorEpsilon=0.01, errorrho=0.01,
             errorcp=0.02, errors=0.01, errorA=0.001, errorkplate=0.01, errorLchar=0.01, errork=0.01,
             errors_paint=0.05, errorcp_paint=0.1, errorlambda_paint=0.0612, errorrho_paint=0.0325,
             errordTdt=0.1, errord2Tdx2=0.1, errord2Tdy2=0.1)


def paper_montecarlo(Thot, Tcold, dTdt, d2x, d2y, hfs, cond, n, seed, Nu_ref):
    """Replica of the sampling loop of the paper version of Calculate_HeatTransfer_Error.m."""
    rng = np.random.default_rng(seed)
    e = ERROR
    s, rho, cp, A, eps = hfs["s"], hfs["rho"], hfs["cp"], hfs["H"] * hfs["W"], hfs["epsilon"]
    kplate = hfs["k"]
    s_p, rho_p, cp_p, lam_p = hfs["s_paint"], hfs["rho_paint"], hfs["cp_paint"], hfs["lambda_paint"]
    L, V, I = cond["L"], cond["V"], cond["I"]
    Tamb_cold, Tamb_hot = cond["Tamb"]
    sigma = STEFAN_BOLTZMANN
    k = air_thermal_conductivity((Thot + Tcold[:, :, None]) / 2.0)  # paper version: (Thot+Tcold)/2
    Num_ref = Nu_ref.mean(axis=2)
    Nuf_ref = np.abs(Nu_ref - Num_ref[:, :, None]).mean()
    out = dict(errorNu=np.empty(n), errorNu_p=np.empty(n), errorNuf=np.empty(n), errorNuf_p=np.empty(n))
    nrm = lambda mu, sig: rng.normal(mu, np.abs(sig))
    for i in range(n):
        t0 = time.time()
        c_Thot = nrm(Thot, e["errorT"]); c_Tcold = nrm(Tcold, e["errorT"])
        c_Tc, c_Th = nrm(Tamb_cold, e["errorTamb"]), nrm(Tamb_hot, e["errorTamb"])
        c_V, c_I = nrm(V, V * e["errorV"]), nrm(I, I * e["errorI"])
        c_eps, c_rho, c_cp, c_s = nrm(eps, eps * e["errorEpsilon"]), nrm(rho, rho * e["errorrho"]), nrm(cp, cp * e["errorcp"]), nrm(s, s * e["errors"])
        c_kx, c_ky = nrm(kplate, kplate * e["errorkplate"]), nrm(kplate, kplate * e["errorkplate"])
        c_A, c_L = nrm(A, A * e["errorA"]), nrm(L, L * e["errorLchar"])
        c_k = nrm(k, k * e["errork"])
        c_cp_p, c_s_p = nrm(cp_p, cp_p * e["errorcp_paint"]), nrm(s_p, s_p * e["errors_paint"])
        c_rho_p, c_lam_p = nrm(rho_p, rho_p * e["errorrho_paint"]), nrm(lam_p, lam_p * e["errorlambda_paint"])
        c_kp = c_lam_p * c_s_p
        ratio = c_Th / c_Tc
        q = c_Thot**4
        q -= c_Th**4
        q *= -sigma * c_eps                                # NO 'sides' factor (paper code)
        q += c_V * c_I / c_A
        tmp = nrm(dTdt, dTdt * e["errordTdt"]); tmp *= (c_rho * c_s * c_cp + c_rho_p * c_s_p * c_cp_p); q -= tmp
        tmp = nrm(d2x, d2x * e["errord2Tdx2"]); tmp *= (c_kx + c_kp); q -= tmp
        tmp = nrm(d2y, d2y * e["errord2Tdy2"]); tmp *= (c_ky + c_kp); q += tmp   # sign error (paper code)
        del tmp
        c_Thot -= (c_Tcold * ratio)[:, :, None]
        np.divide(q, c_Thot, out=q)                        # h
        del c_Thot
        q *= c_L
        q /= c_k                                           # Nu
        Nu = q
        Num = Nu.mean(axis=2)
        out["errorNu"][i] = Num.mean()
        out["errorNuf"][i] = np.abs(Nu - Num[:, :, None]).mean()
        out["errorNu_p"][i] = abs(out["errorNu"][i] - Nu_ref.mean()) / Nu_ref.mean() * 100
        out["errorNuf_p"][i] = abs(out["errorNuf"][i] - Nuf_ref) / Nuf_ref * 100
        print(f"  paper-replica sample {i + 1}/{n} ({time.time() - t0:.1f} s)", flush=True)
    return out


def stats(v):
    v = np.asarray(v, dtype=float).ravel()
    return dict(mean=float(v.mean()), std=float(v.std(ddof=1)), min=float(v.min()), max=float(v.max()), n=int(v.size))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, type=Path)
    ap.add_argument("--ref", required=True, type=Path)
    ap.add_argument("--samples", type=int, default=40)
    ap.add_argument("--seed", type=int, default=2025)
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent / "output" / "montecarlo")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    import h5py
    with h5py.File(args.ref / "error_metrics_1000_f.mat", "r") as f:
        ref = {k: np.asarray(f[k][()]).ravel() for k in ("errorNu", "errorNu_p", "errorNuf", "errorNuf_p")}

    case = io.load_sj_case(args.data)
    hfs, cond = case["HFS"], case["Conditions"]
    Tcold = case["Tcold"].mean(axis=2)
    dx, dy, dt = case["dx"], case["dy"], 1.0 / case["f_acq"]
    t0 = time.time()
    T1, nmod = pod_filter(case["Thot"], criterion="HardThreshold", beta=2000 / (399 * 604), return_nmod=True, verbose=False)
    T2 = gaussian_filter3(T1, (3, 3, 0.1), (9, 3, 1))
    del T1
    Thot, dTdt, d2x, d2y = sgolay32_filter(T2, (5, 5, 3), (dx, dy, dt))
    del T2
    print(f"filtering done (Nmod={nmod}) in {time.time() - t0:.0f} s", flush=True)

    results = {"reference_paper_1000_samples": {k: stats(v) for k, v in ref.items()},
               "samples": args.samples, "seed": args.seed}

    for film in ("adiabatic", "ambient"):
        det = calculate_heat_transfer(Thot, Tcold, hfs, cond, compute=("Nu", "h"), time_der=True, spatial_der=True,
                                      dTdt=dTdt, d2Tdx2=d2x, d2Tdy2=d2y, film_temperature=film, verbose=False)
        results[f"deterministic_{film}"] = {"Nu_mean": float(det["Nu"].mean()),
                                            "Nuf_mean": float(np.abs(det["Nu"] - det["Nu"].mean(axis=2, keepdims=True)).mean())}
        if film == "adiabatic":
            Nu_ref = det["Nu"]
        t0 = time.time()
        err = dict(ERROR, samples=args.samples)
        mc = montecarlo_uncertainty(Thot, Tcold, hfs, cond, err, compute=("Nu",), time_der=True, spatial_der=True,
                                    dTdt=dTdt, d2Tdx2=d2x, d2Tdy2=d2y, reference=det, mode="snapshots",
                                    seed=args.seed, film_temperature=film, verbose=True)
        results[f"python_corrected_{film}"] = {k: stats(v) for k, v in mc.items()}
        print(f"corrected MC ({film}) done in {time.time() - t0:.0f} s: "
              f"errorNu {mc['errorNu'].mean():.3f} +- {mc['errorNu'].std(ddof=1):.3f}, "
              f"errorNu_p {mc['errorNu_p'].mean():.2f} %, errorNuf_p {mc['errorNuf_p'].mean():.2f} %", flush=True)
        (args.out / "montecarlo_check.json").write_text(json.dumps(results, indent=1))

    t0 = time.time()
    rep = paper_montecarlo(Thot, Tcold, dTdt, d2x, d2y, hfs, cond, args.samples, args.seed, Nu_ref)
    results["paper_replica"] = {k: stats(v) for k, v in rep.items()}
    print(f"paper-replica MC done in {time.time() - t0:.0f} s: errorNu {rep['errorNu'].mean():.3f} +- "
          f"{rep['errorNu'].std(ddof=1):.3f}, errorNu_p {rep['errorNu_p'].mean():.2f} %, "
          f"errorNuf_p {rep['errorNuf_p'].mean():.2f} %", flush=True)
    (args.out / "montecarlo_check.json").write_text(json.dumps(results, indent=1))
    print(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
