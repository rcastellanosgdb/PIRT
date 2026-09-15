#!/usr/bin/env python
"""Cross-validation of the Python PIRT implementation against the results
obtained with the MATLAB toolbox for the paper

    I. Robledo, J. Alfaro, C. Sanmiguel Vila, R. Castellanos, "Unsteady
    convective heat transfer of an impinging sweeping jet: A discussion on the
    effect of spatiotemporal filtering", Exp. Therm. Fluid Sci. 169 (2025)
    111526. https://doi.org/10.1016/j.expthermflusci.2025.111526

The reference files (``result_*.mat``, ``IMG_Proc_*.mat``) were produced by the
scripts in ``PIRT_paper_SJ_2025/SJ_processing`` (``SJ_main_3_dT.m``,
``SJ_main_5_filter_effect.m``, ``SJ_main_6_temp_an.m``, ``SJ_main_8_img_filt.m``,
``SJ_main_fderivadas.m``) with the MATLAB code of commit fa07252. Each case
below reproduces the corresponding PIRT configuration with the Python package
and compares the outputs.

Usage
-----
    python compare_with_paper_results.py --data <folder with Thot.mat, Tcold.mat, ...>
                                         --ref <folder with result_*.mat>
                                         [--out output] [--cases FD POD_sgolay ...]
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
import pirt  # noqa: E402
from pirt import io  # noqa: E402
from pirt.filters import gaussian_filter3, mpod_filter, pod_filter, sgolay32_filter, wiener3  # noqa: E402
from pirt.heat_transfer import calculate_heat_transfer  # noqa: E402

# --------------------------------------------------------------------------- #
# Filter specifications used in the paper
# --------------------------------------------------------------------------- #
BETA = 2000 / (399 * 604)
F_POD = ("POD", dict(criterion="HardThreshold", beta=BETA))
F_WIENER = ("wiener3", dict(kernel=(7, 7, 1)))
F_GAUSS = ("gaussian", dict(filter_size=(9, 3, 1), sigma=(3, 3, 0.1)))
F_SGOLAY = ("sgolay32", dict(kernel_size=(5, 5, 3)))          # h filled at run time
F_SGOLAY11 = ("sgolay32", dict(kernel_size=(11, 11, 3)))
F_MPOD = ("mPOD", dict(mode="Peak Removal", threshold=0.99, fpeaks=(1.2, 51, 95), w=1,
                       legacy_transform=True))   # MATLAB back-transform used for the paper
F_MPOD_EXACT = ("mPOD", dict(mode="Peak Removal", threshold=0.99, fpeaks=(1.2, 51, 95), w=1,
                             legacy_transform=False))

# Film temperature used for k_air in Nu. The MATLAB code evaluated k_air at
# (T_w + T_cold)/2 until commit 5a4dd70 (2025-02-21), when it was changed to
# (T_w + T_amb,hot)/2. Reference files dated before that use 'adiabatic'.
FILM_BY_PREFIX = {"TS_": "ambient"}
DEFAULT_FILM = "adiabatic"

# case name -> (filter chain, time_der, spatial_der, reference file, reference kind)
#   kind 'maps'  : Num/Nuf (+ dT) time-averaged maps
#   kind 'snaps' : Nu_snap = first 100 snapshots of Nu - mean(Nu)
#   kind 'series': time_series = mean of Nu' over pixels (199:201, 345:347)
CASES = {
    # SJ_main_5_filter_effect.m (finite differences vs sgolay)
    "FD":                 ((), True, False, "result_FD.mat", "maps"),
    "POD_FD":             ((F_POD,), True, False, "result_POD_FD.mat", "maps"),
    "POD_wiener_FD":      ((F_POD, F_WIENER), True, False, "result_POD_wiener_FD.mat", "maps"),
    "POD_gauss_FD":       ((F_POD, F_GAUSS), True, False, "result_POD_gauss_FD.mat", "maps"),
    # SJ_main_3_dT.m (filter combinations, unsteady term only)
    "POD_sgolay":         ((F_POD, F_SGOLAY), True, False, "result_POD_sgolay.mat", "maps"),
    "wiener_sgolay":      ((F_WIENER, F_SGOLAY), True, False, "result_wiener_sgolay.mat", "maps"),
    "gauss_sgolay":       ((F_GAUSS, F_SGOLAY), True, False, "result_gauss_sgolay.mat", "maps"),
    "POD_wiener_sgolay":  ((F_POD, F_WIENER, F_SGOLAY), True, False, "result_POD_wiener_sgolay.mat", "maps"),
    "POD_gauss_sgolay":   ((F_POD, F_GAUSS, F_SGOLAY), True, False, "result_POD_gauss_sgolay.mat", "maps"),
    "mPOD_sgolay":        ((F_MPOD, F_SGOLAY), True, False, "result_mPOD_sgolay.mat", "maps"),
    "mPOD_wiener_sgolay": ((F_MPOD, F_WIENER, F_SGOLAY), True, False, "result_mPOD_wiener_sgolay.mat", "maps"),
    # same as mPOD_sgolay but with the exact inverse Fourier transform of the correlation
    # matrix (quantifies the effect of the legacy transform of multiscale_POD_filter.m)
    "mPOD_sgolay_exactF": ((F_MPOD_EXACT, F_SGOLAY), True, False, "result_mPOD_sgolay.mat", "maps"),
    # SJ_main_fderivadas.m (effect of the derivative terms, Fig. 5)
    "sgolay_ND":          ((F_POD, F_GAUSS, F_SGOLAY), False, False, "result_sgolay_ND.mat", "maps"),
    "sgolay_TD":          ((F_POD, F_GAUSS, F_SGOLAY), True, False, "result_sgolay_TD.mat", "maps"),
    "sgolay_SD":          ((F_POD, F_GAUSS, F_SGOLAY), False, True, "result_sgolay_SD.mat", "maps"),
    "sgolay_D":           ((F_POD, F_GAUSS, F_SGOLAY), True, True, "result_sgolay_D.mat", "maps"),
    # SJ_main_8_img_filt.m (instantaneous fluctuation snapshots, Fig. 9)
    "IMG_no_filt":        ((), True, False, "IMG_Proc_no_filt.mat", "snaps"),
    "IMG_sgolay":         ((F_SGOLAY,), True, False, "IMG_Proc_sgolay.mat", "snaps"),
    "IMG_POD_sgolay":     ((F_POD, F_SGOLAY), True, False, "IMG_Proc_POD_Sgolay.mat", "snaps"),
    "IMG_POD_gauss_sgolay": ((F_POD, F_GAUSS, F_SGOLAY), True, False, "IMG_Proc_POD_gauss_Sgolay.mat", "snaps"),
    # SJ_main_6_temp_an.m (time series at the turn-around point, Fig. 8)
    "TS_normal":          ((), True, False, "result_normal.mat", "series"),
    "TS_W":               ((F_WIENER,), True, False, "result_W.mat", "series"),
    "TS_gauss":           ((F_GAUSS,), True, False, "result_gauss.mat", "series"),
    "TS_sgolay":          ((F_SGOLAY,), True, False, "result_sgolay.mat", "series"),
    "TS_POD_W":           ((F_POD, F_WIENER), True, False, "result_POD_W.mat", "series"),
    "TS_POD_gauss":       ((F_POD, F_GAUSS), True, False, "result_POD_gauss.mat", "series"),
    "TS_POD_W_sgolay":    ((F_POD, F_WIENER, F_SGOLAY), True, False, "result_POD_W_Sgolay.mat", "series"),
    "TS_POD_gauss_sgolay": ((F_POD, F_GAUSS, F_SGOLAY), True, False, "result_POD_gauss_Sgolay.mat", "series"),
    # SJ_main_10_qjqk.m (relative magnitude of the unsteady and tangential terms)
    "qjqk_POD_gauss_sgolay11": ((F_POD, F_GAUSS, F_SGOLAY11), True, True, "result_qjqk_POD_gauss_sgolay.mat", "qjqk"),
}


def load_reference(path: Path) -> dict:
    """Load a MATLAB v7.3 reference file, transposing to the MATLAB layout."""
    import h5py
    out = {}
    with h5py.File(path, "r") as f:
        for k in f.keys():
            if k == "#refs#":
                continue
            arr = f[k][()]
            out[k] = np.ascontiguousarray(arr.T)
    return out


def compare(name: str, py: np.ndarray, ref: np.ndarray) -> dict:
    py = np.asarray(py, dtype=np.float64).squeeze()
    ref = np.asarray(ref, dtype=np.float64).squeeze()
    if py.shape != ref.shape:
        return {"quantity": name, "shape_py": py.shape, "shape_ref": ref.shape, "error": "shape mismatch"}
    diff = py - ref
    rms_ref = float(np.sqrt(np.mean(ref**2)))
    denom = np.maximum(np.abs(ref), 1e-12)
    corr = float(np.corrcoef(py.ravel(), ref.ravel())[0, 1]) if py.size > 1 else 1.0
    return {
        "quantity": name, "shape": list(py.shape),
        "ref_mean": float(ref.mean()), "py_mean": float(py.mean()),
        "max_abs_diff": float(np.abs(diff).max()),
        "rms_diff": float(np.sqrt(np.mean(diff**2))),
        "rms_diff_rel": float(np.sqrt(np.mean(diff**2)) / rms_ref),
        "max_rel_diff": float(np.max(np.abs(diff) / denom)),
        "median_rel_diff": float(np.median(np.abs(diff) / denom)),
        "corr": corr,
    }


class Runner:
    def __init__(self, data_dir: Path, out_dir: Path, verbose=True):
        self.verbose = verbose
        t0 = time.time()
        self.case = io.load_sj_case(data_dir)
        self.Thot = np.ascontiguousarray(self.case["Thot"], dtype=np.float64)
        Tcold = self.case["Tcold"]
        self.Tcold = Tcold.mean(axis=2) if Tcold.ndim == 3 else Tcold  # as PIRT.go()
        self.h = (self.case["dx"], self.case["dy"], 1.0 / self.case["f_acq"])
        self.out_dir = out_dir
        self._cache: dict = {}
        self._last_sg = (None, None)
        self.nmods: dict = {}
        print(f"Data loaded in {time.time() - t0:.1f} s: Thot {self.Thot.shape} {self.Thot.dtype}, "
              f"dx={self.h[0]:.6g} m, dt={self.h[2]:.6g} s")

    def filtered(self, chain: tuple):
        """Apply a filter chain (memoised by prefix). Returns (Thot, derivs).

        Chains ending in the Savitzky-Golay filter (which produces four arrays)
        are only kept in a single-entry cache to bound the memory footprint.
        """
        key = repr(chain)
        if key in self._cache:
            return self._cache[key]
        if chain and chain[-1][0] == "sgolay32" and self._last_sg[0] == key:
            return self._last_sg[1]
        if not chain:
            res = (self.Thot, {})
        else:
            Tprev, derivs = self.filtered(chain[:-1])
            ftype, params = chain[-1]
            derivs = dict(derivs)
            t0 = time.time()
            if ftype == "POD":
                Tnew, nmod = pod_filter(Tprev, return_nmod=True, verbose=False, **params)
                self.nmods[key] = nmod
                print(f"  POD: Nmod = {nmod} ({time.time() - t0:.1f} s)")
            elif ftype == "mPOD":
                Tnew, nmod = mpod_filter(Tprev, self.case["f_acq"], return_nmod=True, verbose=False, **params)
                self.nmods[key] = nmod
                print(f"  mPOD: Nmod = {nmod} ({time.time() - t0:.1f} s)")
            elif ftype == "wiener3":
                Tnew = wiener3(Tprev, **params)
                print(f"  wiener3 ({time.time() - t0:.1f} s)")
            elif ftype == "gaussian":
                Tnew = gaussian_filter3(Tprev, **params)
                print(f"  gaussian ({time.time() - t0:.1f} s)")
            elif ftype == "sgolay32":
                Tnew, dTdt, d2x, d2y = sgolay32_filter(Tprev, h=self.h, **params)
                derivs = {"dTdt": dTdt, "d2Tdx2": d2x, "d2Tdy2": d2y}
                print(f"  sgolay32 ({time.time() - t0:.1f} s)")
            else:
                raise ValueError(ftype)
            res = (Tnew, derivs)
        if chain and chain[-1][0] == "sgolay32":
            self._last_sg = (key, res)
        else:
            self._cache[key] = res
        return res

    def run_case(self, name: str, ref_dir: Path, film: str | None = None) -> list[dict]:
        chain, time_der, spatial_der, ref_file, kind = CASES[name]
        if film is None:
            film = next((f for pre, f in FILM_BY_PREFIX.items() if name.startswith(pre)), DEFAULT_FILM)
        print(f"\n=== {name}: filters={[c[0] for c in chain]} TimeDer={time_der} SpatialDer={spatial_der} "
              f"film={film}")
        ref_path = ref_dir / ref_file
        if not ref_path.exists():
            print(f"  reference {ref_path} not found, skipping")
            return [{"case": name, "error": "reference not found"}]
        ref = load_reference(ref_path)
        Th, derivs = self.filtered(chain)
        t0 = time.time()
        res = calculate_heat_transfer(
            Th, self.Tcold, self.case["HFS"], self.case["Conditions"], compute=("Nu",),
            time_der=time_der, spatial_der=spatial_der,
            dTdt=derivs.get("dTdt") if time_der else None,
            d2Tdx2=derivs.get("d2Tdx2") if spatial_der else None,
            d2Tdy2=derivs.get("d2Tdy2") if spatial_der else None,
            film_temperature=film, return_terms=(kind == "qjqk"), verbose=False)
        Nu = res["Nu"]
        print(f"  heat transfer ({time.time() - t0:.1f} s): Nu shape {Nu.shape}")
        Num = Nu.mean(axis=2)
        rows = []
        if kind == "maps":
            Nuf = np.abs(Nu - Num[:, :, None]).mean(axis=2)
            rows.append(compare("Num", Num, ref["Num"]))
            rows.append(compare("Nuf", Nuf, ref["Nuf"]))
            if "dT" in ref:
                dTdt = derivs.get("dTdt") if derivs else res.get("dTdt_hot")
                rows.append(compare("dT", np.abs(dTdt).mean(axis=2), ref["dT"]))
            np.savez_compressed(self.out_dir / f"py_{name}.npz", Num=Num, Nuf=Nuf)
        elif kind == "snaps":
            snaps = (Nu - Num[:, :, None])[:, :, :100]
            rows.append(compare("Nu_snap", snaps, ref["Nu_snap"]))
        elif kind == "series":
            Nuf3 = Nu - Num[:, :, None]
            ts = Nuf3[198:201, 344:347, :].mean(axis=(0, 1))
            rows.append(compare("time_series", ts, ref["time_series"].ravel()))
        elif kind == "qjqk":
            Nuf = np.abs(Nu - Num[:, :, None]).mean(axis=2)
            rows.append(compare("Num", Num, ref["Num"]))
            rows.append(compare("Nuf", Nuf, ref["Nuf"]))
            terms = res["terms"]
            qj = terms["q_joule"]
            rows.append(compare("qunsrel", np.abs(terms["q_unsteady"]).mean(axis=2) / qj, ref["qunsrel"]))
            rows.append(compare("qkrel", np.abs(terms["q_tangential"]).mean(axis=2) / qj, ref["qkrel"]))
        for r in rows:
            r["case"] = name
            r["film"] = film
            for k in range(1, len(chain) + 1):
                if repr(chain[:k]) in self.nmods:
                    r["nmod"] = self.nmods[repr(chain[:k])]
            print("  " + json.dumps(r))
        return rows


def write_report(rows: list[dict], path: Path):
    lines = ["# Python vs MATLAB (paper) results", "",
             "Reference: MATLAB PIRT (commit fa07252) results used in Robledo et al., ETFS 169 (2025) 111526, "
             "full dataset 399 x 604 x 2000 snapshots.", "",
             "| case | quantity | film | Nmod | shape | ref mean | py mean | max abs diff | RMS diff | RMS diff / RMS ref | median rel diff | max rel diff | corr |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        if "error" in r:
            lines.append(f"| {r.get('case')} | {r.get('quantity', '')} | | | | | | {r['error']} | | | | | |")
            continue
        lines.append(
            f"| {r['case']} | {r['quantity']} | {r.get('film', '')} | {r.get('nmod', '')} | "
            f"{'x'.join(map(str, r['shape']))} | {r['ref_mean']:.6g} | "
            f"{r['py_mean']:.6g} | {r['max_abs_diff']:.3e} | {r['rms_diff']:.3e} | {r['rms_diff_rel']:.3e} | "
            f"{r['median_rel_diff']:.3e} | {r['max_rel_diff']:.3e} | {r['corr']:.8f} |")
    path.write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, type=Path, help="folder with Thot.mat, Tcold.mat, Resolution.mat, TestConditions.mat, CONF_DATA.STR")
    ap.add_argument("--ref", required=True, type=Path, help="folder with the paper result_*.mat files")
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent / "output")
    ap.add_argument("--cases", nargs="*", default=list(CASES), help="subset of cases to run")
    ap.add_argument("--film", choices=["ambient", "adiabatic"], default=None,
                    help="override the film temperature definition for all cases")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    print(f"pirt {pirt.__version__} | numpy {np.__version__}")
    runner = Runner(args.data, args.out)
    rows = []
    for name in args.cases:
        try:
            rows.extend(runner.run_case(name, args.ref, film=args.film))
        except Exception as exc:  # keep going, report the failure
            import traceback
            traceback.print_exc()
            rows.append({"case": name, "error": repr(exc)})
        (args.out / "comparison.json").write_text(json.dumps(rows, indent=1))
        write_report(rows, args.out / "comparison.md")
    print("\nDone. Report written to", args.out / "comparison.md")


if __name__ == "__main__":
    main()
