#!/usr/bin/env python
"""Build ``REPORT.md`` from the JSON outputs of the validation scripts."""
from __future__ import annotations

import json
import sys
from pathlib import Path

here = Path(__file__).resolve().parent
cmp_path = Path(sys.argv[1]) if len(sys.argv) > 1 else here / "output" / "full" / "comparison.json"
mc_path = Path(sys.argv[2]) if len(sys.argv) > 2 else here / "output" / "montecarlo" / "montecarlo_check.json"

rows = json.loads(cmp_path.read_text())
ok = [r for r in rows if "error" not in r]
maps = [r for r in ok if r["quantity"] in ("Num", "Nuf", "dT", "qunsrel", "qkrel") and r["case"] != "mPOD_sgolay_exactF"]
other = [r for r in ok if r["quantity"] in ("Nu_snap", "time_series") ]
def fmt(x):
    return f"{x:.1e}"
summary = []
if maps:
    summary += [f"* {len(maps)} time-averaged maps compared (Num, Nuf, dT, q ratios) over {len({r['case'] for r in maps})} "
                f"pipelines: RMS relative difference between {fmt(min(r['rms_diff_rel'] for r in maps))} and "
                f"{fmt(max(r['rms_diff_rel'] for r in maps))}; maximum pixel-wise relative difference "
                f"{fmt(max(r['max_rel_diff'] for r in maps))}."]
if other:
    summary += [f"* {len(other)} snapshot-level comparisons (fluctuation snapshots and time series): RMS relative "
                f"difference between {fmt(min(r['rms_diff_rel'] for r in other))} and "
                f"{fmt(max(r['rms_diff_rel'] for r in other))} (the maximum pixel-wise relative differences of these "
                f"zero-mean fields are dominated by values close to zero)."]
nm = sorted({(r['case'], r['nmod']) for r in ok if r.get('nmod') is not None})
if nm:
    groups = {}
    for c, n in nm:
        groups.setdefault(n, []).append(c)
    summary += ["* Number of modes selected by the Python code: " + "; ".join(
        f"{n} modes in {len(cs)} case(s) ({', '.join(cs) if len(cs) <= 4 else ', '.join(cs[:3]) + ', ...'})"
        for n, cs in sorted(groups.items())) + "."]
ex = [r for r in ok if r["case"] == "mPOD_sgolay_exactF"]
if ex:
    summary += ["* Effect of the exact inverse transform in the mPOD filter (case `mPOD_sgolay_exactF`, compared "
                "with the legacy MATLAB result of `mPOD_sgolay`): " + "; ".join(
                    f"{r['quantity']}: field mean {r['ref_mean']:.4g} (legacy) vs {r['py_mean']:.4g} (exact), "
                    f"RMS rel. diff. {fmt(r['rms_diff_rel'])}" for r in ex) + "."]
errs = [r for r in rows if "error" in r]
if errs:
    summary += [f"* {len(errs)} case(s) could not be compared: " + ", ".join(str(r.get('case')) for r in errs) + "."]

lines = ["# Validation report: Python `pirt` vs MATLAB PIRT (paper results)", "",
         "Reference: results of the MATLAB toolbox (commit `fa07252`) used in Robledo et al., ",
         "Exp. Therm. Fluid Sci. 169 (2025) 111526, on the full sweeping-jet dataset ",
         "(399 x 604 pixels x 2000 snapshots). Python: `pirt` 1.0.0 (this repository), ",
         "same raw data, same configuration as the MATLAB scripts (see `compare_with_paper_results.py`).", "",
         "Metrics: `max |py-ref|`, RMS of the difference relative to the RMS of the reference, ",
         "median and maximum of `|py-ref|/|ref|`, and the Pearson correlation. `film` is the film-temperature ",
         "definition used for `k_air` (see README). `Nmod` is the number of POD/mPOD modes selected by the Python code ",
         "(MATLAB logs: 70 for POD with the hard threshold, 48 for mPOD with the Elbow criterion).", "",
         "## Summary", "", *summary, "", "## Detailed comparison", "",
         "| case | quantity | film | Nmod | shape | ref mean | py mean | max abs diff | RMS diff / RMS ref | median rel diff | max rel diff | corr |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|"]
for r in rows:
    if "error" in r:
        lines.append(f"| {r.get('case')} | {r.get('quantity', '')} | | | | | | {r['error']} | | | | |")
        continue
    lines.append(f"| {r['case']} | {r['quantity']} | {r.get('film', '')} | {r.get('nmod', '')} | "
                 f"{'x'.join(map(str, r['shape']))} | {r['ref_mean']:.6g} | {r['py_mean']:.6g} | "
                 f"{r['max_abs_diff']:.2e} | {r['rms_diff_rel']:.2e} | {r['median_rel_diff']:.2e} | "
                 f"{r['max_rel_diff']:.2e} | {r['corr']:.10f} |")

if mc_path.exists():
    mc = json.loads(mc_path.read_text())
    lines += ["", "## Monte Carlo uncertainty (`montecarlo_check.py`)", "",
              f"Samples: {mc.get('samples')} per variant (reference file: 1000 samples). Pipeline POD + Gaussian + ",
              "Savitzky-Golay, Nu with unsteady and tangential terms, uncertainties of `SJ_main_12_error.m`.", "",
              "| variant | quantity | mean | std | min | max |", "|---|---|---|---|---|---|"]
    for variant, block in mc.items():
        if not isinstance(block, dict) or "samples" in variant:
            continue
        if all(isinstance(v, dict) for v in block.values()):
            for q, st in block.items():
                lines.append(f"| {variant} | {q} | {st['mean']:.4g} | {st['std']:.3g} | {st['min']:.4g} | {st['max']:.4g} |")
        else:
            for q, v in block.items():
                lines.append(f"| {variant} | {q} | {v:.4g} | | | |")

(here / "REPORT.md").write_text("\n".join(lines) + "\n")
print("REPORT.md written with", len(rows), "rows")
