# Cross-validation of the Python and MATLAB implementations

The MATLAB toolbox (commit `fa07252`) was used to produce the results of

> I. Robledo, J. Alfaro, C. Sanmiguel Vila, R. Castellanos, *Unsteady convective
> heat transfer of an impinging sweeping jet: A discussion on the effect of
> spatiotemporal filtering*, Exp. Therm. Fluid Sci. 169 (2025) 111526.

The processing scripts of that work (`SJ_main_3_dT.m`, `SJ_main_5_filter_effect.m`,
`SJ_main_6_temp_an.m`, `SJ_main_8_img_filt.m`, `SJ_main_fderivadas.m`,
`SJ_main_10_qjqk.m`) saved, for each filtering pipeline, the time-averaged
Nusselt number `Num`, the mean absolute fluctuation `Nuf`, the mean absolute
time derivative `dT`, the first 100 fluctuation snapshots `Nu_snap` and/or the
fluctuation time series at the turn-around point of the jet. These files are
the reference used here.

`compare_with_paper_results.py` reproduces every configuration with the
Python package on the same raw data (399 × 604 pixels × 2000 snapshots,
`Thot.mat`/`Tcold.mat` in MATLAB v7.3 format) and compares the outputs
(maximum and RMS differences, relative to the reference, and correlation).

```bash
python validation/compare_with_paper_results.py --data <folder with Thot.mat, Tcold.mat, Resolution.mat, TestConditions.mat, CONF_DATA.STR> \
                                                 --ref  <folder with the result_*.mat / IMG_Proc_*.mat files> \
                                                 --out  validation/output
```
The raw data and reference files are not part of the repository (they are
available from the authors upon request). The report obtained on the UC3M
servers is kept in [`REPORT.md`](REPORT.md).

Two aspects must be taken into account when reading the report:

* **Film temperature.** Commit `5a4dd70` (2025‑02‑21) changed the film
  temperature used for `k_air` from `(T_w + T_cold)/2` to `(T_w + T_amb,hot)/2`.
  Reference files written before that date are compared with
  `film_temperature='adiabatic'`, the later ones (time series of Fig. 8) with
  `'ambient'`.
* **mPOD back-transformation.** The MATLAB `multiscale_POD_filter.m` used for
  the paper applies the forward DFT matrix on both sides also when
  transforming the filtered correlation matrix back to the time domain
  (time-reversed temporal modes). The cases `mPOD_*` reproduce this with
  `legacy_transform=True`; the case `mPOD_sgolay_exactF` uses the exact
  inverse (the default since v1.1 in both implementations) and quantifies
  the effect on the final maps.
