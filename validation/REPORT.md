# Validation report: Python `pirt` vs MATLAB PIRT (paper results)

Reference: results of the MATLAB toolbox (commit `fa07252`) used in Robledo et al., 
Exp. Therm. Fluid Sci. 169 (2025) 111526, on the full sweeping-jet dataset 
(399 x 604 pixels x 2000 snapshots). Python: `pirt` 1.0.0 (this repository), 
same raw data, same configuration as the MATLAB scripts (see `compare_with_paper_results.py`).

Metrics: `max |py-ref|`, RMS of the difference relative to the RMS of the reference, 
median and maximum of `|py-ref|/|ref|`, and the Pearson correlation. `film` is the film-temperature 
definition used for `k_air` (see README). `Nmod` is the number of POD/mPOD modes selected by the Python code 
(MATLAB logs: 70 for POD with the hard threshold, 48 for mPOD with the Elbow criterion).

## Summary

* 41 time-averaged maps compared (Num, Nuf, dT, q ratios) over 16 pipelines: RMS relative difference between 1.5e-13 and 4.6e-11; maximum pixel-wise relative difference 1.7e-10.
* 12 snapshot-level comparisons (fluctuation snapshots and time series): RMS relative difference between 2.6e-16 and 2.5e-11 (the maximum pixel-wise relative differences of these zero-mean fields are dominated by values close to zero).
* Number of modes selected by the Python code: 48 modes in 3 case(s) (mPOD_sgolay, mPOD_sgolay_exactF, mPOD_wiener_sgolay); 70 modes in 17 case(s) (IMG_POD_gauss_sgolay, IMG_POD_sgolay, POD_FD, ...).
* Effect of the exact inverse transform in the mPOD filter (case `mPOD_sgolay_exactF`, compared with the legacy MATLAB result of `mPOD_sgolay`): Num: field mean 37.86 (legacy) vs 37.85 (exact), RMS rel. diff. 1.2e-03; Nuf: field mean 10.42 (legacy) vs 12.8 (exact), RMS rel. diff. 2.4e-01; dT: field mean 0.7325 (legacy) vs 0.8977 (exact), RMS rel. diff. 2.4e-01.

## Detailed comparison

| case | quantity | film | Nmod | shape | ref mean | py mean | max abs diff | RMS diff / RMS ref | median rel diff | max rel diff | corr |
|---|---|---|---|---|---|---|---|---|---|---|---|
| FD | Num | adiabatic |  | 399x604 | 37.8574 | 37.8574 | 2.40e-10 | 1.61e-13 | 8.29e-16 | 2.91e-12 | 1.0000000000 |
| FD | Nuf | adiabatic |  | 399x604 | 65.8635 | 65.8635 | 4.04e-10 | 1.59e-13 | 8.56e-16 | 2.92e-12 | 1.0000000000 |
| POD_FD | Num | adiabatic | 70 | 399x604 | 37.8554 | 37.8554 | 4.94e-10 | 1.75e-13 | 1.06e-15 | 6.03e-12 | 1.0000000000 |
| POD_FD | Nuf | adiabatic | 70 | 399x604 | 14.7143 | 14.7143 | 2.05e-10 | 1.88e-13 | 4.15e-14 | 5.85e-12 | 1.0000000000 |
| POD_wiener_FD | Num | adiabatic | 70 | 399x604 | 38.1124 | 38.1124 | 2.39e-10 | 1.71e-13 | 3.66e-14 | 2.91e-12 | 1.0000000000 |
| POD_wiener_FD | Nuf | adiabatic | 70 | 399x604 | 12.3597 | 12.3597 | 1.69e-10 | 1.32e-12 | 8.51e-13 | 1.56e-11 | 1.0000000000 |
| POD_gauss_FD | Num | adiabatic | 70 | 399x604 | 37.8522 | 37.8522 | 2.48e-10 | 1.62e-13 | 1.09e-14 | 3.03e-12 | 1.0000000000 |
| POD_gauss_FD | Nuf | adiabatic | 70 | 399x604 | 12.2807 | 12.2807 | 8.54e-11 | 2.32e-13 | 1.23e-13 | 3.19e-12 | 1.0000000000 |
| IMG_no_filt | Nu_snap | adiabatic |  | 399x604x100 | -0.335168 | -0.335168 | 3.30e-09 | 1.58e-13 | 6.26e-16 | 5.07e-08 | 1.0000000000 |
| TS_normal | time_series | ambient |  | 2000 | -1.42109e-15 | 2.75691e-15 | 1.14e-13 | 2.60e-16 | 2.28e-16 | 5.76e-13 | 1.0000000000 |
| TS_W | time_series | ambient |  | 2000 | 1.7252e-14 | 3.97904e-16 | 4.47e-09 | 1.90e-11 | 1.74e-11 | 1.55e-08 | 1.0000000000 |
| TS_gauss | time_series | ambient |  | 2000 | 3.29692e-14 | 5.68434e-16 | 2.71e-10 | 1.36e-12 | 1.22e-12 | 4.77e-09 | 1.0000000000 |
| TS_POD_W | time_series | ambient | 70 | 2000 | -1.12408e-14 | 1.0516e-15 | 3.85e-09 | 2.54e-11 | 2.25e-11 | 7.55e-08 | 1.0000000000 |
| TS_POD_gauss | time_series | ambient | 70 | 2000 | -1.70246e-14 | -3.78009e-15 | 4.51e-10 | 1.75e-12 | 1.48e-12 | 2.01e-09 | 1.0000000000 |
| TS_sgolay | time_series | ambient |  | 1998 | 4.86498e-15 | 5.69003e-16 | 3.23e-10 | 1.92e-12 | 1.86e-12 | 5.23e-09 | 1.0000000000 |
| IMG_sgolay | Nu_snap | adiabatic |  | 399x604x100 | -0.229482 | -0.229482 | 8.72e-10 | 4.02e-12 | 4.11e-12 | 1.05e-03 | 1.0000000000 |
| POD_sgolay | Num | adiabatic | 70 | 399x604 | 37.8474 | 37.8474 | 2.63e-10 | 1.95e-13 | 7.00e-14 | 3.20e-12 | 1.0000000000 |
| POD_sgolay | Nuf | adiabatic | 70 | 399x604 | 12.382 | 12.382 | 9.36e-11 | 2.69e-13 | 1.55e-13 | 3.13e-12 | 1.0000000000 |
| POD_sgolay | dT | adiabatic | 70 | 399x604 | 0.866461 | 0.866461 | 8.24e-13 | 1.89e-13 | 1.33e-13 | 1.27e-12 | 1.0000000000 |
| IMG_POD_sgolay | Nu_snap | adiabatic | 70 | 399x604x100 | -0.232362 | -0.232362 | 9.41e-10 | 6.58e-12 | 6.83e-12 | 3.22e-04 | 1.0000000000 |
| POD_gauss_sgolay | Num | adiabatic | 70 | 399x604 | 37.8439 | 37.8439 | 2.38e-10 | 1.99e-13 | 7.83e-14 | 2.96e-12 | 1.0000000000 |
| POD_gauss_sgolay | Nuf | adiabatic | 70 | 399x604 | 12.1843 | 12.1843 | 9.84e-11 | 2.78e-13 | 1.65e-13 | 3.06e-12 | 1.0000000000 |
| POD_gauss_sgolay | dT | adiabatic | 70 | 399x604 | 0.85151 | 0.85151 | 8.31e-13 | 1.97e-13 | 1.39e-13 | 1.47e-12 | 1.0000000000 |
| sgolay_ND | Num | adiabatic | 70 | 399x604 | 37.8548 | 37.8548 | 2.37e-10 | 2.07e-13 | 9.65e-14 | 3.02e-12 | 1.0000000000 |
| sgolay_ND | Nuf | adiabatic | 70 | 399x604 | 0.322302 | 0.322302 | 5.93e-12 | 4.54e-13 | 1.88e-13 | 5.65e-12 | 1.0000000000 |
| sgolay_TD | Num | adiabatic | 70 | 399x604 | 37.8439 | 37.8439 | 2.38e-10 | 1.99e-13 | 7.83e-14 | 2.96e-12 | 1.0000000000 |
| sgolay_TD | Nuf | adiabatic | 70 | 399x604 | 12.1843 | 12.1843 | 9.84e-11 | 2.78e-13 | 1.65e-13 | 3.06e-12 | 1.0000000000 |
| sgolay_SD | Num | adiabatic | 70 | 399x604 | 37.8111 | 37.8111 | 6.56e-10 | 4.93e-12 | 5.05e-12 | 8.03e-12 | 1.0000000000 |
| sgolay_SD | Nuf | adiabatic | 70 | 399x604 | 0.349484 | 0.349484 | 1.15e-11 | 3.80e-12 | 3.05e-12 | 1.32e-11 | 1.0000000000 |
| sgolay_D | Num | adiabatic | 70 | 399x604 | 37.8003 | 37.8003 | 6.62e-10 | 4.96e-12 | 5.08e-12 | 8.09e-12 | 1.0000000000 |
| sgolay_D | Nuf | adiabatic | 70 | 399x604 | 12.185 | 12.185 | 9.90e-11 | 2.79e-13 | 1.67e-13 | 3.01e-12 | 1.0000000000 |
| IMG_POD_gauss_sgolay | Nu_snap | adiabatic | 70 | 399x604x100 | -0.23252 | -0.23252 | 1.01e-09 | 6.82e-12 | 7.07e-12 | 8.16e-05 | 1.0000000000 |
| TS_POD_gauss_sgolay | time_series | ambient | 70 | 1998 | 7.53929e-15 | 2.13376e-15 | 4.05e-10 | 2.69e-12 | 2.31e-12 | 9.16e-09 | 1.0000000000 |
| POD_wiener_sgolay | Num | adiabatic | 70 | 399x604 | 37.9939 | 37.9939 | 2.34e-10 | 1.93e-13 | 7.04e-14 | 2.98e-12 | 1.0000000000 |
| POD_wiener_sgolay | Nuf | adiabatic | 70 | 399x604 | 12.1993 | 12.1993 | 1.09e-10 | 7.15e-13 | 4.88e-13 | 5.54e-12 | 1.0000000000 |
| POD_wiener_sgolay | dT | adiabatic | 70 | 399x604 | 0.850312 | 0.850312 | 3.25e-12 | 7.29e-13 | 4.75e-13 | 5.61e-12 | 1.0000000000 |
| TS_POD_W_sgolay | time_series | ambient | 70 | 1998 | -3.88914e-14 | 3.04417e-15 | 2.89e-09 | 1.92e-11 | 1.59e-11 | 7.76e-08 | 1.0000000000 |
| wiener_sgolay | Num | adiabatic |  | 399x604 | 37.9939 | 37.9939 | 2.34e-10 | 1.93e-13 | 7.03e-14 | 2.93e-12 | 1.0000000000 |
| wiener_sgolay | Nuf | adiabatic |  | 399x604 | 15.7631 | 15.7631 | 1.18e-10 | 5.21e-13 | 3.19e-13 | 4.02e-12 | 1.0000000000 |
| wiener_sgolay | dT | adiabatic |  | 399x604 | 1.1209 | 1.1209 | 3.34e-12 | 5.08e-13 | 3.04e-13 | 3.42e-12 | 1.0000000000 |
| gauss_sgolay | Num | adiabatic |  | 399x604 | 37.8439 | 37.8439 | 2.34e-10 | 1.99e-13 | 7.85e-14 | 2.99e-12 | 1.0000000000 |
| gauss_sgolay | Nuf | adiabatic |  | 399x604 | 16.2385 | 16.2385 | 1.13e-10 | 2.49e-13 | 1.29e-13 | 2.99e-12 | 1.0000000000 |
| gauss_sgolay | dT | adiabatic |  | 399x604 | 1.15916 | 1.15916 | 8.74e-13 | 1.48e-13 | 1.01e-13 | 9.20e-13 | 1.0000000000 |
| mPOD_sgolay | Num | adiabatic | 48 | 399x604 | 37.8597 | 37.8597 | 2.64e-10 | 4.10e-13 | 2.99e-13 | 3.26e-12 | 1.0000000000 |
| mPOD_sgolay | Nuf | adiabatic | 48 | 399x604 | 10.4153 | 10.4153 | 1.69e-09 | 3.87e-11 | 2.64e-11 | 1.71e-10 | 1.0000000000 |
| mPOD_sgolay | dT | adiabatic | 48 | 399x604 | 0.732486 | 0.732486 | 1.34e-10 | 4.51e-11 | 2.67e-11 | 1.71e-10 | 1.0000000000 |
| mPOD_wiener_sgolay | Num | adiabatic | 48 | 399x604 | 38.0064 | 38.0064 | 2.53e-10 | 4.09e-13 | 2.99e-13 | 3.22e-12 | 1.0000000000 |
| mPOD_wiener_sgolay | Nuf | adiabatic | 48 | 399x604 | 10.3545 | 10.3545 | 1.65e-09 | 3.92e-11 | 2.67e-11 | 1.65e-10 | 1.0000000000 |
| mPOD_wiener_sgolay | dT | adiabatic | 48 | 399x604 | 0.725649 | 0.725649 | 1.31e-10 | 4.55e-11 | 2.71e-11 | 1.65e-10 | 1.0000000000 |
| mPOD_sgolay_exactF | Num | adiabatic | 48 | 399x604 | 37.8597 | 37.8546 | 1.86e-01 | 1.22e-03 | 8.58e-04 | 7.08e-03 | 0.9999934096 |
| mPOD_sgolay_exactF | Nuf | adiabatic | 48 | 399x604 | 10.4153 | 12.8038 | 8.88e+00 | 2.44e-01 | 2.38e-01 | 1.00e+00 | 0.9880435008 |
| mPOD_sgolay_exactF | dT | adiabatic | 48 | 399x604 | 0.732486 | 0.897714 | 4.16e-01 | 2.43e-01 | 2.38e-01 | 1.00e+00 | 0.9272224129 |
| qjqk_POD_gauss_sgolay11 | Num | adiabatic | 70 | 399x604 | 37.8044 | 37.8044 | 2.60e-10 | 4.65e-13 | 3.67e-13 | 3.33e-12 | 1.0000000000 |
| qjqk_POD_gauss_sgolay11 | Nuf | adiabatic | 70 | 399x604 | 12.0822 | 12.0822 | 1.07e-10 | 5.15e-13 | 3.58e-13 | 4.12e-12 | 1.0000000000 |
| qjqk_POD_gauss_sgolay11 | qunsrel | adiabatic | 70 | 399x604 | 0.276593 | 0.276593 | 6.02e-13 | 4.29e-13 | 3.02e-13 | 3.26e-12 | 1.0000000000 |
| qjqk_POD_gauss_sgolay11 | qkrel | adiabatic | 70 | 399x604 | 0.00827196 | 0.00827196 | 3.39e-13 | 9.09e-12 | 1.41e-11 | 9.51e-11 | 1.0000000000 |
