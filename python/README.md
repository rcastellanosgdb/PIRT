# PIRT - Processing InfraRed Thermography (Python)

`pirt` is the Python implementation of the PIRT toolbox for infrared
thermography measurements of convective heat transfer with the heated-thin-foil
technique: filtering of temperature snapshot sequences (POD, multi-scale POD,
Savitzky–Golay, Wiener, Gaussian, spectral cut-off), energy balance of the foil
(h, Nu, St, with unsteady and tangential-conduction terms) and Monte Carlo
uncertainty estimation. It mirrors the MATLAB toolbox in [`../matlab`](../matlab)
(same inputs, same algorithms) and reproduces its results to machine precision
(see [`../validation`](../validation)).

## Installation

```bash
pip install "git+https://github.com/rcastellanosgdb/PIRT.git#subdirectory=python"
# or, from a clone of the repository
pip install -e python/            # add [io] for h5py (MATLAB v7.3 files), [plot], [test]
```
Requirements: Python ≥ 3.9, NumPy ≥ 1.22, SciPy ≥ 1.8; `h5py` to read MATLAB
v7.3 `.mat` files; `matplotlib` for the example plots.

## Quick start

```python
import numpy as np
from pirt import PIRT, io

case = io.load_sj_case("SJ_processing/")          # example dataset (see below)
Thot, Tcold = case["Thot"], case["Tcold"]         # (ny, nx, nt) and (ny, nx)
dx, dy, dt = case["dx"], case["dy"], 1 / case["f_acq"]

obj = PIRT(
    Thot=Thot, Tcold=Tcold,
    filters=[
        {"type": "POD", "criterion": "HardThreshold", "beta": 2000 / (399 * 604)},
        {"type": "gaussian", "filter_size": (9, 3, 1), "sigma": (3, 3, 0.1)},
        {"type": "sgolay32", "kernel_size": (5, 5, 3), "h": (dx, dy, dt)},
    ],
    heat_transfer=("Nu", "h"), time_der=True, spatial_der=True,
    HFS=case["HFS"], conditions=case["Conditions"],
)
obj.go()
Nu = obj.result["Nu"]                              # (ny, nx, nt - 2)
Nu_mean = Nu.mean(axis=2)
Nu_fluct = np.abs(Nu - Nu_mean[:, :, None]).mean(axis=2)
```
A complete example is given in `examples/Impinging_Sweeping_Jet/sj_main.py`.
The test case (a fraction of the sweeping-jet dataset of Robledo et al., 2025)
can be downloaded from the
[DropBox link](https://www.dropbox.com/scl/fo/c9mxtyxs61hgn8up7iq38/h?rlkey=ic9glj0a8edfne8t4gikis2ae&dl=0).

## Conventions

* Arrays are `(ny, nx)` images or `(ny, nx, nt)` sequences: **rows = y**
  (spacing `dy`), **columns = x** (spacing `dx`), third axis = time (`dt`).
  Vector-valued filter parameters follow the array-axis order
  `(rows, columns, time)`; the Savitzky–Golay spacing is `h = (dx, dy, dt)`.
* Temperatures in K. Fields with a mean below 100 are assumed to be in °C and
  converted (with a warning), as in MATLAB.
* Only `Thot` is filtered; `Tcold` is time-averaged before the balance.
* `crop=((x1, x2), (y1, y2))` uses 0-based, half-open ranges (Python slicing);
  the MATLAB `[x1 x2; y1 y2]` (1-based, inclusive) corresponds to
  `((x1-1, x2), (y1-1, y2))`.
* Filter specifications and the `HFS`/`conditions`/`error` dictionaries use the
  **MATLAB field names**. Filter parameter names are matched ignoring case and
  underscores (`Kernel_size`, `FilterSize`, `f_acquisition`, ...), and the
  MATLAB-like nested form `{"Type": "POD", "Parameters": {...}}` is accepted.

## The `PIRT` class

```python
PIRT(Thot, Tcold=None, *, filters=None, crop=None, heat_transfer=(), time_der=False,
     spatial_der=False, HFS=None, conditions=None, custom_q=None, error=None,
     error_method=None, montecarlo_mode="mean", film_temperature="ambient",
     output_dir=None, verbose=True, seed=None)
obj.go()          # runs filtering -> heat transfer -> uncertainty
obj.result        # dict with the outputs
```
MATLAB-style keywords are also accepted (`Filter=`, `CalculateHeatTransfer=`,
`TimeDer=`, `SpatialDer=`, `Conditions=`, `Error=`, `CustomQ=`, `Crop=`,
`CalculateHeatTransferError=`).

### Filters (`filters=[{...}, ...]`, applied sequentially)

| `type` | Function | Parameters |
|---|---|---|
| `POD` | `pod_filter` | `nmod`, or `criterion` (`Spectrum`, `Elbow`, `HardThreshold`) with `threshold` / `beta` |
| `mPOD` | `mpod_filter` | `type`/`mode` (`Peak Removal`, `Freq Deco`, `Both`), `f_acquisition`, `fpeaks`, `w`, `N_regions`, `f_min`, `f_max`, `threshold`, `nmod`, `legacy_transform` |
| `sgolay32` | `sgolay32_filter` | `kernel_size` (odd, `(rows, cols, time)`), `h` = `(dx, dy, dt)` |
| `wiener3` | `wiener3` | `kernel`, `noise` |
| `gaussian` | `gaussian_filter3` | `sigma`, `filter_size`, `padding`, `filter_domain` |
| `cutoff` | `cutoff_filter_3d` | `spatial={"fl":..., "fh":...}` (normalised to Nyquist), `temporal={"fs":..., "fl":..., "fh":...}` [Hz] |

`fl` is the cut-off of a low-pass (frequencies below are kept) and `fh` the
cut-off of a high-pass; if both are given the band between them is kept.

### Heat transfer (`heat_transfer=("h", "Nu", "St")`)

`HFS` keys: `s`, `rho`, `cp`, `epsilon`, `Area` (or `H` and `W`), `k`
(**conductance** thickness × conductivity, W/K; foil model), `Type`
(`'Foil'`/`'PCB'`), `lambdax`, `lambday` (PCB), `sides` (default 1),
`s_paint`, `rho_paint`, `cp_paint`, `lambda_paint` (defaults with a warning).

`conditions` keys: `L` (characteristic length), `V`, `I`, `Tamb=(T_cold, T_hot)`,
`dt`, `dx`, `dy` (finite-difference derivatives when no `sgolay32` filter is
used), `Uinf`, `rhoinf`, `cpinf` (Stanton number).

`custom_q`: list of extra heat fluxes (scalars, `(ny, nx)` or `(ny, nx, nt)` arrays).

`film_temperature`: `'ambient'` (default; `T_film = (T_w + T_inf)/2`, current
MATLAB behaviour) or `'adiabatic'` (`(T_w + T_aw)/2`, definition written in the
paper and used by the MATLAB runs before 2025‑02‑21).

`output_dir`: if given, the large arrays are saved as `.npy` files and kept in
`result` as memory maps (equivalent of the MATLAB `'Output','file'` option).

### Uncertainty (`error_method="montecarlo"`, `error={...}`)

Keys of `error` (MATLAB names): `errorT`, `errorTamb` (**absolute**, K);
relative: `errorV`, `errorI`, `errorEpsilon`, `errorrho`, `errorcp`, `errors`,
`errorA`, `errorkplate`, `errorLchar`, `errork`, `errors_paint`,
`errorcp_paint`, `errorlambda_paint`, `errorrho_paint`, `errordTdt`,
`errord2Tdx2`, `errord2Tdy2`, `errorUinf`, `errorrhoinf`, `errorcpinf`;
`samples` (default 1000). `montecarlo_mode='mean'` perturbs the time-averaged
maps (MATLAB behaviour); `'snapshots'` perturbs the whole sequence and also
returns the uncertainty of the mean absolute fluctuation (`errorNuf`,
`errorNuf_p`), as used for Table 2 of Robledo et al. (2025).

### Results (`obj.result`)

`Thot_new`, `Nmod_hot`, `dTdt_hot`, `d2Tdx2_hot`, `d2Tdy2_hot`, `noise_hot`,
`h`, `Nu`, `St`, `errorh`, `errorh_p`, `errorNu`, `errorNu_p`, `errorSt`,
`errorSt_p` (+ `errorNuf`, `errorNuf_p`, ... in `'snapshots'` mode).

## Functional API

All the building blocks are available as functions:
`pirt.pod_filter`, `pirt.find_nmod`, `pirt.optimal_svht_coef`,
`pirt.mpod_filter`, `pirt.sgolay32_filter`, `pirt.sgolay32_coef`,
`pirt.wiener3`, `pirt.gaussian_filter3`, `pirt.spatial_cutoff_filter`,
`pirt.temporal_cutoff_filter`, `pirt.cutoff_filter_3d`, `pirt.derivative_fd`,
`pirt.calculate_heat_transfer`, `pirt.montecarlo_uncertainty`, and the I/O
helpers `pirt.io.load_mat` (v5/v7 and v7.3 files), `pirt.io.read_conf_data`,
`pirt.io.load_sj_case`.

## Notes on the equivalence with the MATLAB toolbox

* Filters, derivatives and the energy balance reproduce the MATLAB results to
  round-off (differences ≤ 1e‑10 relative on the paper dataset). Random-based
  results (Monte Carlo) agree statistically.
* `sgolay32_filter` applies the kernels through their exact separable
  decomposition (`method='auto'`); `method='direct'` performs the 3-D
  convolutions as MATLAB does (identical results, slower).
* The temporal cut-off filter uses a Kaiser-window FIR designed with the same
  specifications as MATLAB's `lowpass`/`highpass`/`bandpass` defaults (60 dB
  stopband, steepness 0.85) but not the identical equiripple design.
* `mpod_filter` transforms the filtered correlation matrix back with the
  exact inverse DFT (default). `legacy_transform=True` reproduces the
  back-transformation of the MATLAB code used for the paper (forward DFT
  applied twice, which time-reverses the temporal modes and damps the
  reconstructed fluctuations by ~23 % on the paper dataset).
* `find_nmod` reproduces the MATLAB indexing of the `Elbow` and
  `HardThreshold` criteria (see its docstring).
* The Stanton number uses the fluid properties (`rhoinf`, `cpinf`, `Uinf`).

## Tests

```bash
pip install -e "python/[test]"
pytest python/tests
```
The tests check the filters, derivatives and the energy balance against
analytical and brute-force references. The cross-validation against the MATLAB
results of Robledo et al. (2025) (full 399×604×2000 dataset) is documented in
[`../validation`](../validation).

## License

GNU General Public License v3 — I. Robledo, J. Alfaro, R. Castellanos,
Universidad Carlos III de Madrid.
