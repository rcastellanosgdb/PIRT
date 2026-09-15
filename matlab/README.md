# PIRT - Processing InfraRed Thermography (MATLAB)
<img src="../logo/PIRT_color.png" alt="PIRT" title="PIRT - Processing InfraRed Thermography" width="300">

MATLAB toolbox for filtering infrared thermography snapshots, computing
convective heat-transfer maps with the heated-thin-foil technique and
estimating their uncertainty. A Python implementation with the same interface
and validated against this code is available in [`../python`](../python).

## Installation

Add the source folders to the MATLAB path:
```matlab
addpath('matlab/src'); addpath('matlab/utils');
```
Requirements: MATLAB R2022a or newer (the input parser uses `anynan`). Some filters need toolboxes:
`imgaussfilt3` (Image Processing) for the Gaussian filter, `lowpass`/`highpass`/
`bandpass` (Signal Processing) for the temporal cut-off filter, `smooth`
(Curve Fitting) for the `Elbow` truncation criterion and `normrnd`
(Statistics) for the Monte Carlo uncertainty.

## Examples / test cases

[Link to DropBox](https://www.dropbox.com/scl/fo/c9mxtyxs61hgn8up7iq38/h?rlkey=ic9glj0a8edfne8t4gikis2ae&dl=0)
— a fraction of the impinging sweeping jet dataset of Robledo et al. (2025).
Place the `.mat` files in `examples/Impinging_Sweeping_Jet/SJ_processing/` and
run `examples/Impinging_Sweeping_Jet/SJ_main.m`.

## Introduction

**PIRT** is an object-based toolbox: a `PIRT` object stores the temperature
maps and the configuration introduced with labels (as in MATLAB built-in
functions), and `obj = obj.go()` performs the computations. Two independent
computation blocks are available: the **filter module** and the **heat
transfer module** (with an optional uncertainty estimation).

Conventions used throughout the toolbox:

* Temperature matrices are `ny x nx` images or `ny x nx x nt` sequences
  (rows, columns, snapshots). **x runs along the columns** (spacing `dx`) and
  **y along the rows** (spacing `dy`); the third dimension is time (spacing
  `dt`).
* Temperatures may be given in K or in °C: fields whose mean value is below
  100 are assumed to be in °C and converted to K (a warning is issued).
* Only `Thot` is filtered. `Tcold` is averaged over time (if 3-D) before the
  energy balance.

### Attributes

- `Thot`: 2-D or 3-D matrix of the heated (wall) temperature measurements.
- `Tcold`: 2-D or 3-D matrix of the adiabatic-wall (unheated) temperatures; same image size as `Thot`.
- `Calculate_filter`: boolean flag that controls the utilization of the filter module.
- `filter_params`: struct with the parameters of the filter module.
- `CalculateHeatTransfer`: boolean flag that controls the utilization of the heat transfer module.
- `CalculateHeatTransferError`: boolean flag that controls the uncertainty computation.
- `CalculateHeatTransferErrorMethod`: method used to estimate the uncertainty.
- `HeatTransfer_params`: struct with the parameters of the heat transfer module.
- `cropping_points`: 2x2 matrix with the cropping limits.
- `output`: controls whether results are kept in memory or written to files.
- `result`: struct with the outputs.

## Initialization

### Temperature maps
```matlab
pirt_object = PIRT('Thot', Thot_matrix, 'Tcold', Tcold_matrix)
```
The order of the labels does not matter, as long as the value follows its label.

An optional `'Crop'` input restricts the images to a region: a 2x2 matrix with
the x limits (columns) in the first row and the y limits (rows) in the second,
1-based and inclusive: `PIRT(..., 'Crop', [x1 x2; y1 y2])`.

### Filter module
```matlab
pirt_object = PIRT('Thot', Thot, 'Filter', Filter_info)
```
`Filter_info` is a struct (or struct array, for several filters applied
sequentially) with the fields `Type` and `Parameters`:
```matlab
%% One filter
Filter_info.Type = 'POD';
Filter_info.Parameters.Criterion = 'HardThreshold';
%% Several filters, applied one after the other
Filter(1).Type = 'POD';      Filter(1).Parameters.Criterion = 'HardThreshold';
Filter(2).Type = 'gaussian'; Filter(2).Parameters.FilterSize = [9,3,1]; Filter(2).Parameters.Sigma = [3,3,0.1];
Filter(3).Type = 'sgolay32'; Filter(3).Parameters.Kernel_size = [5,5,3]; Filter(3).Parameters.h = [dx,dy,dt];
pirt_object = PIRT('Thot', Thot, 'Filter', Filter);
```
Available filters and parameters (`Filter.Parameters.<name>`):

| `Type` | Description | Parameters |
|---|---|---|
| `'POD'` | Modal filter based on the proper orthogonal decomposition (`POD_filter.m`). | `Nmod` (number of modes) **or** `Criterion` (`'Spectrum'`, `'Elbow'`, `'HardThreshold'`) with `Threshold` (Spectrum/Elbow) or `beta` (aspect ratio `nt/(ny*nx)`, HardThreshold; defaults to the data aspect ratio). See `findNmod.m`. |
| `'mPOD'` | Multi-scale POD (Mendez et al., JFM 2019), `multiscale_POD_filter.m`. | `Type` (`'Peak Removal'`, `'Freq Deco'` or `'Both'`), `f_acquisition` [Hz], `Threshold` (Elbow criterion, optional). Peak removal: `fpeaks` [Hz], `w` (half window [Hz]). Frequency decoupling: `N_regions`, `f_min`, `f_max` [Hz]. `legacy_transform` (default `false`; `true` reproduces the back-transformation of the code used in the paper, see `CHANGELOG.md`). |
| `'sgolay32'` | Savitzky–Golay smoothing with a quadratic polynomial in a 3-D window; also provides the time derivative and the second spatial derivatives used in the energy balance (`sgolay32_filter.m`). | `Kernel_size` (odd, scalar or `[w_rows w_cols w_time]`), `h` (spacing `[dx dy dt]`). The first and last `floor(w_time/2)` snapshots are removed. |
| `'wiener3'` | 3-D pixel-wise adaptive Wiener filter (`wiener3.m`). | `kernel` (`[rows cols time]`, default 3), `noise` (noise power; estimated if omitted and returned in `result.noise_hot`). |
| `'gaussian'` | 3-D Gaussian filter (`imgaussfilt3`). | `Sigma` (scalar or `[rows cols time]`), `FilterSize` (odd), `Padding`, `FilterDomain`. |
| `'cutoff'` | Spectral low-/high-/band-pass filter in space and/or time (`Cutoff_Filter_3D.m`). | `Spatial`: struct with `fl` (low-pass cut-off) and/or `fh` (high-pass cut-off), normalised to the Nyquist frequency, scalar or `[f_rows f_cols]`. `Temporal`: struct with `fs` [Hz] and `fl` and/or `fh` [Hz]. |

The number of modes selected by the modal filters is stored in `result.Nmod_hot`.

### Heat transfer module

The convective heat-transfer coefficient is computed from the energy balance of the heated thin foil

$$h = \frac{q''_{j} - q''_{r} - q''_{k} - q''_{u} + \sum_{n=1}^{N} q''_{\mathrm{custom},n}}{T_\mathrm{w}-T_\mathrm{aw}}$$

with the Joule heating $q''_{j}=V I/A$, the radiative losses
$q''_{r}=n_\mathrm{sides}\,\sigma \epsilon (T_{w}^4 - T_{\infty}^4)$
(surroundings as a black body at the free-stream temperature), the tangential
conduction $q''_{k} = (k + k_p)\nabla^2 T_w$ for a foil — where `k = s·λ` and
`k_p = s_p·λ_p` are the conductances (thickness × conductivity) of foil and
paint — or $q''_k = (s\lambda_x + k_p)\partial^2_x T_w + (s\lambda_y + k_p)\partial^2_y T_w$
for a PCB ([Torre et al., 2018](https://doi.org/10.1016/j.ijheatmasstransfer.2018.06.106)),
and the unsteady term $q''_{u} = (\rho c_p s + \rho_p c_{p,p} s_p)\,\partial T_w/\partial t$.
$T_\mathrm{aw} = T_\mathrm{cold}\,T_{\infty,hot}/T_{\infty,cold}$.
Nu $= h L/k_{air}(T_{film})$ with $T_{film} = (T_w + T_{\infty,hot})/2$ and
St $= h/(\rho_\infty c_{p,\infty} U_\infty)$.

The label `'CalculateHeatTransfer'` activates the module, followed by any
combination of `'h'`, `'Nu'` and `'St'`:
```matlab
pirt_object = PIRT(..., 'CalculateHeatTransfer', 'h', 'Nu', 'St', 'HFS', HFS, 'Conditions', Conditions)
```
The flags `'TimeDer'` and `'SpatialDer'` include the unsteady and tangential
terms. Their derivatives are taken from the `sgolay32` filter if it is part of
the filter chain; otherwise they are computed with second-order finite
differences (`Derivative_FD.m`), which requires `Conditions.dt` and
`Conditions.dx`, `Conditions.dy`.

`HFS` (heated foil sensor) is a struct with the fields:
- `HFS.s`: thickness [m].
- `HFS.rho`: density [kg/m³].
- `HFS.cp`: specific heat [J/(kg K)].
- `HFS.epsilon`: emissivity of the face oriented towards the camera [-].
- `HFS.Area` [m²], or `HFS.H` and `HFS.W` [m] for a rectangular sensor.
- `HFS.k`: thermal **conductance** of the foil, `s·λ` [W/K] (foil model).
- `HFS.Type`: `'Foil'` (default) or `'PCB'`.
- `HFS.lambdax`, `HFS.lambday`: thermal conductivities [W/(m K)] in x and y (PCB model only).
- `HFS.sides` (optional, default 1): number of faces exchanging heat by radiation.
- `HFS.s_paint`, `HFS.rho_paint`, `HFS.cp_paint`, `HFS.lambda_paint` (optional): paint layer properties; defaults 21.81 µm, 1300 kg/m³, 5000 J/(kg K), 1.38 W/(m K) with a warning.

`Conditions` is a struct with the fields:
- `Conditions.L`: characteristic length [m].
- `Conditions.V`, `Conditions.I`: voltage [V] and current [A] supplied to the foil.
- `Conditions.Tamb`: `[T_cold, T_hot]` free-stream temperatures during the cold and hot acquisitions.
- `Conditions.dt` [s], `Conditions.dx`, `Conditions.dy` [m]: resolutions (finite-difference derivatives).
- `Conditions.Uinf` [m/s], `Conditions.rhoinf` [kg/m³], `Conditions.cpinf` [J/(kg K)]: free-stream velocity, density and specific heat (Stanton number only).

#### Custom heat balance terms
Extra heat fluxes are introduced with the label `'CustomQ'` followed by a cell
(even for a single term). Each element can be a scalar, a 2-D matrix with the
image size or a 3-D matrix with the size of the (filtered) `Thot`.

#### Output to files
`PIRT(..., 'Output', 'file', path)` writes the large results (`Thot_filtered`,
derivatives, `h`, `Nu`, `St`) to `.mat` files in `path` instead of keeping them
in `result`.

### Uncertainty estimation
```matlab
pirt_object = PIRT(..., 'CalculateHeatTransferError', 'Montecarlo', 'Error', Error_struct)
```
A Monte Carlo propagation (Minkina & Dudzik, 2009) perturbs every input of the
balance with a normal distribution and evaluates the spatially averaged `h`,
`Nu` and/or `St` of each sample (`result.errorh`, `result.errorNu`,
`result.errorSt`, vectors with one value per sample) and their percentage
deviation from the deterministic result (`result.error*_p`). The `'Moffat'`
method is not implemented (it falls back to Monte Carlo). Required fields of
`Error`:
- `Error.errorT`, `Error.errorTamb`: **absolute** uncertainties [K] of the temperature maps and of `Tamb`.
- Relative uncertainties: `errorV`, `errorI`, `errorEpsilon`, `errorrho`, `errorcp`, `errors`, `errorA`, `errorkplate`, `errorLchar`, `errork` (air conductivity), `errors_paint`, `errorcp_paint`, `errorlambda_paint`, `errorrho_paint`.
- `errordTdt` (with `'TimeDer'`), `errord2Tdx2`, `errord2Tdy2` (with `'SpatialDer'`).
- `errorUinf`, `errorrhoinf`, `errorcpinf` (Stanton number).
- `Error.samples` (optional, default 1000): number of samples.

## Results
`pirt_object.result` contains (when computed): `Thot_new` (filtered `Thot`),
`Nmod_hot`, `dTdt_hot`, `d2Tdx2_hot`, `d2Tdy2_hot`, `noise_hot`, `h`, `Nu`,
`St`, `errorh`, `errorh_p`, `errorNu`, `errorNu_p`, `errorSt`, `errorSt_p`.

## Tests
`tests/TestPIRTSynthetic.m` contains unit tests on synthetic data
(`runtests('tests/TestPIRTSynthetic.m')` from the `matlab` folder). The results
of this toolbox on the dataset of Robledo et al. (2025) are the reference of
the cross-validation of the Python version, see [`../validation`](../validation).

## Changes with respect to the version used in the paper
See [`../CHANGELOG.md`](../CHANGELOG.md).
