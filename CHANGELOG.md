# Changelog

## v1.1 (September 2026) — repository reorganisation and Python version

The code used in Robledo et al., Exp. Therm. Fluid Sci. 169 (2025) 111526
corresponds to commit `fa07252` (MATLAB only, `src/` and `utils/` at the
repository root). This version reorganises the repository into independent
`matlab/` and `python/` folders and adds a Python implementation of the whole
toolbox, cross-validated against the MATLAB results of that paper
(`validation/`).

### MATLAB toolbox (`matlab/`)

Reorganisation:

* `src/`, `utils/` and `examples/` moved to `matlab/`; the example script uses
  paths relative to its own location.
* `Validation_Polyfilter/` (legacy comparison of the polynomial filter) moved
  out of the class folder `@PIRT` to `matlab/legacy/`, together with the
  unrelated script `gaussianfilter.m`.
* New `matlab/tests/TestPIRTSynthetic.m` (MATLAB unit tests on synthetic data;
  written without access to a MATLAB licence, see the file header).

Fixes (all documented in the corresponding function headers):

* `Derivative_FD.m`: the finite-difference derivatives treated the **rows** as
  the x direction and the **columns** as y, i.e. the opposite convention to
  `sgolay32_filter.m` and to the `Crop` input of PIRT (x = columns). The
  convention is now x = columns (spacing `dx`), y = rows (spacing `dy`) in all
  the toolbox. Results were only affected with `dx ~= dy` or with the
  anisotropic PCB model. The temporal-size check `size(T,1)<3` was corrected to
  `size(T,3)<3`.
* Stanton number: `St = h/(rhoinf*cp*Uinf)` used the specific heat of the
  **foil** (`HFS.cp`); it now uses the fluid properties
  `Conditions.rhoinf`, `Conditions.cpinf` (new required input for St) and
  `Conditions.Uinf`, both in `Calculate_HeatTransfer.m` and in the Monte Carlo
  estimation (`errorcpinf` input).
* Monte Carlo (`Calculate_HeatTransfer_Error.m`): the radiative term did not
  include the number of exposed sides (`HFS.sides`), the y contribution of the
  tangential term had the wrong sign (`- kx*d2Tdx2 + ky*d2Tdy2`), and the film
  temperature for `k_air` was `(Thot+Tcold)/2` instead of the `(Thot+Tamb)/2`
  used by `Calculate_HeatTransfer.m`. All three are now consistent with the
  deterministic balance.
* `multiscale_POD_filter.m`: the filtered correlation matrix was transformed
  back to the time domain with the **forward** DFT matrix applied on both
  sides (`K = real(F*Khat*F)` instead of `real(F'*Khat*F')`). Since `F*F` is
  the index-reversal permutation, this returned the filtered correlation
  matrix with its time index reversed: same eigenvalues, but time-reversed
  temporal modes onto which the data were projected. On the sweeping-jet
  dataset of the paper the exact inverse leaves the mean Nusselt number
  unchanged (0.1 %) but increases the fluctuating Nusselt number and the
  mean |dT/dt| by about 23 % (`validation/REPORT.md`, case
  `mPOD_sgolay_exactF`), i.e. the legacy transform damped the reconstructed
  fluctuations. The exact inverse is now the default; the previous behaviour
  is available with `Filter.Parameters.legacy_transform = true` (MATLAB) or
  `legacy_transform=True` (Python) to reproduce the mPOD results of the paper.
* Spectral cut-off filters: `lowpass_kernel.m`/`highpass_kernel.m` built the
  frequency grid as `(1:m)-m/2`, offset by one bin with respect to the
  zero-frequency location of `fftshift`, and mirrored one half of the mask;
  the resulting mask was not Hermitian-symmetric and the filtered field was
  made real with `abs()`, which destroys the sign of zero-mean (high-pass)
  fields. Both kernels now share `elliptic_mask.m` (grid centred on the DC bin,
  exact symmetry) and `Spatial_Cutoff_Filter.m` takes the real part.

Known behaviours that were **kept** for backward compatibility (documented,
also reproduced by the Python version):

* `findNmod` with the `Elbow` criterion discards the first 5 eigenvalue
  ratios and returns the index relative to the truncated vector (not shifted
  back); with `HardThreshold` the returned `Nmod` is the index of the first
  singular value *below* the threshold (that mode is included).
* `findNmod`/`find_nmod` quirks listed above.
* Only `Thot` is filtered; `Tcold` is time-averaged before the balance (the
  `selection` parameter mentioned in old documentation was never implemented).
* `errorT` and `errorTamb` are **absolute** uncertainties [K]; the remaining
  `error*` inputs are relative.

Note on the film temperature: commit `5a4dd70` (2025-02-21) changed the film
temperature used for `k_air` in the Nusselt number from `(Thot+Tcold)/2` to
`(Thot+Tamb_hot)/2`. Most result files of the paper (January 2025) were
produced with the former definition, which is also the one written in the
paper (`T_film = (T_w + T_aw)/2`); the two differ by about 0.06 % in Nu for
that dataset. The Python version exposes the choice through the
`film_temperature` option (`'ambient'`, default and current MATLAB behaviour,
or `'adiabatic'`).

### Python package (`python/`)

First release (`pirt` 1.0.0): `PIRT` class with the same label-based
interface, filters (`pod_filter`, `mpod_filter`, `sgolay32_filter`, `wiener3`,
`gaussian_filter3`, cut-off filters), `derivative_fd`,
`calculate_heat_transfer`, `montecarlo_uncertainty` (with the additional
`mode='snapshots'` propagation used for the fluctuation uncertainty of the
paper) and I/O helpers for MATLAB `.mat` files and the example test case. The
Savitzky–Golay filter uses an exact separable decomposition of the kernels
(~5x faster than the 3-D convolution, identical results to round-off).
