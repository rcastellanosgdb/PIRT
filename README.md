# PIRT - Processing InfraRed Thermography
<img src="logo/PIRT_color.png" alt="PIRT logo" title="PIRT - Processing InfraRed Thermography" width="300">

**PIRT** is a toolbox for processing infrared thermography (IRT) measurements of
convective heat transfer with the **heated-thin-foil** technique: it filters
sequences of temperature maps (POD, multi-scale POD, Savitzky–Golay, Wiener,
Gaussian and spectral cut-off filters), solves the energy balance of the foil
to obtain the convective heat-transfer coefficient, the Nusselt number and the
Stanton number (including unsteady and tangential-conduction terms), and
estimates their uncertainty with a Monte Carlo method.

The toolbox is available in two self-contained implementations that share the
same interface (label-based inputs), the same algorithms and the same results:

| | MATLAB | Python |
|---|---|---|
| Folder | [`matlab/`](matlab/) | [`python/`](python/) |
| Documentation | [`matlab/README.md`](matlab/README.md) | [`python/README.md`](python/README.md) |
| Requirements | MATLAB R2022a+ (Image Processing, Signal Processing, Statistics, Curve Fitting toolboxes for some filters) | Python ≥ 3.9, NumPy, SciPy (h5py for MATLAB v7.3 files) |
| Entry point | `PIRT(...)` class in `matlab/src/@PIRT` | `pirt.PIRT` class |

The two versions have been cross-validated against the results published in
Robledo et al. (2025) (see [`validation/`](validation/)): the Python
implementation reproduces the MATLAB maps of the Nusselt number and of its
fluctuations to round-off level (RMS relative differences below 5e-11) for all
the filtering pipelines of that work.

## Getting one version only

Each folder is independent and can be used on its own:

* **MATLAB**: download the `matlab/` folder (or clone the repository) and add
  `matlab/src` and `matlab/utils` to the MATLAB path.
* **Python**: install directly from the repository subfolder
  ```bash
  pip install "git+https://github.com/rcastellanosgdb/PIRT.git#subdirectory=python"
  ```
  or clone the repository and run `pip install -e python/`.

To download a single folder without cloning everything you can use a sparse
checkout:
```bash
git clone --filter=blob:none --sparse https://github.com/rcastellanosgdb/PIRT.git
cd PIRT && git sparse-checkout set python   # or matlab
```

## Examples / test cases

[Link to DropBox](https://www.dropbox.com/scl/fo/c9mxtyxs61hgn8up7iq38/h?rlkey=ic9glj0a8edfne8t4gikis2ae&dl=0)

The test case is a fraction of the dataset of the impinging sweeping jet
experiment of Robledo et al. (2025). Example scripts that process it are given
in `matlab/examples/Impinging_Sweeping_Jet/SJ_main.m` and
`python/examples/Impinging_Sweeping_Jet/sj_main.py`.

## Repository layout

```
PIRT/
├── matlab/        MATLAB toolbox (src/, utils/, examples/, tests/, legacy/)
├── python/        Python package `pirt` (pip-installable), examples/, tests/
├── validation/    Cross-validation of both implementations against the paper results
├── logo/
├── CHANGELOG.md   Changes with respect to the code used in the paper (commit fa07252)
└── LICENSE        GNU GPL v3
```

## Physical model

The convective heat-transfer coefficient is obtained from the energy balance of
the heated thin foil (Astarita & Carlomagno, 2012):

$$h = \frac{q''_j - q''_r - q''_k - q''_u + \sum_n q''_{\mathrm{custom},n}}{T_w - T_{aw}}$$

with the Joule heating $q''_j = V I / A$, the radiative losses
$q''_r = n_{\mathrm{sides}}\,\sigma \varepsilon (T_w^4 - T_\infty^4)$, the
tangential conduction $q''_k = (s k_f + s_p k_p)\nabla^2 T_w$ (or
$(s\lambda_x + s_p k_p)\partial^2_x T_w + (s\lambda_y + s_p k_p)\partial^2_y T_w$
for printed circuit boards), the unsteady term
$q''_u = [(\rho c)_f s + (\rho c)_p s_p]\,\partial T_w/\partial t$ and
$T_{aw} = T_{cold}\,T_{\infty,hot}/T_{\infty,cold}$. The Nusselt number is
$\mathrm{Nu} = h L / k_{air}(T_{film})$ and the Stanton number
$\mathrm{St} = h/(\rho_\infty c_{p,\infty} U_\infty)$.

## Reference

If you use PIRT please cite:

> I. Robledo, J. Alfaro, C. Sanmiguel Vila, R. Castellanos, *Unsteady convective
> heat transfer of an impinging sweeping jet: A discussion on the effect of
> spatiotemporal filtering*, Experimental Thermal and Fluid Science 169 (2025)
> 111526. https://doi.org/10.1016/j.expthermflusci.2025.111526

## Authors and license

I. Robledo, J. Alfaro, R. Castellanos — Universidad Carlos III de Madrid,
Experimental Aerodynamics and Propulsion Laboratory. Distributed under the GNU
General Public License v3 (see `LICENSE`).
