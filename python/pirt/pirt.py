"""The :class:`PIRT` driver class (port of the MATLAB ``@PIRT`` class).

Example
-------
>>> from pirt import PIRT
>>> obj = PIRT(Thot=Thot, Tcold=Tcold,
...            filters=[{"type": "POD", "criterion": "HardThreshold", "beta": 2000/(399*604)},
...                     {"type": "gaussian", "filter_size": (9, 3, 1), "sigma": (3, 3, 0.1)},
...                     {"type": "sgolay32", "kernel_size": (5, 5, 3), "h": (dx, dy, dt)}],
...            heat_transfer=("Nu", "h"), time_der=True, spatial_der=True,
...            HFS=hfs, conditions=conditions)
>>> obj.go()
>>> Nu = obj.result["Nu"]

The MATLAB label-based interface is mirrored: filter specifications can be
given either as flat dictionaries (``{"type": ..., <parameters>}``) or as
MATLAB-like ``{"Type": ..., "Parameters": {...}}`` dictionaries; parameter
names are matched case-insensitively and ignoring underscores (e.g.
``Kernel_size``, ``FilterSize``, ``f_acquisition``).
"""
from __future__ import annotations

import os
import warnings
from pathlib import Path

import numpy as np

from .filters import (cutoff_filter_3d, gaussian_filter3, mpod_filter, pod_filter,
                      sgolay32_filter, wiener3)
from .heat_transfer import calculate_heat_transfer
from .uncertainty import montecarlo_uncertainty

__all__ = ["PIRT"]

_FILTER_TYPES = {
    "pod": "POD", "mpod": "mPOD", "sgolay32": "sgolay32", "sgolay": "sgolay32",
    "wiener3": "wiener3", "wiener": "wiener3", "gaussian": "gaussian", "gauss": "gaussian",
    "cutoff": "cutoff", "cutofffilter": "cutoff",
}

# canonical parameter names per filter (keys: normalised = lower-case without underscores)
_PARAM_NAMES = {
    "POD": {"criterion": "criterion", "threshold": "threshold", "nmod": "nmod", "beta": "beta"},
    "mPOD": {"type": "mode", "mode": "mode", "fpeaks": "fpeaks", "fpeak": "fpeaks", "w": "w",
             "nregions": "n_regions", "nf": "n_regions", "fmin": "f_min", "fmax": "f_max",
             "facquisition": "f_acquisition", "facq": "f_acquisition", "fs": "f_acquisition",
             "threshold": "threshold", "nmod": "nmod", "legacytransform": "legacy_transform"},
    "sgolay32": {"kernelsize": "kernel_size", "w": "kernel_size", "framelen": "kernel_size",
                 "h": "h"},
    "wiener3": {"kernel": "kernel", "noise": "noise"},
    "gaussian": {"sigma": "sigma", "filtersize": "filter_size", "padding": "padding",
                 "filterdomain": "filter_domain"},
    "cutoff": {"spatial": "spatial", "temporal": "temporal"},
}

_ERROR_METHODS = {"montecarlo": "Montecarlo", "monte carlo": "Montecarlo", "moffat": "Moffat"}


def _norm_key(key: str) -> str:
    return str(key).replace("_", "").replace(" ", "").lower()


def _normalize_filter_spec(spec) -> dict:
    """Return ``{"type": <canonical>, "params": {<canonical>: value}}``."""
    if not isinstance(spec, dict):
        raise TypeError("Each filter must be a dict with a 'type' entry and its parameters")
    items = {_norm_key(k): (k, v) for k, v in spec.items()}
    if "type" not in items:
        raise ValueError("Each filter specification requires a 'type' entry")
    ftype_raw = items.pop("type")[1]
    ftype_key = _norm_key(ftype_raw)
    if ftype_key not in _FILTER_TYPES:
        raise ValueError(f"PIRT: '{ftype_raw}' is not a valid filter type "
                         f"(POD, mPOD, sgolay32, wiener3, gaussian, cutoff)")
    ftype = _FILTER_TYPES[ftype_key]
    raw_params = {}
    if "parameters" in items:  # MATLAB-like nested struct
        nested = items.pop("parameters")[1] or {}
        raw_params.update(nested)
    raw_params.update({orig: v for orig, v in items.values()})
    names = _PARAM_NAMES[ftype]
    params = {}
    for k, v in raw_params.items():
        nk = _norm_key(k)
        if nk not in names:
            raise ValueError(f"PIRT: unknown parameter '{k}' for the {ftype} filter "
                             f"(valid: {sorted(set(names.values()))})")
        params[names[nk]] = v
    return {"type": ftype, "params": params}


class PIRT:
    """Processing InfraRed Thermography: filtering + heat-transfer computation.

    Parameters
    ----------
    Thot : ndarray, shape (ny, nx) or (ny, nx, nt)
        Wall temperature maps with the foil heated (third dimension = time).
    Tcold : ndarray, shape (ny, nx) or (ny, nx, nt), optional
        Adiabatic-wall temperature maps (time-averaged before use).
    filters : dict or sequence of dict, optional
        Filters applied sequentially to ``Thot``. Each dict has a ``'type'``
        (``'POD'``, ``'mPOD'``, ``'sgolay32'``, ``'wiener3'``, ``'gaussian'``,
        ``'cutoff'``) and the filter parameters (see the filter functions).
    crop : ((x1, x2), (y1, y2)), optional
        Half-open, 0-based pixel ranges ``[x1, x2)`` (columns) and ``[y1, y2)``
        (rows) kept from the images. (MATLAB: ``[x1 x2; y1 y2]`` with 1-based
        inclusive indices, i.e. ``x1_py = x1_ml - 1`` and ``x2_py = x2_ml``.)
    heat_transfer : sequence of {'h', 'Nu', 'St'}, optional
        Quantities to compute. If empty the heat-transfer module is skipped.
    time_der, spatial_der : bool
        Include the unsteady and tangential-conduction terms.
    HFS, conditions : dict
        Sensor and test data (MATLAB field names, see
        :func:`pirt.heat_transfer.prepare_hfs` / ``prepare_conditions``).
    custom_q : sequence, optional
        Extra heat-flux terms [W/m^2] (scalars, 2-D or 3-D arrays).
    error : dict, optional
        Uncertainties for the Monte Carlo estimation (MATLAB field names).
    error_method : {'Montecarlo', 'Moffat'}, optional
        Requests the uncertainty estimation. ``'Moffat'`` is not implemented
        and falls back to ``'Montecarlo'`` with a warning (as in MATLAB).
    montecarlo_mode : {'mean', 'snapshots'}
        See :func:`pirt.uncertainty.montecarlo_uncertainty`.
    film_temperature : {'ambient', 'adiabatic'}
        Film temperature used for the air conductivity in Nu, see
        :func:`pirt.heat_transfer.calculate_heat_transfer`.
    output_dir : str or Path, optional
        If given, the large arrays produced (filtered temperatures,
        derivatives, h/Nu/St) are stored as ``.npy`` files in that folder and
        kept in ``result`` as read-only memory maps.
    verbose : bool
        Print progress messages.
    seed : int, optional
        Seed of the Monte Carlo random generator.

    Attributes
    ----------
    result : dict
        Outputs: ``Thot_new``, ``Nmod_hot``, ``dTdt_hot``, ``d2Tdx2_hot``,
        ``d2Tdy2_hot``, ``noise_hot``, ``h``, ``Nu``, ``St``, ``error*``.
    """

    def __init__(self, Thot=None, Tcold=None, *, filters=None, crop=None,
                 heat_transfer=(), time_der: bool = False, spatial_der: bool = False,
                 HFS: dict | None = None, conditions: dict | None = None, custom_q=None,
                 error: dict | None = None, error_method: str | None = None,
                 montecarlo_mode: str = "mean", film_temperature: str = "ambient",
                 output_dir=None, verbose: bool = True, seed=None, **kwargs):
        # MATLAB-style aliases
        aliases = {"Filter": "filters", "Crop": "crop", "CalculateHeatTransfer": "heat_transfer",
                   "TimeDer": "time_der", "SpatialDer": "spatial_der", "Conditions": "conditions",
                   "CustomQ": "custom_q", "Error": "error", "CalculateHeatTransferError": "error_method"}
        local = dict(filters=filters, crop=crop, heat_transfer=heat_transfer, time_der=time_der,
                     spatial_der=spatial_der, conditions=conditions, custom_q=custom_q, error=error,
                     error_method=error_method)
        for k, v in kwargs.items():
            if k not in aliases:
                raise TypeError(f"PIRT: unexpected keyword argument '{k}'")
            local[aliases[k]] = v

        self.verbose = bool(verbose)
        self.seed = seed
        self.output_dir = Path(output_dir) if output_dir is not None else None

        self.Thot = self._check_temperature(Thot, "Thot")
        self.Tcold = self._check_temperature(Tcold, "Tcold")
        if self.Thot is None and self.Tcold is None:
            raise ValueError("PIRT: input temperature is required. Either Thot or Tcold must be introduced")
        if self.Thot is None:
            raise ValueError("PIRT: Thot is required (Tcold alone cannot be processed)")

        self.cropping_points = self._check_crop(local["crop"])

        # Filter module
        flt = local["filters"]
        if flt is None:
            self.filter_params = []
        else:
            if isinstance(flt, dict):
                flt = [flt]
            self.filter_params = [_normalize_filter_spec(f) for f in flt]
        self.calculate_filter = len(self.filter_params) > 0

        # Heat-transfer module
        ht = local["heat_transfer"]
        if ht is None:
            ht = ()
        if isinstance(ht, str):
            ht = (ht,)
        if isinstance(ht, bool):
            ht = ("h",) if ht else ()
        canon = {"h": "h", "nu": "Nu", "st": "St"}
        self.heat_transfer = []
        for q in ht:
            if str(q).lower() not in canon:
                raise ValueError(f"PIRT: '{q}' is not valid; use 'h', 'Nu' and/or 'St'")
            self.heat_transfer.append(canon[str(q).lower()])
        self.calculate_heat_transfer = len(self.heat_transfer) > 0
        self.time_der = bool(local["time_der"])
        self.spatial_der = bool(local["spatial_der"])
        self.HFS = dict(HFS) if HFS else None
        self.conditions = dict(local["conditions"]) if local["conditions"] else None
        self.custom_q = list(local["custom_q"]) if local["custom_q"] is not None else None
        if self.calculate_heat_transfer:
            if self.Tcold is None:
                raise ValueError("PIRT: in order to compute the heat transfer both Thot and Tcold must be introduced")
            if not self.HFS or not self.conditions:
                raise ValueError("PIRT: both the HFS and Conditions data must be introduced to calculate the heat transfer")
            if self.HFS.get("Type", "Foil") not in ("Foil", "PCB"):
                raise ValueError("PIRT: HFS 'Type' has to be either 'Foil' or 'PCB'")
            if "Type" not in self.HFS:
                warnings.warn("PIRT: the HFS will be treated as a uniform thin foil")
        elif self.verbose:
            print("Heat transfer will not be computed")

        # Uncertainty module
        em = local["error_method"]
        self.calculate_heat_transfer_error = em is not None and em is not False
        self.error = dict(local["error"]) if local["error"] else None
        self.error_method = None
        if self.calculate_heat_transfer_error:
            if not self.calculate_heat_transfer:
                raise ValueError("PIRT: the heat transfer must be computed to estimate its uncertainty")
            method = "Montecarlo" if em is True else _ERROR_METHODS.get(str(em).strip().lower())
            if method is None:
                warnings.warn("PIRT: unknown error method; the Montecarlo method will be used")
                method = "Montecarlo"
            if method == "Moffat":
                warnings.warn("PIRT: the Moffat method is not implemented, Montecarlo will be used")
                method = "Montecarlo"
            self.error_method = method
            if not self.error:
                raise ValueError("PIRT: the 'Error' data must be introduced to estimate the uncertainty")
        self.montecarlo_mode = montecarlo_mode
        if film_temperature not in ("ambient", "adiabatic"):
            raise ValueError("PIRT: film_temperature must be 'ambient' or 'adiabatic'")
        self.film_temperature = film_temperature

        self.result: dict = {}
        if self.verbose:
            print("PIRT: Information saved, run obj.go() to perform calculations")

    # ------------------------------------------------------------------ #
    @staticmethod
    def _check_temperature(T, name):
        if T is None:
            return None
        T = np.asarray(T)
        if not np.issubdtype(T.dtype, np.number):
            raise TypeError(f"PIRT: {name} must be a numeric array")
        if T.ndim not in (2, 3):
            raise ValueError(f"PIRT: {name} must be a 2-D or 3-D array")
        if np.issubdtype(T.dtype, np.integer):
            T = T.astype(np.float64)
        if np.isnan(T).any():
            warnings.warn(f"PIRT: {name} contains NaN elements")
        if np.isinf(T).any():
            warnings.warn(f"PIRT: {name} contains Inf elements")
        return T

    @staticmethod
    def _check_crop(crop):
        if crop is None:
            return None
        arr = np.asarray(crop)
        if arr.size != 4:
            raise ValueError("PIRT: crop must be ((x1, x2), (y1, y2)) with 0-based half-open ranges")
        arr = arr.reshape(2, 2).astype(int)
        if np.any(arr < 0) or arr[0, 1] <= arr[0, 0] or arr[1, 1] <= arr[1, 0]:
            raise ValueError("PIRT: invalid cropping ranges")
        return arr

    def _store(self, name, arr):
        """Keep ``arr`` in ``result`` (optionally offloaded to disk as .npy)."""
        if self.output_dir is not None and isinstance(arr, np.ndarray) and arr.size > 1:
            self.output_dir.mkdir(parents=True, exist_ok=True)
            fname = self.output_dir / f"{name}.npy"
            np.save(fname, np.ascontiguousarray(arr))
            if self.verbose:
                print(f"--> Saving {name} into {os.fspath(fname)}")
            arr = np.load(fname, mmap_mode="r")
        self.result[name] = arr
        return arr

    # ------------------------------------------------------------------ #
    def set_filters(self, filters, append: bool = False):
        """Replace (or append to) the filter chain."""
        if isinstance(filters, dict):
            filters = [filters]
        new = [_normalize_filter_spec(f) for f in filters]
        self.filter_params = self.filter_params + new if append else new
        self.calculate_filter = len(self.filter_params) > 0
        return self

    def set_heat_transfer_data(self, HFS=None, conditions=None, time_der=None, spatial_der=None):
        """Update the heat-transfer inputs."""
        if HFS is not None:
            self.HFS = dict(HFS)
        if conditions is not None:
            self.conditions = dict(conditions)
        if time_der is not None:
            self.time_der = bool(time_der)
        if spatial_der is not None:
            self.spatial_der = bool(spatial_der)
        return self

    # ------------------------------------------------------------------ #
    def _apply_filter(self, idx, spec, Thot):
        ftype, p = spec["type"], dict(spec["params"])
        v = self.verbose
        if ftype in ("POD", "mPOD") and Thot.ndim < 3:
            raise ValueError(f"PIRT: a third dimension is required to apply the {ftype} filter")
        if ftype == "POD":
            if v:
                print("-- POD filter Hot images")
            Thot, nmod = pod_filter(Thot, return_nmod=True, verbose=v, **p)
            self.result["Nmod_hot"][idx] = nmod
        elif ftype == "mPOD":
            if v:
                print("-- Multi-scale POD filter Hot images")
            if "f_acquisition" not in p:
                raise ValueError("PIRT: the acquisition frequency 'f_acquisition' is required for the mPOD filter")
            if "mode" not in p:
                raise ValueError("PIRT: the mPOD 'Type' (peak removal, frequency decoupling or both) is required")
            Thot, nmod = mpod_filter(Thot, return_nmod=True, verbose=v, **p)
            self.result["Nmod_hot"][idx] = nmod
        elif ftype == "sgolay32":
            if v:
                print("-- sgolay32 filter Hot images")
            Thot, dTdt, d2x, d2y = sgolay32_filter(Thot, **p)
            self._store("dTdt_hot", dTdt)
            self._store("d2Tdx2_hot", d2x)
            self._store("d2Tdy2_hot", d2y)
        elif ftype == "wiener3":
            if v:
                print("-- Wiener3 filter Hot images")
            Thot, noise = wiener3(Thot, return_noise=True, **p)
            if "noise" not in p:
                self.result["noise_hot"] = noise
        elif ftype == "gaussian":
            if v:
                print("-- Gaussian filter Hot images")
            Thot = gaussian_filter3(Thot, **p)
        elif ftype == "cutoff":
            if v:
                print("-- Cutoff filter Hot images")
            if "spatial" not in p and "temporal" not in p:
                raise ValueError("PIRT: 'Spatial' or 'Temporal' fields must be introduced in the Cutoff filter")
            Thot = cutoff_filter_3d(Thot, spatial=p.get("spatial"), temporal=p.get("temporal"))
        if v:
            print("--> Thot DONE")
        return Thot

    def go(self):
        """Run the configured modules (filtering, heat transfer, uncertainty)."""
        v = self.verbose
        self.result = {}
        Thot, Tcold = self.Thot, self.Tcold
        if Tcold is not None and Tcold.ndim == 3:
            Tcold = Tcold.mean(axis=2)
        if self.cropping_points is not None:
            (x1, x2), (y1, y2) = self.cropping_points
            Thot = Thot[y1:y2, x1:x2, ...]
            if Tcold is not None:
                Tcold = Tcold[y1:y2, x1:x2]
            else:
                warnings.warn("PIRT: only the hot images were introduced and cropped")

        if self.calculate_filter:
            if v:
                print("*************************************************************")
                print("******************** Temperature Filter *********************")
                print("*************************************************************")
            self.result["Nmod_hot"] = [None] * len(self.filter_params)
            for i, spec in enumerate(self.filter_params):
                Thot = self._apply_filter(i, spec, Thot)
            Thot = self._store("Thot_new", Thot)
        elif v:
            print("******************** No Filter Applied *********************")

        if self.calculate_heat_transfer:
            if v:
                print("*************************************************************")
                print("**************** Heat Transfer Calculation ******************")
                print("*************************************************************")
            dTdt = self.result.get("dTdt_hot") if self.time_der else None
            d2x = self.result.get("d2Tdx2_hot") if self.spatial_der else None
            d2y = self.result.get("d2Tdy2_hot") if self.spatial_der else None
            if self.time_der and dTdt is None and self.calculate_filter:
                warnings.warn("PIRT: unsteady term computed with finite differences (no sgolay32 filter was applied)")
            if self.spatial_der and d2x is None and self.calculate_filter:
                warnings.warn("PIRT: tangential term computed with finite differences (no sgolay32 filter was applied)")
            ht = calculate_heat_transfer(Thot, Tcold, self.HFS, self.conditions,
                                         compute=self.heat_transfer, time_der=self.time_der,
                                         spatial_der=self.spatial_der, dTdt=dTdt, d2Tdx2=d2x,
                                         d2Tdy2=d2y, custom_q=self.custom_q,
                                         film_temperature=self.film_temperature, verbose=v)
            for key, val in ht.items():
                self._store(key, val)

            if self.calculate_heat_transfer_error:
                if v:
                    print("*************************************************************")
                    print("************* Heat Transfer Error Calculation ***************")
                    print("*************************************************************")
                err = montecarlo_uncertainty(
                    Thot, Tcold, self.HFS, self.conditions, self.error,
                    compute=self.heat_transfer, time_der=self.time_der, spatial_der=self.spatial_der,
                    dTdt=self.result.get("dTdt_hot") if self.time_der else None,
                    d2Tdx2=self.result.get("d2Tdx2_hot") if self.spatial_der else None,
                    d2Tdy2=self.result.get("d2Tdy2_hot") if self.spatial_der else None,
                    reference={k: self.result[k] for k in ("h", "Nu", "St") if k in self.result},
                    mode=self.montecarlo_mode, seed=self.seed,
                    film_temperature=self.film_temperature, verbose=v)
                self.result.update(err)
                if v:
                    print("-->DONE")
        return self
