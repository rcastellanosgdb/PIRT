"""Monte Carlo uncertainty estimation of the heat-transfer maps.

Port of the ``montecarlo_method`` of ``Calculate_HeatTransfer_Error.m``
(Minkina & Dudzik, Infrared Thermography: Errors and Uncertainties, Wiley,
2009). Every input of the energy balance is perturbed with a normal
distribution whose standard deviation is the declared uncertainty:

* ``errorT`` and ``errorTamb`` are **absolute** uncertainties [K];
* all the other ``error*`` entries are **relative** uncertainties
  (``sigma = value * error``).

Two propagation modes are available:

* ``mode='mean'`` (default, as in the current MATLAB toolbox): the
  time-averaged temperature and derivative maps are perturbed; each sample
  yields a spatially averaged h / Nu / St.
* ``mode='snapshots'`` (used for the uncertainty figures of Robledo et al.,
  Exp. Therm. Fluid Sci. 169 (2025) 111526): the full snapshot sequence is
  perturbed in each sample, which additionally provides the uncertainty of the
  mean absolute fluctuation ``Nu' = <|Nu - <Nu>_t|>``. It is considerably more
  expensive in memory and time.
"""
from __future__ import annotations

import warnings

import numpy as np

from .heat_transfer import (STEFAN_BOLTZMANN, air_thermal_conductivity, check_celsius,
                            prepare_conditions, prepare_hfs)

__all__ = ["montecarlo_uncertainty", "REQUIRED_ERROR_KEYS"]

REQUIRED_ERROR_KEYS = ("errorT", "errorTamb", "errorV", "errorI", "errorEpsilon", "errorrho",
                       "errorcp", "errors", "errorA", "errorkplate", "errorLchar", "errork",
                       "errors_paint", "errorcp_paint", "errorlambda_paint", "errorrho_paint")


def _validate_error(error: dict, compute, time_der, spatial_der) -> dict:
    if not error:
        raise ValueError("PIRT: the 'Error' data must be introduced to estimate the heat transfer uncertainty")
    err = dict(error)
    if "errorAboard" in err and "errorA" not in err:
        err["errorA"] = err["errorAboard"]
    missing = [k for k in REQUIRED_ERROR_KEYS if k not in err]
    if missing:
        raise ValueError(f"PIRT: missing uncertainty inputs {missing}")
    if "st" in compute:
        for k in ("errorUinf", "errorrhoinf", "errorcpinf"):
            if k not in err:
                raise ValueError(f"PIRT: '{k}' is required to estimate the uncertainty of St")
    if time_der and "errordTdt" not in err:
        raise ValueError("PIRT: 'errordTdt' is required to estimate the uncertainty with the unsteady term")
    if spatial_der and ("errord2Tdx2" not in err or "errord2Tdy2" not in err):
        raise ValueError("PIRT: 'errord2Tdx2' and 'errord2Tdy2' are required with the tangential term")
    if "samples" in err and int(err["samples"]) != err["samples"]:
        raise ValueError("PIRT: the number of Monte Carlo samples must be an integer")
    return err


def montecarlo_uncertainty(Thot, Tcold, hfs: dict, conditions: dict, error: dict, *,
                           compute=("h",), time_der: bool = False, spatial_der: bool = False,
                           dTdt=None, d2Tdx2=None, d2Tdy2=None, reference: dict | None = None,
                           n_samples: int | None = None, mode: str = "mean", seed=None,
                           film_temperature: str = "ambient", verbose: bool = True) -> dict:
    """Monte Carlo propagation of the input uncertainties to h, Nu and St.

    Parameters
    ----------
    Thot, Tcold : ndarray
        Temperature maps used in the heat-transfer computation (i.e. after
        filtering). ``Tcold`` is time-averaged if 3-D.
    hfs, conditions : dict
        As in :func:`pirt.heat_transfer.calculate_heat_transfer`.
    error : dict
        Uncertainties, see module docstring. Optional key ``'samples'``.
    compute : sequence of {'h', 'Nu', 'St'}
    time_der, spatial_der : bool
        Include the unsteady / tangential terms (their derivative maps must be
        passed in ``dTdt`` / ``d2Tdx2``, ``d2Tdy2``; otherwise the terms are
        skipped with a warning, as in MATLAB).
    reference : dict, optional
        Deterministic result (``'h'``, ``'Nu'``, ``'St'`` maps) used to express
        the sampled means as percentage deviations (``error*_p``).
    n_samples : int, optional
        Number of samples (default ``error['samples']`` or 1000).
    mode : {'mean', 'snapshots'}
        Propagation mode, see module docstring.
    seed : int or numpy Generator, optional
        Random seed for reproducibility.

    Returns
    -------
    result : dict
        ``'errorh'``, ``'errorNu'``, ``'errorSt'``: arrays (n_samples,) with the
        spatially averaged quantity of every sample; ``'error*_p'``: absolute
        percentage deviation of each sample from the reference mean. In
        ``mode='snapshots'`` also ``'errorNuf'``/``'errorNuf_p'`` (and the
        equivalents for h and St).
    """
    compute = {c.lower() for c in ([compute] if isinstance(compute, str) else compute)}
    err = _validate_error(error, compute, time_der, spatial_der)
    hfs = prepare_hfs(hfs, spatial_der=spatial_der)
    cond = prepare_conditions(conditions, compute_st="st" in compute)
    mode = mode.lower()
    if mode not in ("mean", "snapshots"):
        raise ValueError("mode must be 'mean' or 'snapshots'")
    rng = seed if isinstance(seed, np.random.Generator) else np.random.default_rng(seed)

    Thot = np.asarray(Thot, dtype=np.float64)
    Tcold = np.asarray(Tcold, dtype=np.float64)
    if Tcold.ndim == 3:
        Tcold = Tcold.mean(axis=2)
    if mode == "mean" and Thot.ndim == 3:
        Thot = Thot.mean(axis=2)
    if mode == "snapshots" and Thot.ndim == 2:
        Thot = Thot[:, :, None]
    Thot = check_celsius(Thot, name="Thot")
    Tcold = check_celsius(Tcold, name="Tcold")
    Tamb_cold = float(check_celsius(cond["Tamb"][0], name="Tamb(cold)"))
    Tamb_hot = float(check_celsius(cond["Tamb"][1], name="Tamb(hot)"))

    def _prep_der(arr, name):
        if arr is None:
            return None
        arr = np.asarray(arr, dtype=np.float64)
        if mode == "mean":
            return arr.mean(axis=2) if arr.ndim == 3 else arr
        if arr.ndim == 2:
            arr = arr[:, :, None]
        if arr.shape != Thot.shape:
            raise ValueError(f"PIRT: {name} must have the same shape as Thot")
        return arr

    dTdt_m = _prep_der(dTdt, "dTdt") if time_der else None
    if time_der and dTdt_m is None:
        warnings.warn("PIRT: the unsteady term is not included in the uncertainty estimation "
                      "because dTdt was not provided")
    d2x_m = _prep_der(d2Tdx2, "d2Tdx2") if spatial_der else None
    d2y_m = _prep_der(d2Tdy2, "d2Tdy2") if spatial_der else None
    if spatial_der and (d2x_m is None or d2y_m is None):
        warnings.warn("PIRT: the tangential term is not included in the uncertainty estimation "
                      "because d2Tdx2/d2Tdy2 were not provided")
        d2x_m = d2y_m = None

    s, rho, cp, A, eps = hfs["s"], hfs["rho"], hfs["cp"], hfs["A"], hfs["epsilon"]
    sides = float(hfs["sides"])
    if hfs["Type"] == "PCB":
        kplatex, kplatey = float(hfs["lambdax"]) * s, float(hfs["lambday"]) * s
    else:
        kplatex = kplatey = hfs["k"]
    s_p, rho_p, cp_p, lam_p = hfs["s_paint"], hfs["rho_paint"], hfs["cp_paint"], hfs["lambda_paint"]
    L_char, V, I = cond["L_char"], cond["V"], cond["I"]
    sigma = STEFAN_BOLTZMANN

    # Nominal air conductivity at the film temperature (same model as the deterministic balance)
    if film_temperature == "ambient":
        k_air = air_thermal_conductivity((Thot + Tamb_hot) / 2.0)
    elif film_temperature == "adiabatic":
        Taw0 = Tcold * (Tamb_hot / Tamb_cold)
        k_air = air_thermal_conductivity((Thot + (Taw0[:, :, None] if Thot.ndim == 3 else Taw0)) / 2.0)
    else:
        raise ValueError("film_temperature must be 'ambient' or 'adiabatic'")

    n = int(n_samples if n_samples is not None else err.get("samples", 1000))
    if n_samples is None and "samples" not in err:
        warnings.warn("PIRT: default number of Monte Carlo samples taken as 1000")

    def nrm(mu, sig):
        return rng.normal(mu, np.abs(sig))

    def rel(mu, e):
        return nrm(mu, mu * e)

    e = err
    want_h, want_nu, want_st = "h" in compute, "nu" in compute, "st" in compute
    out_h, out_nu, out_st = np.empty(n), np.empty(n), np.empty(n)
    outf_h, outf_nu, outf_st = np.empty(n), np.empty(n), np.empty(n)

    for i in range(n):
        if verbose and (i % max(1, n // 10) == 0 or i == n - 1):
            print(f"Sample {i + 1}/{n}", end="\r", flush=True)
        c_Thot = nrm(Thot, e["errorT"])
        c_Tcold = nrm(Tcold, e["errorT"])
        c_Tamb_cold = nrm(Tamb_cold, e["errorTamb"])
        c_Tamb_hot = nrm(Tamb_hot, e["errorTamb"])
        c_V, c_I = rel(V, e["errorV"]), rel(I, e["errorI"])
        c_eps, c_rho, c_cp, c_s = rel(eps, e["errorEpsilon"]), rel(rho, e["errorrho"]), rel(cp, e["errorcp"]), rel(s, e["errors"])
        c_kx, c_ky = rel(kplatex, e["errorkplate"]), rel(kplatey, e["errorkplate"])
        c_A, c_L = rel(A, e["errorA"]), rel(L_char, e["errorLchar"])
        c_k = rel(k_air, e["errork"])
        c_cp_p, c_s_p = rel(cp_p, e["errorcp_paint"]), rel(s_p, e["errors_paint"])
        c_rho_p, c_lam_p = rel(rho_p, e["errorrho_paint"]), rel(lam_p, e["errorlambda_paint"])
        c_kp = c_lam_p * c_s_p
        ratio = c_Tamb_hot / c_Tamb_cold

        # Energy balance, written with in-place operations to limit the peak memory
        q = c_Thot**4
        q -= c_Tamb_hot**4
        q *= -sides * sigma * c_eps
        q += c_V * c_I / c_A
        if dTdt_m is not None:
            c_dTdt = rel(dTdt_m, e["errordTdt"])
            c_dTdt *= c_rho * c_s * c_cp + c_rho_p * c_s_p * c_cp_p
            q -= c_dTdt
            del c_dTdt
        if d2x_m is not None:
            c_d2 = rel(d2x_m, e["errord2Tdx2"])
            c_d2 *= c_kx + c_kp
            q -= c_d2
            c_d2 = rel(d2y_m, e["errord2Tdy2"])
            c_d2 *= c_ky + c_kp
            q -= c_d2
            del c_d2

        Taw = c_Tcold * ratio
        h = c_Thot
        h -= Taw[:, :, None] if c_Thot.ndim == 3 else Taw
        np.divide(q, h, out=h)   # h = q / (Thot - Taw)
        del q, c_Thot

        def _reduce(field, store_mean, store_fluc):
            if field.ndim == 3:
                m = field.mean(axis=2)
                store_fluc[i] = np.abs(field - m[:, :, None]).mean()
                store_mean[i] = m.mean()
            else:
                store_mean[i] = field.mean()

        if want_h:
            _reduce(h, out_h, outf_h)
        if want_nu:
            _reduce(h * (c_L / c_k), out_nu, outf_nu)
        if want_st:
            c_U = rel(cond["Uinf"], e["errorUinf"])
            c_rhoinf = rel(cond["rhoinf"], e["errorrhoinf"])
            c_cpinf = rel(cond["cpinf"], e["errorcpinf"])
            _reduce(h / (c_rhoinf * c_cpinf * c_U), out_st, outf_st)
    if verbose:
        print()

    result = {}

    def _finish(name, samples, fluct):
        result[f"error{name}"] = samples
        ref = None if reference is None else reference.get(name)
        if ref is None:
            warnings.warn(f"PIRT: the percentage error of {name} could not be computed since the "
                          f"reference {name} map was not provided")
        else:
            ref = np.asarray(ref, dtype=np.float64)
            ref_mean = float(np.nanmean(ref))
            result[f"error{name}_p"] = np.abs(samples - ref_mean) / abs(ref_mean) * 100.0
            if mode == "snapshots" and ref.ndim == 3:
                ref_f = float(np.nanmean(np.abs(ref - ref.mean(axis=2, keepdims=True))))
                result[f"error{name}f"] = fluct
                result[f"error{name}f_p"] = np.abs(fluct - ref_f) / abs(ref_f) * 100.0
        if mode == "snapshots" and f"error{name}f" not in result:
            result[f"error{name}f"] = fluct
        if verbose:
            print(f"The mean value of {name} is {np.nanmean(samples):.2f} with a standard "
                  f"deviation {np.nanstd(samples):.2f}")

    if want_h:
        _finish("h", out_h, outf_h)
    if want_nu:
        _finish("Nu", out_nu, outf_nu)
    if want_st:
        _finish("St", out_st, outf_st)
    return result
