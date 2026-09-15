"""Heated-thin-foil energy balance (port of ``Calculate_HeatTransfer.m``).

The convective heat-transfer coefficient is obtained from the energy balance
of the heated thin foil sensor (Astarita & Carlomagno, 2012):

    h = (q_j - q_r - q_k - q_u + sum(q_custom)) / (T_w - T_aw)

with

* ``q_j = V I / A``                                   Joule heating,
* ``q_r = sides * sigma * eps * (T_w^4 - T_inf^4)``   radiation,
* ``q_u = (rho s cp + rho_p s_p cp_p) dT_w/dt``       unsteady (foil + paint),
* ``q_k = (k + k_p) (d2T/dx2 + d2T/dy2)``             tangential conduction,
  or ``(s lambda_x + k_p) d2T/dx2 + (s lambda_y + k_p) d2T/dy2`` for PCBs,

where ``k = s * lambda_foil`` and ``k_p = s_p * lambda_p`` are the thermal
conductances (conductivity x thickness) [W/K] of foil and paint, and
``T_aw = T_cold * T_amb,hot / T_amb,cold``.

Nu = h L / k_air(T_film) with ``T_film = (T_w + T_inf)/2`` and the air thermal
conductivity given by a cubic polynomial fit [W/(m K)]. St = h/(rho_inf cp_inf U_inf).
"""
from __future__ import annotations

import warnings

import numpy as np

from .derivatives import derivative_fd

__all__ = ["calculate_heat_transfer", "air_thermal_conductivity", "check_celsius",
           "STEFAN_BOLTZMANN", "GRAVITY", "PAINT_DEFAULTS", "prepare_hfs", "prepare_conditions"]

STEFAN_BOLTZMANN = 5.67e-8  # [W/(m^2 K^4)]
GRAVITY = 9.80665  # [m/s^2]
CELSIUS_THRESHOLD = 100.0  # [K]; mean temperatures below this are assumed in Celsius
_AIR_K_POLY = (1.5207e-11, -4.8574e-08, 1.0184e-04, -3.9333e-04)

# Default paint properties used when not provided (with a warning), as in MATLAB.
PAINT_DEFAULTS = {
    "s_paint": 21.81e-6,  # [m]   Stafford, Walsh & Egan, Meas. Sci. Technol. 20 (2009)
    "rho_paint": 1300.0,  # [kg/m^3]
    "cp_paint": 5000.0,  # [J/(kg K)]
    "lambda_paint": 1.38,  # [W/(m K)] Raghu & Philip, Meas. Sci. Technol. 17 (2006)
}


def air_thermal_conductivity(T):
    """Thermal conductivity of air [W/(m K)] as a cubic polynomial of T [K]."""
    T = np.asarray(T)
    a, b, c, d = _AIR_K_POLY
    return ((a * T + b) * T + c) * T + d


def check_celsius(T, threshold: float = CELSIUS_THRESHOLD, name: str = "Temperature"):
    """Convert to Kelvin if the mean value is below ``threshold`` (assumed Celsius)."""
    T = np.asarray(T, dtype=float) if np.ndim(T) == 0 else np.asarray(T)
    if np.nanmean(T) < threshold:
        warnings.warn(f"PIRT: {name} was converted to Kelvin. Check your inputs!")
        return T + 273.15
    return T


# --------------------------------------------------------------------------- #
# Input preparation
# --------------------------------------------------------------------------- #
def prepare_hfs(hfs: dict, spatial_der: bool = False) -> dict:
    """Validate and complete the heated-foil-sensor (HFS) dictionary.

    Required keys: ``s`` [m], ``rho`` [kg/m^3], ``cp`` [J/(kg K)], ``epsilon``
    [-], ``Area`` (or ``A``) [m^2] or both ``H`` and ``W`` [m], and ``k``
    [W/K] (conductance ``s*lambda``) for ``Type='Foil'`` or ``lambdax``,
    ``lambday`` [W/(m K)] for ``Type='PCB'``. Optional: ``sides`` (default 1),
    ``s_paint``, ``rho_paint``, ``cp_paint``, ``lambda_paint``.
    """
    if not hfs:
        raise ValueError("PIRT: the HFS data must be introduced to calculate the heat transfer")
    out = dict(hfs)
    htype = str(out.get("Type", "Foil"))
    if htype not in ("Foil", "PCB"):
        raise ValueError("PIRT: HFS 'Type' has to be either 'Foil' or 'PCB'")
    out["Type"] = htype
    for key in ("s", "rho", "cp", "epsilon"):
        if key not in out:
            raise ValueError(f"PIRT: the HFS parameter '{key}' must be introduced to compute the heat transfer")
        out[key] = float(out[key])
    if "A" not in out:
        if "Area" in out:
            out["A"] = float(out["Area"])
        elif "H" in out and "W" in out:
            out["A"] = float(out["H"]) * float(out["W"])
        else:
            raise ValueError("PIRT: either the HFS 'Area' or both 'H' and 'W' must be introduced")
    out["A"] = float(out["A"])
    if htype == "PCB":
        if spatial_der and ("lambdax" not in out or "lambday" not in out):
            raise ValueError("PIRT: 'lambdax' and 'lambday' are required for the tangential term of a PCB")
    else:
        if "k" not in out:
            raise ValueError("PIRT: the HFS conductance 'k' [W/K] (thickness x conductivity) must be introduced")
        out["k"] = float(out["k"])
    if "sides" not in out:
        warnings.warn("PIRT: number of exposed sides was not selected, it will be set to sides=1")
        out["sides"] = 1.0
    for key, default in PAINT_DEFAULTS.items():
        if key not in out:
            warnings.warn(f"PIRT: the paint property '{key}' was not introduced, using {default:g}")
            out[key] = default
        out[key] = float(out[key])
    out["m_plate"] = out["rho"] * out["s"] * out["A"]
    return out


def prepare_conditions(conditions: dict, compute_st: bool = False) -> dict:
    """Validate the test-conditions dictionary.

    Required keys: ``L`` (or ``L_char``) [m], ``V`` [V], ``I`` [A] and
    ``Tamb`` = (T_amb,cold, T_amb,hot). Optional: ``dx``, ``dy`` [m], ``dt``
    [s] (for finite-difference derivatives), ``Uinf`` [m/s], ``rhoinf``
    [kg/m^3] and ``cpinf`` [J/(kg K)] (for the Stanton number).
    """
    if not conditions:
        raise ValueError("PIRT: the Conditions must be introduced to calculate the heat transfer")
    out = dict(conditions)
    if "L_char" not in out:
        if "L" not in out:
            raise ValueError("PIRT: the characteristic length 'L' must be introduced")
        out["L_char"] = out["L"]
    for key in ("V", "I", "L_char"):
        if key not in out:
            raise ValueError(f"PIRT: the condition '{key}' must be introduced to compute the heat transfer")
        out[key] = float(out[key])
    if "Tamb" not in out:
        raise ValueError("PIRT: the ambient temperature 'Tamb' = (T_cold, T_hot) must be introduced")
    tamb = np.atleast_1d(np.asarray(out["Tamb"], dtype=float)).ravel()
    if tamb.size != 2:
        raise ValueError("PIRT: 'Tamb' must have 2 values: the cold and hot ambient temperatures")
    out["Tamb"] = tamb
    if compute_st:
        for key in ("Uinf", "rhoinf", "cpinf"):
            if key not in out:
                raise ValueError(
                    f"PIRT: '{key}' must be introduced to compute the Stanton number "
                    "(St = h / (rhoinf * cpinf * Uinf))"
                )
    return out


# --------------------------------------------------------------------------- #
# Energy balance
# --------------------------------------------------------------------------- #
def _check_custom_q(Thot, q):
    q = np.asarray(q)
    if q.ndim == 0 or q.size == 1:
        return float(q)
    if q.ndim == 2:
        if q.shape != Thot.shape[:2]:
            raise ValueError("PIRT: the sizes of the images and the custom heat matrices must agree")
        return q[:, :, None] if Thot.ndim == 3 else q
    if q.ndim == 3:
        if Thot.ndim != 3 or q.shape != Thot.shape:
            raise ValueError("PIRT: the sizes of the images and the custom heat matrices must agree")
        return q
    raise ValueError("PIRT: error in the custom heat flux dimensions")


def calculate_heat_transfer(Thot, Tcold, hfs: dict, conditions: dict, *,
                            compute=("h",), time_der: bool = False, spatial_der: bool = False,
                            dTdt=None, d2Tdx2=None, d2Tdy2=None, custom_q=None,
                            film_temperature: str = "ambient", return_terms: bool = False,
                            verbose: bool = True) -> dict:
    """Compute h, Nu and/or St maps from the heated-thin-foil energy balance.

    Parameters
    ----------
    Thot : ndarray, shape (ny, nx) or (ny, nx, nt)
        Wall temperature with the foil heated [K] (Celsius is detected and converted).
    Tcold : ndarray, shape (ny, nx) or (ny, nx, nt)
        Adiabatic-wall (unheated) temperature. A 3-D array is time-averaged.
    hfs, conditions : dict
        Sensor and test data, see :func:`prepare_hfs` and :func:`prepare_conditions`.
    compute : sequence of {'h', 'Nu', 'St'}
        Quantities to compute.
    time_der, spatial_der : bool
        Include the unsteady and tangential-conduction terms. The derivatives
        are taken from ``dTdt``/``d2Tdx2``/``d2Tdy2`` when provided (e.g. from
        the Savitzky-Golay filter) and otherwise computed with finite
        differences using ``conditions['dt']``, ``['dx']``, ``['dy']``.
    custom_q : sequence, optional
        Additional heat fluxes [W/m^2] added to the balance (scalars, 2-D or
        3-D arrays matching ``Thot``).
    film_temperature : {'ambient', 'adiabatic'}
        Temperature at which the air conductivity of the Nusselt number is
        evaluated: ``'ambient'`` (default, current MATLAB toolbox) uses
        ``T_film = (T_w + T_amb,hot)/2``; ``'adiabatic'`` uses
        ``T_film = (T_w + T_aw)/2`` with ``T_aw = T_cold * T_amb,hot/T_amb,cold``,
        the definition written in Robledo et al. (2025) and used by the MATLAB
        runs of January 2025 that produced most figures of that paper. The two
        choices differ by ~0.05 % in Nu for the sweeping-jet dataset.
    return_terms : bool
        Also return the individual heat-flux terms (``q_joule``, ``q_rad``,
        ``q_unsteady``, ``q_tangential``) and the derivatives used.

    Returns
    -------
    result : dict
        With keys among ``'h'``, ``'Nu'``, ``'St'`` and, if the derivatives
        were computed here, ``'dTdt_hot'``, ``'d2Tdx2_hot'``, ``'d2Tdy2_hot'``.
    """
    compute = {c.lower() for c in ([compute] if isinstance(compute, str) else compute)}
    unknown = compute - {"h", "nu", "st"}
    if unknown:
        raise ValueError(f"PIRT: unknown heat transfer quantities {sorted(unknown)}; use 'h', 'Nu', 'St'")
    if not compute:
        raise ValueError("PIRT: either 'h', 'Nu' or 'St' must be requested")
    if film_temperature not in ("ambient", "adiabatic"):
        raise ValueError("film_temperature must be 'ambient' or 'adiabatic'")

    Thot = np.asarray(Thot)
    Tcold = np.asarray(Tcold)
    if Thot.ndim not in (2, 3) or Tcold.ndim not in (2, 3):
        raise ValueError("PIRT: Thot and Tcold must be 2-D or 3-D arrays")
    if Tcold.ndim == 3:
        Tcold = Tcold.mean(axis=2)
    if Thot.shape[:2] != Tcold.shape[:2]:
        raise ValueError("PIRT: incompatible image sizes between Thot and Tcold")

    hfs = prepare_hfs(hfs, spatial_der=spatial_der)
    cond = prepare_conditions(conditions, compute_st="st" in compute)

    if verbose:
        print("-- Calculating heat transfer maps")
    Thot = check_celsius(Thot, name="Thot")
    Tcold = check_celsius(Tcold, name="Tcold")
    Tamb_cold = float(check_celsius(cond["Tamb"][0], name="Tamb(cold)"))
    Tamb_hot = float(check_celsius(cond["Tamb"][1], name="Tamb(hot)"))
    ratio = Tamb_hot / Tamb_cold

    s, rho, cp = hfs["s"], hfs["rho"], hfs["cp"]
    s_p, rho_p, cp_p = hfs["s_paint"], hfs["rho_paint"], hfs["cp_paint"]
    kp = s_p * hfs["lambda_paint"]
    A, epsilon, sides = hfs["A"], hfs["epsilon"], float(hfs["sides"])
    L_char, sigma = cond["L_char"], STEFAN_BOLTZMANN

    result: dict = {}
    terms: dict = {}
    work_dtype = np.result_type(Thot.dtype, np.float32)

    # Joule heating [W/m^2]
    q_joule = cond["V"] * cond["I"] / A
    # Radiation [W/m^2]
    q_rad = (sides * sigma * epsilon) * (Thot.astype(work_dtype) ** 4 - Tamb_hot**4)
    q = q_joule - q_rad
    if return_terms:
        terms["q_joule"], terms["q_rad"] = q_joule, q_rad
    del q_rad

    # Unsteady term
    if time_der:
        if dTdt is None:
            if "dt" in cond and Thot.ndim == 3:
                if verbose:
                    print("----Computing unsteady term with finite differences")
                _, _, dTdt = derivative_fd(Thot, temporal=True, dt=cond["dt"])
                result["dTdt_hot"] = dTdt
            else:
                warnings.warn("PIRT: unsteady term could not be computed, the temporal resolution 'dt' is missing "
                              "or Thot has a single snapshot.")
        if dTdt is not None:
            if verbose:
                print("----Adding unsteady term")
            if np.shape(dTdt) != Thot.shape:
                raise ValueError("PIRT: dTdt must have the same shape as Thot")
            q_uns = (rho * s * cp + rho_p * cp_p * s_p) * np.asarray(dTdt)
            q = q - q_uns
            if return_terms:
                terms["q_unsteady"] = q_uns
            del q_uns

    # Tangential conduction term
    if spatial_der:
        if d2Tdx2 is None or d2Tdy2 is None:
            if "dx" in cond and "dy" in cond:
                if verbose:
                    print("----Computing tangential term with finite differences")
                d2Tdx2, d2Tdy2, _ = derivative_fd(Thot, spatial=True, dx=cond["dx"], dy=cond["dy"])
                result["d2Tdx2_hot"], result["d2Tdy2_hot"] = d2Tdx2, d2Tdy2
            else:
                warnings.warn("PIRT: tangential-conduction term could not be computed, the spatial resolution "
                              "'dx', 'dy' is missing.")
        if d2Tdx2 is not None and d2Tdy2 is not None:
            if verbose:
                print("----Adding tangential term")
            if np.shape(d2Tdx2) != Thot.shape or np.shape(d2Tdy2) != Thot.shape:
                raise ValueError("PIRT: d2Tdx2 and d2Tdy2 must have the same shape as Thot")
            if hfs["Type"] == "PCB":
                q_tan = (s * float(hfs["lambdax"]) + kp) * np.asarray(d2Tdx2) \
                    + (s * float(hfs["lambday"]) + kp) * np.asarray(d2Tdy2)
            else:
                q_tan = (hfs["k"] + kp) * (np.asarray(d2Tdx2) + np.asarray(d2Tdy2))
            q = q - q_tan
            if return_terms:
                terms["q_tangential"] = q_tan
            del q_tan

    if custom_q is not None:
        if not isinstance(custom_q, (list, tuple)):
            raise ValueError("PIRT: custom heat terms must be introduced in a list")
        for qc in custom_q:
            q = q + _check_custom_q(Thot, qc)

    # Convective coefficient [W/(m^2 K)]
    Taw = Tcold * ratio
    if Thot.ndim == 3:
        Taw = Taw[:, :, None]
    h = q / (Thot - Taw)
    del q
    if "h" in compute:
        result["h"] = h

    if "nu" in compute:
        if film_temperature == "ambient":
            Tfilm = (Thot + Tamb_hot) / 2.0
        else:
            Tfilm = (Thot + Taw) / 2.0
        kair = air_thermal_conductivity(Tfilm)
        del Tfilm
        result["Nu"] = h * (L_char / kair)
        del kair

    if "st" in compute:
        result["St"] = h / (cond["rhoinf"] * cond["cpinf"] * cond["Uinf"])

    if return_terms:
        result["terms"] = terms
    return result
