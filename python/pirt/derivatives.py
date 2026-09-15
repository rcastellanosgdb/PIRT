"""Finite-difference derivatives of temperature sequences (port of ``Derivative_FD.m``).

Conventions: arrays are ``(ny, nx, nt)``; ``x`` runs along the columns
(axis 1, spacing ``dx``), ``y`` along the rows (axis 0, spacing ``dy``) and
time along axis 2 (spacing ``dt``). This is the same convention as
:func:`pirt.filters.sgolay32_filter`.
"""
from __future__ import annotations

import numpy as np

__all__ = ["derivative_fd", "second_derivative_fd", "time_derivative_fd"]


def _check_step(value, name):
    if value is None:
        raise ValueError(f"The discretisation {name} is required")
    v = np.asarray(value, dtype=float)
    if v.size != 1 or not np.isfinite(v).all() or float(v) <= 0:
        raise ValueError(f"The discretisation {name} must be a single positive number")
    return float(v)


def second_derivative_fd(T, step: float, axis: int) -> np.ndarray:
    """Second derivative along ``axis`` with second-order central differences.

    At the two boundary samples the nearest interior central stencil is
    reused (i.e. the same value as the neighbouring point), as in the MATLAB
    implementation.
    """
    T = np.asarray(T)
    if T.shape[axis] < 4:
        raise ValueError("At least 4 samples are required along the differentiation axis")
    T = np.moveaxis(T, axis, 0)
    d2 = np.empty_like(T, dtype=np.result_type(T.dtype, np.float32))
    inv = 1.0 / step**2
    d2[1:-1] = (T[2:] - 2.0 * T[1:-1] + T[:-2]) * inv
    d2[0] = (T[2] - 2.0 * T[1] + T[0]) * inv
    d2[-1] = (T[-3] - 2.0 * T[-2] + T[-1]) * inv
    return np.moveaxis(d2, 0, axis)


def time_derivative_fd(T, dt: float, axis: int = 2) -> np.ndarray:
    """First derivative along ``axis`` (time): central differences in the
    interior and second-order one-sided formulas at the ends."""
    T = np.asarray(T)
    if T.ndim <= axis or T.shape[axis] < 3:
        raise ValueError("At least 3 snapshots are required to compute the time derivative")
    T = np.moveaxis(T, axis, 0)
    d = np.empty_like(T, dtype=np.result_type(T.dtype, np.float32))
    inv = 1.0 / (2.0 * dt)
    d[1:-1] = (T[2:] - T[:-2]) * inv
    d[0] = (-3.0 * T[0] + 4.0 * T[1] - T[2]) * inv
    d[-1] = (3.0 * T[-1] - 4.0 * T[-2] + T[-3]) * inv
    return np.moveaxis(d, 0, axis)


def derivative_fd(T, spatial: bool = False, temporal: bool = False,
                  dx: float | None = None, dy: float | None = None, dt: float | None = None):
    """Second spatial derivatives and time derivative by finite differences.

    Parameters
    ----------
    T : ndarray, shape (ny, nx) or (ny, nx, nt)
    spatial : bool
        Compute ``d2T/dx2`` (along columns) and ``d2T/dy2`` (along rows).
    temporal : bool
        Compute ``dT/dt`` along the third dimension.
    dx, dy, dt : float
        Grid spacings [m], [m], [s].

    Returns
    -------
    d2Tdx2, d2Tdy2, dTdt : ndarray or None
        ``None`` for the derivatives that were not requested.
    """
    T = np.asarray(T)
    d2Tdx2 = d2Tdy2 = dTdt = None
    if T.ndim not in (2, 3):
        raise ValueError("derivative_fd: T must be a 2-D or 3-D array")
    if spatial:
        dx = _check_step(dx, "dx")
        dy = _check_step(dy, "dy")
        if T.shape[0] < 4 or T.shape[1] < 4:
            raise ValueError("The image is too small for the finite-difference stencils")
        d2Tdx2 = second_derivative_fd(T, dx, axis=1)
        d2Tdy2 = second_derivative_fd(T, dy, axis=0)
    if temporal:
        if T.ndim != 3:
            raise ValueError("No temporal derivative can be computed for a single snapshot")
        dt = _check_step(dt, "dt")
        dTdt = time_derivative_fd(T, dt, axis=2)
    return d2Tdx2, d2Tdy2, dTdt
