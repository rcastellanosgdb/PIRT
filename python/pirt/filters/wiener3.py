"""3-D pixel-wise adaptive Wiener filter (port of ``wiener3.m``).

Reference: J. S. Lim, Two-Dimensional Signal and Image Processing, Prentice
Hall, 1990, p. 548 (3-D extension of MATLAB's ``wiener2``).
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

from ._common import as_triplet, ensure_3d

__all__ = ["wiener3"]


def _box_sum(A: np.ndarray, size: np.ndarray) -> np.ndarray:
    """``convn(A, ones(size), 'same')`` with zero padding.

    ``uniform_filter`` centres even-sized windows one sample to the left of
    MATLAB's ``convn(...,'same')`` convention, hence ``origin=-1`` on even
    axes.
    """
    origin = [-1 if n % 2 == 0 else 0 for n in size]
    out = ndimage.uniform_filter(A, size=tuple(int(n) for n in size), mode="constant",
                                 cval=0.0, origin=origin)
    out *= float(np.prod(size))
    return out


def wiener3(X, kernel=(3, 3, 3), noise: float | None = None, return_noise: bool = False):
    """Pixel-wise adaptive Wiener low-pass filtering of a 3-D array.

    Parameters
    ----------
    X : ndarray, shape (ny, nx, nt)
        Data degraded by constant-power additive noise.
    kernel : int or (3,) sequence
        Neighbourhood size ``(rows, columns, time)`` used to estimate the
        local mean and variance (default ``3``).
    noise : float, optional
        Additive noise power. If omitted it is estimated as the mean of the
        local variances.
    return_noise : bool
        Also return the (estimated or given) noise power.

    Returns
    -------
    f : ndarray
        Filtered data, same shape and dtype as ``X``.
    noise : float
        Only if ``return_noise`` is True.

    Notes
    -----
    Local statistics are computed with zero padding (``convn(...,'same')``),
    exactly as in the MATLAB implementation, so the estimates are biased near
    the borders.
    """
    X = np.asarray(X)
    X3 = ensure_3d(X)
    nhood = as_triplet(kernel, "kernel", dtype=int)
    if np.any(nhood < 1):
        raise ValueError("wiener3: kernel sizes must be positive integers")
    if not np.issubdtype(X3.dtype, np.floating):
        raise TypeError("wiener3: only floating-point input is supported")

    g = X3.astype(np.float64, copy=False)  # im2double
    n_el = float(np.prod(nhood))

    local_mean = _box_sum(g, nhood) / n_el
    local_var = _box_sum(g * g, nhood) / n_el
    local_var -= local_mean**2

    if noise is None:
        noise = float(local_var.mean())

    f = g - local_mean
    gain = np.maximum(local_var - noise, 0.0)
    denom = np.maximum(local_var, noise)
    f /= denom
    f *= gain
    f += local_mean

    f = f.astype(X3.dtype, copy=False)
    if X.ndim == 2:
        f = f[:, :, 0]
    return (f, noise) if return_noise else f
