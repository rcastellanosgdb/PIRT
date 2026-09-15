"""3-D Gaussian smoothing equivalent to MATLAB's ``imgaussfilt3``."""
from __future__ import annotations

import numpy as np
from scipy import ndimage

from ._common import as_triplet, ensure_3d

__all__ = ["gaussian_filter3", "gaussian_kernel_1d"]

_PADDING = {"replicate": "nearest", "circular": "wrap", "symmetric": "reflect"}


def gaussian_kernel_1d(sigma: float, size: int) -> np.ndarray:
    """Normalised 1-D Gaussian kernel of odd ``size`` (as ``imgaussfilt3`` builds it)."""
    size = int(size)
    if size < 1 or size % 2 == 0:
        raise ValueError("Gaussian filter size must be a positive odd integer")
    x = np.arange(size, dtype=float) - (size - 1) / 2.0
    k = np.exp(-(x**2) / (2.0 * float(sigma) ** 2))
    return k / k.sum()


def gaussian_filter3(X, sigma=0.5, filter_size=None, padding="replicate", filter_domain="auto"):
    """Separable 3-D Gaussian filtering (``imgaussfilt3`` semantics).

    Parameters
    ----------
    X : ndarray, shape (ny, nx, nt)
    sigma : float or (3,) sequence
        Standard deviation of the kernel along ``(rows, columns, time)``
        (default 0.5).
    filter_size : int or (3,) sequence, optional
        Odd kernel size per dimension. Default ``2*ceil(2*sigma)+1``.
    padding : {'replicate', 'circular', 'symmetric'} or float
        Boundary handling (MATLAB names). A number pads with that constant.
    filter_domain : {'auto', 'spatial', 'frequency'}
        Accepted for interface compatibility; the filter is always applied in
        the spatial domain (results are identical up to round-off).

    Returns
    -------
    Xf : ndarray
        Filtered array, same shape and dtype as ``X``.
    """
    X = np.asarray(X)
    X3 = ensure_3d(X)
    sig = as_triplet(sigma, "sigma", dtype=float)
    if np.any(sig <= 0):
        raise ValueError("sigma must be positive")
    if filter_size is None:
        fs = 2 * np.ceil(2 * sig).astype(int) + 1
    else:
        fs = as_triplet(filter_size, "filter_size", dtype=int)
    if np.any(fs % 2 == 0) or np.any(fs < 1):
        raise ValueError("filter_size must contain positive odd integers")
    if str(filter_domain).lower() not in ("auto", "spatial", "frequency"):
        raise ValueError("filter_domain must be 'auto', 'spatial' or 'frequency'")

    if isinstance(padding, str):
        if padding.lower() not in _PADDING:
            raise ValueError("padding must be 'replicate', 'circular', 'symmetric' or a number")
        mode, cval = _PADDING[padding.lower()], 0.0
    else:
        mode, cval = "constant", float(padding)

    work_dtype = np.result_type(X3.dtype, np.float32)
    out = X3.astype(work_dtype, copy=True)
    for axis in range(3):
        if fs[axis] > 1:
            k = gaussian_kernel_1d(sig[axis], fs[axis]).astype(work_dtype)
            out = ndimage.correlate1d(out, k, axis=axis, mode=mode, cval=cval)
    out = out.astype(X3.dtype, copy=False)
    if X.ndim == 2:
        out = out[:, :, 0]
    return out
