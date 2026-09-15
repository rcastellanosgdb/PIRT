"""Small helpers shared by the filter implementations."""
from __future__ import annotations

import numpy as np


def as_triplet(value, name: str, dtype=float) -> np.ndarray:
    """Expand a scalar or a 3-element sequence into a length-3 array.

    Mirrors the MATLAB convention used throughout PIRT, where a scalar
    parameter is applied uniformly to the three array dimensions
    (rows, columns, time).
    """
    arr = np.atleast_1d(np.asarray(value, dtype=dtype)).ravel()
    if arr.size == 1:
        arr = np.repeat(arr, 3)
    elif arr.size != 3:
        raise ValueError(
            f"{name} must be a scalar or a 3-element sequence, got {arr.size} values"
        )
    return arr


def ensure_3d(X: np.ndarray, name: str = "X") -> np.ndarray:
    """Return ``X`` as a 3-D array ``(ny, nx, nt)`` without copying when possible."""
    X = np.asarray(X)
    if X.ndim == 2:
        return X[:, :, np.newaxis]
    if X.ndim != 3:
        raise ValueError(f"{name} must be a 2-D or 3-D array, got ndim={X.ndim}")
    return X


def matlab_round(x):
    """MATLAB ``round``: round half away from zero (NumPy rounds half to even)."""
    x = np.asarray(x, dtype=float)
    return np.sign(x) * np.floor(np.abs(x) + 0.5)
