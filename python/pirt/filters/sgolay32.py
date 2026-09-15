"""Savitzky-Golay filter for 3-D arrays with a second-order polynomial.

Python port of ``sgolay32_coef.m`` and ``sgolay32_filter.m`` (R. Castellanos,
S. Discetti, UC3M). A quadratic polynomial in the three coordinates,

    T = a0 + a1 i + a2 j + a3 k + a4 ij + a5 ik + a6 jk + a7 i^2 + a8 j^2 + a9 k^2,

is fitted by least squares in a moving window of size ``w = (w_rows, w_cols,
w_time)``. The coefficient kernels are applied by 3-D convolution, which yields
the smoothed field (``a0``) and, directly from the polynomial fit, the time
derivative (``a3/dt``) and the second spatial derivatives (``2*a8/dx^2`` along
the columns, ``2*a7/dy^2`` along the rows).

Reference: S. J. Orfanidis, Introduction to Signal Processing, Prentice-Hall,
1995, Chapter 8.
"""
from __future__ import annotations

import warnings

import numpy as np
from scipy import ndimage

from ._common import as_triplet, ensure_3d

__all__ = ["sgolay32_coef", "sgolay32_filter", "remove_edges"]


def sgolay32_coef(kernel_size=3, h=1.0) -> np.ndarray:
    """Least-squares coefficient matrix ``C`` of the 3-D quadratic SG filter.

    Parameters
    ----------
    kernel_size : int or (3,) sequence
        Odd window size in each dimension (rows, columns, time).
    h : float or (3,) sequence
        Grid spacing in each dimension (``dx``, ``dy``, ``dz``). The
        polynomial is expressed in scaled coordinates ``i*dx``, ``j*dy``,
        ``k*dz``.

    Returns
    -------
    C : ndarray, shape (10, prod(kernel_size))
        Row ``r`` contains the convolution coefficients of polynomial term
        ``r`` in the order ``[1, i, j, k, ij, ik, jk, ii, jj, kk]``. The
        columns are ordered with ``i`` (rows) varying fastest, then ``j``
        (columns), then ``k`` (time, in *reversed* order so that a convolution
        - which flips the kernel - returns the correct sign of the time
        derivative). Use ``C[r].reshape(kernel_size, order='F')`` to obtain the
        3-D convolution kernel, as ``reshape`` does in MATLAB.
    """
    w = as_triplet(kernel_size, "kernel_size", dtype=int)
    if np.any(w % 2 == 0) or np.any(w < 1):
        raise ValueError("sgolay32_coef: kernel sizes must be odd positive integers")
    hh = as_triplet(h, "h", dtype=float)
    dx, dy, dz = hh

    wc = w // 2
    range_x = np.arange(-wc[0], wc[0] + 1)
    range_y = np.arange(-wc[1], wc[1] + 1)
    range_z = -np.arange(-wc[2], wc[2] + 1)  # reversed, as in the MATLAB code

    rows = []
    for k in range_z:
        for j in range_y:
            for i in range_x:
                rows.append([1.0, i * dx, j * dy, k * dz,
                             i * j * dx * dy, i * k * dx * dz, j * k * dy * dz,
                             i * i * dx**2, j * j * dy**2, k * k * dz**2])
    M = np.asarray(rows)
    # C = (M'M) \ M'. The pseudo-inverse gives the same result when M has full
    # column rank and handles unit windows (a dimension with w=1 makes the
    # corresponding polynomial terms unidentifiable; their kernels become 0).
    if np.any(w == 1):
        warnings.warn("sgolay32_coef: a kernel size of 1 was given in some dimension; the derivative "
                      "along that dimension cannot be estimated and is set to zero")
    C = np.linalg.pinv(M, rcond=1e-12)
    return C


def remove_edges(X: np.ndarray, kernel_size, crop_time: bool = True) -> np.ndarray:
    """Replace the border regions affected by the convolution.

    The ``w//2`` outermost rows/columns are replaced by the first valid inner
    row/column (edge replication) and, if ``crop_time`` is True, the first and
    last ``w_t//2`` snapshots are removed. Port of the ``remove_edges``
    sub-function of ``sgolay32_filter.m``.
    """
    w = as_triplet(kernel_size, "kernel_size", dtype=int)
    win = w // 2
    Xf = X
    if win[0] > 0:
        Xf[: win[0], :, :] = Xf[win[0] : win[0] + 1, :, :]
        Xf[-win[0] :, :, :] = Xf[-win[0] - 1 : -win[0], :, :]
    if win[1] > 0:
        Xf[:, : win[1], :] = Xf[:, win[1] : win[1] + 1, :]
        Xf[:, -win[1] :, :] = Xf[:, -win[1] - 1 : -win[1], :]
    if win[2] > 0:
        Xf[:, :, : win[2]] = Xf[:, :, win[2] : win[2] + 1]
        Xf[:, :, -win[2] :] = Xf[:, :, -win[2] - 1 : -win[2]]
        if crop_time:
            Xf = Xf[:, :, win[2] : Xf.shape[2] - win[2]]
    return Xf


def _separable_decomposition(C: np.ndarray, w: np.ndarray):
    """Express the four SG kernels as sums of rank-1 (separable) kernels.

    For a rectangular window the least-squares kernels of the quadratic fit
    are linear combinations of ``1``, ``i^2``, ``j^2``, ``k^2`` (smoothing and
    second derivatives) and of ``k`` (time derivative). Returns the 1-D kernels
    (in array order, including the reversed time axis of the coefficient
    matrix) and the combination coefficients, or ``None`` if the decomposition
    is not exact (it always is for this filter; the check guards the fallback).
    """
    wc = w // 2
    u_i = np.arange(-wc[0], wc[0] + 1, dtype=float)
    u_j = np.arange(-wc[1], wc[1] + 1, dtype=float)
    u_k = -np.arange(-wc[2], wc[2] + 1, dtype=float)  # reversed, as in sgolay32_coef
    ones = [np.ones(len(u_i)), np.ones(len(u_j)), np.ones(len(u_k))]
    # basis tensors (rows, cols, time)
    basis = {
        "one": (ones[0], ones[1], ones[2]),
        "i2": (u_i**2, ones[1], ones[2]),
        "j2": (ones[0], u_j**2, ones[2]),
        "k2": (ones[0], ones[1], u_k**2),
        "k1": (ones[0], ones[1], u_k),
    }
    names = list(basis)
    B = np.stack([np.einsum("i,j,k->ijk", *basis[n]).ravel() for n in names], axis=1)
    coefs = {}
    for row in (0, 3, 7, 8):
        kern = C[row].reshape(tuple(w), order="F").ravel()
        scale = np.abs(kern).max()
        if scale < 1e-14 * np.abs(C).max():  # identically zero kernel (e.g. unit time window)
            coefs[row] = dict(zip(names, np.zeros(len(names))))
            continue
        c, *_ = np.linalg.lstsq(B, kern, rcond=None)
        if np.abs(B @ c - kern).max() > 1e-10 * scale:
            return None
        c[np.abs(c) < 1e-14 * np.abs(c).max()] = 0.0
        coefs[row] = dict(zip(names, c))
    return basis, coefs


def _conv1d(A, kern, axis):
    return ndimage.convolve1d(A, kern.astype(A.dtype, copy=False), axis=axis, mode="constant", cval=0.0)


def _sgolay_separable(X3, basis, coefs, w, dx, dy, dt):
    """Apply the SG kernels as sums of separable 1-D convolutions (exact)."""
    one_i, one_j, one_k = basis["one"]
    sq_i, sq_j, sq_k, lin_k = basis["i2"][0], basis["j2"][1], basis["k2"][2], basis["k1"][2]
    outputs = {0: None, 7: None, 8: None}  # a0 (smoothing), ii (rows), jj (cols)

    def accumulate(name, M):
        for row in outputs:
            c = coefs[row][name]
            if c != 0.0:
                outputs[row] = c * M if outputs[row] is None else _axpy(outputs[row], c, M)

    def _axpy(acc, c, M):
        acc += c * M
        return acc

    Y2 = _conv1d(X3, one_k, 2)
    Y21 = _conv1d(Y2, one_j, 1)
    accumulate("one", _conv1d(Y21, one_i, 0))
    accumulate("i2", _conv1d(Y21, sq_i, 0))
    del Y21
    accumulate("j2", _conv1d(_conv1d(Y2, sq_j, 1), one_i, 0))
    del Y2
    accumulate("k2", _conv1d(_conv1d(_conv1d(X3, sq_k, 2), one_j, 1), one_i, 0))
    for row in outputs:
        if outputs[row] is None:
            outputs[row] = np.zeros(X3.shape, dtype=X3.dtype)

    Tfilt = remove_edges(outputs[0], w, crop_time=True)
    d2Tdx2 = outputs[8]
    d2Tdx2 *= 2.0 / dx**2
    d2Tdx2 = remove_edges(d2Tdx2, w, crop_time=True)
    d2Tdy2 = outputs[7]
    d2Tdy2 *= 2.0 / dy**2
    d2Tdy2 = remove_edges(d2Tdy2, w, crop_time=True)

    ck = coefs[3]["k1"]
    dTdt = _conv1d(_conv1d(_conv1d(X3, lin_k, 2), one_j, 1), one_i, 0)
    dTdt *= ck / dt
    dTdt = remove_edges(dTdt, w, crop_time=True)
    return Tfilt, dTdt, d2Tdx2, d2Tdy2


def sgolay32_filter(X, kernel_size=3, h=1.0, return_derivatives: bool = True, method: str = "auto"):
    """Savitzky-Golay smoothing of a 3-D array with derivative estimation.

    Parameters
    ----------
    X : ndarray, shape (ny, nx, nt)
        Snapshots (rows, columns, time).
    kernel_size : int or (3,) sequence
        Odd window sizes ``(w_rows, w_cols, w_time)``. Each must not exceed
        the corresponding array dimension.
    h : float or (3,) sequence
        Physical spacing ``(dx, dy, dt)`` used to scale the derivatives:
        ``dx`` along the columns (x), ``dy`` along the rows (y), ``dt`` along
        time. As in the MATLAB implementation, the coefficient kernels are
        computed with unit spacing and the scaling is applied to the outputs.
    return_derivatives : bool
        If False only the smoothed field is returned.
    method : {'auto', 'separable', 'direct'}
        ``'direct'`` performs the four 3-D convolutions (``convn`` in MATLAB);
        ``'separable'`` exploits the exact rank-1 decomposition of the kernels
        (same result to round-off, roughly ten times faster). ``'auto'`` uses
        the separable algorithm whenever the decomposition is exact.

    Returns
    -------
    Tfilt : ndarray, shape (ny, nx, nt - 2*(w_time//2))
        Smoothed field. The first and last ``w_time//2`` snapshots are
        removed; the spatial borders are replaced by edge replication.
    dTdt, d2Tdx2, d2Tdy2 : ndarray
        Time derivative and second spatial derivatives (same shape as
        ``Tfilt``). Only if ``return_derivatives`` is True.
    """
    X = np.asarray(X)
    X3 = ensure_3d(X)
    w = as_triplet(kernel_size, "kernel_size", dtype=int)
    if np.any(w % 2 == 0):
        raise ValueError("sgolay32_filter: kernel sizes must be odd")
    if np.any(w > np.asarray(X3.shape)):
        raise ValueError(
            f"sgolay32_filter: kernel size {tuple(w)} exceeds the array shape {X3.shape}"
        )
    hh = as_triplet(h, "h", dtype=float)
    dx, dy, dt = hh
    if method not in ("auto", "separable", "direct"):
        raise ValueError("method must be 'auto', 'separable' or 'direct'")
    if not np.issubdtype(X3.dtype, np.floating):
        X3 = X3.astype(np.float64)

    C = sgolay32_coef(w)  # unit spacing, as in MATLAB

    decomposition = None if method == "direct" else _separable_decomposition(C, w)
    if method == "separable" and decomposition is None:
        raise RuntimeError("sgolay32_filter: the separable decomposition is not exact for this kernel")

    if decomposition is not None:
        basis, coefs = decomposition
        Tfilt, dTdt, d2Tdx2, d2Tdy2 = _sgolay_separable(X3, basis, coefs, w, dx, dy, dt)
        if X.ndim == 2:
            Tfilt, dTdt, d2Tdx2, d2Tdy2 = (a[:, :, 0] for a in (Tfilt, dTdt, d2Tdx2, d2Tdy2))
        return (Tfilt, dTdt, d2Tdx2, d2Tdy2) if return_derivatives else Tfilt

    def kernel(row):
        return C[row].reshape(tuple(w), order="F")

    def conv(kern):
        # convn(X, kern, 'same') zero-pads; the affected borders are replaced afterwards
        return ndimage.convolve(X3, kern.astype(X3.dtype, copy=False), mode="constant", cval=0.0)

    Tfilt = remove_edges(conv(kernel(0)), w, crop_time=True)
    if not return_derivatives:
        return Tfilt if X.ndim == 3 else Tfilt[:, :, 0]

    dTdt = conv(kernel(3))
    dTdt /= dt
    dTdt = remove_edges(dTdt, w, crop_time=True)

    d2Tdx2 = conv(kernel(8))  # jj term -> columns (x)
    d2Tdx2 *= 2.0 / dx**2
    d2Tdx2 = remove_edges(d2Tdx2, w, crop_time=True)

    d2Tdy2 = conv(kernel(7))  # ii term -> rows (y)
    d2Tdy2 *= 2.0 / dy**2
    d2Tdy2 = remove_edges(d2Tdy2, w, crop_time=True)

    if X.ndim == 2:
        Tfilt, dTdt, d2Tdx2, d2Tdy2 = (a[:, :, 0] for a in (Tfilt, dTdt, d2Tdx2, d2Tdy2))
    return Tfilt, dTdt, d2Tdx2, d2Tdy2
