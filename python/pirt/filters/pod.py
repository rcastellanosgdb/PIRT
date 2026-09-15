"""Proper Orthogonal Decomposition (POD) filtering and mode-truncation criteria.

Python port of ``POD_filter.m``, ``findNmod.m`` and ``optimal_SVHT_coef.m``
from the PIRT MATLAB toolbox.

References
----------
[1] L. Sirovich, "Turbulence and the dynamics of coherent structures. I.
    Coherent structures", Q. Appl. Math. 45 (1987) 561-571.
[2] M. Raiola, S. Discetti, A. Ianiro, "On PIV random error minimization with
    optimal POD-based low-order reconstruction", Exp. Fluids 56 (2015).
[3] M. Gavish, D. L. Donoho, "The optimal hard threshold for singular values is
    4/sqrt(3)", IEEE Trans. Inf. Theory 60 (2014) 5040-5053.
"""
from __future__ import annotations

import warnings

import numpy as np
from scipy import integrate

from ._common import ensure_3d

__all__ = ["pod_filter", "pod_decomposition", "find_nmod", "optimal_svht_coef", "matlab_smooth"]


# --------------------------------------------------------------------------- #
# Optimal singular value hard threshold (Gavish & Donoho, 2014)
# --------------------------------------------------------------------------- #
def _lambda_star(beta: np.ndarray) -> np.ndarray:
    """Optimal hard-threshold coefficient for known noise level."""
    w = (8.0 * beta) / (beta + 1.0 + np.sqrt(beta**2 + 14.0 * beta + 1.0))
    return np.sqrt(2.0 * (beta + 1.0) + w)


def _inc_marcenko_pastur(x0: float, beta: float, gamma: float = 0.0) -> float:
    """Incomplete Marcenko-Pastur integral from ``x0`` to the upper spectrum edge."""
    if beta > 1:
        raise ValueError("beta must be <= 1")
    top = (1.0 + np.sqrt(beta)) ** 2
    bot = (1.0 - np.sqrt(beta)) ** 2

    def density(x):
        q = (top - x) * (x - bot)
        if q <= 0.0:
            return 0.0
        val = np.sqrt(q) / (beta * x) / (2.0 * np.pi)
        return val * x**gamma if gamma != 0 else val

    # quadl in MATLAB uses an absolute tolerance of 1e-6; quad is more accurate.
    val, _ = integrate.quad(density, x0, top, limit=200, epsabs=1e-10, epsrel=1e-10)
    return val


def _median_marcenko_pastur(beta: float) -> float:
    """Median of the Marcenko-Pastur distribution (bisection as in the MATLAB code)."""
    lobnd = (1.0 - np.sqrt(beta)) ** 2
    hibnd = (1.0 + np.sqrt(beta)) ** 2

    def cdf(x):
        return 1.0 - _inc_marcenko_pastur(x, beta, 0.0)

    change = True
    while change and (hibnd - lobnd > 0.001):
        change = False
        x = np.linspace(lobnd, hibnd, 5)
        y = np.array([cdf(xi) for xi in x])
        if np.any(y < 0.5):
            lobnd = np.max(x[y < 0.5])
            change = True
        if np.any(y > 0.5):
            hibnd = np.min(x[y > 0.5])
            change = True
    return 0.5 * (hibnd + lobnd)


def optimal_svht_coef(beta, sigma_known: bool = False):
    """Optimal location of the hard threshold for singular values.

    Parameters
    ----------
    beta : float or array_like
        Aspect ratio ``m/n`` of the data matrix (``0 < beta <= 1``). For a
        snapshot matrix with ``n_pixels`` rows and ``n_snapshots`` columns,
        ``beta = n_snapshots / n_pixels``.
    sigma_known : bool
        ``True`` if the noise level is known (threshold relative to
        ``sigma*sqrt(n)``), ``False`` if unknown (threshold relative to the
        median singular value).

    Returns
    -------
    coef : ndarray
        Coefficient(s) with the same shape as ``beta``.
    """
    beta_arr = np.atleast_1d(np.asarray(beta, dtype=float))
    if np.any(beta_arr <= 0) or np.any(beta_arr > 1):
        raise ValueError("beta must satisfy 0 < beta <= 1")
    coef = _lambda_star(beta_arr)
    if not sigma_known:
        mp_median = np.array([_median_marcenko_pastur(float(b)) for b in beta_arr])
        coef = coef / np.sqrt(mp_median)
    return coef if np.ndim(beta) else float(coef[0])


# --------------------------------------------------------------------------- #
# Truncation criteria
# --------------------------------------------------------------------------- #
def matlab_smooth(y, span: int = 5) -> np.ndarray:
    """Moving average identical to MATLAB ``smooth(y)`` (``'moving'`` method).

    The window shrinks symmetrically near the ends so that the first and last
    elements are unchanged, the second/penultimate ones average three points,
    and so on.
    """
    y = np.asarray(y, dtype=float).ravel()
    n = y.size
    span = int(span)
    if span % 2 == 0:  # MATLAB reduces even spans by one
        span -= 1
    span = max(1, min(span, n))
    half = span // 2
    out = np.empty(n)
    for i in range(n):
        w = min(i, n - 1 - i, half)
        out[i] = y[i - w : i + w + 1].mean()
    return out


def find_nmod(s, criterion: str | None = None, threshold: float | None = None,
              beta: float | None = None, verbose: bool = True) -> int:
    """Number of POD modes to retain, given the singular values.

    Port of ``findNmod.m``.

    Parameters
    ----------
    s : array_like
        Singular values (vector, or a diagonal matrix as returned by ``svd``).
    criterion : {'Elbow', 'Spectrum', 'HardThreshold'}, optional
        Truncation criterion (case-insensitive). Default ``'Elbow'``.

        * ``'Spectrum'``: cumulative energy fraction (default threshold 0.99).
        * ``'Elbow'``: ratio of consecutive eigenvalues, smoothed with a
          5-point moving average, see [2] (default threshold 0.999).
        * ``'HardThreshold'``: optimal singular value hard threshold with
          unknown noise level, see [3]. Requires ``beta``.
    threshold : float, optional
        Threshold for the ``'Spectrum'`` and ``'Elbow'`` criteria.
    beta : float, optional
        Aspect ratio of the data matrix for ``'HardThreshold'``.

    Returns
    -------
    nmod : int
        Suggested number of modes (at least 1, at most ``len(s)``).

    Notes
    -----
    The implementation reproduces the MATLAB code exactly, including two
    quirks that are kept for backward compatibility:

    * ``'Elbow'`` discards the first ``skip=5`` ratios before searching and
      returns the index relative to the truncated vector (i.e. the result is
      not shifted back by ``skip``).
    * ``'HardThreshold'`` returns the index of the *first* singular value that
      falls below the threshold (that mode is therefore included in the
      reconstruction).
    """
    s = np.asarray(s, dtype=float)
    if s.ndim == 2:
        if s.shape[0] == s.shape[1] and np.allclose(s, np.diag(np.diag(s))):
            s = np.diag(s)
        else:
            raise ValueError("s must be a vector of singular values or a diagonal matrix")
    s = s.ravel()
    N = s.size
    lam = s**2

    crit = (criterion or "Elbow").strip().lower()
    if crit == "spectrum":
        if threshold is None:
            threshold = 0.99
        f = np.cumsum(lam) / np.sum(lam)
    elif crit == "elbow":
        if threshold is None:
            threshold = 0.999
        skip = 5
        ratio = lam[1:] / lam[:-1]
        f = matlab_smooth(ratio[skip - 1 : ratio.size - skip])  # MATLAB f(skip:end-skip)
    elif crit in ("hardthreshold", "hard_threshold", "svht"):
        if beta is None:
            raise ValueError("The data aspect ratio 'beta' is required for the HardThreshold criterion")
        coef = optimal_svht_coef(beta, sigma_known=False)
        threshold = -coef * np.median(s)
        f = -s
    else:
        raise ValueError(f"Unknown truncation criterion '{criterion}'")

    if crit in ("spectrum", "elbow") and threshold > 1:
        warnings.warn("findNmod: the threshold cannot exceed 1; using 0.999 instead")
        threshold = 0.999

    idx = np.flatnonzero(f > threshold)
    nmod = int(idx[0]) + 1 if idx.size else N
    nmod = max(1, min(nmod, N))
    if verbose:
        print(f"The suggested number of modes is: {nmod}")
    return nmod


# --------------------------------------------------------------------------- #
# POD filter
# --------------------------------------------------------------------------- #
def pod_decomposition(X):
    """Snapshot POD of a 3-D array ``(ny, nx, nt)``.

    Returns
    -------
    mean : ndarray, shape (ny*nx,)
        Temporal mean of every pixel.
    s : ndarray, shape (nt,) or (ny*nx,)
        Singular values in descending order.
    V : ndarray
        Temporal modes (columns), such that the mean-removed snapshot matrix
        ``D = U @ diag(s) @ V.T``.
    D : ndarray
        The mean-removed snapshot matrix (pixels x snapshots), float64.
    """
    X = ensure_3d(X)
    nt = X.shape[2]
    D = X.reshape(-1, nt).astype(np.float64, copy=True)
    mean = D.mean(axis=1)
    D -= mean[:, None]

    if nt >= D.shape[0]:
        # More snapshots than pixels: eigen-decompose the spatial correlation.
        lam, U = np.linalg.eigh(D @ D.T)
        order = np.argsort(lam)[::-1]
        lam = np.clip(lam[order], 0.0, None)
        U = U[:, order]
        s = np.sqrt(lam)
        with np.errstate(divide="ignore", invalid="ignore"):
            V = (D.T @ U) / s
        V[:, s <= 0] = 0.0
    else:
        # Usual IR case: many pixels, few snapshots -> temporal correlation.
        lam, V = np.linalg.eigh(D.T @ D)
        order = np.argsort(lam)[::-1]
        lam = np.clip(lam[order], 0.0, None)
        V = V[:, order]
        s = np.sqrt(lam)
    return mean, s, V, D


def pod_filter(X, nmod: int | None = None, criterion: str | None = None,
               threshold: float | None = None, beta: float | None = None,
               return_nmod: bool = False, verbose: bool = True):
    """Modal (POD) low-order reconstruction of a sequence of temperature maps.

    Port of ``POD_filter.m``. The snapshot matrix is mean-removed, decomposed
    with the method of snapshots and reconstructed with the first ``nmod``
    modes.

    Parameters
    ----------
    X : ndarray, shape (ny, nx, nt)
        Temperature snapshots; the third dimension is time.
    nmod : int, optional
        Number of modes to retain. If omitted it is selected with
        :func:`find_nmod` using ``criterion``/``threshold``/``beta``.
    criterion, threshold, beta :
        Passed to :func:`find_nmod`. If ``criterion='HardThreshold'`` and
        ``beta`` is omitted, ``beta = nt / (ny*nx)`` is used (with a warning),
        as in the MATLAB implementation.
    return_nmod : bool
        Also return the number of retained modes.

    Returns
    -------
    Xf : ndarray
        Filtered snapshots, same shape and dtype as ``X``.
    nmod : int
        Only if ``return_nmod`` is True.
    """
    X = np.asarray(X)
    X3 = ensure_3d(X)
    if X3.shape[2] < 2:
        raise ValueError("The POD filter requires at least two snapshots (third dimension)")
    ny, nx, nt = X3.shape

    mean, s, V, D = pod_decomposition(X3)

    if nmod is None:
        crit = (criterion or "Elbow")
        if crit.strip().lower() in ("hardthreshold", "hard_threshold", "svht") and beta is None:
            beta = nt / (ny * nx)
            warnings.warn(
                "POD_filter: the data aspect ratio 'beta' was not introduced; "
                f"using the snapshot-matrix aspect ratio nt/(ny*nx) = {beta:.4g}"
            )
        if criterion is None:
            warnings.warn("POD_filter: the Elbow criterion was selected as default")
        nmod = find_nmod(s, criterion=crit, threshold=threshold, beta=beta, verbose=verbose)
    nmod = int(nmod)
    if nmod < 1 or nmod > s.size:
        raise ValueError(f"nmod must be between 1 and {s.size}, got {nmod}")

    # U_k S_k V_k' = D V_k V_k'  (avoids dividing by tiny singular values)
    Vk = V[:, :nmod]
    Xf = (D @ Vk) @ Vk.T
    del D
    Xf += mean[:, None]
    Xf = Xf.reshape(ny, nx, nt).astype(X.dtype, copy=False)
    if X.ndim == 2:
        Xf = Xf[:, :, 0]
    return (Xf, nmod) if return_nmod else Xf
