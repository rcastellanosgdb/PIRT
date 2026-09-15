"""Multi-scale Proper Orthogonal Decomposition (mPOD) filter.

Python port of ``multiscale_POD_filter.m`` (E. Diaz, R. Castellanos et al.,
UC3M), based on Mendez, Balabane & Buchlin, "Multi-scale proper orthogonal
decomposition of complex fluid flows", J. Fluid Mech. 870 (2019) 988-1036.

The temporal correlation matrix ``K = D' D`` is transformed to the frequency
domain with the DFT matrix, multiplied by a transfer function ``H`` (peak
removal, frequency decoupling or both) and transformed back. The filtered
correlation matrix is then decomposed and the data reconstructed with the
first ``nmod`` modes.
"""
from __future__ import annotations

import warnings

import numpy as np

from ._common import ensure_3d, matlab_round
from .pod import find_nmod

__all__ = ["mpod_filter", "mpod_transfer_function"]

_MODE_ALIASES = {
    "peakrem": 1, "peak removal": 1, "peakremoval": 1, "peak_removal": 1,
    "freq deco": 2, "freqdec": 2, "frec decoupling": 2, "freq decoupling": 2,
    "frequency_decoupling": 2, "frequencydecoupling": 2, "frequency decoupling": 2,
    "both": 3,
}


def _mode_code(mode) -> int:
    if isinstance(mode, (int, np.integer)):
        if int(mode) in (1, 2, 3):
            return int(mode)
        raise ValueError("mode must be 1 (peak removal), 2 (frequency decoupling) or 3 (both)")
    key = str(mode).strip().lower()
    if key in _MODE_ALIASES:
        return _MODE_ALIASES[key]
    raise ValueError(
        f"Unknown mPOD mode '{mode}'. Use 'peak_removal', 'frequency_decoupling' or 'both'."
    )


def _dft_matrix(nt: int) -> np.ndarray:
    """Unitary, symmetric DFT matrix ``F_jk = exp(-2*pi*i*j*k/nt)/sqrt(nt)``.

    Identical to ``conj(w_F.^((0:nt-1)'*(0:nt-1)))/sqrt(nt)`` with
    ``w_F = exp(2*pi*i/nt)`` in the MATLAB code (the exponent is reduced
    modulo ``nt`` for accuracy).
    """
    k = np.arange(nt)
    expo = np.outer(k, k) % nt
    return np.exp(-2j * np.pi * expo / nt) / np.sqrt(nt)


def mpod_transfer_function(nt: int, f_acquisition: float, mode, fpeaks=None, w=None,
                           n_regions=None, f_min=None, f_max=None) -> np.ndarray:
    """Build the ``nt x nt`` transfer function ``H`` used to filter ``K_hat``.

    The index arithmetic reproduces the MATLAB implementation exactly
    (1-based indices, ``round`` half away from zero).
    """
    code = _mode_code(mode)
    df = f_acquisition / nt

    def peak_removal(H):
        if fpeaks is None or w is None:
            raise ValueError("Peak removal requires 'fpeaks' and 'w'")
        for f in np.atleast_1d(np.asarray(fpeaks, dtype=float)):
            idx1 = int(matlab_round((f - w) / df))
            idx2 = int(matlab_round((f + w) / df))
            if idx1 < 1 or idx2 > nt - 1 or idx2 < idx1:
                raise ValueError(
                    f"Peak removal window [{f - w:g}, {f + w:g}] Hz is out of the valid "
                    f"frequency range (df={df:g} Hz, nt={nt})"
                )
            rows = np.r_[idx1 - 1 : idx2, nt - idx2 - 1 : nt - idx1]  # MATLAB [idx1:idx2 nt-idx2:nt-idx1]
            H[rows, :] = 0.0
            H[:, rows] = 0.0
        return H

    def frequency_decoupling(H):
        if n_regions is None or f_min is None or f_max is None:
            raise ValueError("Frequency decoupling requires 'n_regions', 'f_min' and 'f_max'")
        nf = int(n_regions)
        idx0 = int(matlab_round(f_min / df))
        delta = int(matlab_round(((f_max - f_min) / nf) / df))
        for i in range(1, nf + 1):
            fr0 = idx0 + delta * (i - 1)
            fr1 = idx0 + delta * i
            # MATLAB: [frange(1)+1:frange(2)  nt-frange(2)+1:nt-frange(1)]
            rows = np.r_[fr0:fr1, nt - fr1 : nt - fr0]
            rows = rows[(rows >= 0) & (rows < nt)]
            H[np.ix_(rows, rows)] = 1.0
        return H

    if code == 1:
        H = peak_removal(np.ones((nt, nt)))
    elif code == 2:
        H = frequency_decoupling(np.zeros((nt, nt)))
    else:
        H = peak_removal(frequency_decoupling(np.zeros((nt, nt))))
    return H


def mpod_filter(X, f_acquisition: float, mode="peak_removal", fpeaks=None, w=None,
                n_regions=None, f_min=None, f_max=None, threshold: float | None = None,
                nmod: int | None = None, legacy_transform: bool = False,
                return_nmod: bool = True, verbose: bool = True):
    """mPOD filtering of a sequence of temperature maps.

    Parameters
    ----------
    X : ndarray, shape (ny, nx, nt)
        Snapshots; the third dimension is time.
    f_acquisition : float
        Sampling frequency [Hz].
    mode : {'peak_removal', 'frequency_decoupling', 'both'}
        How the transfer function ``H`` is built (MATLAB spellings such as
        ``'Peak Removal'`` or ``'Freq Deco'`` are accepted).
    fpeaks, w : array_like, float
        Frequencies [Hz] to remove and half-width [Hz] of the window around
        each of them (peak removal).
    n_regions, f_min, f_max :
        Number of frequency bands and range [Hz] retained (frequency
        decoupling).
    threshold : float, optional
        Threshold of the ``'Elbow'`` criterion used to select the number of
        modes (default 0.999).
    nmod : int, optional
        Number of modes to retain; if given, the Elbow criterion is skipped.
    legacy_transform : bool
        ``False`` (default) transforms the filtered correlation matrix back
        with the inverse DFT, ``K = real(conj(F) @ (K_hat * H) @ conj(F))``.
        ``True`` reproduces the MATLAB implementation used for Robledo et al.
        (2025), which applies the *forward* DFT matrix on both sides also when
        transforming back (``K = real(F * (K_hat .* H) * F)``). Since ``F @ F``
        is the index-reversal permutation, that returns the filtered
        correlation matrix with its time index reversed
        (``K_legacy[m, n] = K_exact[-m mod nt, -n mod nt]``): the eigenvalues
        are identical but the temporal modes are time-reversed and the data are
        projected onto them, which damps the reconstructed fluctuations (about
        -23 % in the fluctuating Nusselt number of the sweeping-jet dataset,
        see ``validation/REPORT.md``).
    return_nmod : bool
        Also return the number of retained modes.

    Returns
    -------
    Xf : ndarray
        Filtered snapshots (same shape/dtype as ``X``).
    nmod : int
        Number of retained modes (if ``return_nmod``).
    """
    X = np.asarray(X)
    X3 = ensure_3d(X)
    ny, nx, nt = X3.shape
    if nt < 4:
        raise ValueError("The mPOD filter requires a temporal dimension with several snapshots")

    D = X3.reshape(-1, nt).astype(np.float64, copy=True)
    D_mean = D.mean(axis=1)
    D -= D_mean[:, None]

    K = D.T @ D
    F = _dft_matrix(nt)
    K_hat = F @ K @ F
    del K

    H = mpod_transfer_function(nt, f_acquisition, mode, fpeaks=fpeaks, w=w,
                               n_regions=n_regions, f_min=f_min, f_max=f_max)
    K_hat *= H
    del H

    if legacy_transform:
        Kf = (F @ K_hat @ F).real
    else:
        Fc = F.conj()
        Kf = (Fc @ K_hat @ Fc).real
    del K_hat, F

    # svd(K,'econ') of a symmetric PSD matrix: singular values = eigenvalues
    U, S, _ = np.linalg.svd(Kf, hermitian=True)
    del Kf

    if nmod is None:
        if threshold is None:
            warnings.warn("mPOD_filter: the Elbow criterion with default threshold 0.999 is used")
        nmod = find_nmod(S, criterion="Elbow", threshold=threshold, verbose=verbose)
    nmod = int(nmod)
    if nmod < 1 or nmod > nt:
        raise ValueError(f"nmod must be between 1 and {nt}")

    # PHI = D U S^-1/2 ; D_new = PHI_k S_k^1/2 U_k' = D U_k U_k'
    Uk = U[:, :nmod]
    Xf = (D @ Uk) @ Uk.T
    del D
    Xf += D_mean[:, None]
    Xf = Xf.reshape(ny, nx, nt).astype(X.dtype, copy=False)
    if X.ndim == 2:
        Xf = Xf[:, :, 0]
    return (Xf, nmod) if return_nmod else Xf
