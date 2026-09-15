"""Spatial and temporal cut-off (low-/high-/band-pass) filters.

Ports of ``lowpass_kernel.m``, ``highpass_kernel.m``,
``Spatial_Cutoff_Filter.m``, ``Temporal_Cutoff_Filter.m`` and
``Cutoff_Filter_3D.m`` (I. Robledo, UC3M).

Conventions
-----------
* Arrays are ``(ny, nx, nt)`` = (rows, columns, time).
* Spatial cut-off frequencies are normalised to the Nyquist frequency of each
  direction (``0 < f <= 1``). A scalar applies to both directions (circular
  mask); a pair ``(f_rows, f_cols)`` defines an elliptical mask whose semi-axes
  are given in the array-axis order, like the MATLAB 2-element vector.
* Temporal cut-off frequencies are in Hz and require the sampling frequency.
* ``fc_low`` is the cut-off of a low-pass (frequencies below it are kept);
  ``fc_high`` is the cut-off of a high-pass (frequencies above it are kept).
  If both are given the result is a band-pass keeping the frequencies between
  the two values, whichever order they are given in.
"""
from __future__ import annotations

import warnings

import numpy as np
from scipy import signal

from ._common import ensure_3d

__all__ = ["lowpass_kernel", "highpass_kernel", "spatial_cutoff_filter",
           "temporal_cutoff_filter", "cutoff_filter_3d"]


def _elliptic_mask(n: int, m: int, fn: float, fm: float, lowpass: bool) -> np.ndarray:
    """Elliptical low-/high-pass mask for the fft-shifted 2-D spectrum.

    ``fn`` is the normalised cut-off along the columns (``m``) and ``fm`` along
    the rows (``n``), as in the MATLAB kernels. The semi-axes in frequency bins
    are ``floor(f * size/2)``. The mask is centred on the DC bin of
    ``fftshift`` (index ``size//2``) and therefore Hermitian-symmetric, so the
    filtered field is real.

    Note: the MATLAB kernels (``lowpass_kernel.m``/``highpass_kernel.m``) built
    the frequency grid as ``(1:m) - m/2``, which is offset by one bin with
    respect to the DC location; that offset is corrected here (and in the
    MATLAB toolbox from v1.1).
    """
    if fn > 1 or fm > 1:
        raise ValueError("Normalised cut-off frequencies must not exceed 1 (Nyquist)")
    if fn < 0 or fm < 0:
        raise ValueError("Cut-off frequencies must be non-negative")
    fn_b = np.floor(fn * m * 0.5)
    fm_b = np.floor(fm * n * 0.5)
    if fn_b < 1 or fm_b < 1:
        warnings.warn(
            f"cutoff: the normalised cut-off ({fm:g} rows, {fn:g} cols) is below the frequency "
            f"resolution of the {n}x{m} image; the mask degenerates (no filtering in that direction)"
        )
    kx = (np.arange(m) - m // 2)[None, :].astype(float)
    ky = (np.arange(n) - n // 2)[:, None].astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        rho2 = (kx / fn_b) ** 2 + (ky / fm_b) ** 2  # inf where a semi-axis is 0
        rho2 = np.where(np.isnan(rho2), 0.0, rho2)  # 0/0 at the centre
    inside = rho2 <= 1.0
    mask = inside if lowpass else ~inside | (rho2 == 1.0)
    return mask.astype(np.float32)


def lowpass_kernel(n: int, m: int, fn: float, fm: float) -> np.ndarray:
    """Elliptical low-pass mask (``n`` rows, ``m`` columns), centred spectrum."""
    return _elliptic_mask(n, m, fn, fm, lowpass=True)


def highpass_kernel(n: int, m: int, fn: float, fm: float) -> np.ndarray:
    """Elliptical high-pass mask (``n`` rows, ``m`` columns), centred spectrum."""
    return _elliptic_mask(n, m, fn, fm, lowpass=False)


def _pair(f, name):
    f = np.atleast_1d(np.asarray(f, dtype=float)).ravel()
    if f.size == 1:
        return float(f[0]), float(f[0])  # (rows, cols)
    if f.size == 2:
        return float(f[0]), float(f[1])
    raise ValueError(f"{name} must be a scalar or a 2-element sequence (f_rows, f_cols)")


def spatial_cutoff_filter(X, fc_low=None, fc_high=None):
    """2-D spectral cut-off filtering of every snapshot.

    Parameters
    ----------
    X : ndarray, shape (ny, nx) or (ny, nx, nt)
    fc_low : float or (f_rows, f_cols), optional
        Low-pass cut-off (normalised to Nyquist).
    fc_high : float or (f_rows, f_cols), optional
        High-pass cut-off (normalised to Nyquist).

    If neither is given, a low-pass with cut-off at half the Nyquist frequency
    is applied in both directions. If both are given, a band-pass keeping the
    frequencies between them is applied.

    Returns
    -------
    Xf : ndarray
        Filtered array (real part of the inverse FFT), same shape/dtype as X.
    """
    X = np.asarray(X)
    X3 = ensure_3d(X)
    n, m = X3.shape[:2]

    if fc_low is None and fc_high is None:
        H = lowpass_kernel(n, m, 0.5, 0.5)
    elif fc_low is not None and fc_high is None:
        fr, fc = _pair(fc_low, "fc_low")
        H = lowpass_kernel(n, m, fc, fr)
    elif fc_low is None and fc_high is not None:
        fr, fc = _pair(fc_high, "fc_high")
        H = highpass_kernel(n, m, fc, fr)
    else:
        lo_r, lo_c = _pair(fc_low, "fc_low")
        hi_r, hi_c = _pair(fc_high, "fc_high")
        # keep the band between the two cut-offs regardless of their order
        lp = (max(lo_r, hi_r), max(lo_c, hi_c))
        hp = (min(lo_r, hi_r), min(lo_c, hi_c))
        H = lowpass_kernel(n, m, lp[1], lp[0]) * highpass_kernel(n, m, hp[1], hp[0])

    Freq = np.fft.fftshift(np.fft.fft2(X3, axes=(0, 1)), axes=(0, 1))
    Freq *= H[:, :, None]
    Xf = np.fft.ifft2(np.fft.ifftshift(Freq, axes=(0, 1)), axes=(0, 1)).real
    Xf = Xf.astype(X3.dtype, copy=False)
    if X.ndim == 2:
        Xf = Xf[:, :, 0]
    return Xf


def _fir_filter_along_time(X3, fs, fc_low=None, fc_high=None, steepness=0.85, attenuation_db=60.0):
    """Zero-phase FIR filtering along the last axis.

    The design follows the defaults of MATLAB's ``lowpass``/``highpass``/
    ``bandpass``: minimum-order FIR with 60 dB stopband attenuation and a
    transition width set by the ``steepness`` parameter (0.85), delay
    compensated. The exact MATLAB design (equiripple) is not reproduced; the
    Kaiser-window design used here meets the same specifications.
    """
    fnyq = fs / 2.0
    if fc_low is not None and fc_high is not None:
        f1, f2 = sorted((float(fc_high), float(fc_low)))
        if not (0 < f1 < f2 < fnyq):
            raise ValueError("Band-pass cut-offs must satisfy 0 < f_low_edge < f_high_edge < fs/2")
        tw = min((1 - steepness) * f1, (1 - steepness) * (fnyq - f2))
        cutoff = [f1, f2]
        pass_zero = False
    elif fc_low is not None:
        fc = float(fc_low)
        if not (0 < fc < fnyq):
            raise ValueError("The low-pass cut-off must satisfy 0 < fc < fs/2")
        tw = (1 - steepness) * (fnyq - fc)
        cutoff = fc
        pass_zero = True
    else:
        fc = float(fc_high)
        if not (0 < fc < fnyq):
            raise ValueError("The high-pass cut-off must satisfy 0 < fc < fs/2")
        tw = (1 - steepness) * fc
        cutoff = fc
        pass_zero = False

    numtaps, beta = signal.kaiserord(attenuation_db, tw / fnyq)
    numtaps = int(numtaps) | 1  # odd length -> integer group delay
    nt = X3.shape[2]
    if numtaps > nt:
        raise ValueError(
            f"The temporal cut-off filter needs {numtaps} taps but only {nt} snapshots are available"
        )
    taps = signal.firwin(numtaps, cutoff, window=("kaiser", beta), pass_zero=pass_zero, fs=fs)
    Xf = signal.fftconvolve(X3, taps[None, None, :], mode="same", axes=2)
    return Xf


def temporal_cutoff_filter(X, fs: float, fc_low=None, fc_high=None):
    """Low-, high- or band-pass filtering along the time axis (third dimension).

    Parameters
    ----------
    X : ndarray, shape (ny, nx, nt)
    fs : float
        Sampling frequency [Hz].
    fc_low : float, optional
        Low-pass cut-off [Hz] (frequencies below are kept).
    fc_high : float, optional
        High-pass cut-off [Hz] (frequencies above are kept). If both cut-offs
        are given the band between them is kept.
    """
    X = np.asarray(X)
    if X.ndim != 3:
        raise ValueError("temporal_cutoff_filter: only valid for 3-D arrays")
    if fc_low is None and fc_high is None:
        raise ValueError("temporal_cutoff_filter: a cut-off frequency is required")
    Xf = _fir_filter_along_time(X.astype(np.float64, copy=False), float(fs), fc_low, fc_high)
    return Xf.astype(X.dtype, copy=False)


def cutoff_filter_3d(X, spatial: dict | None = None, temporal: dict | None = None):
    """Combined spatial and temporal cut-off filtering (``Cutoff_Filter_3D``).

    Parameters
    ----------
    spatial : dict, optional
        Keys ``'fl'`` (low-pass cut-off) and/or ``'fh'`` (high-pass cut-off),
        normalised to Nyquist. An empty dict applies the default low-pass.
    temporal : dict, optional
        Keys ``'fs'`` (required, Hz), ``'fl'`` and/or ``'fh'`` (Hz).
    """
    if spatial is None and temporal is None:
        raise ValueError("cutoff_filter_3d: neither spatial nor temporal filtering was selected")
    Xf = np.asarray(X)
    if spatial is not None:
        Xf = spatial_cutoff_filter(Xf, fc_low=spatial.get("fl"), fc_high=spatial.get("fh"))
    if temporal is not None:
        if "fs" not in temporal:
            raise ValueError("The sampling frequency 'fs' is required for temporal cut-off filtering")
        Xf = temporal_cutoff_filter(Xf, temporal["fs"], fc_low=temporal.get("fl"),
                                    fc_high=temporal.get("fh"))
    return Xf
