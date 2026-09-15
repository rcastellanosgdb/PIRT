"""Unit tests of the PIRT filters against analytical / brute-force references."""
import numpy as np
import pytest
from scipy import ndimage

from pirt.filters import (find_nmod, gaussian_filter3, gaussian_kernel_1d, highpass_kernel,
                          lowpass_kernel, matlab_smooth, mpod_filter, mpod_transfer_function,
                          optimal_svht_coef, pod_filter, sgolay32_coef, sgolay32_filter,
                          spatial_cutoff_filter, temporal_cutoff_filter, wiener3)


# --------------------------------------------------------------------------- #
# Savitzky-Golay
# --------------------------------------------------------------------------- #
def _brute_force_sg(X, i, j, k, w):
    """Least-squares quadratic fit around (i, j, k) with unit spacing."""
    wc = np.array(w) // 2
    rows, vals = [], []
    for kk in range(-wc[2], wc[2] + 1):
        for jj in range(-wc[1], wc[1] + 1):
            for ii in range(-wc[0], wc[0] + 1):
                rows.append([1, ii, jj, kk, ii * jj, ii * kk, jj * kk, ii * ii, jj * jj, kk * kk])
                vals.append(X[i + ii, j + jj, k + kk])
    return np.linalg.lstsq(np.array(rows), np.array(vals), rcond=None)[0]


@pytest.mark.parametrize("method", ["direct", "separable"])
def test_sgolay32_matches_local_least_squares(rng, method):
    X = rng.standard_normal((9, 11, 7))
    w, h = (5, 3, 3), (0.2, 0.3, 0.05)
    Tf, dTdt, d2x, d2y = sgolay32_filter(X, w, h, method=method)
    c = _brute_force_sg(X, 4, 5, 3, w)
    kt = 3 - 1  # one snapshot cropped at the start
    assert Tf.shape == (9, 11, 5)
    assert Tf[4, 5, kt] == pytest.approx(c[0], rel=1e-10)
    assert dTdt[4, 5, kt] == pytest.approx(c[3] / h[2], rel=1e-10)
    assert d2x[4, 5, kt] == pytest.approx(2 * c[8] / h[0] ** 2, rel=1e-10)  # jj -> columns
    assert d2y[4, 5, kt] == pytest.approx(2 * c[7] / h[1] ** 2, rel=1e-10)  # ii -> rows


def test_sgolay32_exact_on_quadratic_field():
    dx, dy, dt = 0.2, 0.3, 0.05
    y = np.arange(20)[:, None, None] * dy
    x = np.arange(30)[None, :, None] * dx
    t = np.arange(10)[None, None, :] * dt
    T = 3 + 2 * x - y + 0.5 * t + x * y + 4 * x**2 - 2 * y**2 + 3 * t**2 + x * t
    Tf, dTdt, d2x, d2y = sgolay32_filter(T, (5, 5, 3), (dx, dy, dt))
    inner = (slice(2, -2), slice(2, -2), slice(None))
    np.testing.assert_allclose(Tf[inner], T[:, :, 1:-1][inner], atol=1e-10)
    np.testing.assert_allclose(d2x[inner], 8.0, atol=1e-9)
    np.testing.assert_allclose(d2y[inner], -4.0, atol=1e-9)
    expected_dTdt = np.broadcast_to(0.5 + 6 * t + x, T.shape)[:, :, 1:-1]
    np.testing.assert_allclose(dTdt[inner], expected_dTdt[inner], atol=1e-9)


def test_sgolay32_separable_equals_direct(rng):
    X = rng.standard_normal((14, 21, 9)) * 5 + 300
    for w in [(5, 5, 3), (11, 11, 3), (3, 7, 5)]:
        A = sgolay32_filter(X, w, (0.2, 0.3, 0.05), method="direct")
        B = sgolay32_filter(X, w, (0.2, 0.3, 0.05), method="separable")
        for a, b in zip(A, B):
            np.testing.assert_allclose(a, b, rtol=1e-11, atol=1e-9)


def test_sgolay32_edges_replicated_and_time_cropped(rng):
    X = rng.standard_normal((8, 9, 6))
    Tf = sgolay32_filter(X, (5, 3, 3), return_derivatives=False)
    assert Tf.shape == (8, 9, 4)
    np.testing.assert_array_equal(Tf[0], Tf[2])
    np.testing.assert_array_equal(Tf[-1], Tf[-3])
    np.testing.assert_array_equal(Tf[:, 0], Tf[:, 1])


def test_sgolay32_rejects_even_or_large_kernels(rng):
    X = rng.standard_normal((6, 6, 4))
    with pytest.raises(ValueError):
        sgolay32_filter(X, (4, 3, 3))
    with pytest.raises(ValueError):
        sgolay32_filter(X, (7, 3, 3))


def test_sgolay32_coef_shape_and_constant_preservation():
    C = sgolay32_coef((5, 5, 3))
    assert C.shape == (10, 75)
    assert C[0].sum() == pytest.approx(1.0)  # a0 kernel reproduces constants
    assert abs(C[3].sum()) < 1e-12  # derivative kernels kill constants


# --------------------------------------------------------------------------- #
# POD / truncation criteria
# --------------------------------------------------------------------------- #
def test_pod_filter_exact_for_low_rank_data(rng):
    A = rng.standard_normal((30, 40, 3))
    B = rng.standard_normal((3, 50))
    X = np.einsum("ijk,kt->ijt", A, B) + 5.0
    Xf, nmod = pod_filter(X, nmod=3, return_nmod=True, verbose=False)
    assert nmod == 3
    np.testing.assert_allclose(Xf, X, atol=1e-10)
    Xf1 = pod_filter(X, nmod=1, verbose=False)
    assert np.abs(Xf1 - X).max() > 1e-3  # truncation actually removes content


def test_pod_filter_preserves_mean_and_dtype(rng):
    X = (rng.standard_normal((10, 12, 20)) + 300).astype(np.float32)
    Xf = pod_filter(X, nmod=2, verbose=False)
    assert Xf.dtype == np.float32
    np.testing.assert_allclose(Xf.mean(axis=2), X.mean(axis=2), rtol=1e-6)


def test_pod_filter_more_snapshots_than_pixels(rng):
    X = rng.standard_normal((3, 4, 40))
    Xf = pod_filter(X, nmod=12, verbose=False)  # rank is 11 (12 pixels, mean removed)
    np.testing.assert_allclose(Xf, X, atol=1e-10)


def test_optimal_svht_coef_reference_values():
    # Gavish & Donoho (2014): lambda*(1) = 4/sqrt(3); omega(1) ~ 2.858 (Table / eq. (5)).
    assert optimal_svht_coef(1.0, sigma_known=True) == pytest.approx(4 / np.sqrt(3), rel=1e-12)
    assert optimal_svht_coef(1.0, sigma_known=False) == pytest.approx(2.858, abs=2e-3)
    # polynomial approximation of omega(beta), accurate to ~1%
    for b in (0.5, 0.1, 0.01):
        approx = 0.56 * b**3 - 0.95 * b**2 + 1.82 * b + 1.43
        assert optimal_svht_coef(b, sigma_known=False) == pytest.approx(approx, rel=0.015)
    with pytest.raises(ValueError):
        optimal_svht_coef(1.5)


def test_find_nmod_criteria():
    s = np.array([10.0, 5.0, 2.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
    lam = s**2
    energy = np.cumsum(lam) / lam.sum()
    assert find_nmod(s, "Spectrum", threshold=0.99, verbose=False) == int(np.flatnonzero(energy > 0.99)[0]) + 1
    # HardThreshold: first singular value below coef*median(s)
    coef = optimal_svht_coef(0.1, False)
    expected = int(np.flatnonzero(s < coef * np.median(s))[0]) + 1
    assert find_nmod(s, "HardThreshold", beta=0.1, verbose=False) == expected
    # a diagonal matrix is accepted as in MATLAB
    assert find_nmod(np.diag(s), "Spectrum", threshold=0.99, verbose=False) == find_nmod(s, "Spectrum", 0.99, verbose=False)
    with pytest.raises(ValueError):
        find_nmod(s, "HardThreshold", verbose=False)
    with pytest.raises(ValueError):
        find_nmod(s, "unknown", verbose=False)


def test_find_nmod_elbow_replicates_matlab_indexing():
    lam = np.r_[np.geomspace(1e4, 10, 15), np.full(30, 1.0)]  # flat tail after mode 15
    s = np.sqrt(lam)
    ratio = lam[1:] / lam[:-1]
    f = matlab_smooth(ratio[4:-5])  # MATLAB: smooth(f(skip:end-skip)) with skip=5
    expected = int(np.flatnonzero(f > 0.999)[0]) + 1
    assert find_nmod(s, "Elbow", threshold=0.999, verbose=False) == expected


def test_matlab_smooth_end_handling():
    y = np.array([1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0])
    out = matlab_smooth(y, 5)
    assert out[0] == 1.0
    assert out[1] == pytest.approx((1 + 2 + 4) / 3)
    assert out[2] == pytest.approx((1 + 2 + 4 + 8 + 16) / 5)
    assert out[-1] == 64.0
    assert out[-2] == pytest.approx((16 + 32 + 64) / 3)


# --------------------------------------------------------------------------- #
# Wiener
# --------------------------------------------------------------------------- #
def _box_sum_matlab(A, nh):
    """convn(A, ones(nh), 'same') with zero padding (MATLAB centring for even sizes)."""
    out = np.zeros_like(A)
    rngs = [range(-(L // 2) + (1 if L % 2 == 0 else 0), L // 2 + 1) for L in nh]
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            for k in range(A.shape[2]):
                s = 0.0
                for di in rngs[0]:
                    for dj in rngs[1]:
                        for dk in rngs[2]:
                            ii, jj, kk = i + di, j + dj, k + dk
                            if 0 <= ii < A.shape[0] and 0 <= jj < A.shape[1] and 0 <= kk < A.shape[2]:
                                s += A[ii, jj, kk]
                out[i, j, k] = s
    return out


@pytest.mark.parametrize("nh", [(3, 3, 3), (3, 4, 1), (7, 7, 1)])
def test_wiener3_matches_brute_force(rng, nh):
    X = rng.standard_normal((8, 9, 5))
    f = wiener3(X, nh, noise=0.5)
    n = np.prod(nh)
    lm = _box_sum_matlab(X, nh) / n
    lv = _box_sum_matlab(X * X, nh) / n - lm**2
    expected = lm + np.maximum(lv - 0.5, 0) / np.maximum(lv, 0.5) * (X - lm)
    np.testing.assert_allclose(f, expected, atol=1e-12)


def test_wiener3_reduces_noise_and_estimates_power(rng):
    clean = np.full((60, 60, 8), 300.0)
    noisy = clean + 0.3 * rng.standard_normal(clean.shape)
    f, noise = wiener3(noisy, (7, 7, 1), return_noise=True)
    inner = (slice(10, -10), slice(10, -10), slice(None))
    assert np.std(f[inner]) < 0.5 * np.std(noisy[inner])
    assert noise > 0


def test_wiener3_scalar_kernel_and_dtype(rng):
    X = rng.standard_normal((6, 6, 6)).astype(np.float32)
    f = wiener3(X, 3)
    assert f.dtype == np.float32 and f.shape == X.shape


# --------------------------------------------------------------------------- #
# Gaussian
# --------------------------------------------------------------------------- #
def test_gaussian_filter3_equals_direct_correlation(rng):
    X = rng.standard_normal((8, 9, 4))
    g = gaussian_filter3(X, (3, 3, 0.1), (9, 3, 1))
    K = gaussian_kernel_1d(3, 9)[:, None, None] * gaussian_kernel_1d(3, 3)[None, :, None]
    np.testing.assert_allclose(g, ndimage.correlate(X, K, mode="nearest"), atol=1e-13)


def test_gaussian_filter3_defaults_and_constant(rng):
    X = np.full((8, 9, 4), 2.5)
    np.testing.assert_allclose(gaussian_filter3(X), 2.5)
    k = gaussian_kernel_1d(0.5, 3)
    assert k.sum() == pytest.approx(1.0)
    with pytest.raises(ValueError):
        gaussian_filter3(X, sigma=1.0, filter_size=(4, 3, 3))


def test_gaussian_filter3_padding_modes(rng):
    X = rng.standard_normal((10, 10, 3))
    for pad in ("replicate", "symmetric", "circular", 0.0):
        assert gaussian_filter3(X, 1.0, 5, padding=pad).shape == X.shape


# --------------------------------------------------------------------------- #
# Cut-off filters
# --------------------------------------------------------------------------- #
def test_spatial_cutoff_low_high_pass():
    ny, nx, nt = 128, 96, 3
    yy, xx = np.mgrid[0:ny, 0:nx]
    low = np.cos(2 * np.pi * 2 * xx / nx)     # 2 bins  -> normalised 2/48
    high = np.cos(2 * np.pi * 40 * xx / nx)   # 40 bins -> normalised 40/48
    X = (10 + low + high)[:, :, None] * np.ones((1, 1, nt))
    lp = spatial_cutoff_filter(X, fc_low=0.5)
    hp = spatial_cutoff_filter(X, fc_high=0.5)
    bp = spatial_cutoff_filter(X, fc_low=0.5, fc_high=0.03)  # high-pass semi-axes of 1 bin: removes DC only
    np.testing.assert_allclose(lp, np.broadcast_to((10 + low)[:, :, None], X.shape), atol=1e-10)
    np.testing.assert_allclose(hp, np.broadcast_to(high[:, :, None], X.shape), atol=1e-10)
    np.testing.assert_allclose(bp, np.broadcast_to(low[:, :, None], X.shape), atol=1e-10)
    with pytest.warns(UserWarning):  # cut-off below the frequency resolution
        spatial_cutoff_filter(X, fc_high=0.005)


@pytest.mark.parametrize("n, m", [(20, 30), (21, 31), (20, 31)])
def test_cutoff_kernels_are_complementary_and_symmetric(n, m):
    H_lo = lowpass_kernel(n, m, 0.4, 0.6)
    H_hi = highpass_kernel(n, m, 0.4, 0.6)
    assert H_lo.shape == (n, m)
    # every frequency belongs to at least one of them (boundary shared)
    assert np.all(H_lo + H_hi >= 1)
    assert 0 < H_lo.sum() < n * m
    # Hermitian symmetry H(k) = H(-k) about the DC bin of fftshift (index size//2)
    for H in (H_lo, H_hi):
        Hs = np.roll(H[::-1, ::-1], (1 if n % 2 == 0 else 0, 1 if m % 2 == 0 else 0), axis=(0, 1))
        np.testing.assert_array_equal(Hs, H)
    # DC is kept by the low-pass and removed by the high-pass
    assert H_lo[n // 2, m // 2] == 1 and H_hi[n // 2, m // 2] == 0
    with pytest.raises(ValueError):
        lowpass_kernel(n, m, 1.5, 0.5)


def test_temporal_cutoff_filters():
    fs = 100.0
    t = np.arange(600) / fs
    X = (np.sin(2 * np.pi * 1 * t) + np.sin(2 * np.pi * 30 * t))[None, None, :] * np.ones((2, 2, 1))
    lp = temporal_cutoff_filter(X, fs, fc_low=5.0)
    hp = temporal_cutoff_filter(X, fs, fc_high=10.0)
    bp = temporal_cutoff_filter(X, fs, fc_low=40.0, fc_high=10.0)
    sl = slice(80, -80)
    assert np.abs(lp[0, 0, sl] - np.sin(2 * np.pi * t[sl])).max() < 5e-3
    assert np.abs(hp[0, 0, sl] - np.sin(2 * np.pi * 30 * t[sl])).max() < 5e-3
    assert np.abs(bp[0, 0, sl] - np.sin(2 * np.pi * 30 * t[sl])).max() < 5e-3
    with pytest.raises(ValueError):
        temporal_cutoff_filter(X, fs)


# --------------------------------------------------------------------------- #
# mPOD
# --------------------------------------------------------------------------- #
def test_mpod_transfer_function_indices():
    nt, fs = 400, 80.0  # df = 0.2 Hz
    H = mpod_transfer_function(nt, fs, "Peak Removal", fpeaks=[30], w=2)
    # MATLAB: idx1 = round(28/0.2) = 140, idx2 = round(32/0.2) = 160 -> rows 140:160 and 240:260 (1-based)
    zero_rows = np.flatnonzero(H[:, 0] == 0)
    np.testing.assert_array_equal(zero_rows, np.r_[139:160, 239:260])
    Hd = mpod_transfer_function(nt, fs, "frequency_decoupling", n_regions=2, f_min=1, f_max=5)
    # idx0 = 5, delta = 10 -> bands (1-based) 6:15, 16:25 and their mirrors
    kept = np.flatnonzero(Hd.sum(axis=1) > 0)
    np.testing.assert_array_equal(kept, np.r_[5:25, 375:395])
    with pytest.raises(ValueError):
        mpod_transfer_function(nt, fs, "Peak Removal", fpeaks=[0.1], w=1)


@pytest.mark.parametrize("legacy", [True, False])
def test_mpod_peak_removal_suppresses_spurious_frequency(rng, legacy):
    nt, fs = 400, 80.0
    t = np.arange(nt) / fs
    phi1, phi2 = rng.standard_normal((10, 12)), rng.standard_normal((10, 12))
    signal = 300 + phi1[:, :, None] * np.sin(2 * np.pi * 5 * t)
    X = signal + 0.3 * phi2[:, :, None] * np.sin(2 * np.pi * 30 * t)
    Xf, nmod = mpod_filter(X, fs, "Peak Removal", fpeaks=[30], w=2, nmod=1, legacy_transform=legacy,
                           verbose=False)
    assert nmod == 1
    assert np.abs(Xf - signal).max() < 0.1 * np.abs(X - signal).max()


def test_mpod_accepts_matlab_mode_names_and_requires_inputs(rng):
    X = rng.standard_normal((4, 5, 64)) + 300
    for name in ("PeakRem", "Peak Removal", "peak_removal"):
        mpod_filter(X, 64.0, name, fpeaks=[10], w=1, nmod=2, verbose=False)
    for name in ("Freq Deco", "FreqDec", "frequency_decoupling"):
        mpod_filter(X, 64.0, name, n_regions=2, f_min=1, f_max=20, nmod=2, verbose=False)
    with pytest.raises(ValueError):
        mpod_filter(X, 64.0, "Peak Removal", nmod=2, verbose=False)
    with pytest.raises(ValueError):
        mpod_filter(X, 64.0, "nonsense", nmod=2, verbose=False)
