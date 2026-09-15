"""Filtering module of PIRT (Python).

Available filters:

* :func:`pod_filter` - POD low-order reconstruction (with :func:`find_nmod`).
* :func:`mpod_filter` - multi-scale POD (Mendez et al., 2019).
* :func:`sgolay32_filter` - 3-D Savitzky-Golay smoothing with derivatives.
* :func:`wiener3` - 3-D adaptive Wiener filter.
* :func:`gaussian_filter3` - 3-D Gaussian filter (``imgaussfilt3``).
* :func:`spatial_cutoff_filter`, :func:`temporal_cutoff_filter`,
  :func:`cutoff_filter_3d` - spectral cut-off filters.
"""
from .pod import pod_filter, pod_decomposition, find_nmod, optimal_svht_coef, matlab_smooth
from .mpod import mpod_filter, mpod_transfer_function
from .sgolay32 import sgolay32_coef, sgolay32_filter, remove_edges
from .wiener3 import wiener3
from .gaussian import gaussian_filter3, gaussian_kernel_1d
from .cutoff import (lowpass_kernel, highpass_kernel, spatial_cutoff_filter,
                     temporal_cutoff_filter, cutoff_filter_3d)

__all__ = [
    "pod_filter", "pod_decomposition", "find_nmod", "optimal_svht_coef", "matlab_smooth",
    "mpod_filter", "mpod_transfer_function",
    "sgolay32_coef", "sgolay32_filter", "remove_edges",
    "wiener3",
    "gaussian_filter3", "gaussian_kernel_1d",
    "lowpass_kernel", "highpass_kernel", "spatial_cutoff_filter",
    "temporal_cutoff_filter", "cutoff_filter_3d",
]
