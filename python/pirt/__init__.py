"""PIRT - Processing InfraRed Thermography (Python version).

Toolbox for filtering infrared temperature snapshots and computing convective
heat-transfer maps (h, Nu, St) with the heated-thin-foil technique, including
Monte Carlo uncertainty estimation.

Universidad Carlos III de Madrid - Experimental Aerodynamics and Propulsion Lab.
"""
from .pirt import PIRT
from .filters import (pod_filter, find_nmod, optimal_svht_coef, mpod_filter, sgolay32_filter,
                      sgolay32_coef, wiener3, gaussian_filter3, spatial_cutoff_filter,
                      temporal_cutoff_filter, cutoff_filter_3d)
from .derivatives import derivative_fd
from .heat_transfer import calculate_heat_transfer, air_thermal_conductivity
from .uncertainty import montecarlo_uncertainty
from . import io

__version__ = "1.0.0"
__all__ = [
    "PIRT", "pod_filter", "find_nmod", "optimal_svht_coef", "mpod_filter", "sgolay32_filter",
    "sgolay32_coef", "wiener3", "gaussian_filter3", "spatial_cutoff_filter",
    "temporal_cutoff_filter", "cutoff_filter_3d", "derivative_fd", "calculate_heat_transfer",
    "air_thermal_conductivity", "montecarlo_uncertainty", "io",
]
