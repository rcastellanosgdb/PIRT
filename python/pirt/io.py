"""Input/output helpers: MATLAB ``.mat`` files and PIRT test-case folders."""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np

__all__ = ["load_mat", "read_conf_data", "load_sj_case", "save_results"]


def load_mat(path, variables=None, squeeze: bool = True) -> dict:
    """Load numeric variables from a MATLAB ``.mat`` file (v5/v7 or v7.3).

    Arrays are returned with the MATLAB shape (v7.3 files store the data
    transposed and are transposed back here).

    Parameters
    ----------
    path : str or Path
    variables : sequence of str, optional
        Names to load (default: all numeric variables).
    squeeze : bool
        Convert 1x1 arrays to Python floats.
    """
    path = Path(path)
    out = {}
    is_hdf5 = False
    try:
        import h5py
        is_hdf5 = h5py.is_hdf5(path)
    except ImportError:
        h5py = None
    if is_hdf5:  # MATLAB v7.3 files are HDF5 containers
        with h5py.File(path, "r") as f:
            names = list(variables) if variables else [k for k in f.keys() if k != "#refs#"]
            for k in names:
                obj = f[k]
                if isinstance(obj, h5py.Dataset):
                    arr = obj[()]
                    if isinstance(arr, np.ndarray) and arr.dtype.kind in "biuf":
                        out[k] = np.ascontiguousarray(arr.T)
    else:
        import scipy.io as sio
        try:
            data = sio.loadmat(path, variable_names=list(variables) if variables else None)
        except NotImplementedError as exc:  # v7.3 file but h5py not installed
            raise ImportError("Reading MATLAB v7.3 files requires h5py (pip install h5py)") from exc
        for k, v in data.items():
            if k.startswith("__"):
                continue
            out[k] = np.asarray(v)
    if squeeze:
        for k, v in out.items():
            if isinstance(v, np.ndarray) and v.size == 1:
                out[k] = float(v.ravel()[0])
    return out


def read_conf_data(path) -> dict:
    """Parse a ``CONF_DATA.STR`` file (heated-thin-foil sensor data).

    The fixed line layout of the original MATLAB ``load_case`` is used: after
    five header lines the characteristic length ``L`` is read, then one header
    line, then ``s``, ``rho``, ``cp``, ``lambda``, ``H``, ``W`` and ``epsilon``
    (one value per line, the first token of each line).
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()

    def first_number(line):
        return float(line.split()[0])

    keys = ["s", "rho", "cp", "lambda", "H", "W", "epsilon"]
    conf = {"L": first_number(lines[5])}
    for i, k in enumerate(keys):
        conf[k] = first_number(lines[7 + i])
    return conf


def load_sj_case(path, n_snapshots: int | None = None, dtype=None) -> dict:
    """Load an impinging-sweeping-jet test case in the PIRT example format.

    The folder must contain ``Resolution.mat`` (``dx``, ``dy`` in px/mm),
    ``TestConditions.mat`` (``V``, ``I``, ``Uinf``, ``Tamb`` [C], ``dt``),
    ``Thot.mat`` (``Timage_hot``), ``Tcold.mat`` (``Timage_cold``) and
    ``CONF_DATA.STR``. This reproduces the ``load_case`` function of the
    MATLAB example ``SJ_main.m`` (including the paint properties of the
    experiment of Robledo et al., 2025).

    Returns
    -------
    case : dict
        Keys ``Thot``, ``Tcold``, ``HFS``, ``Conditions``, ``f_acq`` [Hz],
        ``dx``, ``dy`` [m/px].
    """
    path = Path(path)
    res = load_mat(path / "Resolution.mat")
    tc = load_mat(path / "TestConditions.mat")
    dt = float(tc["dt"])
    conditions = {
        "V": float(tc["V"]), "I": float(tc["I"]), "Uinf": float(tc["Uinf"]),
        "Tamb": np.asarray(tc["Tamb"], dtype=float).ravel() + 273.15,
        "dt": dt,
    }
    Thot = load_mat(path / "Thot.mat", ["Timage_hot"])["Timage_hot"]
    Tcold = load_mat(path / "Tcold.mat", ["Timage_cold"])["Timage_cold"]
    if n_snapshots is not None:
        Thot = Thot[:, :, :n_snapshots]
        if Tcold.ndim == 3:
            Tcold = Tcold[:, :, :n_snapshots]
    if dtype is not None:
        Thot = Thot.astype(dtype, copy=False)
        Tcold = Tcold.astype(dtype, copy=False)

    conf = read_conf_data(path / "CONF_DATA.STR")
    conditions["L"] = conf["L"]
    hfs = {
        "s": conf["s"], "rho": conf["rho"], "cp": conf["cp"],
        "k": conf["lambda"] * conf["s"],  # conductance [W/K]
        "H": conf["H"], "W": conf["W"], "epsilon": conf["epsilon"],
        "sides": 2,
        "s_paint": 42e-6, "cp_paint": 3061.5, "rho_paint": 1261.175, "lambda_paint": 1.38,
    }
    dx = 1.0 / float(res["dx"]) / 1000.0  # px/mm -> m/px
    dy = 1.0 / float(res["dy"]) / 1000.0
    conditions["dx"], conditions["dy"] = dx, dy
    print(f"Data loaded for case {os.fspath(path)}")
    return {"Thot": Thot, "Tcold": Tcold, "HFS": hfs, "Conditions": conditions,
            "f_acq": 1.0 / dt, "dx": dx, "dy": dy}


def save_results(path, **arrays):
    """Save arrays to a compressed ``.npz`` file (``np.savez_compressed``)."""
    np.savez_compressed(path, **arrays)
