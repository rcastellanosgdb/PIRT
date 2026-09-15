#!/usr/bin/env python
"""Impinging sweeping jet example (Python counterpart of ``SJ_main.m``).

Processes the sample dataset of the sweeping-jet experiment of Robledo et al.,
Exp. Therm. Fluid Sci. 169 (2025) 111526, with the filtering pipeline selected
in that work (POD + rectangular Gaussian + Savitzky-Golay) and computes the
Nusselt number including the unsteady and tangential-conduction terms.

Data
----
Download the test case from the link given in the repository README and place
``Thot.mat``, ``Tcold.mat``, ``Resolution.mat`` and ``TestConditions.mat`` next
to ``CONF_DATA.STR`` in a folder (default ``./SJ_processing/``), or pass the
folder as first argument.

    python sj_main.py [path/to/SJ_processing]
"""
import sys
import time
from pathlib import Path

import numpy as np

from pirt import PIRT, io

path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent / "SJ_processing"

case = io.load_sj_case(path)
Thot, Tcold = case["Thot"], case["Tcold"]
HFS, Conditions = case["HFS"], case["Conditions"]
dx, dy, dt = case["dx"], case["dy"], 1.0 / case["f_acq"]

for name, val in HFS.items():
    print(f"HFS.{name} = {val}")
for name, val in Conditions.items():
    print(f"Conditions.{name} = {val}")

# --- POD + Gaussian + Savitzky-Golay ---------------------------------------
filters = [
    {"type": "POD", "criterion": "HardThreshold", "beta": 2000 / (399 * 604)},
    {"type": "gaussian", "filter_size": (9, 3, 1), "sigma": (3, 3, 0.1)},
    {"type": "sgolay32", "kernel_size": (5, 5, 3), "h": (dx, dy, dt)},
]

t0 = time.time()
output = PIRT(Thot=Thot, Tcold=Tcold,                 # temperature maps
              heat_transfer=("Nu", "h"),              # quantities to compute
              time_der=True, spatial_der=True,        # terms of the energy balance
              HFS=HFS, conditions=Conditions,         # sensor and test data
              filters=filters)                        # filtering pipeline
output.go()
print(f"Elapsed time: {time.time() - t0:.1f} s")

# --- Extract results ----------------------------------------------------------
Nu = output.result["Nu"]
Num = Nu.mean(axis=2)                       # time-averaged Nusselt number
Nuf = np.abs(Nu - Num[:, :, None]).mean(axis=2)  # mean absolute fluctuation Nu'
print(f"Nu: mean {Num.mean():.2f}, max {Num.max():.2f}; Nu': mean {Nuf.mean():.2f}, max {Nuf.max():.2f}")

# --- Plot ---------------------------------------------------------------------
try:
    import matplotlib.pyplot as plt
except ImportError:  # plotting is optional
    sys.exit(0)

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
for ax, field, title in zip(axes, (Num, Nuf), ("Mean Nu map", "Fluctuating Nu map")):
    im = ax.imshow(field, origin="upper", cmap="viridis")
    ax.set_title(title)
    ax.set_axis_off()
    fig.colorbar(im, ax=ax, shrink=0.8)
fig.tight_layout()
plt.show()
