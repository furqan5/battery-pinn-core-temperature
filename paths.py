"""Repository and dataset locations, resolved rather than hardcoded.

Added for the public release. The scripts originally carried absolute paths from
the machine they were developed on, which made them unrunnable elsewhere.

Layout expected by default:

    <repo>/                 this file
    <repo>/Data_Sets/       datasets, unpacked as described in DATA.md (gitignored)

Override the dataset location with the ``BATTERY_PINN_DATA`` environment
variable if the datasets live elsewhere -- point it at the directory that
*contains* ``kxsbr4x3j2-2``, ``Li-ion Battery Dataset from NASA PCoE`` and
``EKF-Battery-Impedance-Temperature-master``.

    set BATTERY_PINN_DATA=D:\\datasets          (Windows)
    export BATTERY_PINN_DATA=/data/battery      (POSIX)

No measurement data is redistributed with this repository. See DATA.md.
"""

from __future__ import annotations

import os

#: Repository root -- the directory holding this file.
REPO = os.path.dirname(os.path.abspath(__file__))

#: Directory containing the unpacked datasets.
DATA_ROOT = os.environ.get("BATTERY_PINN_DATA", os.path.join(REPO, "Data_Sets"))

#: NASA Ames Prognostics battery ageing set (Saha & Goebel 2007).
NASA_ROOT = os.path.join(DATA_ROOT, "Li-ion Battery Dataset from NASA PCoE")

#: Catenaro & Onori galvanostatic discharge set (doi:10.17632/kxsbr4x3j2.1).
MULTIRATE_ROOT = os.path.join(
    DATA_ROOT, "kxsbr4x3j2-2", "galvanostatic_discharge_test"
)

#: Richardson & Howey A123 26650 records, distributed with their EKF code.
RICHARDSON_ROOT = os.path.join(
    DATA_ROOT,
    "EKF-Battery-Impedance-Temperature-master",
    "EKF-Battery-Impedance-Temperature-master",
)

#: Prior-work notebook package, if present (not tracked in this branch).
COREFIELD_ROOT = os.path.join(
    REPO, "CoreField_Battery_Package_v1.0", "CoreField_Battery_Package_v1.0"
)


def require(path: str, what: str) -> str:
    """Return ``path``, or raise with a pointer to DATA.md if it is missing.

    Failing early with an actionable message beats a FileNotFoundError from
    somewhere three calls deeper.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{what} not found at:\n  {path}\n\n"
            "No measurement data ships with this repository. See DATA.md for the "
            "source and the expected layout, or set BATTERY_PINN_DATA to the "
            "directory containing the datasets."
        )
    return path
