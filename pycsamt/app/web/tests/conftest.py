# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Shared fixtures for pycsamt.app.web tests.

This directory is a sibling of ``pycsamt/app/desktop/tests`` and
``pycsamt/app/mapview/tests`` (previously all three lived together, flat,
under ``pycsamt/app/tests``), so pytest does not pick up any conftest.py
from those directories automatically. Unlike desktop, this suite is pure
Dash -- no test here imports PySide6 -- so none of the Qt offscreen/qapp/
Shiboken-teardown fixtures are needed or duplicated; only the Dash app
fixture and the EDI fixture some callback tests use.
"""

from __future__ import annotations

import os

import pytest

# pycsamt.app.web.callbacks.inversion pulls in torch at module scope, and a
# "Traditional" inversion test calls scipy.optimize.least_squares for real
# in the same process. On Windows/conda that combination can abort the
# whole interpreter ("OMP: Error #15: Initializing libiomp5md.dll, but
# found ... already initialized") -- duplicate OpenMP runtime registration,
# not a bug in either library. Must be set before either one is imported.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")


@pytest.fixture(scope="session")
def web_app():
    """Single fully-wired Dash app (all callbacks registered) for reuse."""
    from pycsamt.app.web.app import create_app

    return create_app()


@pytest.fixture(scope="session")
def simulated_edi(tmp_path_factory):
    """Return a path to a minimal synthetic EDI file."""
    import numpy as np

    tmp = tmp_path_factory.mktemp("edi")
    edi_path = tmp / "SIM001.edi"

    # Build the minimal >INFO / >HEAD / >FREQ / >ZXXR … blocks
    nfreq = 8
    freqs = np.logspace(2, -1, nfreq)
    z_real = np.ones(nfreq) * 10.0
    np.ones(nfreq) * 8.0

    lines = [
        ">HEAD",
        " DATAID=SIM001",
        " LAT=48:30:0.0",
        " LONG=7:45:0.0",
        " ELEV=200.0",
        ">INFO",
        " MAXINFO=999",
        ">DEFINEMEAS",
        " MAXCHAN=7",
        " MAXRUN=999",
        " MAXMEAS=9999",
        " UNITS=M",
        " REFTYPE=CART",
        " REFLAT=48:30:0.0",
        " REFLONG=7:45:0.0",
        " REFELEV=200.0",
        ">=MTSECT",
        " SECTID=SIM001",
        f" NFREQ={nfreq}",
        " HX=1001.001",
        " HY=1002.001",
        " HZ=1003.001",
        " EX=1004.001",
        " EY=1005.001",
        f">FREQ // {nfreq}",
    ]
    lines.append("  " + "  ".join(f"{f:.6E}" for f in freqs))
    for comp in (
        "ZXXR",
        "ZXXI",
        "ZXYR",
        "ZXYI",
        "ZYXR",
        "ZYXI",
        "ZYYR",
        "ZYYI",
    ):
        lines.append(f">{comp} // {nfreq}")
        lines.append("  " + "  ".join(f"{v:.6E}" for v in z_real))
    lines.append(">END")

    edi_path.write_text("\n".join(lines), encoding="utf-8")
    return edi_path
