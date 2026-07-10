# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Shared fixtures for pycsamt.app tests.

Qt is initialised once per session with the offscreen platform so
tests can instantiate Qt widgets without a real display.
"""

from __future__ import annotations

import os

import pytest

# ── Qt / offscreen setup ───────────────────────────────────────────────────

@pytest.fixture(scope="session", autouse=True)
def qt_offscreen():
    """Force Qt to use the offscreen (headless) platform for all tests."""
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    yield


@pytest.fixture(scope="session")
def qapp():
    """Single QApplication shared across the whole test session."""
    from PySide6.QtWidgets import QApplication

    existing = QApplication.instance()
    if existing is not None:
        yield existing
        return

    app = QApplication(["pytest", "-platform", "offscreen"])
    app.setApplicationName("pycsamt-test")
    yield app


# ── EDI fixtures ───────────────────────────────────────────────────────────

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
        ">MTSECT",
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
    for comp in ("ZXXR", "ZXXI", "ZXYR", "ZXYI", "ZYXR", "ZYXI", "ZYYR", "ZYYI"):
        lines.append(f">{comp} // {nfreq}")
        lines.append("  " + "  ".join(f"{v:.6E}" for v in z_real))
    lines.append(">END")

    edi_path.write_text("\n".join(lines), encoding="utf-8")
    return edi_path
