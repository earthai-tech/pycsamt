# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Shared time-series sample for the IoT gallery examples.

The real KAP103 recording (``data/MT/TS/kap103as.ts``, ~27 MB) is not
bundled with the repository, so on a clean or CI checkout the IoT
examples that need raw magnetotelluric time series fall back to a small
**synthetic** record with the same channels, sampling interval and
station metadata. The pipeline being demonstrated is identical — only the
numbers differ — and each example notes when the fallback is in use.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np

from pycsamt.ts import TSData, read_ts

STATION = "kap103"
DT = 5.0  # seconds (0.2 Hz), matching the real record
CHANNELS = ["HX", "HY", "HZ", "EX", "EY"]
LAT, LON = -32.1388893, 20.4675007


def repo_root() -> Path:
    """Repository root during docs builds and local runs."""
    env_root = os.environ.get("PYCSAMT_DOCS_REPO_ROOT")
    if env_root:
        return Path(env_root)
    return Path(__file__).resolve().parents[3]


def real_ts_path() -> Path:
    return repo_root() / "data" / "MT" / "TS" / "kap103as.ts" / "kap103as.ts"


def _synthetic(n: int = 65536, seed: int = 103) -> TSData:
    """A small MT-plausible five-channel record.

    Magnetic channels are red (Brownian) noise; the electric channels are
    a resistive transfer function of the magnetic field plus noise, so the
    impedance recovered downstream by ``ts_to_z`` is well-defined rather
    than degenerate. Enough samples for several ``nfft=8192`` windows.
    """
    rng = np.random.default_rng(seed)

    def red(scale: float) -> np.ndarray:
        x = np.cumsum(rng.standard_normal(n))
        x -= x.mean()
        x /= x.std() + 1e-12
        return (x * scale).astype(float)

    hx, hy = red(2.0), red(2.0)
    hz = 0.15 * hx + red(0.3)
    ex = 40.0 * hy + red(1.5)  # Zxy-like
    ey = -40.0 * hx + red(1.5)  # Zyx-like
    data = {"HX": hx, "HY": hy, "HZ": hz, "EX": ex, "EY": ey}
    return TSData(
        data=data,
        dt=DT,
        station=f"{STATION}-synthetic",
        lat=LAT,
        lon=LON,
    )


def load_ts_sample() -> tuple[TSData, bool]:
    """Return ``(record, is_real)``.

    The real KAP103 record when it is present on disk, otherwise a
    synthetic mirror so the example still runs in a clean checkout.
    """
    path = real_ts_path()
    if path.is_file():
        return read_ts(str(path)), True
    return _synthetic(), False
