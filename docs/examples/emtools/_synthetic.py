"""Synthetic site generator for ``pycsamt.emtools`` gallery scripts that
need site-like data but have no suitable bundled real dataset (or want to
show an idealized, noise-free response alongside real data).

Not an example itself: see the note at the top of ``_datasets.py``.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np


class SyntheticTipper:
    """``.tipper`` is (n_freq, 2) complex — matches ``pycsamt``'s Tipper."""

    def __init__(self, tipper: np.ndarray, freq: np.ndarray):
        self.tipper = tipper
        self.freq = freq


class SyntheticImpedance:
    """A flat, isotropic Z — just enough for ``ensure_sites`` to recognize
    the object as a genuine site (it looks for a Z/Tipper pair)."""

    def __init__(self, freq: np.ndarray, rho: float = 100.0):
        z = np.zeros((freq.size, 2, 2), dtype=complex)
        amp = np.sqrt(5.0 * freq * rho)
        z[:, 0, 1] = amp * (1 + 1j) / np.sqrt(2)
        z[:, 1, 0] = -amp * (1 + 1j) / np.sqrt(2)
        self.z = z
        self.freq = freq


class SyntheticSite:
    """Minimal site-like object understood by ``pycsamt.emtools``."""

    def __init__(
        self,
        station: str,
        east: float,
        north: float,
        freq: np.ndarray,
        tipper: np.ndarray,
        *,
        rho: float = 100.0,
    ):
        self.station = station
        self.east = float(east)
        self.north = float(north)
        self.freq = freq
        self.Tipper = SyntheticTipper(tipper, freq)
        self.Z = SyntheticImpedance(freq, rho=rho)

    def get_section(self, *_, **__):
        return None


def synthetic_line(
    n: int = 12, spacing: float = 200.0
) -> tuple[list, np.ndarray, np.ndarray, np.ndarray]:
    """A fully made-up straight station line: names, east, north, freq."""
    names = [f"S{i:02d}" for i in range(n)]
    east = np.arange(n) * float(spacing)
    north = np.zeros(n)
    freq = np.logspace(1, 4, 32)
    return names, east, north, freq


def conductor_tipper_sites(
    names: Sequence[str],
    east: np.ndarray,
    north: np.ndarray,
    freq: np.ndarray,
    *,
    conductor_xy: tuple[float, float] | None = None,
    seed: int = 0,
) -> list:
    """Synthetic tipper-bearing sites shaped like one buried conductor.

    Induction vectors point toward ``conductor_xy`` (profile midpoint by
    default), growing with period and decaying with distance — a standard
    qualitative teaching pattern for induction-vector diagnostics. Station
    geometry (``east``/``north``/``freq``) can come from a real survey
    (see ``_datasets.py``) or from :func:`synthetic_line`.
    """
    east = np.asarray(east, float)
    north = np.asarray(north, float)
    xy = np.column_stack([east, north])
    if conductor_xy is None:
        conductor_xy = tuple(xy.mean(axis=0))
    c_xy = np.asarray(conductor_xy, float)

    to_conductor = c_xy - xy
    dist = np.linalg.norm(to_conductor, axis=1) + 1e-6
    direction = to_conductor / dist[:, None]

    period = 1.0 / freq
    period_gain = np.clip(
        np.log10(period / period.min())
        / np.log10(period.max() / period.min()),
        0.05,
        1.0,
    )
    spatial_gain = np.exp(-((dist / dist.max()) ** 2) * 2.0)

    rng = np.random.default_rng(seed)
    sites = []
    for k, name in enumerate(names):
        mag = 0.35 * period_gain * spatial_gain[k]
        phase = np.exp(1j * 0.3 * np.pi * period_gain)
        tx = mag * direction[k, 0] * phase
        ty = mag * direction[k, 1] * phase
        noise = 0.01 * (
            rng.standard_normal((freq.size, 2))
            + 1j * rng.standard_normal((freq.size, 2))
        )
        tipper = np.column_stack([tx, ty]) + noise
        sites.append(SyntheticSite(name, east[k], north[k], freq, tipper))
    return sites
