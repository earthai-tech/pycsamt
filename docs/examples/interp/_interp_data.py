"""Shared synthetic model for the ``pycsamt.interp`` example gallery.

Not an example itself: sphinx-gallery skips files whose name starts with
``_`` (see ``sphinx_gallery_conf["ignore_pattern"]``). Every ``plot_*.py``
script in this folder imports the same synthetic resistivity section from
here, so all interpretation deliverables trace back to one model.

The section is a classic hydrogeological sequence — a resistive dry
overburden over a conductive sand aquifer, a very conductive clay aquitard,
and a resistive basement — with gentle lateral structure. ``wet`` lowers the
aquifer resistivity (more water) for time-lapse examples.
"""

from __future__ import annotations

import numpy as np

from pycsamt.interp import ResistivityModel
from pycsamt.interp.borehole import Borehole, Interval


def demo_model(wet: float = 1.0, seed: int = 0) -> ResistivityModel:
    """Return the shared synthetic 2-D resistivity section.

    Parameters
    ----------
    wet : float, default 1.0
        Aquifer wetness factor; >1 lowers the aquifer resistivity (more
        pore water) for the time-lapse examples.
    seed : int, default 0
        Seed for the small cell-to-cell resistivity scatter.
    """
    rng = np.random.default_rng(seed)
    x = np.linspace(0.0, 2000.0, 44)          # profile distance, m
    z = np.linspace(5.0, 300.0, 34)           # depth, m
    ZZ, XX = np.meshgrid(z, x, indexing="ij")

    rho = np.full_like(ZZ, 200.0)             # dry overburden
    aquifer_top = 20.0 + 8.0 * np.sin(XX / 600.0)
    aquifer_bot = 80.0 + 15.0 * np.sin(XX / 500.0 + 1.0)
    clay_bot = 120.0 + 10.0 * np.cos(XX / 700.0)
    rho[(ZZ >= aquifer_top) & (ZZ < aquifer_bot)] = 40.0 / wet   # sand aquifer
    rho[(ZZ >= aquifer_bot) & (ZZ < clay_bot)] = 12.0            # clay aquitard
    rho[ZZ >= clay_bot] = 1500.0                                 # basement
    rho *= 1.0 + 0.05 * rng.standard_normal(rho.shape)

    return ResistivityModel.from_array(
        np.log10(np.clip(rho, 1.0, 1e5)),
        x_centers=x, z_centers=z,
        station_x=x[::4],
        station_names=[f"S{i:02d}" for i in range(len(x[::4]))],
        method="synthetic",
    )


def demo_boreholes() -> list:
    """Two ground-truth boreholes matching the synthetic sequence."""
    seq = [
        Interval(top=0.0, bottom=20.0, lithology="overburden"),
        Interval(top=20.0, bottom=80.0, lithology="sand aquifer"),
        Interval(top=80.0, bottom=120.0, lithology="clay"),
        Interval(top=120.0, bottom=300.0, lithology="granite basement"),
    ]
    return [
        Borehole("BH-1", intervals=seq, x=500.0),
        Borehole("BH-2", intervals=seq, x=1500.0),
    ]
