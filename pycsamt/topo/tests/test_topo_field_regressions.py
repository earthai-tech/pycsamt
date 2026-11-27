# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Regressions from a field survey: station order and markers."""

from __future__ import annotations

import glob
import os
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pytest  # noqa: E402
from matplotlib.collections import PathCollection  # noqa: E402

from pycsamt.seg.edi import EDIFile  # noqa: E402
from pycsamt.topo import plot_topo_section  # noqa: E402
from pycsamt.topo.extract import extract_chainage  # noqa: E402
from pycsamt.topo.section import _match_by_name  # noqa: E402

_L18 = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                    "data", "AMT", "WILLY_DATA", "L18PLT")


def _l18():
    paths = sorted(glob.glob(os.path.join(_L18, "*.edi")))
    if len(paths) < 4:
        pytest.skip("bundled L18 EDI files not found")
    return [EDIFile(p) for p in paths]


def _grid(n=6):
    x = np.linspace(0.0, 1.0, n)
    z = np.linspace(0.05, 1.0, 10)
    return x, z, np.full((10, n), 2.0)


def test_show_stations_false_draws_no_markers():
    x, z, rho = _grid()
    ax = plot_topo_section((x, z, rho), elevation=np.full(6, 0.3),
                           chainage=x, show_stations=False)
    pts = [c for c in ax.collections if isinstance(c, PathCollection)]
    assert all(len(c.get_offsets()) == 0 for c in pts)


def test_default_markers_follow_inversion_preset():
    from pycsamt.api.station import PYCSAMT_STATION_RENDERING

    x, z, rho = _grid()
    ax = plot_topo_section((x, z, rho), elevation=np.full(6, 0.3),
                           chainage=x)
    pts = [c for c in ax.collections if isinstance(c, PathCollection)
           and len(c.get_offsets())]
    face = matplotlib.colors.to_rgba(
        PYCSAMT_STATION_RENDERING.inversion.marker.facecolor)
    assert tuple(pts[-1].get_facecolor()[0]) == pytest.approx(face)


def test_match_by_name_aligns_unsorted_sites():
    grid = SimpleNamespace(station_names=["L1-0", "L1-30", "L1-150"],
                           station_x=np.array([0.0, 30.0, 150.0]),
                           unit="m")
    # file-name order and DATAID extensions, as from a folder
    names = ["L1-0.atimpd", "L1-150.atimpd", "L1-30.atimpd"]
    elev = np.array([374.0, 326.0, 367.0])
    x_km, el = _match_by_name(grid, names, elev)
    np.testing.assert_allclose(x_km, [0.0, 0.03, 0.15])
    np.testing.assert_allclose(el, [374.0, 367.0, 326.0])
    grid.station_names = ["L1-0", "L1-999", "L1-150"]
    assert _match_by_name(grid, names, elev) is None


def test_chainage_auto_is_order_independent():
    sites = _l18()
    ref = extract_chainage(sites, order="auto")
    perm = np.random.default_rng(3).permutation(len(sites))
    shuffled = extract_chainage([sites[i] for i in perm],
                                order="auto")
    np.testing.assert_allclose(shuffled, ref[perm], atol=1e-9)


def test_chainage_order_validation():
    with pytest.raises(ValueError, match="order"):
        extract_chainage([], order="sorted")
