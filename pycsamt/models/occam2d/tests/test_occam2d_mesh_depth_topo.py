"""Tests for OccamMesh's max_depth cap and OccamData's topo-package
elevation extraction.
"""

import numpy as np
import pytest

from pycsamt.models.occam2d.config import OccamConfig
from pycsamt.models.occam2d.data import OccamData
from pycsamt.models.occam2d.mesh import OccamMesh


class _FakeSite:
    """Minimal duck-typed site accepted by OccamData.from_edi, optionally
    carrying a direct ``.elev`` attribute picked up by
    :func:`pycsamt.topo.extract_elevation`.
    """

    def __init__(self, name, lat, lon, freqs, rho, phs, elev=None):
        self.name = name
        self.coords = (lat, lon, 0.0)
        n = len(freqs)
        self.freq = np.asarray(freqs, dtype=float)
        rho_arr = np.zeros((n, 2, 2), dtype=float)
        rho_arr[:, 0, 1] = rho
        rho_arr[:, 1, 0] = rho
        self.rho = rho_arr
        phs_arr = np.zeros((n, 2, 2), dtype=float)
        phs_arr[:, 0, 1] = phs
        phs_arr[:, 1, 0] = phs - 180.0
        self.phase = phs_arr
        if elev is not None:
            self.elev = float(elev)


def _make_sites(n_sites=4, n_freqs=8, elevs=None):
    freqs = np.logspace(3, -1, n_freqs)
    sites = []
    for i in range(n_sites):
        rho = np.full(n_freqs, 100.0)
        phs = np.full(n_freqs, 45.0)
        elev = elevs[i] if elevs is not None else None
        sites.append(
            _FakeSite(
                name=f"S{i:02d}",
                lat=0.0,
                lon=0.0 + i * 0.01,
                freqs=freqs,
                rho=rho,
                phs=phs,
                elev=elev,
            )
        )
    return sites


def _simple_data(n_sites=4):
    return OccamData.from_edi(_make_sites(n_sites))


# ---------------------------------------------------------------------------
# OccamConfig.max_depth / OccamMesh.from_data depth cap
# ---------------------------------------------------------------------------


def test_default_max_depth_is_1500m():
    assert OccamConfig().max_depth == pytest.approx(1500.0)


def test_mesh_from_data_default_reaches_max_depth_exactly():
    d = _simple_data()
    m = OccamMesh.from_data(d)
    assert m.z_nodes[-1] == pytest.approx(1500.0)
    # fewer layers than the old unbounded n_layers=30 default
    assert m.n_zcells < OccamConfig().n_layers


def test_mesh_from_data_custom_max_depth():
    d = _simple_data()
    cfg = OccamConfig(max_depth=800.0)
    m = OccamMesh.from_data(d, config=cfg)
    assert m.z_nodes[-1] == pytest.approx(800.0)


def test_mesh_from_data_n_layers_still_caps_when_shallower_than_max_depth():
    d = _simple_data()
    cfg = OccamConfig(n_layers=5, max_depth=1500.0)
    m = OccamMesh.from_data(d, config=cfg)
    # geometric series with defaults stays well under 1500 m at 5 layers,
    # so n_layers (not max_depth) is the binding constraint here.
    assert m.n_zcells == 5
    assert m.z_nodes[-1] < 1500.0


def test_mesh_from_data_zwidths_positive_after_truncation():
    d = _simple_data()
    m = OccamMesh.from_data(d)
    assert np.all(m.z_widths > 0)


def test_mesh_from_data_zero_max_depth_disables_cap():
    d = _simple_data()
    cfg = OccamConfig(n_layers=6, max_depth=0.0)
    m = OccamMesh.from_data(d, config=cfg)
    assert m.n_zcells == 6


# ---------------------------------------------------------------------------
# OccamData.from_edi topography extraction (pycsamt.topo)
# ---------------------------------------------------------------------------


def test_from_edi_no_elevation_source_has_no_topography():
    d = _simple_data()
    assert d.elevations.shape == (4,)
    assert not d.has_topography
    assert d.station_elevations() == {}


def test_from_edi_extracts_real_elevation():
    elevs = [1200.0, 1195.0, 1190.0, 1185.0]
    sites = _make_sites(4, elevs=elevs)
    d = OccamData.from_edi(sites)
    assert d.has_topography
    # sites are already west-to-east in ascending order, matching offsets
    assert d.elevations.tolist() == pytest.approx(elevs)
    assert d.station_elevations() == {
        name: elev for name, elev in zip(d.sites, elevs)
    }


def test_from_edi_elevation_follows_chainage_sort():
    # site order given east-to-west; offsets get sorted ascending, and
    # elevation must be reordered along with names/offsets.
    sites = _make_sites(4, elevs=[10.0, 20.0, 30.0, 40.0])
    sites = list(reversed(sites))  # now S03..S00, west-to-east becomes east-to-west
    d = OccamData.from_edi(sites)
    assert d.sites == ["S00", "S01", "S02", "S03"]
    assert d.elevations.tolist() == pytest.approx([10.0, 20.0, 30.0, 40.0])
