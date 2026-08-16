# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.interp._base — ResistivityModel."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.interp._base import ResistivityModel


def _model():
    x = np.array([0.0, 250.0, 500.0])
    z = np.array([5.0, 15.0, 30.0])
    rho = np.log10(
        np.array(
            [[100.0, 120.0, 90.0], [50.0, 60.0, 40.0], [200.0, 210.0, 190.0]]
        )
    )
    return ResistivityModel.from_array(
        rho, x, z, station_names=["S0", "S1", "S2"], method="demo"
    )


def test_from_array_shapes():
    m = _model()
    assert m.n_x == 3
    assert m.n_z == 3
    assert m.station_x.tolist() == m.x_centers.tolist()


def test_from_array_default_station_names():
    x = np.array([0.0, 1.0])
    z = np.array([1.0, 2.0])
    rho = np.zeros((2, 2))
    m = ResistivityModel.from_array(rho, x, z)
    assert m.station_names == ["S000", "S001"]


def test_column_access():
    m = _model()
    col = m.station_column("S1")
    assert col.shape == (3,)
    with pytest.raises(KeyError):
        m.station_column("nope")


def test_repr_summarizes_arrays_instead_of_dumping_them():
    # A dataclass's default __repr__ would print every element of rho_2d;
    # PyCSAMTObject's summarized repr must be used instead (repr=False).
    m = _model()
    text = repr(m)
    assert "ndarray(shape=" in text
    assert "100.0" not in text  # no raw array contents leaked into repr


def test_metadata_mixin_is_lazy_and_mutable():
    m = _model()
    assert m.metadata == {}
    m.update_metadata(crs="EPSG:32633")
    assert m.metadata == {"crs": "EPSG:32633"}


def test_clone_creates_independent_copy():
    m = _model()
    m2 = m.clone(method="cloned")
    assert m2.method == "cloned"
    assert m.method == "demo"


def test_clip_to_stations_drops_padding_keeps_core():
    x = np.linspace(-5000.0, 5000.0, 100)
    z = np.array([10.0, 20.0])
    rho = np.zeros((2, 100))
    m = ResistivityModel.from_array(
        rho, x, z, station_x=np.array([0.0, 500.0, 1000.0]),
        station_names=["S0", "S1", "S2"],
    )
    clipped = m.clip_to_stations(pad_m=100.0)
    assert clipped.x_centers.min() >= -100.0
    assert clipped.x_centers.max() <= 1100.0
    assert clipped.n_x < m.n_x
    # station metadata untouched
    assert clipped.station_x.tolist() == m.station_x.tolist()
    assert clipped.station_names == m.station_names
    assert clipped.method == m.method


def test_clip_to_stations_requires_station_x():
    x = np.array([0.0, 1.0])
    z = np.array([1.0])
    rho = np.zeros((1, 2))
    m = ResistivityModel.from_array(rho, x, z, station_x=np.array([]))
    with pytest.raises(ValueError):
        m.clip_to_stations()


def test_clip_to_depth_drops_deep_cells_keeps_shallow():
    x = np.array([0.0, 100.0])
    z = np.array([10.0, 500.0, 9000.0])
    rho = np.arange(6, dtype=float).reshape(3, 2)
    m = ResistivityModel.from_array(
        rho, x, z, station_x=np.array([0.0]), station_names=["S0"],
    )
    clipped = m.clip_to_depth(1000.0)
    assert clipped.z_centers.tolist() == [10.0, 500.0]
    assert clipped.rho_2d.shape == (2, 2)
    np.testing.assert_array_equal(clipped.rho_2d, rho[:2, :])
    # x/station metadata untouched
    assert clipped.x_centers.tolist() == m.x_centers.tolist()
    assert clipped.station_x.tolist() == m.station_x.tolist()


def test_clip_to_depth_raises_when_nothing_survives():
    x = np.array([0.0])
    z = np.array([500.0, 9000.0])
    rho = np.zeros((2, 1))
    m = ResistivityModel.from_array(rho, x, z)
    with pytest.raises(ValueError):
        m.clip_to_depth(10.0)


class _FakeMesh:
    """Minimal OccamMesh-shaped object: 2 left-pad + 4 core + 2 right-pad
    cells, widths [200, 100 | 50, 50, 50, 50 | 100, 200], matching the
    left/right-symmetric padding pycsamt.models.occam2d.OccamMesh.from_data
    actually builds."""

    def __init__(self):
        x_widths = np.array([200.0, 100.0, 50.0, 50.0, 50.0, 50.0, 100.0, 200.0])
        self.x_nodes = np.concatenate([[0.0], np.cumsum(x_widths)])
        z_widths = np.array([40.0, 40.0])
        self.z_nodes = np.concatenate([[0.0], np.cumsum(z_widths)])


class _FakeData:
    def __init__(self):
        # station chainage in the *data* frame -- starts at 0, not at the
        # mesh's own x=0 (which is the outer edge of the left padding).
        self.offsets = np.array([0.0, 50.0, 100.0, 150.0])
        self.sites = ["S0", "S1", "S2", "S3"]


class _FakeOccamResult:
    def __init__(self):
        self.mesh = _FakeMesh()
        self.data = _FakeData()
        self.rho_2d = np.zeros((2, 8))
        self.final_rms = 1.23


def test_from_occam2d_aligns_x_centers_to_station_offsets():
    # Regression test: x_centers previously stayed in the mesh's own
    # coordinate frame (origin at the outer edge of the padding) while
    # station_x was real chainage (origin at the first station), so every
    # station's nearest x_centers column silently collapsed onto the same
    # outermost padding cell. Padding is built symmetrically on both sides,
    # so the station-carrying core is centred in the mesh; from_occam2d
    # must recover and undo that offset.
    model = ResistivityModel.from_occam2d(_FakeOccamResult())
    for x_station in model.station_x:
        ix = int(np.argmin(np.abs(model.x_centers - x_station)))
        # nearest core cell should be within half a core cell width (25 m)
        assert abs(model.x_centers[ix] - x_station) <= 25.0
    # and every station must resolve to a *different* column, not all the
    # same one
    nearest = {
        int(np.argmin(np.abs(model.x_centers - x))) for x in model.station_x
    }
    assert len(nearest) == len(model.station_x)


# ---------------------------------------------------------------------------
# from_any — the "accept any 2-D inversion result" adapter
# ---------------------------------------------------------------------------


def test_from_any_passes_through_existing_model():
    m = _model()
    assert ResistivityModel.from_any(m) is m


def test_from_any_raw_array_triple():
    x = np.array([0.0, 100.0, 200.0])
    z = np.array([10.0, 50.0])
    rho_log10 = np.array([[2.0, 2.1, 2.2], [2.5, 2.6, 2.7]])
    m = ResistivityModel.from_any((x, z, rho_log10))
    assert isinstance(m, ResistivityModel)
    assert m.method == "array"
    assert m.n_x == 3 and m.n_z == 2
    np.testing.assert_array_equal(m.rho_2d, rho_log10)
    np.testing.assert_array_equal(m.x_centers, x)


def test_from_any_native_occam2d_result():
    m = ResistivityModel.from_any(_FakeOccamResult())
    assert m.method == "occam2d"
    assert m.rms == pytest.approx(1.23)
    assert m.n_z == 2 and m.n_x == 8


def test_from_any_station_overrides_applied():
    x = np.array([0.0, 100.0])
    z = np.array([10.0])
    rho_log10 = np.array([[2.0, 2.1]])
    m = ResistivityModel.from_any(
        (x, z, rho_log10),
        station_x=np.array([5.0, 95.0]),
        station_names=["A", "B"],
    )
    np.testing.assert_array_equal(m.station_x, [5.0, 95.0])
    assert m.station_names == ["A", "B"]


def test_from_any_km_unit_converted_to_metres():
    x_km = np.array([0.0, 0.1, 0.2])
    z_km = np.array([0.01, 0.05])
    rho_log10 = np.array([[2.0, 2.1, 2.2], [2.5, 2.6, 2.7]])
    m = ResistivityModel.from_any((x_km, z_km, rho_log10), unit="km")
    np.testing.assert_allclose(m.x_centers, [0.0, 100.0, 200.0])
    np.testing.assert_allclose(m.z_centers, [10.0, 50.0])


def test_from_any_rejects_unsupported_type():
    with pytest.raises(TypeError):
        ResistivityModel.from_any(object())
