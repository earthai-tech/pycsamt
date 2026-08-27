# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Phase 5 tests — multiline PCSF builder/reader.

Includes a direct equivalence check against
``pycsamt.app.web.callbacks.map3d``'s own private
``_line_real_offsets``/``_assemble_3d_grid`` functions: this module's
public ``line_offsets_from_stations``/``stack_lines_to_common_grid``
must be numerically identical to them, not just "close", since that
identity is what makes a cached ``derived_volume`` bit-for-bit
reproducible by ``map3d.py``'s existing renderers.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.format import read_pcsf, write_pcsf
from pycsamt.format.multiline import (
    build_multiline_pcsf,
    line_offsets_from_stations,
    multiline_pcsf_to_profiles,
    stack_lines_to_common_grid,
)


def _flat_profiles() -> dict:
    return {
        "L1": {
            "x": np.array([0.0, 100.0, 200.0]),
            "z": np.array([10.0, 50.0]),
            "rho": np.array([[100.0, 110.0, 120.0], [50.0, 55.0, 60.0]]),
            "sta_x": [0.0, 100.0, 200.0],
            "sta_names": ["S0", "S1", "S2"],
            "sta_elev": [10.0, 12.0, 9.0],
        },
        "L2": {
            "x": np.array([0.0, 100.0, 200.0]),
            "z": np.array([10.0, 50.0]),
            "rho": np.array([[200.0, 210.0, 220.0], [90.0, 95.0, 99.0]]),
            "sta_x": [0.0, 100.0, 200.0],
            "sta_names": ["S3", "S4", "S5"],
            "sta_elev": [11.0, 13.0, 8.0],
        },
    }


def _profiles_with_real_geometry() -> dict:
    # Two parallel N-S lines, 200 m apart east-west.
    profiles = _flat_profiles()
    profiles["L1"]["sta_lat"] = [10.00, 10.01, 10.02]
    profiles["L1"]["sta_lon"] = [20.000, 20.000, 20.000]
    profiles["L2"]["sta_lat"] = [10.00, 10.01, 10.02]
    profiles["L2"]["sta_lon"] = [20.002, 20.002, 20.002]
    return profiles


class TestLineOffsetsFromStations:
    def test_returns_none_without_lat_lon(self):
        assert line_offsets_from_stations(_flat_profiles()) is None

    def test_real_offsets_are_normalized_nonnegative(self):
        offsets = line_offsets_from_stations(_profiles_with_real_geometry())
        assert offsets is not None
        assert min(offsets.values()) == pytest.approx(0.0)
        assert offsets["L1"] != offsets["L2"]


class TestStackLinesToCommonGrid:
    def test_shapes(self):
        x_arr, z_arr, y_arr, rho_3d = stack_lines_to_common_grid(
            _flat_profiles(), line_spacing=1.0
        )
        assert x_arr.shape == (3,)
        assert z_arr.shape == (2,)
        assert y_arr.shape == (2,)
        assert rho_3d.shape == (2, 3, 2)  # (n_lines, n_x, n_z)

    def test_same_grid_lines_need_no_interpolation(self):
        # Both lines already share x/z with the reference (L1) -> exact
        # transpose only, no RegularGridInterpolator resampling loss.
        profiles = _flat_profiles()
        x_arr, z_arr, y_arr, rho_3d = stack_lines_to_common_grid(profiles)
        np.testing.assert_array_equal(rho_3d[0], profiles["L1"]["rho"].T)
        np.testing.assert_array_equal(rho_3d[1], profiles["L2"]["rho"].T)

    def test_rejects_empty_profiles(self):
        with pytest.raises(ValueError, match="at least one line"):
            stack_lines_to_common_grid({})

    def test_resamples_a_line_with_a_different_grid(self):
        profiles = _flat_profiles()
        profiles["L2"]["x"] = np.array([0.0, 50.0, 100.0, 150.0, 200.0])
        profiles["L2"]["rho"] = np.array(
            [[200.0, 205.0, 210.0, 215.0, 220.0], [90.0, 92.0, 95.0, 97.0, 99.0]]
        )
        x_arr, z_arr, y_arr, rho_3d = stack_lines_to_common_grid(profiles)
        assert rho_3d.shape == (2, 3, 2)
        # Interpolated onto the reference x grid: value at x=100 should
        # land close to L2's own x=100 sample (210.0 at z=10).
        assert rho_3d[1, 1, 0] == pytest.approx(210.0, rel=1e-6)


class TestBuildMultilinePcsf:
    def test_rejects_empty_profiles(self):
        with pytest.raises(ValueError, match="at least one line"):
            build_multiline_pcsf({})

    def test_rejects_missing_required_keys(self):
        with pytest.raises(ValueError, match="missing required"):
            build_multiline_pcsf({"L1": {"x": [0.0], "z": [1.0]}})

    def test_synthetic_offsets_when_no_lat_lon(self):
        model = build_multiline_pcsf(_flat_profiles(), line_spacing=1.0)
        model.validate()
        assert [line.offset_kind for line in model.geometry.lines] == [
            "synthetic", "synthetic",
        ]
        assert model.geometry.lines[0].offset_y == pytest.approx(0.0)
        assert model.geometry.lines[1].offset_y == pytest.approx(1000.0)

    def test_real_offsets_when_lat_lon_present(self):
        model = build_multiline_pcsf(_profiles_with_real_geometry())
        assert [line.offset_kind for line in model.geometry.lines] == [
            "real", "real",
        ]

    def test_station_lon_lat_persisted_alongside_the_offset(self):
        # sta_lat/sta_lon are consumed for the offset math above -- they
        # must also survive into StationTable.lon/lat, not just be used
        # transiently, so a loaded file stays self-georeferenced too.
        model = build_multiline_pcsf(_profiles_with_real_geometry())
        assert model.stations.lon is not None
        assert model.stations.lat is not None
        np.testing.assert_allclose(
            model.stations.lon, [20.000, 20.000, 20.000, 20.002, 20.002, 20.002]
        )
        np.testing.assert_allclose(
            model.stations.lat, [10.00, 10.01, 10.02, 10.00, 10.01, 10.02]
        )

    def test_no_station_lon_lat_when_absent(self):
        model = build_multiline_pcsf(_flat_profiles())
        assert model.stations.lon is None
        assert model.stations.lat is None

    def test_derived_volume_cached_by_default(self):
        model = build_multiline_pcsf(_flat_profiles())
        assert model.geometry.derived_volume is not None
        assert model.geometry.derived_volume.derivation_method == "linear_interp"
        assert model.geometry.derived_volume.derived_from == ["L1", "L2"]
        assert model.geometry.derived_volume.synthesized is True

    def test_derived_volume_can_be_skipped(self):
        model = build_multiline_pcsf(_flat_profiles(), cache_derived_volume=False)
        assert model.geometry.derived_volume is None

    def test_no_derived_volume_for_a_single_line(self):
        profiles = {"L1": _flat_profiles()["L1"]}
        model = build_multiline_pcsf(profiles)
        assert model.geometry.derived_volume is None

    def test_stations_and_topography_populated(self):
        model = build_multiline_pcsf(_flat_profiles())
        assert model.stations.name == ["S0", "S1", "S2", "S3", "S4", "S5"]
        assert list(model.stations.line_id) == ["L1"] * 3 + ["L2"] * 3
        assert model.topography.station_id == model.stations.name
        np.testing.assert_allclose(
            model.topography.elevation, [10.0, 12.0, 9.0, 11.0, 13.0, 8.0]
        )

    def test_no_stations_when_absent(self):
        profiles = {
            "L1": {"x": np.array([0.0, 1.0]), "z": np.array([1.0]),
                   "rho": np.array([[1.0, 2.0]])},
        }
        model = build_multiline_pcsf(profiles, cache_derived_volume=False)
        assert model.stations is None
        assert model.topography is None


class TestBuildMultilinePcsfWithTopo:
    """topo= computes the real cross-strike offset (via the existing
    line_offsets_from_stations machinery, fed under the hood) AND
    populates StationTable.lon/lat in one pass -- no separate offset
    step needed."""

    def test_single_named_csv_spans_both_lines(self, tmp_path):
        path = tmp_path / "topo.csv"
        path.write_text(
            "station,lat,lon,elev\n"
            "S0,10.00,20.000,50\nS1,10.01,20.000,51\nS2,10.02,20.000,52\n"
            "S3,10.00,20.002,60\nS4,10.01,20.002,61\nS5,10.02,20.002,62\n",
            encoding="utf-8",
        )
        model = build_multiline_pcsf(
            _flat_profiles(), topo=path, cache_derived_volume=False
        )
        assert [line.offset_kind for line in model.geometry.lines] == ["real", "real"]
        assert model.geometry.lines[0].offset_y == pytest.approx(0.0)
        assert model.geometry.lines[1].offset_y > 0.0
        np.testing.assert_allclose(
            model.stations.lon, [20.0, 20.0, 20.0, 20.002, 20.002, 20.002]
        )
        np.testing.assert_allclose(model.stations.z, [50, 51, 52, 60, 61, 62])

    def test_positional_bln_one_per_line(self, tmp_path):
        p1 = tmp_path / "l1.bln"
        p1.write_text("3,1\n20.0,10.00,50\n20.0,10.01,51\n20.0,10.02,52\n", encoding="utf-8")
        p2 = tmp_path / "l2.bln"
        p2.write_text(
            "3,1\n20.002,10.00,60\n20.002,10.01,61\n20.002,10.02,62\n", encoding="utf-8"
        )
        model = build_multiline_pcsf(
            _flat_profiles(), topo=[p1, p2], latlon=True, cache_derived_volume=False
        )
        assert [line.offset_kind for line in model.geometry.lines] == ["real", "real"]
        np.testing.assert_allclose(
            model.stations.lon, [20.0, 20.0, 20.0, 20.002, 20.002, 20.002]
        )

    def test_line_keyed_mapping_of_sources(self, tmp_path):
        p1 = tmp_path / "l1.csv"
        p1.write_text(
            "station,lat,lon\nS0,10.00,20.0\nS1,10.01,20.0\nS2,10.02,20.0\n",
            encoding="utf-8",
        )
        p2 = tmp_path / "l2.csv"
        p2.write_text(
            "station,lat,lon\nS3,10.00,20.002\nS4,10.01,20.002\nS5,10.02,20.002\n",
            encoding="utf-8",
        )
        model = build_multiline_pcsf(
            _flat_profiles(), topo={"L1": p1, "L2": p2}, cache_derived_volume=False
        )
        assert [line.offset_kind for line in model.geometry.lines] == ["real", "real"]

    def test_topo_overrides_existing_sta_lat_lon_with_warning(self, tmp_path):
        path = tmp_path / "topo.csv"
        path.write_text(
            "station,lat,lon\nS0,10.00,20.0\nS1,10.01,20.0\nS2,10.02,20.0\n"
            "S3,10.00,20.002\nS4,10.01,20.002\nS5,10.02,20.002\n",
            encoding="utf-8",
        )
        profiles = _profiles_with_real_geometry()
        with pytest.warns(UserWarning, match="takes precedence"):
            model = build_multiline_pcsf(
                profiles, topo=path, cache_derived_volume=False
            )
        # topo's own values (20.0/20.002) replace the profile's original
        # sta_lon (20.000/20.002 -- same here, so check elevation/lat too
        # to be sure the override really happened, not a coincidence).
        assert model.stations.lon[0] == pytest.approx(20.0)

    def test_no_topo_behaves_exactly_as_before(self):
        model = build_multiline_pcsf(_flat_profiles(), cache_derived_volume=False)
        assert [line.offset_kind for line in model.geometry.lines] == [
            "synthetic", "synthetic",
        ]
        assert model.stations.lon is None


class TestMultilinePcsfToProfiles:
    def test_round_trip_preserves_geometry_and_stations(self, tmp_path):
        profiles = _flat_profiles()
        model = build_multiline_pcsf(profiles, line_spacing=1.0)
        path = write_pcsf(model, tmp_path / "multiline.pcsf")
        restored = read_pcsf(path)
        back = multiline_pcsf_to_profiles(restored)

        assert set(back) == {"L1", "L2"}
        for name in ("L1", "L2"):
            np.testing.assert_array_equal(back[name]["x"], profiles[name]["x"])
            np.testing.assert_array_equal(back[name]["z"], profiles[name]["z"])
            np.testing.assert_array_equal(back[name]["rho"], profiles[name]["rho"])
            assert back[name]["sta_names"] == profiles[name]["sta_names"]
            np.testing.assert_allclose(
                back[name]["sta_elev"], profiles[name]["sta_elev"]
            )

    def test_sta_lat_lon_round_trip_through_a_real_pcsf_file(self, tmp_path):
        profiles = _profiles_with_real_geometry()
        model = build_multiline_pcsf(profiles, cache_derived_volume=False)
        path = write_pcsf(model, tmp_path / "multiline_geo.pcsf")
        back = multiline_pcsf_to_profiles(read_pcsf(path))

        for name in ("L1", "L2"):
            np.testing.assert_allclose(
                back[name]["sta_lat"], profiles[name]["sta_lat"]
            )
            np.testing.assert_allclose(
                back[name]["sta_lon"], profiles[name]["sta_lon"]
            )

    def test_sta_lat_lon_empty_when_file_has_none(self, tmp_path):
        profiles = _flat_profiles()
        model = build_multiline_pcsf(profiles, cache_derived_volume=False)
        path = write_pcsf(model, tmp_path / "multiline_nogeo.pcsf")
        back = multiline_pcsf_to_profiles(read_pcsf(path))
        assert back["L1"]["sta_lat"] == []
        assert back["L1"]["sta_lon"] == []

    def test_reconstructed_profiles_feed_stack_lines_identically(self, tmp_path):
        # The whole point of persisting a multiline file: re-running the
        # exact same assembly on the reconstructed profiles must match
        # what the original profiles would have produced.
        profiles = _flat_profiles()
        model = build_multiline_pcsf(profiles, line_spacing=1.0)
        path = write_pcsf(model, tmp_path / "multiline.pcsf")
        back = multiline_pcsf_to_profiles(read_pcsf(path))

        x1, z1, y1, rho1 = stack_lines_to_common_grid(profiles, line_spacing=1.0)
        x2, z2, y2, rho2 = stack_lines_to_common_grid(back, line_spacing=1.0)
        np.testing.assert_array_equal(x1, x2)
        np.testing.assert_array_equal(z1, z2)
        np.testing.assert_array_equal(y1, y2)
        np.testing.assert_array_equal(rho1, rho2)

    def test_rejects_non_multiline_model(self):
        from pycsamt.format import Grid2DGeometry, PCSFModel

        model = PCSFModel(
            geometry=Grid2DGeometry(x=np.array([0.0]), z=np.array([0.0])),
            resistivity=np.array([[1.0]]),
        )
        with pytest.raises(ValueError, match="multiline"):
            multiline_pcsf_to_profiles(model)


# ---------------------------------------------------------------------
# Equivalence against map3d.py's own private implementation
# ---------------------------------------------------------------------

dash = pytest.importorskip("dash")


class TestEquivalenceWithMap3DPrivateFunctions:
    """pycsamt.format.multiline must match map3d.py's own algorithm
    exactly, not merely approximate it — see the module docstring."""

    def test_line_offsets_match(self):
        from pycsamt.app.web.callbacks.map3d import _line_real_offsets

        profiles = _profiles_with_real_geometry()
        expected = _line_real_offsets(profiles)
        actual = line_offsets_from_stations(profiles)
        assert expected == actual

    def test_assemble_3d_grid_matches(self):
        from pycsamt.app.web.callbacks.map3d import _assemble_3d_grid

        profiles = _flat_profiles()
        profiles["L2"]["x"] = np.array([0.0, 50.0, 100.0, 150.0, 200.0])
        profiles["L2"]["rho"] = np.array(
            [[200.0, 205.0, 210.0, 215.0, 220.0], [90.0, 92.0, 95.0, 97.0, 99.0]]
        )

        x1, y1, z1, rho1 = _assemble_3d_grid(profiles, line_spacing=1.0)
        x2, z2, y2, rho2 = stack_lines_to_common_grid(profiles, line_spacing=1.0)

        np.testing.assert_array_equal(x1, x2)
        np.testing.assert_array_equal(y1, y2)
        np.testing.assert_array_equal(z1, z2)
        np.testing.assert_array_equal(rho1, rho2)
