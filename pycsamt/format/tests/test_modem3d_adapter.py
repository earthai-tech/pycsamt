# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Phase 3 tests — the ModEM 3-D adapter, against the real bundled
``data/modem/willy_27freq_watex_line02_sample`` dataset (a 41x50x288
real 3-D inversion volume, 125 stations, 74 iterations) -- the first
genuinely new persisted 3-D volume artifact in the project.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.format import read_pcsf, write_pcsf
from pycsamt.format.adapters.modem3d import modem3d_to_pcsf
from pycsamt.models.modem.results import InversionResult

_DATA_DIR = (
    Path(__file__).parents[3] / "data" / "modem" / "willy_27freq_watex_line02_sample"
)
_SKIP = pytest.mark.skipif(
    not _DATA_DIR.exists(), reason=f"bundled ModEM data not found: {_DATA_DIR}"
)


@pytest.fixture(scope="module")
def result() -> InversionResult:
    return InversionResult(workdir=_DATA_DIR, load_data=True)


@pytest.fixture(scope="module")
def model(result):
    return modem3d_to_pcsf(
        result,
        created_by="pytest",
        crs="local",
        description="Willy L18 line02 27-freq ModEM 3D",
    )


@_SKIP
class TestModEm3DAdapter:
    def test_geometry_kind_and_shape(self, result, model):
        assert model.kind == "grid3d"
        assert model.resistivity.shape == result.model_final.shape
        model.validate()

    def test_resistivity_matches_source_model_exactly(self, result, model):
        # HDF5 float64 storage should be bit-exact from the in-memory
        # array (unlike the ASCII .rho format's own ~5-sig-fig
        # round-trip precision, which is unrelated to PCSF).
        np.testing.assert_array_equal(
            model.resistivity, result.model_final.rho_linear
        )
        np.testing.assert_array_equal(
            model.resistivity_native, result.model_final.rho_loge
        )
        assert model.resistivity_native_encoding == "ln"

    def test_axis_order_is_z_y_x(self, result, model):
        nz, ny, nx = (
            result.model_final.nz,
            result.model_final.ny,
            result.model_final.nx,
        )
        assert model.resistivity.shape == (nz, ny, nx)
        assert model.geometry.z.shape == (nz,)
        assert model.geometry.y.shape == (ny,)
        assert model.geometry.x.shape == (nx,)

    def test_real_grid_origin_and_rotation_preserved(self, result, model):
        # Regression guard for the exact bug found while building this
        # adapter: ModEmModel3D used to silently discard the trailing
        # centre-coordinate/rotation line a real ModEM writer appends.
        np.testing.assert_allclose(
            model.geometry.origin, [-4509.828, -6752.725, 0.0]
        )
        assert model.geometry.rotation_deg == pytest.approx(0.0)
        assert model.geometry.n_air == result.model_final.n_air

    def test_stations_from_real_data_file(self, result, model):
        assert model.stations is not None
        assert len(model.stations.name) == result.data_obs.n_sites
        first = model.stations.name[0]
        expected_x, expected_y, expected_z = result.data_obs.site_coords[first]
        assert model.stations.x[0] == pytest.approx(expected_x)
        assert model.stations.y[0] == pytest.approx(expected_y)
        assert model.stations.z[0] == pytest.approx(expected_z)

    def test_stations_lon_lat_from_real_dat_gg_columns(self, result, model):
        # data/modem/willy_27freq_watex_line02_sample's .dat file carries
        # real GG_Lat/GG_Lon columns -- the adapter must pass them
        # straight through to StationTable.lon/lat, no known_stations
        # match required, so a grid3d PCSF file this adapter writes is
        # self-sufficiently geo-referenced.
        assert result.data_obs.site_lonlat  # real data has them
        assert model.stations.lon is not None
        assert model.stations.lat is not None
        first = model.stations.name[0]
        expected_lon, expected_lat = result.data_obs.site_lonlat[first]
        assert model.stations.lon[0] == pytest.approx(expected_lon)
        assert model.stations.lat[0] == pytest.approx(expected_lat)

    def test_stations_lon_lat_survive_pcsf_round_trip(self, model, tmp_path):
        path = write_pcsf(model, tmp_path / "modem3d_lonlat.pcsf")
        restored = read_pcsf(path)
        np.testing.assert_allclose(restored.stations.lon, model.stations.lon)
        np.testing.assert_allclose(restored.stations.lat, model.stations.lat)

    def test_station_elevations_override_and_populate_topography(self, result):
        from pycsamt.map import load_lines

        willy_dir = (
            Path(__file__).parents[3] / "data" / "AMT" / "WILLY_DATA"
        )
        if not willy_dir.exists():
            pytest.skip(f"bundled WILLY_DATA not found: {willy_dir}")
        willy = load_lines(willy_dir, detect="folder")
        willy_elev = {
            s.id: s.elevation for s in willy.stations if s.elevation is not None
        }
        # Real ModEM station names carry a "23-" survey-year prefix
        # WILLY_DATA's own ids don't -- these are the same physical
        # stations under two naming conventions.
        elev_map = {
            name: willy_elev[name.split("-", 1)[1]]
            for name in result.data_obs.site_names
            if name.startswith("23-") and name.split("-", 1)[1] in willy_elev
        }
        assert len(elev_map) > 100  # sanity: most of the survey matches

        model = modem3d_to_pcsf(result, station_elevations=elev_map)
        model.validate()
        assert model.topography is not None
        assert len(model.topography.station_id) == len(elev_map)
        by_name = dict(zip(model.stations.name, model.stations.z))
        first_name, first_elev = next(iter(elev_map.items()))
        assert by_name[first_name] == pytest.approx(first_elev)
        # An unmatched station keeps ModEM's own recorded z (0.0 is
        # real here, not "unknown").
        unmatched = [n for n in result.data_obs.site_names if n not in elev_map]
        if unmatched:
            assert by_name[unmatched[0]] == pytest.approx(0.0)

    def _renamed_willy_sites(self):
        """Wrap real WILLY_DATA StationRecords under ModEM's own
        "23-"-prefixed site names -- exercises topo_from_sites's
        duck-typed extraction with a real Sites/MapData object, the
        same real-data cross-reference
        test_station_elevations_override_and_populate_topography
        already establishes for station_elevations."""
        from types import SimpleNamespace

        from pycsamt.map import load_lines

        willy_dir = Path(__file__).parents[3] / "data" / "AMT" / "WILLY_DATA"
        if not willy_dir.exists():
            pytest.skip(f"bundled WILLY_DATA not found: {willy_dir}")
        willy = load_lines(willy_dir, detect="folder")
        willy_by_id = {s.id: s for s in willy.stations}
        modem_names = list(InversionResult(workdir=_DATA_DIR, load_data=True).data_obs.site_names)
        renamed = [
            SimpleNamespace(
                id=name,
                longitude=willy_by_id[name.split("-", 1)[1]].longitude,
                latitude=willy_by_id[name.split("-", 1)[1]].latitude,
                elevation=willy_by_id[name.split("-", 1)[1]].elevation,
            )
            for name in modem_names
            if name.startswith("23-") and name.split("-", 1)[1] in willy_by_id
        ]
        assert len(renamed) > 100  # sanity: most of the survey matches
        return renamed

    def test_topo_from_real_sites_object_overrides_dat_lonlat(self, result):
        renamed = self._renamed_willy_sites()
        with pytest.warns(UserWarning, match="takes precedence"):
            model = modem3d_to_pcsf(result, topo=renamed)
        model.validate()
        by_name = dict(zip(model.stations.name, model.stations.lon))
        first = renamed[0]
        assert by_name[first.id] == pytest.approx(first.longitude)
        # every station still ends up with *some* lon/lat: topo covers
        # the "23-" matched majority, the .dat file's own GG_Lat/GG_Lon
        # fills in the rest.
        assert np.all(np.isfinite(model.stations.lon))

    def test_topo_elevation_populates_topography_group(self, result):
        renamed = self._renamed_willy_sites()
        with pytest.warns(UserWarning):
            model = modem3d_to_pcsf(result, topo=renamed)
        assert model.topography is not None
        assert len(model.topography.station_id) == len(renamed)
        by_name = dict(zip(model.stations.name, model.stations.z))
        assert by_name[renamed[0].id] == pytest.approx(renamed[0].elevation)

    def test_no_topo_behaves_exactly_as_before(self, result, model):
        # Regression guard: passing no topo must not change any
        # existing behaviour (same fixture 'model' as the rest of this
        # test class, built with no topo/station_elevations at all).
        assert model.stations.lon is not None  # from the .dat file itself
        assert model.topography is None

    def test_history_from_log(self, result, model):
        assert set(model.history) == {
            "iteration", "rms", "objective", "model_norm", "lagrange", "alpha",
        }
        n = result.log.iterations.size
        for arr in model.history.values():
            assert arr.shape == (n,)
        np.testing.assert_allclose(model.history["rms"][-1], result.final_rms)

    def test_metadata_carries_provenance(self, result, model):
        assert model.metadata["mode"] == "3d"
        assert model.metadata["n_iter"] == result.n_iter
        assert model.metadata["final_rms"] == pytest.approx(result.final_rms)
        assert "iter_0073" in model.metadata["model_keys"]

    def test_full_round_trip_through_pcsf_file_is_bit_exact(
        self, result, model, tmp_path
    ):
        path = write_pcsf(model, tmp_path / "modem3d_real.pcsf")
        restored = read_pcsf(path)

        assert restored.source_backend == "modem3d"
        np.testing.assert_array_equal(restored.resistivity, model.resistivity)
        np.testing.assert_array_equal(
            restored.geometry.origin, model.geometry.origin
        )
        np.testing.assert_array_equal(
            restored.geometry.x_nodes, model.geometry.x_nodes
        )
        assert restored.stations.name == model.stations.name
        np.testing.assert_allclose(restored.history["rms"], model.history["rms"])
        assert restored.metadata["n_iter"] == result.n_iter

    def test_raises_for_non_3d_result(self, monkeypatch, result):
        monkeypatch.setattr(result, "mode", "2d", raising=False)
        with pytest.raises(ValueError, match="3-D ModEM result"):
            modem3d_to_pcsf(result)

    def test_raises_without_any_model(self):
        empty = InversionResult.__new__(InversionResult)
        empty.mode = "3d"
        empty.model_final = None
        empty.model_initial = None
        with pytest.raises(ValueError, match="model_final/model_initial"):
            modem3d_to_pcsf(empty)
