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
        # Sanity: most of the survey should cross-match. Only L18PLT and
        # L22PLT are bundled in the repo (the rest of WILLY_DATA is
        # gitignored), so a thin checkout matches far fewer -- skip rather
        # than fail when the full survey is not present.
        if len(elev_map) <= 100:
            pytest.skip(
                f"only {len(elev_map)} ModEM/WILLY stations cross-match; "
                "the full WILLY_DATA survey is not bundled in this checkout"
            )

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
        # Only L18PLT/L22PLT are bundled (rest of WILLY_DATA is
        # gitignored); skip when the full survey is not present.
        if len(renamed) <= 100:
            pytest.skip(
                f"only {len(renamed)} ModEM/WILLY stations cross-match; "
                "the full WILLY_DATA survey is not bundled in this checkout"
            )
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


def _synthetic_model(*, air_layers: int = 3, air_ohm_m: float = 1e12):
    """A tiny ModEmModel3D whose top *air_layers* z-layers are air fill.

    Mimics a real student ModEM run written with ``n_air == 0`` that
    still leaves the model cells above topography at 1e10-1e13 ohm.m
    (see fig08 / the Baohuashan discrepancy this feature fixes).
    """
    from pycsamt.models.modem.model3d import ModEmModel3D

    nz, ny, nx = 6, 4, 5
    m = ModEmModel3D()
    m.x_widths = np.full(nx, 100.0)
    m.y_widths = np.full(ny, 100.0)
    m.z_widths = np.full(nz, 50.0)
    rho = np.full((nz, ny, nx), 200.0)  # 200 ohm.m earth everywhere
    rho[:air_layers, :, :] = air_ohm_m  # complete air layers on top
    rho[air_layers, 0, 0] = air_ohm_m  # + one ragged air cell below them
    m.rho_loge = np.log(rho)
    m.n_air = 0
    m.log_type = "LOGE"
    m.origin = np.zeros(3)
    m.rotation = 0.0
    return m


def _fake_result(model):
    from types import SimpleNamespace

    return SimpleNamespace(
        mode="3d",
        model_final=model,
        model_initial=None,
        data_obs=None,
        data_pred=None,
        log=None,
        workdir="synthetic",
        final_rms=1.5,
        n_iter=10,
        models=[],
    )


class TestAirFillMasking:
    def test_air_fill_masked_by_default(self):
        model = _synthetic_model(air_layers=3)
        pcsf = modem3d_to_pcsf(_fake_result(model))
        pcsf.validate()

        rho = pcsf.resistivity
        # Every air-fill cell (the 3 top layers + the ragged cell) is NaN.
        assert np.isnan(rho[:3]).all()
        assert np.isnan(rho[3, 0, 0])
        # Resolved earth is untouched (bar exp/log float round-off).
        np.testing.assert_allclose(rho[3:][~np.isnan(rho[3:])], 200.0)
        # Native array keeps the raw values for provenance.
        np.testing.assert_array_equal(
            pcsf.resistivity_native, model.rho_loge
        )
        # Complete leading air layers are recorded as n_air.
        assert pcsf.geometry.n_air == 3
        meta = pcsf.metadata["air_mask"]
        assert meta["threshold_ohm_m"] == 1e8
        assert meta["n_cells_masked"] == 3 * 4 * 5 + 1
        assert meta["n_air_detected"] == 3

    def test_air_threshold_none_keeps_exact(self):
        model = _synthetic_model(air_layers=3)
        pcsf = modem3d_to_pcsf(_fake_result(model), air_threshold_ohm_m=None)
        np.testing.assert_array_equal(
            pcsf.resistivity, model.rho_linear
        )
        assert np.isfinite(pcsf.resistivity).all()
        assert pcsf.geometry.n_air == 0
        assert "air_mask" not in pcsf.metadata

    def test_clean_model_is_untouched(self):
        # No cell above the threshold -> a pure passthrough, no metadata.
        model = _synthetic_model(air_layers=0, air_ohm_m=200.0)
        pcsf = modem3d_to_pcsf(_fake_result(model))
        np.testing.assert_array_equal(pcsf.resistivity, model.rho_linear)
        assert "air_mask" not in pcsf.metadata
        assert pcsf.geometry.n_air == 0


def _fake_data(z_values, *, comment="", names=None, lonlat=None):
    """A minimal ModEmData-like for _stations_from_modem_data."""
    from types import SimpleNamespace

    names = names or [f"18-{i + 1:03d}" for i in range(len(z_values))]
    coords = {
        n: (float(i) * 100.0, 0.0, float(z))
        for i, (n, z) in enumerate(zip(names, z_values))
    }
    return SimpleNamespace(
        site_names=list(names),
        site_coords=coords,
        site_lonlat=dict(lonlat or {}),
        comment=comment,
    )


class TestLonLatQuantizationWarning:
    def _pcsf(self, lonlat):
        model = _synthetic_model(air_layers=0, air_ohm_m=200.0)
        res = _fake_result(model)
        res.data_obs = _fake_data(
            [0.0] * len(lonlat), names=list(lonlat), lonlat=lonlat
        )
        return res

    def test_warns_on_3dp_quantized_lonlat(self):
        # a ModEM '-R' rewrite echo: every value is exactly 3 dp
        ll = {
            "18-001": (119.127, 32.118),
            "18-002": (119.127, 32.118),
            "18-003": (119.126, 32.119),
            "18-004": (119.127, 32.120),
        }
        with pytest.warns(UserWarning, match="quantized to 3 decimal"):
            modem3d_to_pcsf(self._pcsf(ll))

    def test_no_warning_on_full_precision_lonlat(self):
        ll = {
            "18-001": (119.1269, 32.1179),
            "18-002": (119.1265, 32.1188),
            "18-003": (119.1266, 32.1197),
            "18-004": (119.1264, 32.1206),
        }
        import warnings as _w

        with _w.catch_warnings():
            _w.simplefilter("error")  # any UserWarning -> test failure
            modem3d_to_pcsf(self._pcsf(ll))


class TestStationZConvention:
    _COMMENT = "Baohuashan. Z(m) is depth below model top (top = 224 m a.s.l.)"

    def _pcsf(self, data, **kw):
        model = _synthetic_model(air_layers=0, air_ohm_m=200.0)
        res = _fake_result(model)
        res.data_obs = data
        return modem3d_to_pcsf(res, **kw)

    def test_depth_down_comment_flips_to_absolute_elevation(self):
        data = _fake_data([125.0, 114.0, 143.0], comment=self._COMMENT)
        with pytest.warns(UserWarning, match="positive-down depth"):
            pcsf = self._pcsf(data)
        # elevation = 224 - Z
        np.testing.assert_allclose(pcsf.stations.z, [99.0, 110.0, 81.0])
        assert pcsf.metadata["station_z"] == {
            "convention": "auto",
            "flipped": True,
            "datum_masl": 224.0,
        }

    def test_auto_leaves_plain_elevation_alone(self):
        # No depth-down marker in the comment -> trust the column.
        data = _fake_data([99.0, 110.0, 81.0], comment="a normal survey")
        pcsf = self._pcsf(data)
        np.testing.assert_allclose(pcsf.stations.z, [99.0, 110.0, 81.0])
        assert pcsf.metadata["station_z"]["flipped"] is False

    def test_flat_zero_placeholder_untouched(self):
        data = _fake_data([0.0, 0.0, 0.0], comment=self._COMMENT)
        pcsf = self._pcsf(data)
        np.testing.assert_array_equal(pcsf.stations.z, [0.0, 0.0, 0.0])
        assert pcsf.metadata["station_z"]["flipped"] is False

    def test_depth_down_without_datum_is_relative(self):
        data = _fake_data(
            [125.0, 114.0, 143.0], comment="Z is depth below datum"
        )
        with pytest.warns(UserWarning, match="relative elevation"):
            pcsf = self._pcsf(data)
        # deepest station (143) -> 0, shallower ones positive
        np.testing.assert_allclose(pcsf.stations.z, [18.0, 29.0, 0.0])
        assert pcsf.metadata["station_z"]["datum_masl"] is None

    def test_force_elevation_keeps_verbatim(self):
        data = _fake_data([125.0, 114.0, 143.0], comment=self._COMMENT)
        pcsf = self._pcsf(data, station_z_convention="elevation")
        np.testing.assert_allclose(pcsf.stations.z, [125.0, 114.0, 143.0])
        assert pcsf.metadata["station_z"]["flipped"] is False

    def test_force_depth_down_flips_even_without_comment(self):
        data = _fake_data([10.0, 20.0], comment="")
        pcsf = self._pcsf(data, station_z_convention="depth_down")
        np.testing.assert_allclose(pcsf.stations.z, [10.0, 0.0])
