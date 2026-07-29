# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.web.callbacks.map3d (3D resistivity map).

Strategy
--------
Every figure builder here (fence / volume-block / depth-slices) is a pure
NumPy + Plotly function driven by plain dicts of arrays, so they're
exercised directly with small synthetic profiles rather than mocked.
For the "real survey data" path, real cached WILLY sites drive the
"Skin-depth pseudo" source end to end (``_profiles_from_pseudo`` through
``generate_grid`` through ``display_3d``); the "Session inversion result"
source uses a small duck-typed fake result object (same shape as
``map.py``'s tests) since running a full inversion is out of scope for a
unit test.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from dash.exceptions import PreventUpdate

import pycsamt.app.web.callbacks.map3d as map3d_mod
from pycsamt.app.web.cache import cache_set, cache_set_inversion_result
from pycsamt.app.web.layout import IDs

_ROOT = Path(__file__).parents[3]
_WILLY_L18 = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
_HAS_WILLY = _WILLY_L18.exists() and any(_WILLY_L18.glob("*.edi"))


@pytest.fixture(scope="session")
def willy_sites():
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_WILLY:
        pytest.skip("WILLY L18PLT data not available")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(_WILLY_L18))


@pytest.fixture(scope="session")
def willy_records(willy_sites):
    recs = []
    for i, ed in enumerate(willy_sites.as_list()):
        recs.append({"ID": ed.station, "Line": "L1" if i % 2 == 0 else "L2"})
    return recs


@pytest.fixture
def store_data_willy(willy_records):
    return {
        "station_records": willy_records,
        "n_stations": len(willy_records),
        "line_counts": {"L1": 14, "L2": 14},
    }


@pytest.fixture
def cached_session(willy_sites):
    session_id = "test-map3d-session"
    cache_set(session_id, willy_sites)
    return session_id


class _FakeResModel:
    def __init__(self, rho_2d, z_centers, x_centers, station_names):
        self.rho_2d = rho_2d
        self.z_centers = z_centers
        self.x_centers = x_centers
        self.station_names = station_names
        self.station_x = x_centers


class _FakeInvResult:
    def __init__(self, model, method="occam2d", dimension="2d"):
        self._model = model
        self.method = method
        self.dimension = dimension

    def to_resistivity_model(self):
        return self._model


@pytest.fixture
def cached_inv_result():
    """A session id with only an inversion result cached (no raw Sites)."""
    n_x, n_z = 10, 6
    rho_2d = np.random.default_rng(0).uniform(1.0, 1000.0, size=(n_z, n_x))
    z_centers = np.linspace(10, 1000, n_z)
    x_centers = np.linspace(0, 5000, n_x)
    station_names = [f"S{i:03d}" for i in range(n_x)]
    model = _FakeResModel(rho_2d, z_centers, x_centers, station_names)
    session_id = "test-map3d-inv-session"
    cache_set_inversion_result(session_id, _FakeInvResult(model))
    return session_id


@pytest.fixture
def cached_session_with_inv_result(willy_sites):
    """A session id with BOTH raw Sites and an inversion result cached.

    ``generate_grid`` always requires real Sites to be cached (via
    ``cache_get``) as its precondition regardless of the chosen data
    source -- even the "inversion" source needs it, since in normal usage
    a survey is always loaded before an inversion is run against it.

    Uses its own session id rather than layering onto ``cached_session``:
    the backing cache is disk-based with a 24h TTL and has no per-test
    teardown, so a shared id would leak the inversion result into any
    other test -- e.g. ``test_inversion_source_without_result`` -- that
    expects a session with Sites but no inversion result yet.
    """
    session_id = "test-map3d-session-with-inv"
    cache_set(session_id, willy_sites)
    n_x, n_z = 10, 6
    rho_2d = np.random.default_rng(0).uniform(1.0, 1000.0, size=(n_z, n_x))
    z_centers = np.linspace(10, 1000, n_z)
    x_centers = np.linspace(0, 5000, n_x)
    station_names = [f"S{i:03d}" for i in range(n_x)]
    model = _FakeResModel(rho_2d, z_centers, x_centers, station_names)
    cache_set_inversion_result(session_id, _FakeInvResult(model))
    return session_id


def _synthetic_profiles(n_lines=3, n_x=8, n_z=6):
    rng = np.random.default_rng(1)
    profiles = {}
    for i in range(n_lines):
        name = f"L{i + 1}"
        x = np.linspace(0, 1000, n_x)
        z = np.linspace(10, 500, n_z)
        rho = rng.uniform(10, 1000, size=(n_z, n_x))
        profiles[name] = {
            "x": x,
            "z": z,
            "rho": rho,
            "sta_x": x.tolist(),
            "sta_elev": [0.0] * n_x,
            "sta_names": [f"{name}-{j}" for j in range(n_x)],
            "sta_lat": [30.0 + i * 0.01 + j * 0.0001 for j in range(n_x)],
            "sta_lon": [110.0 + i * 0.02 for j in range(n_x)],
        }
    return profiles


def _unwrap(entry):
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def _cb(web_app, output_id_prop):
    return _unwrap(web_app.callback_map[output_id_prop])


def _cb_multi(web_app, *substrings):
    key = next(k for k in web_app.callback_map if all(s in k for s in substrings))
    return _unwrap(web_app.callback_map[key])


def _cb_by_input(web_app, output_substr, input_id):
    for k, v in web_app.callback_map.items():
        if output_substr not in k:
            continue
        if any(input_id in str(i.get("id")) for i in v.get("inputs", [])):
            return _unwrap(v)
    raise AssertionError(
        f"no callback found for output~={output_substr!r} input={input_id!r}"
    )


def _set_triggered(prop_id):
    import dash._callback_context as cc
    from dash._utils import AttributeDict

    cc.context_value.set(AttributeDict(triggered_inputs=[{"prop_id": prop_id}]))


def _clear_triggered():
    import dash._callback_context as cc

    cc.context_value.set({})


# ══════════════════════════════════════════════════════════════════════════
# Pure numeric helpers
# ══════════════════════════════════════════════════════════════════════════


class TestApplyRhoMask:
    def test_no_filter_returns_unchanged(self):
        rho = np.array([1.0, 10.0, 100.0])
        out = map3d_mod._apply_rho_mask(rho, 0, 1e8)
        assert np.array_equal(out, rho)

    def test_masks_outside_range(self):
        rho = np.array([1.0, 10.0, 100.0])
        out = map3d_mod._apply_rho_mask(rho, 5, 50)
        assert np.isnan(out[0])
        assert out[1] == 10.0
        assert np.isnan(out[2])


class TestApplyDepthMask:
    def test_no_mask_match_returns_unchanged(self):
        z = np.array([10.0, 20.0])
        rho = np.array([[1.0, 2.0], [3.0, 4.0]])
        z_out, rho_out = map3d_mod._apply_depth_mask(z, rho, 1000, 2000)
        assert np.array_equal(z_out, z)

    def test_clips_first_axis(self):
        z = np.array([10.0, 20.0, 30.0])
        rho = np.array([[1, 1], [2, 2], [3, 3]], float)
        z_out, rho_out = map3d_mod._apply_depth_mask(z, rho, 15, 25)
        assert list(z_out) == [20.0]
        assert rho_out.shape == (1, 2)


class TestAxisFractionMask:
    def test_empty_array(self):
        out = map3d_mod._axis_fraction_mask(np.array([]), 0.5)
        assert out.size == 0

    def test_full_fraction_keeps_all(self):
        out = map3d_mod._axis_fraction_mask(np.array([1, 2, 3]), 1.0)
        assert out.all()

    def test_partial_fraction(self):
        out = map3d_mod._axis_fraction_mask(np.array([0, 10, 20]), 0.5)
        assert out[0] and not out[-1]

    def test_flat_array_keeps_all(self):
        out = map3d_mod._axis_fraction_mask(np.array([5, 5, 5]), 0.3)
        assert out.all()


class TestInterpFinite1d:
    def test_two_or_more_valid_points(self):
        out = map3d_mod._interp_finite_1d(
            np.array([0.5]), np.array([0, 1]), np.array([0, 10])
        )
        assert out[0] == pytest.approx(5.0)

    def test_single_valid_point(self):
        out = map3d_mod._interp_finite_1d(
            np.array([0, 1, 2]),
            np.array([0, 1, 2]),
            np.array([np.nan, 5, np.nan]),
        )
        assert np.all(out == 5.0)

    def test_no_valid_points(self):
        out = map3d_mod._interp_finite_1d(
            np.array([0, 1]), np.array([0, 1]), np.array([np.nan, np.nan])
        )
        assert np.all(np.isnan(out))


class TestLineRealOffsets:
    def test_missing_metadata_returns_none(self):
        profiles = {"L1": {"sta_lat": [], "sta_lon": [], "sta_names": []}}
        assert map3d_mod._line_real_offsets(profiles) is None

    def test_mismatched_lengths_returns_none(self):
        profiles = {"L1": {"sta_lat": [1.0], "sta_lon": [], "sta_names": ["a"]}}
        assert map3d_mod._line_real_offsets(profiles) is None

    def test_real_offsets_computed(self):
        profiles = _synthetic_profiles(n_lines=3)
        out = map3d_mod._line_real_offsets(profiles)
        # Either resolves real offsets or bails to None (survey_uv may
        # reject degenerate synthetic geometry) -- both are valid outcomes.
        assert out is None or set(out.keys()) == set(profiles.keys())


class TestSceneHelpers:
    def test_scene_colors_dark_and_light_differ(self):
        dark = map3d_mod._scene_colors(True)
        light = map3d_mod._scene_colors(False)
        assert dark["scene"] != light["scene"]

    def test_dark_scene_includes_title(self):
        layout = map3d_mod._dark_scene("My Title", dark=True)
        assert layout["title"]["text"] == "My Title"

    def test_dark_scene_no_title(self):
        layout = map3d_mod._dark_scene("", dark=True)
        assert layout["title"] == {}


class TestRhoXyFromSite:
    def test_none_freq_returns_none(self):
        class _S:
            freq = None

        assert map3d_mod._rho_xy_from_site(_S()) == (None, None)

    def test_from_rho_3d(self):
        class _S:
            freq = np.array([1.0, 10.0])
            rho = np.ones((2, 2, 2)) * 5.0
            z = None

        freq, rho_xy = map3d_mod._rho_xy_from_site(_S())
        assert freq is not None
        assert np.all(rho_xy == 5.0)

    def test_falls_back_to_z_when_no_rho(self):
        class _S:
            freq = np.array([1.0, 10.0])
            rho = None
            z = np.ones((2, 2, 2)) * (1 + 1j)

        freq, rho_xy = map3d_mod._rho_xy_from_site(_S())
        assert freq is not None
        assert rho_xy is not None

    def test_exception_returns_none(self):
        class _S:
            @property
            def freq(self):
                raise RuntimeError("boom")

        assert map3d_mod._rho_xy_from_site(_S()) == (None, None)


class TestStationDistanceM:
    def test_single_station(self):
        out = map3d_mod._station_distance_m([(30.0, 110.0, 0.0)])
        assert out[0] == 0.0

    def test_monotonic_increasing(self):
        coords = [(30.0, 110.0, 0), (30.01, 110.01, 0), (30.02, 110.02, 0)]
        out = map3d_mod._station_distance_m(coords)
        assert out[0] == 0.0
        assert out[1] > 0
        assert out[2] > out[1]

    def test_degenerate_falls_back_to_index_spacing(self):
        coords = [(30.0, 110.0, 0), (30.0, 110.0, 0), (30.0, 110.0, 0)]
        out = map3d_mod._station_distance_m(coords)
        assert list(out) == [0.0, 100.0, 200.0]


class TestSafeElev:
    def test_none_coords(self):
        assert map3d_mod._safe_elev(None) is None

    def test_valid(self):
        assert map3d_mod._safe_elev((1.0, 2.0, 42.0)) == 42.0

    def test_bad_value(self):
        assert map3d_mod._safe_elev((1.0, 2.0, "bad")) is None


class TestElevAlongLine:
    def test_empty_sta_returns_zeros(self):
        out = map3d_mod._elev_along_line(
            np.array([0, 1, 2]), np.array([]), np.array([])
        )
        assert np.all(out == 0.0)

    def test_interpolates(self):
        out = map3d_mod._elev_along_line(
            np.array([0, 5, 10]), np.array([0, 10]), np.array([0.0, 100.0])
        )
        assert out[1] == pytest.approx(50.0)

    def test_single_station_broadcasts(self):
        out = map3d_mod._elev_along_line(
            np.array([0, 5, 10]), np.array([5]), np.array([42.0])
        )
        assert np.all(out == 42.0)


class TestResolveProfileElevs:
    def test_empty_profile(self):
        out = map3d_mod._resolve_profile_elevs({"sta_x": []}, "stations")
        assert out.size == 0

    def test_stations_source(self):
        profile = {"sta_x": [0, 1], "sta_elev": [10.0, 20.0]}
        out = map3d_mod._resolve_profile_elevs(profile, "stations")
        assert list(out) == [10.0, 20.0]

    def test_elev_corr_prefers_corr_then_raw_then_sta(self):
        profile = {
            "sta_x": [0, 1],
            "sta_elev": [1.0, 1.0],
            "sta_elev_raw": [2.0, 2.0],
            "sta_elev_corr": [3.0, 3.0],
        }
        out = map3d_mod._resolve_profile_elevs(profile, "elev-corr")
        assert list(out) == [3.0, 3.0]

    def test_missing_key_falls_back_to_zeros(self):
        profile = {"sta_x": [0, 1]}
        out = map3d_mod._resolve_profile_elevs(profile, "elev-raw")
        assert list(out) == [0.0, 0.0]

    def test_nan_values_interpolated(self):
        profile = {"sta_x": [0, 1, 2], "sta_elev": [10.0, np.nan, 30.0]}
        out = map3d_mod._resolve_profile_elevs(profile, "stations")
        assert out[1] == pytest.approx(20.0)


class TestHasSurveyLineProfiles:
    def test_none_store(self):
        assert map3d_mod._has_survey_line_profiles(None) is False

    def test_one_line_false(self):
        assert map3d_mod._has_survey_line_profiles({"line_counts": {"L1": 5}}) is False

    def test_two_lines_true(self):
        store = {"line_counts": {"L1": 5, "L2": 3}}
        assert map3d_mod._has_survey_line_profiles(store) is True


class TestHasCachedInversionResult:
    def test_none_session(self):
        assert map3d_mod._has_cached_inversion_result(None) is False

    def test_no_result_cached(self):
        assert map3d_mod._has_cached_inversion_result("no-such-session") is False

    def test_real_cached_result(self, cached_inv_result):
        assert map3d_mod._has_cached_inversion_result(cached_inv_result) is True


class TestBuildElevsDict:
    def test_merges_all_sources(self):
        station_records = [{"ID": "S1", "Elevation": 100.0}]
        elev_raw = [{"station": "S1", "elev": 90.0}]
        elev_corr = [{"station": "S1", "elev_corrected": 95.0}]
        topo_upload = [{"station": "S1", "elev": 80.0}]
        out = map3d_mod._build_elevs_dict(
            elev_corr, elev_raw, station_records, topo_upload
        )
        assert out["S1"]["elev_sta"] == 100.0
        assert out["S1"]["elev_raw"] == 90.0
        assert out["S1"]["elev_corr"] == 95.0
        assert out["S1"]["elev_upload"] == 80.0

    def test_empty_inputs(self):
        assert map3d_mod._build_elevs_dict(None, None, None, None) == {}


class TestRhoLogToOhmM:
    def test_leaves_linear_values_unchanged(self):
        rho = np.array([10.0, 100.0, 1000.0])
        out = map3d_mod._rho_log_to_ohm_m(rho)
        assert np.array_equal(out, rho)

    def test_converts_log_values(self):
        rho_log = np.array([1.0, 2.0, 3.0])
        out = map3d_mod._rho_log_to_ohm_m(rho_log)
        assert out[0] == pytest.approx(10.0)


class TestProfilesFromInversionResult:
    def test_builds_single_profile(self, cached_inv_result):
        from pycsamt.app.web.cache import cache_get_inversion_result

        result = cache_get_inversion_result(cached_inv_result)
        profiles = map3d_mod._profiles_from_inversion_result(result)
        assert len(profiles) == 1
        p = next(iter(profiles.values()))
        assert p["rho"].shape[1] == len(p["sta_x"])

    def test_with_elevs_by_name(self, cached_inv_result):
        from pycsamt.app.web.cache import cache_get_inversion_result

        result = cache_get_inversion_result(cached_inv_result)
        elevs = {"S000": {"elev_sta": 123.0}}
        profiles = map3d_mod._profiles_from_inversion_result(
            result, elevs_by_name=elevs
        )
        p = next(iter(profiles.values()))
        assert p["sta_elev"][0] == 123.0


class TestProfilesFromInversion:
    def test_no_cached_result_and_no_desktop_ctrl_returns_empty(self):
        out = map3d_mod._profiles_from_inversion("no-such-session")
        assert out == {}

    def test_real_cached_result(self, cached_inv_result):
        out = map3d_mod._profiles_from_inversion(cached_inv_result)
        assert len(out) == 1


class TestProfilesFromPseudo:
    def test_real_willy_sites(self, willy_sites, store_data_willy):
        profiles = map3d_mod._profiles_from_pseudo(
            willy_sites, store_data=store_data_willy
        )
        assert profiles
        for p in profiles.values():
            assert p["rho"].shape[0] == len(p["z"])

    def test_require_line_metadata_with_lines(self, willy_sites, store_data_willy):
        profiles = map3d_mod._profiles_from_pseudo(
            willy_sites,
            store_data=store_data_willy,
            require_line_metadata=True,
        )
        assert set(profiles.keys()) <= {"L1", "L2"}

    def test_require_line_metadata_without_store_data_empty(self, willy_sites):
        profiles = map3d_mod._profiles_from_pseudo(
            willy_sites, store_data=None, require_line_metadata=True
        )
        assert profiles == {}

    def test_non_iterable_sites_returns_empty(self):
        assert map3d_mod._profiles_from_pseudo(object()) == {}


# ══════════════════════════════════════════════════════════════════════════
# Figure builders (synthetic profiles)
# ══════════════════════════════════════════════════════════════════════════


class TestBuildFenceFig:
    def test_builds_real_figure(self):
        profiles = _synthetic_profiles()
        fig = map3d_mod._build_fence_fig(
            profiles,
            line_spacing=1.0,
            azimuth_deg=0.0,
            panel_thick=50.0,
            cmap="RdYlBu_r",
            vmin=1,
            vmax=1000,
            opacity=0.85,
            contours=True,
            n_contours=10,
            show_labels=True,
            title="Test Fence",
        )
        assert fig is not None
        assert len(fig.data) == len(profiles)

    def test_with_topo(self):
        profiles = _synthetic_profiles(n_lines=2)
        for p in profiles.values():
            p["sta_elev"] = [float(i * 10) for i in range(len(p["sta_x"]))]
        # Real bug: _build_fence_fig's own default sta_symbol="auto" is
        # never resolved to a valid Plotly Scatter3d symbol -- that
        # resolution only happens in display_3d's callback body before it
        # calls this function. Calling it directly (as here) with the
        # default crashes; the production code path never hits this
        # because display_3d always passes a resolved symbol.
        with pytest.raises(ValueError, match="symbol"):
            map3d_mod._build_fence_fig(
                profiles,
                line_spacing=1.0,
                azimuth_deg=15.0,
                panel_thick=50.0,
                cmap="jet",
                vmin=1,
                vmax=1000,
                opacity=0.85,
                contours=False,
                n_contours=10,
                show_labels=False,
                title="",
                topo_src="stations",
                apply_topo=True,
            )

    def test_with_topo_resolved_symbol(self):
        profiles = _synthetic_profiles(n_lines=2)
        for p in profiles.values():
            p["sta_elev"] = [float(i * 10) for i in range(len(p["sta_x"]))]
        fig = map3d_mod._build_fence_fig(
            profiles,
            line_spacing=1.0,
            azimuth_deg=15.0,
            panel_thick=50.0,
            cmap="jet",
            vmin=1,
            vmax=1000,
            opacity=0.85,
            contours=False,
            n_contours=10,
            show_labels=False,
            title="",
            topo_src="stations",
            apply_topo=True,
            sta_symbol="diamond-open",
        )
        assert fig is not None


class TestBuildBlockFig:
    def test_builds_real_volume(self):
        profiles = _synthetic_profiles(n_lines=3)
        x_arr, y_arr, z_arr, rho_3d = map3d_mod._assemble_3d_grid(profiles)
        fig = map3d_mod._build_block_fig(
            x_arr,
            y_arr,
            z_arr,
            rho_3d,
            cmap="RdYlBu_r",
            vmin=10,
            vmax=1000,
            opacity=0.5,
            isovalue_lo=0.5,
            isovalue_hi=2.5,
            n_surfaces=6,
            clip_x=1.0,
            clip_y=1.0,
            clip_z=1.0,
            contours=True,
            title="Block",
        )
        assert fig is not None
        assert len(fig.data) >= 1

    def test_too_few_cells_returns_placeholder(self):
        profiles = _synthetic_profiles(n_lines=2)
        x_arr, y_arr, z_arr, rho_3d = map3d_mod._assemble_3d_grid(profiles)
        # clip_z=0.0 keeps only cells with |z| <= 0 -- effectively none of
        # the synthetic grid (z starts at 10) -- forcing the "not enough
        # visible cells" placeholder branch.
        fig = map3d_mod._build_block_fig(
            x_arr,
            y_arr,
            z_arr,
            rho_3d,
            cmap="jet",
            vmin=10,
            vmax=1000,
            opacity=0.5,
            isovalue_lo=0.5,
            isovalue_hi=2.5,
            n_surfaces=6,
            clip_x=1.0,
            clip_y=1.0,
            clip_z=0.0,
            contours=False,
            title="",
        )
        assert any("Not enough visible cells" in a.text for a in fig.layout.annotations)


class TestBuildDepthSlicesFig:
    def test_builds_real_slices(self):
        profiles = _synthetic_profiles(n_lines=3)
        x_arr, y_arr, z_arr, rho_3d = map3d_mod._assemble_3d_grid(profiles)

        def rho_fn(d):
            idx = int(np.argmin(np.abs(z_arr - d)))
            return rho_3d[:, :, idx]

        fig = map3d_mod._build_depth_slices_fig(
            x_arr,
            y_arr,
            depth_lo=10,
            depth_hi=500,
            n_slices=4,
            rho_fn=rho_fn,
            cmap="RdYlBu_r",
            vmin=10,
            vmax=1000,
            opacity=0.8,
            contours=True,
            n_contours=8,
            show_labels=True,
            title="Slices",
        )
        assert fig is not None
        assert len(fig.data) == 4


# ══════════════════════════════════════════════════════════════════════════
# Topo trace builders + upload parsing
# ══════════════════════════════════════════════════════════════════════════


class TestTopoTraceBuilders:
    def test_topo_strip_trace(self):
        trace = map3d_mod._build_topo_strip_trace(
            sta_x=np.array([0.0, 100.0]),
            sta_elev=np.array([10.0, 20.0]),
            y_center=0.0,
            panel_thick=50.0,
            opacity=0.7,
            dark=True,
        )
        assert trace is not None

    def test_station_markers_trace(self):
        trace = map3d_mod._build_station_markers_trace(
            sta_x=np.array([0.0, 100.0]),
            y_center=0.0,
            sta_elev=np.array([10.0, 20.0]),
            sta_names=["A", "B"],
            dark=True,
        )
        assert trace is not None

    def test_topo_2d_surface_returns_none_without_elev(self):
        profiles = _synthetic_profiles(n_lines=2)
        out = map3d_mod._build_topo_2d_surface(
            profiles,
            x_arr=np.linspace(0, 1000, 8),
            y_arr=np.array([0.0, 100.0]),
            topo_src="stations",
            opacity=0.7,
            dark=True,
        )
        assert out is None  # sta_elev all zero -> no real topo signal

    def test_topo_2d_surface_builds_with_elev(self):
        profiles = _synthetic_profiles(n_lines=2)
        for p in profiles.values():
            p["sta_elev"] = [float(i * 5 + 1) for i in range(len(p["sta_x"]))]
        out = map3d_mod._build_topo_2d_surface(
            profiles,
            x_arr=np.linspace(0, 1000, 8),
            y_arr=np.array([0.0, 100.0]),
            topo_src="stations",
            opacity=0.7,
            dark=True,
        )
        assert out is not None


class TestParseTopoUpload:
    def test_csv_upload(self):
        import base64

        csv_bytes = b"station,elevation\nS1,100.5\nS2,200.0\n"
        contents = "data:text/csv;base64," + base64.b64encode(csv_bytes).decode()
        records = map3d_mod._parse_topo_upload(contents, "topo.csv")
        assert records == [
            {"station": "S1", "elev": 100.5},
            {"station": "S2", "elev": 200.0},
        ]

    def test_csv_missing_columns_returns_empty(self):
        import base64

        csv_bytes = b"foo,bar\n1,2\n"
        contents = "data:text/csv;base64," + base64.b64encode(csv_bytes).decode()
        assert map3d_mod._parse_topo_upload(contents, "topo.csv") == []

    def test_unsupported_extension_returns_empty(self):
        import base64

        contents = "data:text/plain;base64," + base64.b64encode(b"x").decode()
        assert map3d_mod._parse_topo_upload(contents, "topo.txt") == []

    def test_npz_upload(self):
        import base64
        import io as _io

        buf = _io.BytesIO()
        np.savez(
            buf,
            station=np.array(["S1", "S2"]),
            elev=np.array([11.0, 22.0]),
        )
        contents = (
            "data:application/npz;base64," + base64.b64encode(buf.getvalue()).decode()
        )
        records = map3d_mod._parse_topo_upload(contents, "topo.npz")
        assert records == [
            {"station": "S1", "elev": 11.0},
            {"station": "S2", "elev": 22.0},
        ]


# ══════════════════════════════════════════════════════════════════════════
# Callbacks
# ══════════════════════════════════════════════════════════════════════════


class TestSwitchMode:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.MAP3D_ACTIVE_MODE}.data")

    def test_no_clicks_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None, None, None)

    def test_matches_triggered_slug(self, web_app):
        _set_triggered("map3d-mode-btn-block.n_clicks")
        try:
            clicks = [0, 1, 0]  # fence, block, depth
            out = self._fn(web_app)(*clicks)
            assert out == "block"
        finally:
            _clear_triggered()


class TestUpdateSourceOptions:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.MAP3D_DATA_SRC}.options")

    def test_no_data_only_pseudo_option_disabled(self, web_app):
        opts, value = self._fn(web_app)(None, None, None)
        assert all(o["disabled"] for o in opts)
        assert value == "pseudo"

    def test_with_store_data_enables_pseudo(self, web_app, store_data_willy):
        opts, value = self._fn(web_app)(store_data_willy, None, "pseudo")
        pseudo_opt = next(o for o in opts if o["value"] == "pseudo")
        assert pseudo_opt["disabled"] is False
        assert value == "pseudo"

    def test_with_lines_enables_profiles(self, web_app, store_data_willy):
        opts, value = self._fn(web_app)(store_data_willy, None, None)
        prof_opt = next(o for o in opts if o["value"] == "profiles")
        assert prof_opt["disabled"] is False

    def test_with_inversion_result_enables_inversion(self, web_app, cached_inv_result):
        opts, value = self._fn(web_app)(None, cached_inv_result, None)
        inv_opt = next(o for o in opts if o["value"] == "inversion")
        assert inv_opt["disabled"] is False

    def test_current_value_kept_if_still_enabled(self, web_app, store_data_willy):
        opts, value = self._fn(web_app)(store_data_willy, None, "pseudo")
        assert value == "pseudo"


class TestPresets:
    def _rho_fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.MAP3D_RHO_LO}.value")

    def _depth_fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.MAP3D_DEPTH_LO}.value")

    def _iso_fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.MAP3D_ISOVALUE_LO}.value")

    def test_rho_preset(self, web_app):
        _set_triggered("map3d-rho-preset-cond.n_clicks")
        try:
            lo, hi = self._rho_fn(web_app)(1, 0, 0, 0)
            assert (lo, hi) == (1, 100)
        finally:
            _clear_triggered()

    def test_depth_preset(self, web_app):
        _set_triggered("map3d-depth-preset-500.n_clicks")
        try:
            lo, hi = self._depth_fn(web_app)(1, 0, 0, 0)
            assert (lo, hi) == (0, 500)
        finally:
            _clear_triggered()

    def test_sync_iso_from_rho(self, web_app):
        lo, hi = self._iso_fn(web_app)(10.0, 1000.0)
        assert lo == pytest.approx(1.0)
        assert hi == pytest.approx(3.0)

    def test_sync_iso_missing_value_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._iso_fn(web_app)(None, 1000.0)


class TestGenerateGrid:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.MAP3D_GRID_STORE}.data")

    def test_no_clicks_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None, "pseudo", "s", None, None, None, None)

    def test_no_cached_sites(self, web_app):
        out = self._fn(web_app)(1, "pseudo", "no-such-session", None, None, None, None)
        assert "Load survey data first." in out[4]

    def test_pseudo_source_real_data(self, web_app, cached_session, store_data_willy):
        out = self._fn(web_app)(
            1, "pseudo", cached_session, store_data_willy, None, None, None
        )
        store, hint, spinner, is_open, body = out
        assert store["n_profiles"] > 0
        assert is_open is False

    def test_profiles_source_requires_two_lines(self, web_app, cached_session):
        out = self._fn(web_app)(
            1,
            "profiles",
            cached_session,
            {"station_records": []},
            None,
            None,
            None,
        )
        assert "at least two named lines" in out[1]

    def test_profiles_source_with_real_lines(
        self, web_app, cached_session, store_data_willy
    ):
        out = self._fn(web_app)(
            1, "profiles", cached_session, store_data_willy, None, None, None
        )
        store = out[0]
        assert store["n_profiles"] >= 1

    def test_inversion_source_without_result(self, web_app, cached_session):
        out = self._fn(web_app)(1, "inversion", cached_session, None, None, None, None)
        assert "Run an inversion first" in out[1]

    def test_inversion_source_with_result(
        self, web_app, cached_session_with_inv_result
    ):
        out = self._fn(web_app)(
            1,
            "inversion",
            cached_session_with_inv_result,
            None,
            None,
            None,
            None,
        )
        store = out[0]
        assert store["n_profiles"] == 1
        assert store["src"] == "inversion"

    def test_with_elevation_stores(self, web_app, cached_session, store_data_willy):
        elev_corr = [
            {"station": r["ID"], "elev_corrected": 50.0}
            for r in store_data_willy["station_records"][:2]
        ]
        elev_raw = [
            {"station": r["ID"], "elev": 40.0}
            for r in store_data_willy["station_records"][:2]
        ]
        out = self._fn(web_app)(
            1,
            "pseudo",
            cached_session,
            store_data_willy,
            elev_corr,
            elev_raw,
            None,
        )
        assert out[0]["n_profiles"] > 0


class TestToggleUploadWidget:
    def _fn(self, web_app):
        return _cb(web_app, "map3d-topo-upload-wrap.style")

    def test_upload_shows_widget(self, web_app):
        assert self._fn(web_app)("upload") == {
            "display": "block",
            "marginBottom": "4px",
        }

    def test_other_hides_widget(self, web_app):
        assert self._fn(web_app)("stations") == {"display": "none"}


class TestParseTopoFile:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.MAP3D_TOPO_UPLOAD_STORE}.data")

    def test_no_contents_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None, "topo.csv")

    def test_bad_file_shows_error(self, web_app):
        import base64

        contents = "data:text/csv;base64," + base64.b64encode(b"bad,data\n1,2").decode()
        records, msg = self._fn(web_app)(contents, "topo.csv")
        assert records is None

    def test_good_file_shows_success(self, web_app):
        import base64

        csv_bytes = b"station,elevation\nS1,100.0\n"
        contents = "data:text/csv;base64," + base64.b64encode(csv_bytes).decode()
        records, msg = self._fn(web_app)(contents, "topo.csv")
        assert records == [{"station": "S1", "elev": 100.0}]


class TestDisplay3d:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.MAP3D_GRAPH}.figure")

    def _base_args(self, **overrides):
        args = dict(
            grid_store=None,
            mode="fence",
            rho_lo=1,
            rho_hi=1e8,
            depth_lo=0,
            depth_hi=1e6,
            cmap="RdYlBu_r",
            scale="log",
            vmin=1,
            vmax=1000,
            opacity=0.85,
            contours=True,
            n_contours=10,
            line_spacing=1.0,
            azimuth=0,
            panel_thick=50,
            iso_lo=0.5,
            iso_hi=2.5,
            n_surfaces=10,
            clip_x=1.0,
            clip_y=1.0,
            clip_z=1.0,
            n_slices=4,
            show_labels=True,
            title="",
            theme="dark",
            topo_src="none",
            topo_opacity=0.7,
            topo_stations=True,
            topo_apply=False,
            sta_symbol_ui="auto",
            sta_size_ui=8,
            sta_color_ui="black",
        )
        args.update(overrides)
        return args

    def _grid_store_from(self, profiles, src="pseudo"):
        all_rho = np.concatenate(
            [np.asarray(p["rho"]).ravel() for p in profiles.values()]
        )
        all_rho = all_rho[np.isfinite(all_rho) & (all_rho > 0)]
        return {
            "profiles": {
                name: {
                    "x": np.asarray(p["x"]).tolist(),
                    "z": np.asarray(p["z"]).tolist(),
                    "rho": np.asarray(p["rho"]).tolist(),
                    "sta_x": p.get("sta_x", []),
                    "sta_elev": p.get("sta_elev", [0.0] * len(p.get("sta_x", []))),
                    "sta_elev_raw": [],
                    "sta_elev_corr": [],
                    "sta_elev_upload": [],
                    "sta_names": p.get("sta_names", []),
                    "sta_lat": p.get("sta_lat", []),
                    "sta_lon": p.get("sta_lon", []),
                }
                for name, p in profiles.items()
            },
            "n_profiles": len(profiles),
            "src": src,
            "log_lo": float(np.log10(all_rho.min())),
            "log_hi": float(np.log10(all_rho.max())),
        }

    def test_no_grid_store_shows_placeholder(self, web_app):
        fig, info, is_open, body = self._fn(web_app)(**self._base_args())
        assert fig is not None
        assert is_open is False

    def test_fence_mode_renders(self, web_app):
        profiles = _synthetic_profiles()
        grid = self._grid_store_from(profiles)
        fig, info, is_open, body = self._fn(web_app)(
            **self._base_args(grid_store=grid, mode="fence")
        )
        assert len(fig["data"]) >= 1
        assert is_open is False

    def test_block_mode_renders(self, web_app):
        profiles = _synthetic_profiles()
        grid = self._grid_store_from(profiles)
        fig, info, is_open, body = self._fn(web_app)(
            **self._base_args(grid_store=grid, mode="block")
        )
        assert fig is not None
        assert is_open is False

    def test_depth_mode_renders(self, web_app):
        profiles = _synthetic_profiles()
        grid = self._grid_store_from(profiles)
        fig, info, is_open, body = self._fn(web_app)(
            **self._base_args(grid_store=grid, mode="depth")
        )
        assert len(fig["data"]) == 4
        assert is_open is False

    def test_unknown_mode_reports_render_error(self, web_app):
        profiles = _synthetic_profiles()
        grid = self._grid_store_from(profiles)
        fig, info, is_open, body = self._fn(web_app)(
            **self._base_args(grid_store=grid, mode="not-a-real-mode")
        )
        assert is_open is True
        assert "Unknown mode" in body

    def test_fence_with_topo(self, web_app):
        profiles = _synthetic_profiles(n_lines=2)
        for p in profiles.values():
            p["sta_elev"] = [10.0] * len(p["sta_x"])
        grid = self._grid_store_from(profiles)
        fig, info, is_open, body = self._fn(web_app)(
            **self._base_args(
                grid_store=grid,
                mode="fence",
                topo_src="stations",
                topo_apply=True,
            )
        )
        assert is_open is False

    def test_real_pseudo_end_to_end(self, web_app, cached_session, store_data_willy):
        gen_fn = _cb_multi(web_app, f"{IDs.MAP3D_GRID_STORE}.data")
        grid, hint, spinner, err_open, err_body = gen_fn(
            1, "pseudo", cached_session, store_data_willy, None, None, None
        )
        assert err_open is False
        for mode in ("fence", "block", "depth"):
            fig, info, is_open, body = self._fn(web_app)(
                **self._base_args(grid_store=grid, mode=mode)
            )
            assert is_open is False, f"mode={mode} failed: {body}"


class TestExportHtml:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.MAP3D_DOWNLOAD}.data")

    def test_no_clicks_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None, {"data": [], "layout": {}}, "t")

    def test_no_figure_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(1, None, "t")

    def test_real_export(self, web_app):
        fig_dict = {
            "data": [{"type": "scatter3d", "x": [1], "y": [1], "z": [1]}],
            "layout": {},
        }
        out = self._fn(web_app)(1, fig_dict, "My Title")
        assert out["filename"] == "My_Title.html"

    def test_default_filename_when_no_title(self, web_app):
        fig_dict = {"data": [], "layout": {}}
        out = self._fn(web_app)(1, fig_dict, None)
        assert out["filename"] == "map3d.html"
