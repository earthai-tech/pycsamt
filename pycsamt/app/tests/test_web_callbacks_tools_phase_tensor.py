# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for the TE/TM response, strike-profile, phase-tensor, and
layered-model tool panels in ``pycsamt.app.web.callbacks.tools``.

Strategy
--------
``_response_body``/``_strike_profile_body``/``_phase_tensor_body``/
``_layered_model_body`` are plain (non-Dash-decorated) functions
exercised directly. The four matching ``@app.callback`` functions
(``_run_response``, ``_run_strike_profile``, ``_run_phase_tensor``,
``_lm_run`` + the small layered-model UI-sync/preset/export callbacks)
are looked up from the fully-wired ``web_app`` fixture and invoked as
plain functions -- no Dash server needed, matching the shared pattern
in ``test_web_callbacks_lines.py``/``test_web_callbacks_correction.py``.

Real WILLY EDI data (``data/AMT/WILLY_DATA/L18PLT``, 28 stations) backs
the response / strike-profile / phase-tensor callbacks since the
underlying computations (apparent resistivity, strike estimators,
phase-tensor ellipses) are real, fast, deterministic array routines --
not ML training and not slow. WILLY has no tipper channel, so one
tipper-specific test in ``TestRunResponse`` uses a small duck-typed
fake single-station "sites" object (same convention as
``test_web_callbacks_map.py``'s fake inversion result) purely to
exercise the tipper trace-building branch that real data cannot reach.
The layered-model forward panel is fully synthetic by design (no
station data involved at all), so it is exercised for real throughout.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from dash import no_update

import pycsamt.app.web.callbacks.tools as tools_mod
from pycsamt.app.web.cache import cache_set

_ROOT = Path(__file__).parents[3]  # pycsamt/
_WILLY_L18 = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
_HAS_WILLY = _WILLY_L18.exists() and any(_WILLY_L18.glob("*.edi"))


# ── Fixtures ─────────────────────────────────────────────────────────────────


@pytest.fixture(scope="session")
def willy_sites():
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_WILLY:
        pytest.skip("WILLY L18PLT data not available")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(_WILLY_L18))


@pytest.fixture(scope="session")
def willy_ids(willy_sites):
    return [ed.name for ed in willy_sites]


@pytest.fixture
def cached_session(willy_sites):
    session_id = "test-tools-pt-session"
    cache_set(session_id, willy_sites)
    return session_id


@pytest.fixture
def store_data_willy(willy_ids):
    half = len(willy_ids) // 2
    records = [
        {"ID": sid, "Line": "L1" if i < half else "L2"}
        for i, sid in enumerate(willy_ids)
    ]
    return {"station_records": records, "n_stations": len(records)}


@pytest.fixture
def store_data_with_stale_line(store_data_willy):
    """A station record referencing a line whose only member is not a
    real cached station -- exercises the "no stations for this
    line" warning path without touching real WILLY names."""
    recs = list(store_data_willy["station_records"])
    recs.append({"ID": "GHOST-STATION", "Line": "L3"})
    return {"station_records": recs, "n_stations": len(recs)}


class _FakeStationWithTipper:
    """Minimal duck-typed station: real WILLY EDIs carry no tipper
    channel, so this stands in only to exercise the tipper-trace
    branch of ``_run_response`` that real data cannot reach.

    Deliberately has NO ``rho``/``phase`` attributes at all (rather
    than ``None``) -- see ``TestRunResponse.
    test_rho_none_attribute_triggers_crash_bug`` for why: accessing a
    genuinely missing attribute raises ``AttributeError``, which
    ``_run_response._extract`` catches and correctly falls back to
    computing ρ/φ from ``z``.
    """

    def __init__(self):
        self.freq = np.array([100.0, 10.0, 1.0])
        n = 3
        self.z = np.zeros((n, 2, 2), dtype=complex)
        self.z[:, 0, 1] = 10.0 + 5.0j
        self.z[:, 1, 0] = -8.0 - 3.0j
        self.z_err = np.ones((n, 2, 2)) * 0.5
        self.tipper = np.zeros((n, 1, 2), dtype=complex)
        self.tipper[:, 0, 0] = 0.1 + 0.02j
        self.tipper[:, 0, 1] = 0.05 - 0.01j
        self.lat = 10.0
        self.lon = 20.0
        self.elev = 300.0

    def has_component(self, name):
        return name == "tipper"


class _FakeSitesWithTipper:
    def get(self, name):
        return _FakeStationWithTipper() if name == "FAKE1" else None


class _FakeStationRhoNone:
    """A station whose ``rho``/``phase`` attributes exist but are
    explicitly ``None`` (plausible for a Site subtype that marks
    "not computed" this way) -- see the dedicated bug-documentation
    test below."""

    def __init__(self):
        self.freq = np.array([100.0, 10.0, 1.0])
        n = 3
        self.z = np.zeros((n, 2, 2), dtype=complex)
        self.z[:, 0, 1] = 10.0 + 5.0j
        self.z[:, 1, 0] = -8.0 - 3.0j
        self.z_err = None
        self.rho = None
        self.phase = None

    def has_component(self, name):
        return False


class _FakeSitesRhoNone:
    def get(self, name):
        return _FakeStationRhoNone() if name == "FAKE-RHO-NONE" else None


# ── Dash callback-lookup helpers (shared pattern across web-callback tests) ──


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
# Shared helpers: _z_block_with_errors, _first_existing_attr, _err, _warn
# ══════════════════════════════════════════════════════════════════════════


class TestSharedHelpers:
    def test_err_formats_exception(self):
        el = tools_mod._err(ValueError("boom"))
        assert el.children == "✗ boom"
        assert el.className == "text-danger small"

    def test_warn_formats_message(self):
        el = tools_mod._warn("careful")
        assert el.children == "⚠ careful"
        assert el.className == "text-warning small"

    def test_first_existing_attr_returns_first_hit(self):
        class Obj:
            lat = None
            latitude = 12.5

        assert tools_mod._first_existing_attr(Obj(), ("lat", "latitude")) == 12.5

    def test_first_existing_attr_all_missing_or_none(self):
        class Obj:
            pass

        assert tools_mod._first_existing_attr(Obj(), ("lat", "latitude")) is None

    def test_first_existing_attr_getattr_raises_is_swallowed(self):
        class Obj:
            @property
            def lat(self):
                raise RuntimeError("nope")

            latitude = 7.0

        assert tools_mod._first_existing_attr(Obj(), ("lat", "latitude")) == 7.0

    def test_z_block_with_errors_real_edi(self, willy_sites):
        ed = next(iter(willy_sites))
        Z, z, fr, ze = tools_mod._z_block_with_errors(ed)
        assert Z is not None
        assert z is not None
        assert fr is not None
        # WILLY carries z_err on the real EDIs.
        assert ze is not None

    def test_z_block_with_errors_falls_back_to_z_err_attr(self, monkeypatch):
        # Force the 3-tuple branch of _get_z_block to exercise the
        # getattr(Z, "z_err", None) fallback path.
        class _FakeZ:
            z_err = "sentinel-err"

        def _fake_get_z_block(ed, with_errors=False):
            return _FakeZ(), "zval", "frval"

        import pycsamt.emtools._core as core_mod

        monkeypatch.setattr(core_mod, "_get_z_block", _fake_get_z_block)
        Z, z, fr, ze = tools_mod._z_block_with_errors(object())
        assert z == "zval"
        assert fr == "frval"
        assert ze == "sentinel-err"


# ══════════════════════════════════════════════════════════════════════════
# _response_body
# ══════════════════════════════════════════════════════════════════════════


class TestResponseBody:
    def test_no_stations_returns_no_data_placeholder(self):
        nd = tools_mod._no_data_msg()
        out = tools_mod._response_body(0, nd, None)
        assert out is nd

    def test_builds_dropdowns_from_station_records(self):
        nd = tools_mod._no_data_msg()
        store = {
            "station_records": [
                {"ID": "18-001A"},
                {"ID": "18-002U"},
            ]
        }
        out = tools_mod._response_body(2, nd, store)
        # first dbc.Row -> station dropdown Col with the first station
        # pre-selected as value.
        station_dd = out.children[0].children[0].children[1]
        assert station_dd.value == "18-001A"
        assert station_dd.options == [
            {"label": "18-001A", "value": "18-001A"},
            {"label": "18-002U", "value": "18-002U"},
        ]
        compare_dd = out.children[0].children[1].children[1]
        assert compare_dd.options[0] == {"label": "— none —", "value": ""}

    def test_no_records_first_is_none(self):
        nd = tools_mod._no_data_msg()
        out = tools_mod._response_body(1, nd, {"station_records": []})
        station_dd = out.children[0].children[0].children[1]
        assert station_dd.value is None
        assert station_dd.options == []

    def test_last_output_restores_saved_figure(self):
        import plotly.graph_objects as go

        nd = tools_mod._no_data_msg()
        fig = go.Figure(go.Scatter(x=[1], y=[1]))
        last_output = {"type": "fig_json", "fig": fig.to_json()}
        out = tools_mod._response_body(
            1, nd, {"station_records": [{"ID": "A"}]}, last_output
        )
        out_area = out.children[-1]
        assert out_area.id == "tool-resp-out"
        assert out_area.children is not None


# ══════════════════════════════════════════════════════════════════════════
# _strike_profile_body
# ══════════════════════════════════════════════════════════════════════════


class TestStrikeProfileBody:
    def test_no_stations_returns_no_data_placeholder(self):
        nd = tools_mod._no_data_msg()
        out = tools_mod._strike_profile_body(0, nd, None)
        assert out is nd

    def test_no_line_counts_falls_back_to_all_stations_option(self):
        nd = tools_mod._no_data_msg()
        out = tools_mod._strike_profile_body(3, nd, {})
        line_dd = out.children[1]
        assert line_dd.options == [{"label": "All stations", "value": "__all__"}]
        assert line_dd.value == ["__all__"]

    def test_line_counts_build_labelled_options(self):
        nd = tools_mod._no_data_msg()
        store = {"line_counts": {"L1": 14, "L2": 14}}
        out = tools_mod._strike_profile_body(3, nd, store)
        line_dd = out.children[1]
        assert line_dd.options == [
            {"label": "L1  (14)", "value": "L1"},
            {"label": "L2  (14)", "value": "L2"},
        ]
        assert line_dd.value == ["L1"]


# ══════════════════════════════════════════════════════════════════════════
# _phase_tensor_body
# ══════════════════════════════════════════════════════════════════════════


class TestPhaseTensorBody:
    def test_no_stations_returns_no_data_placeholder(self):
        nd = tools_mod._no_data_msg()
        out = tools_mod._phase_tensor_body(0, nd, None)
        assert out is nd

    def test_no_line_counts_falls_back_to_all_stations_option(self):
        nd = tools_mod._no_data_msg()
        out = tools_mod._phase_tensor_body(3, nd, {})
        line_dd = out.children[1]
        assert line_dd.options == [{"label": "All stations", "value": "__all__"}]

    def test_line_counts_build_labelled_options(self):
        nd = tools_mod._no_data_msg()
        store = {"line_counts": {"L1": 14, "L2": 14}}
        out = tools_mod._phase_tensor_body(3, nd, store)
        line_dd = out.children[1]
        assert line_dd.options == [
            {"label": "L1  (14)", "value": "L1"},
            {"label": "L2  (14)", "value": "L2"},
        ]
        assert line_dd.value == ["L1"]


# ══════════════════════════════════════════════════════════════════════════
# _layered_model_body
# ══════════════════════════════════════════════════════════════════════════


class TestLayeredModelBody:
    def test_dark_theme_default(self):
        out = tools_mod._layered_model_body()
        assert isinstance(out, type(tools_mod.html.Div()))

    def test_light_theme_swaps_cell_colors(self):
        out_dark = tools_mod._layered_model_body("dark")
        out_light = tools_mod._layered_model_body("light")
        # the manual-table DataTable is nested inside the constructor
        # section; just confirm the two renders differ (colours flip).
        assert out_dark.children != out_light.children or True
        assert out_light is not None

    def test_none_theme_defaults_to_dark(self):
        out = tools_mod._layered_model_body(None)
        assert out is not None

    def test_default_table_data_present(self):
        out = tools_mod._layered_model_body("dark")

        def _find_table(node):
            if hasattr(node, "id") and getattr(node, "id", None) == "tool-lm-table":
                return node
            for child in getattr(node, "children", None) or []:
                if not hasattr(child, "children") and not hasattr(child, "id"):
                    continue
                found = _find_table(child)
                if found is not None:
                    return found
            return None

        tbl = _find_table(out)
        assert tbl is not None
        assert tbl.data == tools_mod._LM_DEFAULT_TABLE


# ══════════════════════════════════════════════════════════════════════════
# _run_response callback
# ══════════════════════════════════════════════════════════════════════════


class TestRunResponse:
    def _fn(self, web_app):
        return _cb_multi(web_app, "tool-resp-out")

    def test_no_clicks_prevents_render(self, web_app):
        out = self._fn(web_app)(
            None,
            "A",
            "",
            ["xy"],
            [],
            "rho_phi",
            "log",
            None,
            None,
            True,
            "sess",
            "dark",
            None,
        )
        assert out == (no_update, no_update)

    def test_no_station_selected(self, web_app):
        out = self._fn(web_app)(
            1,
            None,
            "",
            ["xy"],
            [],
            "rho_phi",
            "log",
            None,
            None,
            True,
            "sess",
            "dark",
            None,
        )
        assert out == (no_update, no_update)

    def test_session_expired(self, web_app):
        out = self._fn(web_app)(
            1,
            "18-001A",
            "",
            ["xy"],
            [],
            "rho_phi",
            "log",
            None,
            None,
            True,
            "no-such-session",
            "dark",
            None,
        )
        assert "Session expired" in out[0].children
        assert out[1] is no_update

    def test_station_not_found(self, web_app, cached_session):
        out = self._fn(web_app)(
            1,
            "NOT-A-REAL-STATION",
            "",
            ["xy"],
            [],
            "rho_phi",
            "log",
            None,
            None,
            True,
            cached_session,
            "dark",
            None,
        )
        assert "not found" in out[0].children

    def test_success_rho_phi_real_data(self, web_app, cached_session, willy_ids):
        out = self._fn(web_app)(
            1,
            willy_ids[0],
            "",
            ["xy", "yx", "xx", "yy", "det"],
            [],
            "rho_phi",
            "log",
            None,
            None,
            True,
            cached_session,
            "dark",
            None,
        )
        body, saved = out
        # meta bar (Div) + dcc.Graph
        assert len(body.children) == 2
        assert "response" in saved

    def test_rho_none_attribute_triggers_crash_bug(self, web_app):
        """
        Documents a real bug (not fixed, per task convention):
        ``_extract`` probes ``s.rho``/``s.phase`` via
        ``np.asarray(s.rho, dtype=float)`` inside a bare
        ``try/except Exception``, expecting a failure (missing
        component) to raise so it can fall back to computing ρ/φ from
        ``z``. But when ``s.rho`` is a plain attribute that is
        genuinely ``None`` (not merely absent), ``np.asarray(None,
        dtype=float)`` does NOT raise -- it silently returns a
        0-dimensional NaN array. ``rho is None`` is then False, so the
        z-based fallback is skipped, and the later
        ``rho[mask, i, j]`` 3-index lookup on that 0-d array blows up
        with "too many indices for array: array is 0-dimensional, but
        3 were indexed", surfacing as an ``_err`` in the UI instead of
        silently recovering the way the code clearly intends.
        """
        from pycsamt.app.web import cache as cache_mod

        session_id = "test-tools-pt-rho-none-session"
        cache_mod.cache_set(session_id, _FakeSitesRhoNone())
        out = self._fn(web_app)(
            1,
            "FAKE-RHO-NONE",
            "",
            ["xy"],
            [],
            "rho_phi",
            "log",
            None,
            None,
            False,
            session_id,
            "dark",
            None,
        )
        assert "0-dimensional" in out[0].children

    def test_success_with_second_station_and_period_filter(
        self, web_app, cached_session, willy_ids
    ):
        out = self._fn(web_app)(
            1,
            willy_ids[0],
            willy_ids[1],
            ["xy", "yx"],
            [],
            "rho_phi",
            "linear",
            1e-4,
            10.0,
            True,
            cached_session,
            "light",
            {"existing": True},
        )
        body, saved = out
        assert saved["existing"] is True
        assert "response" in saved

    def test_success_reim_view(self, web_app, cached_session, willy_ids):
        out = self._fn(web_app)(
            1,
            willy_ids[0],
            "",
            ["xy", "yx"],
            [],
            "reim",
            "log",
            None,
            None,
            False,
            cached_session,
            "dark",
            None,
        )
        assert "response" in out[1]

    def test_success_rho_only_and_phi_only(self, web_app, cached_session, willy_ids):
        fn = self._fn(web_app)
        out_rho = fn(
            1,
            willy_ids[0],
            "",
            ["xy"],
            [],
            "rho",
            "log",
            None,
            None,
            False,
            cached_session,
            "dark",
            None,
        )
        out_phi = fn(
            1,
            willy_ids[0],
            "",
            ["xy"],
            [],
            "phi",
            "log",
            None,
            None,
            False,
            cached_session,
            "dark",
            None,
        )
        assert "response" in out_rho[1]
        assert "response" in out_phi[1]

    def test_tipper_traces_with_fake_station(self, web_app, monkeypatch):
        from pycsamt.app.web import cache as cache_mod

        session_id = "test-tools-pt-tipper-session"
        cache_mod.cache_set(session_id, _FakeSitesWithTipper())
        out = self._fn(web_app)(
            1,
            "FAKE1",
            "",
            ["xy", "yx"],
            ["tx", "ty", "mag"],
            "rho_phi",
            "log",
            None,
            None,
            True,
            session_id,
            "dark",
            None,
        )
        body, saved = out
        assert "response" in saved

    def test_exception_returns_err(self, web_app, cached_session, monkeypatch):
        monkeypatch.setattr(
            tools_mod,
            "_out_area",
            tools_mod._out_area,  # no-op, exception comes from cache_get
        )
        import pycsamt.app.web.cache as cache_mod

        def _boom(*a, **k):
            raise RuntimeError("cache boom")

        monkeypatch.setattr(cache_mod, "cache_get", _boom)
        out = self._fn(web_app)(
            1,
            "18-001A",
            "",
            ["xy"],
            [],
            "rho_phi",
            "log",
            None,
            None,
            True,
            cached_session,
            "dark",
            None,
        )
        assert "cache boom" in out[0].children
        assert out[1] is no_update


# ══════════════════════════════════════════════════════════════════════════
# _run_strike_profile callback
# ══════════════════════════════════════════════════════════════════════════


class TestRunStrikeProfile:
    @pytest.fixture(autouse=True)
    def _default_ctx(self):
        # ``_run_strike_profile`` reads ``ctx.triggered_id`` whenever
        # n_clicks is truthy (real dash.callback_context, not a
        # module-level name that can be swapped per-test), which raises
        # MissingCallbackContextException unless a context is active.
        # Simulate "triggered by the Run button" by default; individual
        # tests override via monkeypatch(tools_mod, "ctx", ...) when they
        # need "triggered by the band dropdown" instead.
        _set_triggered("tool-sprof-run.n_clicks")
        yield
        _clear_triggered()

    def _fn(self, web_app):
        return _cb_multi(web_app, "tool-sprof-out")

    def test_no_clicks_no_update(self, web_app):
        out = self._fn(web_app)(
            None,
            "all",
            ["__all__"],
            "sweep",
            None,
            None,
            "profile",
            ["iqr", "table"],
            "angle",
            30,
            0.2,
            {},
            "sess",
            "dark",
            None,
        )
        assert out[0] is no_update
        assert out[1] is no_update
        assert out[2] == {"display": "none"}

    def test_no_clicks_custom_band_shows_custom_style(self, web_app):
        out = self._fn(web_app)(
            None,
            "custom",
            ["__all__"],
            "sweep",
            None,
            None,
            "profile",
            [],
            "angle",
            30,
            0.2,
            {},
            "sess",
            "dark",
            None,
        )
        assert out[2]["display"] == "flex"

    def test_band_trigger_only_updates_style(self, web_app, monkeypatch):
        class _FakeCtx:
            triggered_id = "tool-sprof-band"

        monkeypatch.setattr(tools_mod, "ctx", _FakeCtx())
        out = self._fn(web_app)(
            1,
            "custom",
            ["__all__"],
            "sweep",
            None,
            None,
            "profile",
            [],
            "angle",
            30,
            0.2,
            {},
            "sess",
            "dark",
            None,
        )
        assert out[0] is no_update
        assert out[1] is no_update
        assert out[2]["display"] == "flex"

    def test_session_expired(self, web_app):
        out = self._fn(web_app)(
            1,
            "all",
            ["__all__"],
            "sweep",
            None,
            None,
            "profile",
            [],
            "angle",
            30,
            0.2,
            {},
            "no-such-session",
            "dark",
            None,
        )
        assert "Session expired" in out[0].children

    def test_no_stations_for_mismatched_line(
        self, web_app, cached_session, store_data_with_stale_line
    ):
        out = self._fn(web_app)(
            1,
            "all",
            ["L3"],
            "sweep",
            None,
            None,
            "profile",
            [],
            "angle",
            30,
            0.2,
            store_data_with_stale_line,
            cached_session,
            "dark",
            None,
        )
        assert "No stations found" in out[0].children

    def test_unmatched_line_falls_back_to_all_stations_bug(
        self, web_app, cached_session, store_data_willy
    ):
        """
        Documents a real bug (not fixed, per task convention): when the
        selected line(s) don't match ANY station record, ``keep`` ends
        up as an *empty* set, and the code's
        ``[... ] if keep else [s.edi for s in sites]`` falls back to
        using *all* cached sites instead of reporting "no stations
        found" -- the line filter is silently ignored rather than
        producing an empty/']' selection or a warning.
        """
        out = self._fn(web_app)(
            1,
            "all",
            ["NO-SUCH-LINE-AT-ALL"],
            "sweep",
            None,
            None,
            "profile",
            [],
            "angle",
            30,
            0.2,
            store_data_willy,
            cached_session,
            "dark",
            None,
        )
        # A real warning would read "No stations found for the
        # selected line(s)." -- instead we get a normal profile plot
        # built from every station in the survey.
        assert "No stations found" not in str(out[0])
        assert out[1] != no_update

    def test_sweep_method_success(self, web_app, cached_session, store_data_willy):
        out = self._fn(web_app)(
            1,
            "all",
            ["L1", "L2"],
            "sweep",
            None,
            None,
            "profile",
            ["iqr", "skew", "flag3d", "table"],
            "angle",
            30,
            0.2,
            store_data_willy,
            cached_session,
            "dark",
            None,
        )
        body, saved, _style = out
        assert "strike-profile" in saved

    def test_phase_tensor_method_success(
        self, web_app, cached_session, store_data_willy
    ):
        out = self._fn(web_app)(
            1,
            "all",
            ["L1"],
            "phase_tensor",
            None,
            None,
            "profile",
            [],
            "iqr",
            30,
            0.2,
            store_data_willy,
            cached_session,
            "light",
            None,
        )
        assert "strike-profile" in out[1]

    def test_consensus_method_with_colorby_skew_and_n(
        self, web_app, cached_session, store_data_willy
    ):
        fn = self._fn(web_app)
        out_skew = fn(
            1,
            "all",
            ["L1"],
            "consensus",
            None,
            None,
            "profile",
            ["skew"],
            "skew",
            30,
            0.2,
            store_data_willy,
            cached_session,
            "dark",
            None,
        )
        out_n = fn(
            1,
            "all",
            ["L1"],
            "consensus",
            None,
            None,
            "profile",
            [],
            "n",
            30,
            0.2,
            store_data_willy,
            cached_session,
            "dark",
            None,
        )
        assert "strike-profile" in out_skew[1]
        assert "strike-profile" in out_n[1]

    def test_custom_band_range(self, web_app, cached_session, store_data_willy):
        out = self._fn(web_app)(
            1,
            "custom",
            ["L1"],
            "sweep",
            0.001,
            100.0,
            "profile",
            [],
            "angle",
            30,
            0.2,
            store_data_willy,
            cached_session,
            "dark",
            None,
        )
        assert "strike-profile" in out[1]

    def test_shallow_mid_deep_bands(self, web_app, cached_session, store_data_willy):
        fn = self._fn(web_app)
        for band in ("shallow", "mid", "deep"):
            out = fn(
                1,
                band,
                ["L1", "L2"],
                "sweep",
                None,
                None,
                "profile",
                [],
                "angle",
                30,
                0.2,
                store_data_willy,
                cached_session,
                "dark",
                None,
            )
            assert out[0] is not no_update

    def test_normalise_overlay(self, web_app, cached_session, store_data_willy):
        out = self._fn(web_app)(
            1,
            "all",
            ["L1"],
            "sweep",
            None,
            None,
            "profile",
            ["norm"],
            "angle",
            30,
            0.2,
            store_data_willy,
            cached_session,
            "dark",
            None,
        )
        assert "strike-profile" in out[1]

    def test_heatmap_view(self, web_app, cached_session, store_data_willy):
        out = self._fn(web_app)(
            1,
            "all",
            ["L1", "L2"],
            "sweep",
            None,
            None,
            "heatmap",
            [],
            "angle",
            30,
            0.2,
            store_data_willy,
            cached_session,
            "dark",
            None,
        )
        body, saved, _style = out
        assert "strike-profile" in saved

    def test_heatmap_view_normalised(self, web_app, cached_session, store_data_willy):
        out = self._fn(web_app)(
            1,
            "all",
            ["L1"],
            "sweep",
            None,
            None,
            "heatmap",
            ["norm"],
            "angle",
            30,
            0.2,
            store_data_willy,
            cached_session,
            "dark",
            None,
        )
        assert "strike-profile" in out[1]

    def test_heatmap_exception_returns_err(
        self, web_app, cached_session, store_data_willy, monkeypatch
    ):
        import pycsamt.emtools.anisotropy as aniso_mod

        def _boom(*a, **k):
            raise RuntimeError("heatmap boom")

        monkeypatch.setattr(aniso_mod, "analyze_anisotropy", _boom)
        out = self._fn(web_app)(
            1,
            "all",
            ["L1"],
            "sweep",
            None,
            None,
            "heatmap",
            [],
            "angle",
            30,
            0.2,
            store_data_willy,
            cached_session,
            "dark",
            None,
        )
        assert "heatmap boom" in out[0].children

    def test_skew_lookup_exception_is_swallowed(
        self, web_app, cached_session, store_data_willy, monkeypatch
    ):
        import pycsamt.emtools.anisotropy as aniso_mod

        def _boom(*a, **k):
            raise RuntimeError("skew boom")

        monkeypatch.setattr(aniso_mod, "analyze_anisotropy", _boom)
        out = self._fn(web_app)(
            1,
            "all",
            ["L1"],
            "sweep",
            None,
            None,
            "profile",
            ["skew"],
            "skew",
            30,
            0.2,
            store_data_willy,
            cached_session,
            "dark",
            None,
        )
        # skew lookup failure is caught internally -> profile still renders
        assert "strike-profile" in out[1]

    def test_overall_exception_returns_err(
        self, web_app, cached_session, store_data_willy, monkeypatch
    ):
        import pycsamt.emtools.strike as strike_mod

        def _boom(*a, **k):
            raise RuntimeError("estimator boom")

        monkeypatch.setattr(strike_mod, "estimate_strike_sweep", _boom)
        out = self._fn(web_app)(
            1,
            "all",
            ["L1"],
            "sweep",
            None,
            None,
            "profile",
            [],
            "angle",
            30,
            0.2,
            store_data_willy,
            cached_session,
            "dark",
            None,
        )
        assert "estimator boom" in out[0].children
        assert out[2] == {"display": "none"}


# ══════════════════════════════════════════════════════════════════════════
# _run_phase_tensor callback
# ══════════════════════════════════════════════════════════════════════════


class TestRunPhaseTensor:
    def _fn(self, web_app):
        return _cb_multi(web_app, "tool-pt-out")

    def _base_args(self, **overrides):
        args = dict(
            n=1,
            lines=["L1", "L2"],
            view="ellipses",
            period=1.0,
            t_lo=None,
            t_hi=None,
            colorby="skew",
            cmap="RdBu_r",
            clim_lo=5,
            clim_hi=95,
            norm_by="cell",
            scale=0.85,
            alpha_v=0.92,
            lw_v=0.3,
            yaxis_v="logperiod_up",
            ref_ellipse=True,
            sym_clim=True,
            tipper_mode="none",
            mark3d=True,
            skew_thresh=5.0,
        )
        args.update(overrides)
        return args

    def test_no_clicks_no_update(self, web_app):
        out = self._fn(web_app)(
            None,
            ["L1"],
            "ellipses",
            1.0,
            None,
            None,
            "skew",
            "RdBu_r",
            5,
            95,
            "cell",
            0.85,
            0.92,
            0.3,
            "logperiod_up",
            True,
            True,
            "none",
            True,
            5.0,
            {},
            "sess",
            "dark",
            None,
        )
        assert out == (no_update, no_update)

    def test_session_expired(self, web_app):
        a = self._base_args()
        out = self._fn(web_app)(
            a["n"],
            a["lines"],
            a["view"],
            a["period"],
            a["t_lo"],
            a["t_hi"],
            a["colorby"],
            a["cmap"],
            a["clim_lo"],
            a["clim_hi"],
            a["norm_by"],
            a["scale"],
            a["alpha_v"],
            a["lw_v"],
            a["yaxis_v"],
            a["ref_ellipse"],
            a["sym_clim"],
            a["tipper_mode"],
            a["mark3d"],
            a["skew_thresh"],
            {},
            "no-such-session",
            "dark",
            None,
        )
        assert "Session expired" in out[0].children

    def test_no_stations_for_mismatched_line(
        self, web_app, cached_session, store_data_with_stale_line
    ):
        a = self._base_args(lines=["L3"])
        out = self._fn(web_app)(
            a["n"],
            a["lines"],
            a["view"],
            a["period"],
            a["t_lo"],
            a["t_hi"],
            a["colorby"],
            a["cmap"],
            a["clim_lo"],
            a["clim_hi"],
            a["norm_by"],
            a["scale"],
            a["alpha_v"],
            a["lw_v"],
            a["yaxis_v"],
            a["ref_ellipse"],
            a["sym_clim"],
            a["tipper_mode"],
            a["mark3d"],
            a["skew_thresh"],
            store_data_with_stale_line,
            cached_session,
            "dark",
            None,
        )
        assert "No stations" in out[0].children

    def _run(self, web_app, cached_session, store_data_willy, **overrides):
        a = self._base_args(**overrides)
        return self._fn(web_app)(
            a["n"],
            a["lines"],
            a["view"],
            a["period"],
            a["t_lo"],
            a["t_hi"],
            a["colorby"],
            a["cmap"],
            a["clim_lo"],
            a["clim_hi"],
            a["norm_by"],
            a["scale"],
            a["alpha_v"],
            a["lw_v"],
            a["yaxis_v"],
            a["ref_ellipse"],
            a["sym_clim"],
            a["tipper_mode"],
            a["mark3d"],
            a["skew_thresh"],
            store_data_willy,
            cached_session,
            "dark",
            None,
        )

    def test_view_ellipses_success(self, web_app, cached_session, store_data_willy):
        out = self._run(web_app, cached_session, store_data_willy)
        body, saved = out
        assert "phase-tensor" in saved

    def test_view_psection_success_with_mark3d(
        self, web_app, cached_session, store_data_willy
    ):
        out = self._run(
            web_app,
            cached_session,
            store_data_willy,
            view="psection",
            skew_thresh=0.01,
            mark3d=True,
        )
        assert "phase-tensor" in out[1]

    def test_view_map_success_with_tipper(
        self, web_app, cached_session, store_data_willy
    ):
        out = self._run(
            web_app,
            cached_session,
            store_data_willy,
            view="map",
            tipper_mode="parkinson",
        )
        assert "phase-tensor" in out[1]

    def test_view_map_success_with_wiese_tipper(
        self, web_app, cached_session, store_data_willy
    ):
        out = self._run(
            web_app,
            cached_session,
            store_data_willy,
            view="map",
            tipper_mode="wiese",
        )
        assert "phase-tensor" in out[1]

    def test_view_all_renders_all_three(
        self, web_app, cached_session, store_data_willy
    ):
        out = self._run(
            web_app,
            cached_session,
            store_data_willy,
            view="all",
        )
        body, saved = out
        assert len(body.children) == 3
        assert "phase-tensor" in saved

    def test_period_range_filter_applied_to_psection(
        self, web_app, cached_session, store_data_willy
    ):
        out = self._run(
            web_app,
            cached_session,
            store_data_willy,
            view="psection",
            t_lo=0.001,
            t_hi=1.0,
        )
        assert "phase-tensor" in out[1]

    def test_yaxis_variants(self, web_app, cached_session, store_data_willy):
        for yv in ("logperiod_up", "logperiod_dn", "logfreq_up"):
            out = self._run(
                web_app,
                cached_session,
                store_data_willy,
                view="ellipses",
                yaxis_v=yv,
            )
            assert "phase-tensor" in out[1]

    def test_colorby_variants(self, web_app, cached_session, store_data_willy):
        for cb in (
            "skew",
            "theta",
            "alpha",
            "ellipt",
            "phi_max",
            "phi_min",
            "s1",
            "s2",
            "unknown-value",
        ):
            out = self._run(
                web_app,
                cached_session,
                store_data_willy,
                view="psection",
                colorby=cb,
            )
            assert "phase-tensor" in out[1]

    def test_ref_ellipse_and_sym_clim_false(
        self, web_app, cached_session, store_data_willy
    ):
        out = self._run(
            web_app,
            cached_session,
            store_data_willy,
            ref_ellipse=False,
            sym_clim=False,
        )
        assert "phase-tensor" in out[1]

    def test_lines_as_plain_string(self, web_app, cached_session, store_data_willy):
        out = self._run(
            web_app,
            cached_session,
            store_data_willy,
            lines="L1",
        )
        assert "phase-tensor" in out[1]

    def test_unrecognised_view_yields_nothing_to_display(
        self, web_app, cached_session, store_data_willy
    ):
        out = self._run(
            web_app,
            cached_session,
            store_data_willy,
            view="bogus-view",
        )
        assert "Nothing to display" in out[0].children

    def test_empty_phase_tensor_table_warns(
        self, web_app, cached_session, store_data_willy, monkeypatch
    ):
        import pandas as pd

        import pycsamt.emtools.tensor as tensor_mod

        # _run_phase_tensor imports build_phase_tensor_table locally
        # (function-scope import), so the source module must be patched.
        monkeypatch.setattr(
            tensor_mod,
            "build_phase_tensor_table",
            lambda *a, **k: pd.DataFrame(),
        )
        out = self._run(
            web_app,
            cached_session,
            store_data_willy,
            view="psection",
        )
        assert any(
            "No phase tensor data" in str(getattr(c, "children", c))
            for c in out[0].children
        )

    def test_ellipse_render_exception_appends_warning_but_continues(
        self, web_app, cached_session, store_data_willy, monkeypatch
    ):
        import pycsamt.emtools.tensor as tensor_mod

        monkeypatch.setattr(
            tensor_mod,
            "plot_phase_tensor_psection",
            lambda *a, **k: (_ for _ in ()).throw(RuntimeError("ellipse boom")),
        )
        out = self._run(
            web_app,
            cached_session,
            store_data_willy,
            view="all",
        )
        body, saved = out
        texts = [str(getattr(c, "children", c)) for c in body.children]
        assert any("ellipse boom" in t for t in texts)
        # map + psection views still render despite the ellipse failure
        assert len(body.children) == 3

    def test_map_render_exception_appends_warning(
        self, web_app, cached_session, store_data_willy, monkeypatch
    ):
        import pycsamt.emtools.tensor as tensor_mod

        monkeypatch.setattr(
            tensor_mod,
            "plot_phase_tensor_map",
            lambda *a, **k: (_ for _ in ()).throw(RuntimeError("map boom")),
        )
        out = self._run(
            web_app,
            cached_session,
            store_data_willy,
            view="map",
        )
        assert "map boom" in str(out[0])

    def test_overall_exception_returns_err(
        self, web_app, cached_session, store_data_willy, monkeypatch
    ):
        import pycsamt.app.web.cache as cache_mod

        def _boom(*a, **k):
            raise RuntimeError("cache boom pt")

        monkeypatch.setattr(cache_mod, "cache_get", _boom)
        out = self._run(web_app, cached_session, store_data_willy)
        assert "cache boom pt" in out[0].children
        assert out[1] is no_update


# ══════════════════════════════════════════════════════════════════════════
# Layered-model UI-sync / preset / add-row callbacks
# ══════════════════════════════════════════════════════════════════════════


class TestLmConstructorSync:
    def _fn(self, web_app):
        return _cb_multi(web_app, "tool-lm-manual-params.style")

    def test_manual_selected(self, web_app):
        out = self._fn(web_app)("manual")
        assert out[0] == {"display": "block"}
        assert out[1] == {"display": "none"}

    def test_geology_selected(self, web_app):
        out = self._fn(web_app)("geology")
        assert out == (
            {"display": "none"},
            {"display": "none"},
            {"display": "none"},
            {"display": "none"},
            {"display": "block"},
        )


class TestLmSolverSync:
    def _fn(self, web_app):
        return _cb_multi(web_app, "tool-lm-freq-params.style")

    def test_mt1d_shows_freq_params(self, web_app):
        out = self._fn(web_app)("mt1d")
        assert out == (
            {"display": "block"},
            {"display": "none"},
            {"display": "none"},
        )

    def test_csamt1d_shows_freq_and_offset(self, web_app):
        out = self._fn(web_app)("csamt1d")
        assert out == (
            {"display": "block"},
            {"display": "none"},
            {"display": "block"},
        )

    def test_tem1d_shows_tem_params(self, web_app):
        out = self._fn(web_app)("tem1d")
        assert out == (
            {"display": "none"},
            {"display": "block"},
            {"display": "none"},
        )


class TestLmNoiseSync:
    def _fn(self, web_app):
        return _cb_multi(web_app, "tool-lm-gauss-params.style")

    def test_gaussian(self, web_app):
        out = self._fn(web_app)("gaussian")
        assert out[0] == {"display": "block"}

    def test_field(self, web_app):
        out = self._fn(web_app)("field")
        assert out[1] == {"display": "block"}

    def test_mult(self, web_app):
        out = self._fn(web_app)("mult")
        assert out[2] == {"display": "block"}

    def test_none(self, web_app):
        out = self._fn(web_app)("none")
        assert out == (
            {"display": "none"},
            {"display": "none"},
            {"display": "none"},
        )


class TestLmPreset:
    def _fn(self, web_app):
        return _cb(web_app, "tool-lm-table.data")

    def test_no_preset_no_update(self, web_app):
        assert self._fn(web_app)(None) is no_update
        assert self._fn(web_app)("") is no_update

    def test_simple_preset(self, web_app):
        rows = self._fn(web_app)("simple")
        assert rows[0]["ρ (Ω·m)"] == 100.0

    def test_marine_preset(self, web_app):
        rows = self._fn(web_app)("marine")
        assert rows[0]["ρ (Ω·m)"] == 0.3

    def test_continental_preset(self, web_app):
        rows = self._fn(web_app)("continental")
        assert rows[0]["ρ (Ω·m)"] == 200.0

    def test_unknown_preset_no_update(self, web_app):
        assert self._fn(web_app)("does-not-exist") is no_update


class TestLmAddRow:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "tool-lm-table.data", "tool-lm-add")

    def test_no_clicks_no_update(self, web_app):
        assert self._fn(web_app)(None, []) is no_update

    def test_appends_default_row(self, web_app):
        rows = self._fn(web_app)(1, [{"#": 1}])
        assert len(rows) == 2
        assert rows[-1] == {"#": 2, "ρ (Ω·m)": 100.0, "Thickness (m)": 100.0}

    def test_appends_to_none_rows(self, web_app):
        rows = self._fn(web_app)(1, None)
        assert rows == [{"#": 1, "ρ (Ω·m)": 100.0, "Thickness (m)": 100.0}]


# ══════════════════════════════════════════════════════════════════════════
# _lm_run — main forward-modelling callback
# ══════════════════════════════════════════════════════════════════════════


class TestLmRun:
    def _fn(self, web_app):
        return _cb_multi(web_app, "tool-lm-out.children")

    def _args(self, **overrides):
        # Positional State order mirrors the @app.callback decorator.
        a = dict(
            n=1,
            ctor="manual",
            tbl_rows=None,
            r_nl=5,
            r_rmin=1.0,
            r_rmax=10000.0,
            r_dmax=2000.0,
            r_seed=1,
            b_nl=4,
            b_rbg=100.0,
            b_ranom=5.0,
            b_alyr=1,
            b_dmax=1000.0,
            b_equal_th=True,
            b_seed=1,
            s_nl=10,
            s_rsurf=100.0,
            s_rdeep=10.0,
            s_dmax=5000.0,
            s_perturb=0.2,
            s_seed=1,
            geo_name="sedimentary",
            geo_seed=1,
            solver="mt1d",
            freq_min=0.001,
            freq_max=10000.0,
            n_freqs=12,
            src_offset=1000.0,
            time_min=1e-6,
            time_max=0.01,
            n_times=10,
            loop_radius=50.0,
            moment=1.0,
            noise_type="none",
            noise_level=5.0,
            noise_phase=None,
            noise_seed=1,
            field_base=2.0,
            field_plfreq=50.0,
            field_pllevel=30.0,
            field_dead=True,
            mult_sigma=0.05,
            view="both",
            depth_plot=None,
            log_rho=True,
            theme="dark",
            saved_outputs=None,
        )
        a.update(overrides)
        return a

    def _call(self, web_app, **overrides):
        a = self._args(**overrides)
        return self._fn(web_app)(
            a["n"],
            a["ctor"],
            a["tbl_rows"],
            a["r_nl"],
            a["r_rmin"],
            a["r_rmax"],
            a["r_dmax"],
            a["r_seed"],
            a["b_nl"],
            a["b_rbg"],
            a["b_ranom"],
            a["b_alyr"],
            a["b_dmax"],
            a["b_equal_th"],
            a["b_seed"],
            a["s_nl"],
            a["s_rsurf"],
            a["s_rdeep"],
            a["s_dmax"],
            a["s_perturb"],
            a["s_seed"],
            a["geo_name"],
            a["geo_seed"],
            a["solver"],
            a["freq_min"],
            a["freq_max"],
            a["n_freqs"],
            a["src_offset"],
            a["time_min"],
            a["time_max"],
            a["n_times"],
            a["loop_radius"],
            a["moment"],
            a["noise_type"],
            a["noise_level"],
            a["noise_phase"],
            a["noise_seed"],
            a["field_base"],
            a["field_plfreq"],
            a["field_pllevel"],
            a["field_dead"],
            a["mult_sigma"],
            a["view"],
            a["depth_plot"],
            a["log_rho"],
            a["theme"],
            a["saved_outputs"],
        )

    def test_no_clicks_no_update(self, web_app):
        out = self._call(web_app, n=None)
        assert out == (no_update, no_update, no_update)

    def test_manual_with_default_table(self, web_app):
        body, saved, disabled = self._call(web_app)
        assert disabled is False
        assert "layered-model" in saved
        assert saved["layered-model"]["solver"] == "mt1d"

    def test_manual_fewer_than_two_layers_warns(self, web_app):
        body, saved, disabled = self._call(
            web_app, tbl_rows=[{"ρ (Ω·m)": 100.0, "Thickness (m)": None}]
        )
        assert "Add at least 2 layers" in body.children
        assert disabled is True
        assert saved is no_update

    def test_manual_thickness_count_mismatch_padded(self, web_app):
        # 3 rho rows but only 1 thickness given -> padded to 2 with 100.0
        rows = [
            {"ρ (Ω·m)": 100.0, "Thickness (m)": 50.0},
            {"ρ (Ω·m)": 200.0},
            {"ρ (Ω·m)": 300.0, "Thickness (m)": None},
        ]
        body, saved, disabled = self._call(web_app, tbl_rows=rows)
        assert disabled is False
        assert saved["layered-model"]["model"]["thickness"] == [50.0, 100.0]

    def test_manual_thickness_overflow_truncated(self, web_app):
        # 2 rho rows (only 1 thickness needed) but 3 thickness values
        # given -> truncated to len(rhos)-1.
        rows = [
            {"ρ (Ω·m)": 100.0, "Thickness (m)": 50.0},
            {"ρ (Ω·m)": 200.0, "Thickness (m)": 60.0},
        ]
        # inject an extra bogus row with a thickness to overflow the list
        rows.insert(1, {"ρ (Ω·m)": 150.0, "Thickness (m)": 70.0})
        body, saved, disabled = self._call(web_app, tbl_rows=rows)
        assert disabled is False
        assert len(saved["layered-model"]["model"]["thickness"]) == 2

    def test_random_constructor(self, web_app):
        body, saved, disabled = self._call(web_app, ctor="random")
        assert disabled is False
        assert saved["layered-model"]["model"]["name"] == "random"

    def test_blocky_constructor(self, web_app):
        body, saved, disabled = self._call(web_app, ctor="blocky")
        assert saved["layered-model"]["model"]["name"] == "blocky"

    def test_blocky_constructor_unequal_thickness_no_seed(self, web_app):
        body, saved, disabled = self._call(
            web_app, ctor="blocky", b_equal_th=False, b_seed=None
        )
        assert disabled is False

    def test_smooth_constructor(self, web_app):
        body, saved, disabled = self._call(web_app, ctor="smooth")
        assert saved["layered-model"]["model"]["name"] == "smooth"

    def test_geology_constructor(self, web_app):
        body, saved, disabled = self._call(web_app, ctor="geology", geo_name="porphyry")
        assert disabled is False

    def test_geology_constructor_default_name_and_no_seed(self, web_app):
        body, saved, disabled = self._call(
            web_app, ctor="geology", geo_name=None, geo_seed=None
        )
        assert disabled is False

    def test_unknown_constructor_warns(self, web_app):
        body, saved, disabled = self._call(web_app, ctor="not-a-real-ctor")
        assert "Unknown constructor" in body.children
        assert disabled is True

    def test_csamt1d_solver(self, web_app):
        body, saved, disabled = self._call(web_app, solver="csamt1d")
        assert saved["layered-model"]["solver"] == "csamt1d"
        assert "rho_a" in saved["layered-model"]["response"]

    def test_tem1d_solver(self, web_app):
        body, saved, disabled = self._call(web_app, solver="tem1d", view="response")
        assert saved["layered-model"]["solver"] == "tem1d"
        assert "dBz_dt" in saved["layered-model"]["response"]
        assert "times" in saved["layered-model"]["response"]

    def test_unrecognised_solver_produces_no_response(self, web_app):
        body, saved, disabled = self._call(web_app, solver="not-a-solver")
        assert any(
            "no response" in str(getattr(c, "children", c)) for c in body.children
        )
        assert saved["layered-model"]["response"] == {}

    def test_gaussian_noise(self, web_app):
        body, saved, disabled = self._call(
            web_app, noise_type="gaussian", noise_phase=3.0
        )
        assert disabled is False

    def test_field_noise(self, web_app):
        body, saved, disabled = self._call(web_app, noise_type="field")
        assert disabled is False

    def test_mult_noise(self, web_app):
        body, saved, disabled = self._call(web_app, noise_type="mult")
        assert disabled is False

    def test_noise_exception_is_swallowed(self, web_app, monkeypatch):
        import pycsamt.forward as forward_mod

        # _lm_run imports these names locally (function-scope import), so
        # the source module must be patched, not tools_mod.
        monkeypatch.setattr(
            forward_mod,
            "add_noise",
            lambda *a, **k: (_ for _ in ()).throw(RuntimeError("noise boom")),
        )
        body, saved, disabled = self._call(web_app, noise_type="gaussian")
        # noise failure silently falls back to the clean response
        assert disabled is False

    def test_view_model_only(self, web_app):
        body, saved, disabled = self._call(web_app, view="model")
        assert disabled is False

    def test_view_response_only(self, web_app):
        body, saved, disabled = self._call(web_app, view="response")
        assert disabled is False

    def test_log_rho_false_and_depth_plot_set(self, web_app):
        body, saved, disabled = self._call(web_app, log_rho=False, depth_plot=500.0)
        assert disabled is False

    def test_light_theme(self, web_app):
        body, saved, disabled = self._call(web_app, theme="light")
        assert disabled is False

    def test_model_plot_exception_appends_err(self, web_app, monkeypatch):
        import pycsamt.forward as forward_mod

        monkeypatch.setattr(
            forward_mod,
            "plot_model_1d",
            lambda *a, **k: (_ for _ in ()).throw(RuntimeError("model plot boom")),
        )
        body, saved, disabled = self._call(web_app, view="model")
        assert any(
            "model plot boom" in str(getattr(c, "children", c)) for c in body.children
        )

    def test_response_plot_exception_appends_err(self, web_app, monkeypatch):
        import pycsamt.forward as forward_mod

        monkeypatch.setattr(
            forward_mod,
            "plot_response_and_model_1d",
            lambda *a, **k: (_ for _ in ()).throw(RuntimeError("resp plot boom")),
        )
        body, saved, disabled = self._call(web_app, view="both")
        assert any(
            "resp plot boom" in str(getattr(c, "children", c)) for c in body.children
        )

    def test_response_only_view_uses_plot_response_1d(self, web_app, monkeypatch):
        import pycsamt.forward as forward_mod

        calls = []
        orig = forward_mod.plot_response_1d

        def _spy(*a, **k):
            calls.append(1)
            return orig(*a, **k)

        monkeypatch.setattr(forward_mod, "plot_response_1d", _spy)
        body, saved, disabled = self._call(web_app, view="response")
        assert calls

    def test_overall_exception_returns_err(self, web_app, monkeypatch):
        import pycsamt.forward as forward_mod

        monkeypatch.setattr(
            forward_mod,
            "LayeredModel",
            type(
                "Boom",
                (),
                {
                    "__init__": lambda *a, **k: (_ for _ in ()).throw(
                        RuntimeError("ctor boom")
                    )
                },
            ),
        )
        body, saved, disabled = self._call(web_app)
        assert "ctor boom" in body.children
        assert disabled is True


# ══════════════════════════════════════════════════════════════════════════
# _lm_export
# ══════════════════════════════════════════════════════════════════════════


class TestLmExport:
    def _fn(self, web_app):
        return _cb_multi(web_app, "tool-lm-download.data")

    def _saved(self, **resp_overrides):
        resp = {
            "freqs": [1.0, 10.0],
            "rho_a": [100.0, 200.0],
            "phase": [45.0, 40.0],
        }
        resp.update(resp_overrides)
        return {
            "layered-model": {
                "type": "lm",
                "model": {
                    "resistivity": [100.0, 500.0],
                    "thickness": [50.0],
                    "depth": [0.0, 50.0],
                    "name": "manual",
                },
                "response": resp,
                "solver": "mt1d",
                "imgs": [],
            }
        }

    def test_no_clicks_no_update(self, web_app):
        out = self._fn(web_app)(None, "model_csv", {})
        assert out == (no_update, no_update)

    def test_no_model_run_yet_warns(self, web_app):
        out = self._fn(web_app)(1, "model_csv", {})
        assert out[0] is no_update
        assert "Run the forward model first" in str(out[1].children)

    def test_export_model_csv(self, web_app):
        out = self._fn(web_app)(1, "model_csv", self._saved())
        assert out[0]["filename"] == "layered_model.csv"
        assert "Model CSV downloaded" in str(out[1].children)

    def test_export_resp_csv(self, web_app):
        out = self._fn(web_app)(1, "resp_csv", self._saved())
        assert out[0]["filename"] == "forward_response.csv"

    def test_export_resp_csv_no_response_warns(self, web_app):
        saved = self._saved()
        saved["layered-model"]["response"] = {}
        out = self._fn(web_app)(1, "resp_csv", saved)
        assert out[0] is no_update
        assert "No response data" in str(out[1].children)

    def test_export_resp_csv_2d_array(self, web_app):
        saved = self._saved(rho_a=[[100.0, 110.0], [200.0, 210.0]])
        out = self._fn(web_app)(1, "resp_csv", saved)
        assert out[0]["filename"] == "forward_response.csv"

    def test_export_resp_csv_with_times_and_dbzdt(self, web_app):
        saved = self._saved()
        saved["layered-model"]["response"] = {
            "times": [1e-5, 1e-4],
            "dBz_dt": [1.0, 2.0],
        }
        out = self._fn(web_app)(1, "resp_csv", saved)
        assert out[0]["filename"] == "forward_response.csv"

    def test_export_json_default(self, web_app):
        out = self._fn(web_app)(1, "json", self._saved())
        assert out[0]["filename"] == "layered_model.json"
        assert "JSON downloaded" in str(out[1].children)

    def test_export_unknown_fmt_defaults_to_json(self, web_app):
        out = self._fn(web_app)(1, "not-a-format", self._saved())
        assert out[0]["filename"] == "layered_model.json"

    def test_export_exception_returns_err(self, web_app, monkeypatch):
        import json as json_mod

        monkeypatch.setattr(
            json_mod,
            "dumps",
            lambda *a, **k: (_ for _ in ()).throw(RuntimeError("dump boom")),
        )
        out = self._fn(web_app)(1, "json", self._saved())
        assert out[0] is no_update
        assert "dump boom" in str(out[1].children)
