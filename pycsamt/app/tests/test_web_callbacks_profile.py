# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.web.callbacks.profile (Profile View page).

Real data
---------
A small (first 6 stations) real-EDI subset of data/AMT/WILLY_DATA/L18PLT/
drives ``PlotController``/``build_multi_pseudosection`` for real through the
Dash callback (matplotlib PNG + Plotly pseudosection tabs alike). Nothing
here is genuinely expensive, so nothing is mocked for the render paths
themselves -- only ``cache_get``/``fig_to_src`` are monkeypatched in the two
tests that specifically target the outer exception handler.

Pre-existing gap (not a regression, just worth knowing)
---------------------------------------------------------
``PlotController`` has no ``draw_2d_section`` method at all, so the
"section" tab's ``_draw_section()`` helper *always* takes its
``except AttributeError`` branch and shows "2D section not available
(run inversion first)" instead of a real section -- there is currently no
way to reach a real 2-D section render through this tab.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from dash import no_update
from dash.exceptions import PreventUpdate

from pycsamt.app.web.callbacks.profile import _SKIP, _comp_bar, _draw_section
from pycsamt.app.web.layout import IDs

_ROOT = Path(__file__).parents[3]
_WILLY_L18 = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
_HAS_WILLY = _WILLY_L18.exists() and any(_WILLY_L18.glob("*.edi"))


# ── Dash callback-lookup helpers (shared pattern across web-callback tests) ──


def _unwrap(entry):
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def _cb(web_app, output_id_prop):
    return _unwrap(web_app.callback_map[output_id_prop])


def _cb_multi(web_app, *substrings):
    key = next(
        k for k in web_app.callback_map if all(s in k for s in substrings)
    )
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


# ── Fixtures ─────────────────────────────────────────────────────────────────


@pytest.fixture(scope="session")
def willy_sites():
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_WILLY:
        pytest.skip("WILLY L18PLT data not available")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(_WILLY_L18))


@pytest.fixture(scope="session")
def willy_subset(willy_sites):
    """Small (6-station) real Sites subset -- keeps per-render calls fast."""
    from pycsamt.emtools import ensure_sites

    return ensure_sites(willy_sites.as_list()[:6])


@pytest.fixture(scope="session")
def willy_station_id(willy_subset):
    return willy_subset.as_list()[0].station


@pytest.fixture
def willy_subset_records(willy_subset):
    return [
        {"ID": ed.station, "Line": "L1" if i % 2 == 0 else "L2"}
        for i, ed in enumerate(willy_subset.as_list())
    ]


@pytest.fixture
def cached_session(willy_subset):
    from pycsamt.app.web.cache import cache_set

    session_id = "test-profile-session"
    cache_set(session_id, willy_subset)
    return session_id


# ── 1. _comp_bar helper ──────────────────────────────────────────────────────


class TestCompBar:
    def test_rho_prefix_and_pills(self):
        children = _comp_bar(("xy", "yx"), "rho")
        assert children[0].children == "ρₐ sections:"
        assert len(children) == 3
        assert children[1].children == "Z_XY"
        assert "te" in children[1].className
        assert children[2].children == "Z_YX"
        assert "tm" in children[2].className

    def test_phi_prefix(self):
        children = _comp_bar(("yx",), "phi")
        assert children[0].children == "φ sections:"
        assert children[1].children == "Z_YX"

    def test_unknown_component_falls_back_to_upper_no_class(self):
        children = _comp_bar(("zz",), "rho")
        assert children[1].children == "ZZ"
        assert children[1].className == "prof-ps-pill"


# ── 2. _draw_section helper ──────────────────────────────────────────────────


class TestDrawSection:
    def test_missing_draw_2d_section_shows_placeholder(self, willy_subset):
        from pycsamt.app.desktop.controllers.plot_controller import (
            PlotController,
        )

        ctrl = PlotController()
        ctrl.set_sites(willy_subset)
        src = _draw_section(ctrl)
        assert src.startswith("data:image/png;base64,")

    def test_generic_exception_also_shows_a_valid_placeholder(
        self, willy_subset, monkeypatch
    ):
        from pycsamt.app.desktop.controllers.plot_controller import (
            PlotController,
        )

        monkeypatch.setattr(
            PlotController,
            "draw_2d_section",
            lambda self, ax: (_ for _ in ()).throw(ValueError("weird")),
            raising=False,
        )
        ctrl = PlotController()
        ctrl.set_sites(willy_subset)
        src = _draw_section(ctrl)
        assert src.startswith("data:image/png;base64,")


# ── 3. switch_tab callback ───────────────────────────────────────────────────


class TestSwitchTab:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.PROF_ACTIVE_TAB}.data")

    def _ctx(self, monkeypatch, triggered_id):
        import pycsamt.app.web.callbacks.profile as profile_mod

        monkeypatch.setattr(
            profile_mod, "ctx", type("C", (), {"triggered_id": triggered_id})()
        )

    def test_no_trigger_no_update(self, web_app, monkeypatch):
        self._ctx(monkeypatch, None)
        assert self._fn(web_app)(*([1] * 7)) is no_update

    def test_extracts_tab_name_from_button_id(self, web_app, monkeypatch):
        self._ctx(monkeypatch, "prof-tab-rho-ps")
        out = self._fn(web_app)(*([1] * 7))
        assert out == "rho-ps"


# ── 4. pub_back callback ─────────────────────────────────────────────────────


class TestPubBack:
    def _fn(self, web_app):
        return _cb_by_input(
            web_app, f"{IDs.PROF_ACTIVE_TAB}.data", "prof-pub-back-btn"
        )

    def test_no_clicks_no_update(self, web_app):
        assert self._fn(web_app)(None) is no_update

    def test_clicks_returns_rho_phi(self, web_app):
        assert self._fn(web_app)(3) == "rho-phi"


# ── 5. update_profile callback ───────────────────────────────────────────────


class TestUpdateProfile:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.IMG_RHO_PHI}.src")

    def test_missing_store_data_returns_skip(self, web_app):
        out = self._fn(web_app)(
            "rho-phi", {}, None, "dark", False, ["xy", "yx"], 0, "auto", None,
            "sid", None,
        )
        assert out == _SKIP

    def test_sites_not_cached_returns_skip(self, web_app):
        out = self._fn(web_app)(
            "rho-phi",
            {},
            {"station_records": []},
            "dark",
            False,
            ["xy", "yx"],
            0,
            "auto",
            None,
            "no-such-session",
            None,
        )
        assert out == _SKIP

    def test_all_lines_muted_returns_warning_images(
        self, web_app, willy_subset, cached_session, willy_subset_records
    ):
        store_data = {"station_records": willy_subset_records}
        als = {"active": ["GhostLine"], "all": ["GhostLine", "L1"]}
        out = self._fn(web_app)(
            "rho-phi",
            {},
            store_data,
            "dark",
            False,
            ["xy", "yx"],
            0,
            "auto",
            None,
            cached_session,
            als,
        )
        (
            rho_phi,
            tipper,
            pt,
            rho_ps_fig,
            phi_ps_fig,
            section,
            pub,
            rho_hdr,
            phi_hdr,
            toast,
            body,
        ) = out
        assert rho_phi.startswith("data:image/png;base64,")
        assert rho_phi == tipper == pt == section == pub
        assert rho_ps_fig is no_update
        assert phi_ps_fig is no_update
        assert rho_hdr is no_update
        assert phi_hdr is no_update
        assert toast is False

    def test_rho_ps_tab_returns_plotly_figure_and_header(
        self, web_app, willy_subset, cached_session
    ):
        import plotly.graph_objects as go

        out = self._fn(web_app)(
            "rho-ps",
            {},
            {"station_records": []},
            "dark",
            False,
            ["xy", "yx"],
            0,
            "auto",
            None,
            cached_session,
            None,
        )
        rho_ps_fig = out[3]
        rho_hdr = out[7]
        assert isinstance(rho_ps_fig, go.Figure)
        assert rho_hdr[0].children == "ρₐ sections:"
        assert out[-2] is False  # toast_error

    def test_phi_ps_tab_returns_plotly_figure_and_header(
        self, web_app, willy_subset, cached_session
    ):
        import plotly.graph_objects as go

        out = self._fn(web_app)(
            "phi-ps",
            {},
            {"station_records": []},
            "light",
            False,
            ["xy"],
            0,
            "auto",
            None,
            cached_session,
            None,
        )
        phi_ps_fig = out[4]
        phi_hdr = out[8]
        assert isinstance(phi_ps_fig, go.Figure)
        assert phi_hdr[0].children == "φ sections:"

    def test_rho_phi_tab_renders_real_png(
        self, web_app, willy_subset, cached_session, willy_station_id
    ):
        out = self._fn(web_app)(
            "rho-phi",
            {"station_id": willy_station_id},
            {"station_records": []},
            "dark",
            True,
            ["xy", "yx"],
            0,
            "auto",
            [-2.0, 2.0],
            cached_session,
            None,
        )
        rho_phi_src = out[0]
        assert rho_phi_src.startswith("data:image/png;base64,")
        assert len(rho_phi_src) > 2000
        assert out[1] is no_update  # tipper untouched
        assert out[-2] is False

    def test_pt_tab_renders_real_png(
        self, web_app, willy_subset, cached_session, willy_station_id
    ):
        out = self._fn(web_app)(
            "pt",
            {"station_id": willy_station_id},
            {"station_records": []},
            "dark",
            False,
            ["xy", "yx"],
            0,
            "0_90",
            None,
            cached_session,
            None,
        )
        pt_src = out[2]
        assert pt_src.startswith("data:image/png;base64,")

    def test_tipper_tab_renders_real_png(
        self, web_app, willy_subset, cached_session, willy_station_id
    ):
        out = self._fn(web_app)(
            "tipper",
            {"station_id": willy_station_id},
            {"station_records": []},
            "light",
            False,
            ["xy", "yx"],
            0,
            "-90_90",
            None,
            cached_session,
            None,
        )
        tipper_src = out[1]
        assert tipper_src.startswith("data:image/png;base64,")

    def test_section_tab_renders_placeholder_png(
        self, web_app, willy_subset, cached_session, willy_station_id
    ):
        out = self._fn(web_app)(
            "section",
            {"station_id": willy_station_id},
            {"station_records": []},
            "dark",
            False,
            ["xy", "yx"],
            0,
            "auto",
            None,
            cached_session,
            None,
        )
        section_src = out[5]
        assert section_src.startswith("data:image/png;base64,")

    def test_pub_tab_renders_real_png(
        self, web_app, willy_subset, cached_session, willy_station_id
    ):
        out = self._fn(web_app)(
            "pub",
            {"station_id": willy_station_id},
            {"station_records": []},
            "dark",
            False,
            ["xy", "yx"],
            0,
            "-360_360",
            None,
            cached_session,
            None,
        )
        pub_src = out[6]
        assert pub_src.startswith("data:image/png;base64,")

    def test_unknown_active_tab_falls_back_to_rho_phi(
        self, web_app, willy_subset, cached_session
    ):
        out = self._fn(web_app)(
            "not-a-real-tab",
            {},
            {"station_records": []},
            "dark",
            False,
            None,
            0,
            "auto",
            None,
            cached_session,
            None,
        )
        assert out[0].startswith("data:image/png;base64,")

    def test_exception_inside_render_returns_error_toast(
        self, web_app, willy_subset, cached_session, monkeypatch
    ):
        import pycsamt.app.web.callbacks.profile as profile_mod

        monkeypatch.setattr(
            profile_mod,
            "fig_to_src",
            lambda fig: (_ for _ in ()).throw(RuntimeError("boom")),
        )
        out = self._fn(web_app)(
            "rho-phi",
            {},
            {"station_records": []},
            "dark",
            False,
            ["xy", "yx"],
            0,
            "auto",
            None,
            cached_session,
            None,
        )
        toast, body = out[-2], out[-1]
        assert toast is True
        assert "boom" in body
        assert "Profile render failed (rho-phi)" in body


# ── 6. navigate_station callback ─────────────────────────────────────────────


class TestNavigateStation:
    def _fn(self, web_app):
        return _cb_by_input(
            web_app, "store-selection.data", IDs.PROFILE_PAGE_PREV
        )

    def _ctx(self, monkeypatch, triggered_id):
        import pycsamt.app.web.callbacks.profile as profile_mod

        monkeypatch.setattr(
            profile_mod, "ctx", type("C", (), {"triggered_id": triggered_id})()
        )

    def test_no_options_no_update(self, web_app, monkeypatch):
        self._ctx(monkeypatch, IDs.PROFILE_PAGE_PREV)
        assert self._fn(web_app)(1, None, "A", None) is no_update
        assert self._fn(web_app)(1, None, "A", []) is no_update

    def test_no_trigger_no_update(self, web_app, monkeypatch):
        self._ctx(monkeypatch, None)
        out = self._fn(web_app)(
            None, None, "A", [{"value": "A"}, {"value": "B"}]
        )
        assert out is no_update

    def test_current_id_not_in_options_no_update(self, web_app, monkeypatch):
        self._ctx(monkeypatch, IDs.PROFILE_PAGE_PREV)
        out = self._fn(web_app)(
            1, None, "ZZZ", [{"value": "A"}, {"value": "B"}]
        )
        assert out is no_update

    def test_prev_wraps_around_to_last(self, web_app, monkeypatch):
        self._ctx(monkeypatch, IDs.PROFILE_PAGE_PREV)
        out = self._fn(web_app)(
            1, None, "A", [{"value": "A"}, {"value": "B"}]
        )
        assert out == {"station_id": "B"}

    def test_next_advances(self, web_app, monkeypatch):
        self._ctx(monkeypatch, IDs.PROFILE_PAGE_NEXT)
        out = self._fn(web_app)(
            None, 1, "A", [{"value": "A"}, {"value": "B"}]
        )
        assert out == {"station_id": "B"}


# ── 7. export_profile_tab callback ───────────────────────────────────────────


class TestExportProfileTab:
    def _fn(self, web_app):
        return _cb(web_app, "profile-page-download.data")

    def test_no_clicks_no_update(self, web_app):
        assert self._fn(web_app)(None, "rho-phi", "src", "", "", "", "") is no_update

    def test_unknown_tab_no_update(self, web_app):
        assert self._fn(web_app)(1, "rho-ps", "x", "", "", "", "") is no_update

    def test_missing_or_non_image_src_no_update(self, web_app):
        assert self._fn(web_app)(1, "rho-phi", None, "", "", "", "") is no_update
        assert (
            self._fn(web_app)(1, "rho-phi", "not-a-data-uri", "", "", "", "")
            is no_update
        )

    def test_valid_export_returns_download_dict(self, web_app):
        out = self._fn(web_app)(
            1, "rho-phi", "data:image/png;base64,ABC123", "", "", "", ""
        )
        assert out == {
            "content": "ABC123",
            "filename": "pycsamt_rho_phi.png",
            "base64": True,
            "type": "image/png",
        }

    def test_valid_export_for_each_tab_uses_the_right_src_slot(self, web_app):
        out = self._fn(web_app)(
            1,
            "section",
            "src-rho-phi",
            "src-tipper",
            "src-pt",
            "data:image/png;base64,SECTIONDATA",
            "src-pub",
        )
        assert out["content"] == "SECTIONDATA"
        assert out["filename"] == "pycsamt_2d_section.png"
