# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.mapview.callbacks.controls."""

from __future__ import annotations

import types

import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")

import pycsamt.app.mapview.callbacks.controls as controls_mod


def _capture():
    """Capture every ``app.callback``-decorated closure by function name.

    ``_register_group_visibility`` uses ``app.clientside_callback`` instead
    (pure-JS callbacks with no Python body to unit-test), so the fake
    registrar just no-ops that method.
    """
    captured: dict = {}

    class _App:
        def callback(self, *a, **k):
            def deco(fn):
                captured[fn.__name__] = fn
                return fn

            return deco

        def clientside_callback(self, *a, **k):
            return None

    controls_mod.register_controls(_App())
    return captured


def _set_trigger(monkeypatch, triggered_id):
    """``controls.py`` binds ``ctx`` at module-import time (``from dash
    import ... ctx``), so the module's own bound name must be patched
    directly rather than ``dash.ctx``."""
    monkeypatch.setattr(
        controls_mod, "ctx", types.SimpleNamespace(triggered_id=triggered_id)
    )


class TestFmtFreq:
    def test_none_returns_dash(self):
        from pycsamt.app.mapview.callbacks.controls import _fmt_freq

        assert _fmt_freq(None) == "—"

    def test_khz_range(self):
        from pycsamt.app.mapview.callbacks.controls import _fmt_freq

        assert _fmt_freq(1500) == "1.5 kHz"

    def test_hz_range(self):
        from pycsamt.app.mapview.callbacks.controls import _fmt_freq

        assert _fmt_freq(50) == "50 Hz"

    def test_sub_hz_range_includes_period(self):
        from pycsamt.app.mapview.callbacks.controls import _fmt_freq

        result = _fmt_freq(0.1)
        assert "Hz" in result and "s)" in result


class TestPresets:
    def test_depth_presets_cover_full_and_bands(self):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.callbacks.controls import _DEPTH_PRESETS

        assert _DEPTH_PRESETS[IDs.BTN_DEPTH_FULL] == (None, None)
        assert _DEPTH_PRESETS[IDs.BTN_DEPTH_500] == (0, 500)
        assert _DEPTH_PRESETS[IDs.BTN_DEPTH_2K] == (0, 2000)

    def test_toolbar_depth_select_presets_cover_full_and_bands(self):
        from pycsamt.app.mapview.callbacks.controls import (
            _DEPTH_SELECT_PRESETS,
        )

        assert _DEPTH_SELECT_PRESETS["full"] == (None, None)
        assert _DEPTH_SELECT_PRESETS["500"] == (0, 500)
        assert _DEPTH_SELECT_PRESETS["1000"] == (0, 1000)
        assert _DEPTH_SELECT_PRESETS["2000"] == (0, 2000)

    def test_rho_presets_are_ordered_bands(self):
        from pycsamt.app.mapview.callbacks.controls import _RHO_PRESETS

        for lo, hi in _RHO_PRESETS.values():
            assert lo < hi


class TestRegisterControls:
    def test_register_controls_is_callable(self):
        from pycsamt.app.mapview.callbacks.controls import register_controls

        assert callable(register_controls)

    def test_expected_outputs_wired(self):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.app import create_app

        app = create_app()
        cb_outputs = str(app.callback_map)
        assert IDs.STORE_CONTROLS in cb_outputs
        assert IDs.CTL_DEPTH_LO in cb_outputs
        assert IDs.CTL_RHO_LO in cb_outputs

    def test_capture_finds_all_python_callbacks(self):
        captured = _capture()
        assert {
            "sync",
            "apply_preset",
            "apply_toolbar_preset",
            "apply_rho_preset",
            "populate",
            "gather",
        } <= set(captured)


class TestNiceStep:
    @pytest.mark.parametrize(
        "raw, expected",
        [
            (0.5, 1.0),
            (1.5, 2.0),
            (3.0, 5.0),
            (7.0, 10.0),
            (600.0, 1000.0),
        ],
    )
    def test_snaps_up_to_next_nice_value(self, raw, expected):
        assert controls_mod._nice_step(raw) == expected


class TestSplitNames:
    def test_blank_returns_none(self):
        assert controls_mod._split_names("") is None
        assert controls_mod._split_names("   ") is None
        assert controls_mod._split_names(None) is None

    def test_splits_on_comma_semicolon_and_whitespace(self):
        assert controls_mod._split_names("A, B; C  D") == ("A", "B", "C", "D")

    def test_all_separators_only_returns_none(self):
        assert controls_mod._split_names(" ,; ") is None


class TestAsFraction:
    def test_clamps_to_0_1_range(self):
        assert controls_mod._as_fraction(2.0) == 1.0
        assert controls_mod._as_fraction(-1.0) == 0.0
        assert controls_mod._as_fraction(0.5) == 0.5

    def test_invalid_falls_back_to_1(self):
        assert controls_mod._as_fraction(None) == 1.0
        assert controls_mod._as_fraction("nope") == 1.0


class TestAsFloat:
    def test_valid_string_parses(self):
        assert controls_mod._as_float("2.5", 0.0) == 2.5

    def test_invalid_falls_back_to_default(self):
        assert controls_mod._as_float(None, 3.0) == 3.0
        assert controls_mod._as_float("nope", -1.0) == -1.0


class TestSyncMapOverlay:
    def test_edi_survey_defaults(self):
        captured = _capture()
        sync = captured["sync"]
        opts, value, lo, hi, step, depth, marks = sync(
            {"is_inversion": False}, None, None
        )
        assert opts == controls_mod._EDI_OVERLAYS
        assert value == "index"
        assert lo == 0
        assert hi == 1000
        assert marks[0] == "0"

    def test_inversion_survey_uses_inv_overlays_and_depth_range(self):
        captured = _capture()
        sync = captured["sync"]
        opts, value, lo, hi, step, depth, marks = sync(
            {"is_inversion": True, "depth_range": [10.0, 510.0]},
            "depth_rho",
            None,
        )
        assert opts == controls_mod._INV_OVERLAYS
        assert value == "depth_rho"
        assert lo == 10
        assert hi == 510

    def test_current_overlay_invalid_for_mode_falls_back(self):
        captured = _capture()
        sync = captured["sync"]
        _, value, *_ = sync({"is_inversion": True}, "skin_depth", None)
        assert value == controls_mod._INV_OVERLAYS[0]["value"]

    def test_current_depth_within_range_is_preserved(self):
        captured = _capture()
        sync = captured["sync"]
        _, _, lo, hi, _, depth, _ = sync(
            {"is_inversion": True, "depth_range": [0.0, 100.0]}, None, 40
        )
        assert depth == 40

    def test_current_depth_out_of_range_recentered(self):
        captured = _capture()
        sync = captured["sync"]
        _, _, lo, hi, _, depth, _ = sync(
            {"is_inversion": True, "depth_range": [0.0, 100.0]}, None, 999
        )
        assert lo <= depth <= hi

    def test_none_store_uses_edi_defaults(self):
        captured = _capture()
        sync = captured["sync"]
        opts, *_ = sync(None, None, None)
        assert opts == controls_mod._EDI_OVERLAYS


class TestDepthAndRhoPresets:
    def test_apply_preset_dispatches_on_triggered_id(self, monkeypatch):
        from pycsamt.app.mapview._ids import IDs

        captured = _capture()
        _set_trigger(monkeypatch, IDs.BTN_DEPTH_500)
        assert captured["apply_preset"](1, None, None, None) == (0, 500)

    def test_apply_preset_unknown_trigger_returns_none_none(self, monkeypatch):
        captured = _capture()
        _set_trigger(monkeypatch, "unknown-id")
        assert captured["apply_preset"](1, None, None, None) == (None, None)

    def test_apply_toolbar_preset_maps_select_values(self):
        captured = _capture()
        apply_toolbar_preset = captured["apply_toolbar_preset"]
        assert apply_toolbar_preset("1000") == (0, 1000)
        assert apply_toolbar_preset("full") == (None, None)
        assert apply_toolbar_preset("bogus") == (None, None)

    def test_apply_rho_preset_dispatches_on_triggered_id(self, monkeypatch):
        from pycsamt.app.mapview._ids import IDs

        captured = _capture()
        _set_trigger(monkeypatch, IDs.BTN_RHO_COND)
        assert captured["apply_rho_preset"](1, None, None, None) == (1, 100)

    def test_apply_rho_preset_unknown_trigger(self, monkeypatch):
        captured = _capture()
        _set_trigger(monkeypatch, None)
        assert captured["apply_rho_preset"](None, None, None, None) == (
            None,
            None,
        )


class TestFreqSlider:
    def test_populate_no_frequencies_returns_zeros(self):
        captured = _capture()
        populate = captured["populate"]
        assert populate(None) == (0, 0, None, 0)
        assert populate({}) == (0, 0, None, 0)

    def test_populate_builds_sparse_marks(self):
        captured = _capture()
        populate = captured["populate"]
        freqs = [100.0, 50.0, 10.0, 1.0, 0.1]
        lo, hi, marks, value = populate({"frequencies": freqs})
        assert lo == 0
        assert hi == len(freqs) - 1
        assert value == 0
        assert set(marks) == {0, 2, 4}

    def test_populate_two_frequencies_dedupes_marks(self):
        captured = _capture()
        populate = captured["populate"]
        lo, hi, marks, value = populate({"frequencies": [10.0, 1.0]})
        # n=2 -> {0, n//2=1, n-1=1} dedupes to {0, 1}
        assert set(marks) == {0, 1}


class TestGather:
    def _base_kwargs(self, **overrides):
        kwargs = dict(
            overlay=None,
            quantity=None,
            component=None,
            freq_idx=None,
            cmap=None,
            log=None,
            mode3d=None,
            labels=None,
            opacity=None,
            azimuth=None,
            spacing=None,
            depth_lo=None,
            depth_hi=None,
            n_slices=None,
            surfaces=None,
            contours=None,
            scale=None,
            vmin=None,
            vmax=None,
            crange_plo=None,
            crange_phi=None,
            rho_cutoff=None,
            rho_lo=None,
            rho_hi=None,
            topo=None,
            terrain=None,
            basemap=None,
            contour_levels=None,
            contour_mode=None,
            show_sta=None,
            sta_labels=None,
            sta_label_angle=None,
            sta_label_density=None,
            sta_label_names=None,
            sta_max=None,
            sta_symbol=None,
            sta_size=None,
            sta_color=None,
            contour_enable=None,
            marker_size=None,
            map_opacity=None,
            profiles=None,
            map_stations=None,
            map_depth=None,
            crs_mode=None,
            utm_zone=None,
            utm_hem=None,
            epsg=None,
            contour_interp=None,
            contour_smooth=None,
            contour_res=None,
            aspect=None,
            x_unit=None,
            depth_unit=None,
            smooth_sections=None,
            section_res=None,
            vol_smooth=None,
            geo_legend_visible=None,
            geo_legend_style=None,
            geo_fill=None,
            ve=None,
            store=None,
        )
        kwargs.update(overrides)
        return kwargs

    def test_all_defaults_no_store(self):
        captured = _capture()
        gather = captured["gather"]
        controls, label = gather(**self._base_kwargs())
        assert controls["vertical_exaggeration"] is None  # off
        assert controls["overlay"] == "index"
        assert controls["quantity"] == "rho"
        assert controls["frequency"] is None
        assert controls["basemap"] == "esri-satellite"
        assert controls["station_label_names"] is None
        assert controls["geology_legend"] is True
        assert label == "—"

    def test_frequency_selected_from_store(self):
        captured = _capture()
        gather = captured["gather"]
        controls, label = gather(
            **self._base_kwargs(
                freq_idx=1, store={"frequencies": [100.0, 10.0, 1.0]}
            )
        )
        assert controls["frequency"] == 10.0
        assert "Hz" in label

    def test_frequency_index_clamped_to_valid_range(self):
        captured = _capture()
        gather = captured["gather"]
        controls, _ = gather(
            **self._base_kwargs(
                freq_idx=99, store={"frequencies": [100.0, 10.0]}
            )
        )
        assert controls["frequency"] == 10.0
        controls, _ = gather(
            **self._base_kwargs(
                freq_idx=-5, store={"frequencies": [100.0, 10.0]}
            )
        )
        assert controls["frequency"] == 100.0

    def test_truthy_overrides_pass_through(self):
        captured = _capture()
        gather = captured["gather"]
        controls, _ = gather(
            **self._base_kwargs(
                overlay="rho",
                quantity="phase",
                component="yx",
                cmap="viridis",
                log=True,
                mode3d="volume",
                labels=True,
                opacity=0.5,
                azimuth=45.0,
                spacing=2.0,
                depth_lo=10,
                depth_hi=200,
                n_slices=4,
                surfaces=20,
                contours=True,
                scale="linear",
                vmin=1.0,
                vmax=100.0,
                topo=True,
                terrain=True,
                basemap="carto-positron",
                contour_levels=8,
                contour_mode="lines",
                show_sta=True,
                sta_labels=True,
                sta_label_angle=30,
                sta_label_density="0.5",
                sta_label_names="A, B",
                sta_max=5,
                sta_symbol="circle",
                sta_size=8,
                sta_color="#ff0000",
                contour_enable=True,
                marker_size=15,
                map_opacity=50,
                profiles=True,
                map_stations=True,
                map_depth=250,
                crs_mode="utm",
                utm_zone="30",
                utm_hem="S",
                epsg="32630",
                contour_interp="linear",
                contour_smooth=2.0,
                contour_res=200,
                aspect="equal",
                x_unit="km",
                depth_unit="ft",
                smooth_sections=True,
                section_res=50,
                vol_smooth=0.3,
                geo_legend_visible=False,
                geo_legend_style="line",
                geo_fill="hatch",
            )
        )
        assert controls["overlay"] == "rho"
        assert controls["component"] == "yx"
        assert controls["log"] is True
        assert controls["mode3d"] == "volume"
        assert controls["station_label_fraction"] == 0.5
        assert controls["station_label_names"] == ("A", "B")
        assert controls["station_max"] == 5
        assert controls["map_depth"] == 250.0
        assert controls["crs_mode"] == "utm"
        assert controls["geology_legend"] is False
        assert controls["geology_legend_style"] == "line"
        assert controls["geology_fill"] == "hatch"
        assert controls["volume_smoothing"] == 0.3
        assert controls["contour_res"] == 200
        assert controls["section_res"] == 50


def test_gather_vertical_exaggeration_values():
    gather = _capture()["gather"]
    base = TestGather()._base_kwargs
    assert gather(**base(ve=0))[0]["vertical_exaggeration"] == 0.0  # auto
    assert gather(**base(ve="2.5"))[0]["vertical_exaggeration"] == 2.5
    assert gather(**base(ve=""))[0]["vertical_exaggeration"] is None
