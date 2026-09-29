# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.web.callbacks.lines (survey-line management panel).

Pure dict/store manipulation -- no matplotlib, no cached Sites required.
``detect_lines_from_station_ids``/``apply_line_renames`` (pycsamt.site.lines)
are real, fast, deterministic string-parsing helpers, so they are exercised
directly rather than mocked.
"""

from __future__ import annotations

from dash import no_update
from dash.exceptions import PreventUpdate

from pycsamt.app.web.layout import IDs

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


def _records(pairs):
    """[(id, line), ...] -> station_records list."""
    return [{"ID": sid, "Line": ln} for sid, ln in pairs]


STORE_TWO_LINES = {
    "station_records": _records(
        [("A1", "L1"), ("A2", "L1"), ("B1", "L2"), ("B2", "L2")]
    )
}


# ── 1. init_active_lines ─────────────────────────────────────────────────────


class TestInitActiveLines:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.STORE_ACTIVE_LINES}.data")

    def test_no_store_data(self, web_app):
        out = self._fn(web_app)(None, None, "replace")
        assert out == {"active": [], "all": []}

    def test_first_load_defaults_all_active(self, web_app):
        out = self._fn(web_app)(STORE_TWO_LINES, None, "replace")
        assert out == {"active": ["L1", "L2"], "all": ["L1", "L2"]}

    def test_ignores_dash_placeholder_line(self, web_app):
        store = {"station_records": _records([("A1", "—"), ("A2", "L1")])}
        out = self._fn(web_app)(store, None, "replace")
        assert out == {"active": ["L1"], "all": ["L1"]}

    def test_replace_mode_resets_all_active(self, web_app):
        prev = {"active": ["L1"], "all": ["L1", "L2"]}
        out = self._fn(web_app)(STORE_TWO_LINES, prev, "replace")
        assert out == {"active": ["L1", "L2"], "all": ["L1", "L2"]}

    def test_append_mode_preserves_prior_choices_new_lines_active(self, web_app):
        prev = {"active": ["L1"], "all": ["L1", "L2"]}
        store = {
            "station_records": _records(
                [
                    ("A1", "L1"),
                    ("B1", "L2"),
                    ("C1", "L3"),
                ]
            )
        }
        out = self._fn(web_app)(store, prev, "append")
        # L2 was previously inactive -> stays inactive; L3 is new -> active
        assert set(out["active"]) == {"L1", "L3"}
        assert set(out["all"]) == {"L1", "L2", "L3"}

    def test_append_mode_without_prior_all_falls_back_to_replace(self, web_app):
        out = self._fn(web_app)(STORE_TWO_LINES, {"active": [], "all": []}, "append")
        assert out == {"active": ["L1", "L2"], "all": ["L1", "L2"]}


# ── 2. toggle_lines_offcanvas ────────────────────────────────────────────────


class TestToggleOffcanvas:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.LINES_OFFCANVAS}.is_open")

    def test_no_clicks_no_update(self, web_app):
        assert self._fn(web_app)(None, False) is no_update

    def test_toggles_open(self, web_app):
        assert self._fn(web_app)(1, False) is True

    def test_toggles_closed(self, web_app):
        assert self._fn(web_app)(1, True) is False


# ── 3. sync_mode_ui ──────────────────────────────────────────────────────────


class TestSyncModeUi:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.BTN_LINES_DETECT}.style")

    def test_default_folder_mode(self, web_app):
        detect_style, rename_style, hint = self._fn(web_app)(None)
        assert detect_style == {"display": "none"}
        assert rename_style == {"display": "none"}
        assert "subfolder" in hint

    def test_auto_mode_shows_detect_btn(self, web_app):
        detect_style, rename_style, hint = self._fn(web_app)("auto")
        assert detect_style == {"display": "block"}
        assert rename_style == {"display": "none"}

    def test_edit_mode_shows_rename_btn(self, web_app):
        detect_style, rename_style, hint = self._fn(web_app)("edit")
        assert detect_style == {"display": "none"}
        assert rename_style == {"display": "block"}


# ── 4. detect_lines_auto ─────────────────────────────────────────────────────


class TestDetectLinesAuto:
    def _fn(self, web_app):
        return _cb_by_input(web_app, f"{IDs.STORE_DATA}.data", IDs.BTN_LINES_DETECT)

    def test_no_clicks_prevents_update(self, web_app):
        with pytest_raises_prevent_update():
            self._fn(web_app)(None, STORE_TWO_LINES)

    def test_no_store_data_prevents_update(self, web_app):
        with pytest_raises_prevent_update():
            self._fn(web_app)(1, None)

    def test_no_records_prevents_update(self, web_app):
        with pytest_raises_prevent_update():
            self._fn(web_app)(1, {"station_records": []})

    def test_detects_lines_from_id_prefixes(self, web_app):
        store = {
            "station_records": _records(
                [("18-001A", ""), ("18-002B", ""), ("22-010A", "")]
            )
        }
        new_store = self._fn(web_app)(1, store)
        lines = {r["ID"]: r["Line"] for r in new_store["station_records"]}
        assert lines["18-001A"] == "L18"
        assert lines["18-002B"] == "L18"
        assert lines["22-010A"] == "L22"
        assert new_store["n_lines"] == 2
        assert new_store["line_counts"] == {"L18": 2, "L22": 1}


def pytest_raises_prevent_update():
    import pytest

    return pytest.raises(PreventUpdate)


# ── 5. apply_renames ─────────────────────────────────────────────────────────


class TestApplyRenames:
    def _fn(self, web_app):
        return _cb_by_input(
            web_app, f"{IDs.STORE_DATA}.data", IDs.BTN_LINES_RENAME_APPLY
        )

    def test_no_clicks_prevents_update(self, web_app):
        with pytest_raises_prevent_update():
            self._fn(web_app)(None, ["X"], [{"index": "L1"}], STORE_TWO_LINES)

    def test_no_store_data_prevents_update(self, web_app):
        with pytest_raises_prevent_update():
            self._fn(web_app)(1, ["X"], [{"index": "L1"}], None)

    def test_no_input_ids_prevents_update(self, web_app):
        with pytest_raises_prevent_update():
            self._fn(web_app)(1, [], [], STORE_TWO_LINES)

    def test_no_actual_renames_prevents_update(self, web_app):
        # new name equals original -> no rename recorded -> PreventUpdate
        with pytest_raises_prevent_update():
            self._fn(web_app)(
                1,
                ["L1", "L2"],
                [{"index": "L1"}, {"index": "L2"}],
                STORE_TWO_LINES,
            )

    def test_applies_rename_and_rebuilds_counts(self, web_app):
        new_store = self._fn(web_app)(
            1,
            ["West Line", "L2"],
            [{"index": "L1"}, {"index": "L2"}],
            STORE_TWO_LINES,
        )
        lines = {r["ID"]: r["Line"] for r in new_store["station_records"]}
        assert lines["A1"] == "West Line"
        assert lines["B1"] == "L2"
        assert new_store["line_counts"] == {"West Line": 2, "L2": 2}
        assert new_store["n_lines"] == 2


# ── 6. build_lines_panel ─────────────────────────────────────────────────────


class TestBuildLinesPanel:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.LINES_PANEL_BODY}.children")

    def test_closed_no_update(self, web_app):
        assert (
            self._fn(web_app)(False, STORE_TWO_LINES, "folder", None, "dark")
            is no_update
        )

    def test_no_store_data_shows_hint(self, web_app):
        out = self._fn(web_app)(True, None, "folder", None, "dark")
        assert "Load survey data first." in out.children

    def test_no_lines_shows_mode_specific_empty_hint(self, web_app):
        store = {"station_records": _records([("A1", "—")])}
        out = self._fn(web_app)(True, store, "auto", None, "dark")
        assert "Detect Lines" in out.children

    def test_builds_rows_for_each_line(self, web_app):
        out = self._fn(web_app)(True, STORE_TWO_LINES, "folder", None, "dark")
        rows_container = out[0]
        assert len(rows_container.children) == 2

    def test_edit_mode_shows_rename_inputs_and_footer(self, web_app):
        out = self._fn(web_app)(
            True,
            STORE_TWO_LINES,
            "edit",
            {"active": ["L1", "L2"]},
            "dark",
        )
        assert len(out) == 2  # rows container + footer hint
        assert "Apply Renames" in out[1].children

    def test_auto_mode_shows_detected_footer(self, web_app):
        out = self._fn(web_app)(
            True,
            STORE_TWO_LINES,
            "auto",
            {"active": ["L1", "L2"]},
            "dark",
        )
        assert "detected from station ID prefixes" in out[1].children


# ── 7. sync_active_lines ─────────────────────────────────────────────────────


class TestSyncActiveLines:
    def _fn(self, web_app):
        return _cb_by_input(web_app, f"{IDs.STORE_ACTIVE_LINES}.data", "line-switch")

    def test_no_trigger_no_update(self, monkeypatch, web_app):
        import pycsamt.app.web.callbacks.lines as lines_mod

        class _FakeCtx:
            inputs_list = [[]]

        monkeypatch.setattr(lines_mod, "ctx", _FakeCtx())
        out = self._fn(web_app)([], {"active": [], "all": []})
        assert out is no_update

    def test_builds_active_list_from_switch_values(self, monkeypatch, web_app):
        import pycsamt.app.web.callbacks.lines as lines_mod

        class _FakeCtx:
            inputs_list = [
                [
                    {
                        "id": {"type": "line-switch", "index": "L1"},
                        "value": True,
                    },
                    {
                        "id": {"type": "line-switch", "index": "L2"},
                        "value": False,
                    },
                ]
            ]

        monkeypatch.setattr(lines_mod, "ctx", _FakeCtx())
        out = self._fn(web_app)(
            [True, False], {"active": ["L1", "L2"], "all": ["L1", "L2"]}
        )
        assert out == {"active": ["L1"], "all": ["L1", "L2"]}


# ── 8. all_or_none ───────────────────────────────────────────────────────────


class TestAllOrNone:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "line-switch", IDs.BTN_LINES_ALL)

    def test_no_switch_ids_no_update(self, monkeypatch, web_app):
        import pycsamt.app.web.callbacks.lines as lines_mod

        monkeypatch.setattr(
            lines_mod,
            "ctx",
            type("C", (), {"triggered_id": IDs.BTN_LINES_ALL})(),
        )
        out = self._fn(web_app)(1, None, [])
        assert out is no_update

    def test_all_button_sets_all_true(self, monkeypatch, web_app):
        import pycsamt.app.web.callbacks.lines as lines_mod

        monkeypatch.setattr(
            lines_mod,
            "ctx",
            type("C", (), {"triggered_id": IDs.BTN_LINES_ALL})(),
        )
        out = self._fn(web_app)(
            1,
            None,
            [
                {"type": "line-switch", "index": "L1"},
                {"type": "line-switch", "index": "L2"},
            ],
        )
        assert out == [True, True]

    def test_none_button_sets_all_false(self, monkeypatch, web_app):
        import pycsamt.app.web.callbacks.lines as lines_mod

        monkeypatch.setattr(
            lines_mod,
            "ctx",
            type("C", (), {"triggered_id": IDs.BTN_LINES_NONE})(),
        )
        out = self._fn(web_app)(None, 1, [{"type": "line-switch", "index": "L1"}])
        assert out == [False]

    def test_unrelated_trigger_no_update(self, monkeypatch, web_app):
        import pycsamt.app.web.callbacks.lines as lines_mod

        monkeypatch.setattr(
            lines_mod,
            "ctx",
            type("C", (), {"triggered_id": "something-else"})(),
        )
        out = self._fn(web_app)(1, 1, [{"type": "line-switch", "index": "L1"}])
        assert out is no_update
