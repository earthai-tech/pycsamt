# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Additional tests for pycsamt.app.agent_master.callbacks.edi pure helpers
(complementing the basic coverage already in test_app.py)."""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")


class TestColor:
    def test_wraps_around_palette(self):
        from pycsamt.app.agent_master.callbacks.edi import (
            _LINE_COLORS,
            _color,
        )

        assert _color(0) == _LINE_COLORS[0]
        assert _color(len(_LINE_COLORS)) == _LINE_COLORS[0]


class TestFolderToLines:
    def test_groups_by_subfolder(self, tmp_path):
        from pycsamt.app.agent_master.callbacks.edi import _folder_to_lines

        (tmp_path / "L1").mkdir()
        (tmp_path / "L2").mkdir()
        (tmp_path / "L1" / "a.edi").write_text("x")
        (tmp_path / "L1" / "b.edi").write_text("x")
        (tmp_path / "L2" / "c.edi").write_text("x")
        groups = _folder_to_lines(str(tmp_path))
        assert groups["L1"] == sorted(groups["L1"])
        assert len(groups["L1"]) == 2
        assert len(groups["L2"]) == 1

    def test_uppercase_extension_is_found(self, tmp_path):
        from pycsamt.app.agent_master.callbacks.edi import _folder_to_lines

        (tmp_path / "L1").mkdir()
        (tmp_path / "L1" / "a.EDI").write_text("x")
        groups = _folder_to_lines(str(tmp_path))
        assert len(groups["L1"]) == 1

    def test_empty_folder_returns_empty_dict(self, tmp_path):
        from pycsamt.app.agent_master.callbacks.edi import _folder_to_lines

        assert _folder_to_lines(str(tmp_path)) == {}


class TestDetectLinesToFiles:
    def test_groups_by_station_prefix(self, tmp_path):
        from pycsamt.app.agent_master.callbacks.edi import (
            _detect_lines_to_files,
        )

        (tmp_path / "22-001.edi").write_text("x")
        (tmp_path / "22-002.edi").write_text("x")
        (tmp_path / "26-001.edi").write_text("x")
        groups = _detect_lines_to_files(str(tmp_path))
        assert set(groups) == {"L22", "L26"}
        assert len(groups["L22"]) == 2

    def test_empty_folder_returns_empty_dict(self, tmp_path):
        from pycsamt.app.agent_master.callbacks.edi import (
            _detect_lines_to_files,
        )

        assert _detect_lines_to_files(str(tmp_path)) == {}


class TestBuildLinesPanel:
    def test_empty_groups_shows_placeholder(self):
        from pycsamt.app.agent_master.callbacks.edi import (
            _build_lines_panel,
        )

        panel = _build_lines_panel({})
        assert "No EDI files found" in str(panel[0])

    def test_non_editable_shows_span_names(self):
        from pycsamt.app.agent_master.callbacks.edi import (
            _build_lines_panel,
        )

        panel = _build_lines_panel({"L1": ["a.edi", "b.edi"]})
        text = str(panel)
        assert "L1" in text
        assert "2 EDI" in text
        assert "1 line(s) detected" in text

    def test_editable_shows_rename_inputs(self):
        from pycsamt.app.agent_master.callbacks.edi import (
            _build_lines_panel,
        )

        panel = _build_lines_panel({"L1": ["a.edi"]}, editable=True)
        text = str(panel)
        assert "am-line-rename" in text


class TestErrDiv:
    def test_contains_message_and_icon(self):
        from pycsamt.app.agent_master.callbacks.edi import _err_div

        text = str(_err_div("Something broke"))
        assert "Something broke" in text
        assert "bi-exclamation-triangle" in text


class TestRegisterEdi:
    def test_register_edi_is_callable(self):
        from pycsamt.app.agent_master.callbacks.edi import register_edi

        assert callable(register_edi)

    def test_expected_outputs_wired(self):
        from pycsamt.app.agent_master._ids import IDs
        from pycsamt.app.agent_master.app import create_app

        app = create_app()
        cb_outputs = str(app.callback_map)
        assert IDs.STORE_EDI in cb_outputs
        assert IDs.LINES_PANEL in cb_outputs


class TestDetectFromIds:
    def test_groups_by_alpha_numeric_prefix(self):
        from pycsamt.app.agent_master.callbacks.edi import _detect_from_ids

        groups = _detect_from_ids(["L18-001", "L18-002", "L22-001"])
        assert set(groups) == {"L18", "L22"}
        assert groups["L18"] == ["L18-001", "L18-002"]

    def test_numeric_only_prefix_gets_l_prefix(self):
        from pycsamt.app.agent_master.callbacks.edi import _detect_from_ids

        groups = _detect_from_ids(["22-001", "22-002"])
        assert set(groups) == {"L22"}

    def test_no_prefix_match_falls_back_to_first_four_chars(self):
        from pycsamt.app.agent_master.callbacks.edi import _detect_from_ids

        groups = _detect_from_ids(["####station"])
        assert list(groups.keys()) == ["####"]


def _unwrap(entry):
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def _find(agent_app, input_id, output_hint):
    matches = [
        k
        for k, entry in agent_app.callback_map.items()
        if entry["inputs"]
        and entry["inputs"][0]["id"] == input_id
        and output_hint in k
    ]
    assert len(matches) == 1, (input_id, output_hint, matches)
    return _unwrap(agent_app.callback_map[matches[0]])


class TestToggleEdiCanvas:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        matches = [
            k
            for k, entry in agent_app.callback_map.items()
            if entry["inputs"]
            and entry["inputs"][0]["id"] == IDs.BTN_LOAD_EDI
            and k == f"{IDs.CANVAS_EDI}.is_open"
        ]
        assert len(matches) == 1, matches
        return _unwrap(agent_app.callback_map[matches[0]])

    def test_flips_open_state(self, agent_app):
        fn = self._fn(agent_app)
        assert fn(1, None, False) is True
        assert fn(1, None, True) is False


class TestHandleFolderStore:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(agent_app, IDs.FOLDER_STORE, IDs.EDI_PATH_INPUT)

    def test_prevent_update_without_store_data(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(None)

    def test_prevent_update_without_filenames(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn({"filenames": [], "contents": []})

    def test_writes_nested_files_to_temp_dir(self, agent_app):
        import base64
        from pathlib import Path

        content = "data:text/plain;base64," + base64.b64encode(
            b"edi-content"
        ).decode()
        fn = self._fn(agent_app)
        tmp = fn(
            {
                "filenames": ["L18PLT/18-001A.edi"],
                "contents": [content],
            }
        )
        out = Path(tmp) / "L18PLT" / "18-001A.edi"
        assert out.read_bytes() == b"edi-content"


class TestToggleModeBtns:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(agent_app, IDs.LINES_MODE, IDs.BTN_DETECT_LINES)

    def test_auto_mode_shows_detect_button(self, agent_app):
        fn = self._fn(agent_app)
        det, ren = fn("auto")
        assert det == {"display": "block"}
        assert ren == {"display": "none"}

    def test_edit_mode_shows_rename_button(self, agent_app):
        fn = self._fn(agent_app)
        det, ren = fn("edit")
        assert det == {"display": "none"}
        assert ren == {"display": "block"}

    def test_folder_mode_hides_both(self, agent_app):
        fn = self._fn(agent_app)
        det, ren = fn("folder")
        assert det == {"display": "none"}
        assert ren == {"display": "none"}


class TestUpdateLinesPanel:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(agent_app, IDs.EDI_PATH_INPUT, IDs.LINES_PANEL)

    def test_prevent_update_without_path(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(None, None, "folder")

    def test_nonexistent_path_shows_error(self, agent_app, tmp_path):
        from dash import no_update

        fn = self._fn(agent_app)
        panel, status = fn(str(tmp_path / "ghost"), None, "folder")
        assert panel is no_update
        assert "Path not found" in str(status)

    def test_single_file_path_groups_as_default(
        self, agent_app, tmp_path
    ):
        f = tmp_path / "station.edi"
        f.write_text("x")
        fn = self._fn(agent_app)
        panel, status = fn(str(f), None, "folder")
        assert "1 EDI file(s) in 1 line(s)" in str(status)

    def test_folder_path_groups_by_subfolder(self, agent_app, tmp_path):
        sub = tmp_path / "L1"
        sub.mkdir()
        (sub / "a.edi").write_text("x")
        (sub / "b.edi").write_text("x")
        fn = self._fn(agent_app)
        panel, status = fn(str(tmp_path), None, "folder")
        assert "2 EDI file(s) in 1 line(s)" in str(status)

    def test_auto_mode_detect_click_regroups_by_station_id(
        self, agent_app, tmp_path
    ):
        import dash._callback_context as cc
        from dash._utils import AttributeDict
        from pycsamt.app.agent_master._ids import IDs

        sub = tmp_path / "somefolder"
        sub.mkdir()
        (sub / "22-001.edi").write_text("x")
        (sub / "22-002.edi").write_text("x")

        cc.context_value.set(
            AttributeDict(
                triggered_inputs=[
                    {
                        "prop_id": f"{IDs.BTN_DETECT_LINES}.n_clicks",
                        "value": 1,
                    }
                ]
            )
        )
        try:
            fn = self._fn(agent_app)
            panel, status = fn(str(tmp_path), 1, "auto")
            assert "L22" in str(panel)
        finally:
            cc.context_value.set({})


class TestHandleUpload:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(agent_app, IDs.UPLOAD_EDI, IDs.EDI_PATH_INPUT)

    def test_prevent_update_without_contents(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(None, "a.edi")

    def test_writes_single_upload_to_temp_dir(self, agent_app):
        import base64
        from pathlib import Path

        content = "data:text/plain;base64," + base64.b64encode(
            b"hello"
        ).decode()
        fn = self._fn(agent_app)
        tmp = fn(content, "a.edi")
        assert (Path(tmp) / "a.edi").read_bytes() == b"hello"

    def test_writes_multiple_uploads(self, agent_app):
        import base64
        from pathlib import Path

        c1 = "data:text/plain;base64," + base64.b64encode(b"one").decode()
        c2 = "data:text/plain;base64," + base64.b64encode(b"two").decode()
        fn = self._fn(agent_app)
        tmp = fn([c1, c2], ["a.edi", "b.edi"])
        assert (Path(tmp) / "a.edi").read_bytes() == b"one"
        assert (Path(tmp) / "b.edi").read_bytes() == b"two"


class TestConfirmLoad:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(agent_app, IDs.BTN_LOAD_CONFIRM, IDs.STORE_EDI)

    def test_prevent_update_without_click(self, agent_app, tmp_path):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(None, str(tmp_path), "folder", [])

    def test_prevent_update_without_path(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(1, "", "folder", [])

    def test_prevent_update_for_nonexistent_path(self, agent_app, tmp_path):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(1, str(tmp_path / "ghost"), "folder", [])

    def test_single_file_load(self, agent_app, tmp_path):
        f = tmp_path / "one.edi"
        f.write_text("x")
        fn = self._fn(agent_app)
        store, cls, text, is_open = fn(1, str(f), "folder", [])
        assert store["groups"] == {"Default": [str(f)]}
        assert is_open is False
        assert "1 EDI" in text

    def test_folder_mode_load(self, agent_app, tmp_path):
        sub = tmp_path / "L1"
        sub.mkdir()
        (sub / "a.edi").write_text("x")
        fn = self._fn(agent_app)
        store, _cls, _text, _is_open = fn(1, str(tmp_path), "folder", [])
        assert store["mode"] == "folder"
        assert "L1" in store["groups"]

    def test_auto_mode_load_groups_by_station_prefix(
        self, agent_app, tmp_path
    ):
        sub = tmp_path / "any"
        sub.mkdir()
        (sub / "22-001.edi").write_text("x")
        fn = self._fn(agent_app)
        store, _cls, _text, _is_open = fn(1, str(tmp_path), "auto", [])
        assert "L22" in store["groups"]

    def test_edit_mode_applies_rename_values(self, agent_app, tmp_path):
        sub = tmp_path / "OrigName"
        sub.mkdir()
        (sub / "a.edi").write_text("x")
        fn = self._fn(agent_app)
        store, _cls, _text, _is_open = fn(
            1, str(tmp_path), "edit", ["RenamedLine"]
        )
        assert "RenamedLine" in store["groups"]
        assert "OrigName" not in store["groups"]
