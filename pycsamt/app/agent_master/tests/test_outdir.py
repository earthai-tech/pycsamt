# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.agent_master.callbacks.outdir."""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")


class TestListDirs:
    def test_lists_only_subdirectories(self, tmp_path):
        from pycsamt.app.agent_master.callbacks.outdir import _list_dirs

        (tmp_path / "b_dir").mkdir()
        (tmp_path / "a_dir").mkdir()
        (tmp_path / "file.txt").write_text("x")
        assert _list_dirs(str(tmp_path)) == ["a_dir", "b_dir"]

    def test_hidden_dirs_are_excluded(self, tmp_path):
        from pycsamt.app.agent_master.callbacks.outdir import _list_dirs

        (tmp_path / ".hidden").mkdir()
        (tmp_path / "visible").mkdir()
        assert _list_dirs(str(tmp_path)) == ["visible"]

    def test_nonexistent_path_raises_file_not_found(self, tmp_path):
        # only PermissionError is swallowed; a missing path propagates.
        from pycsamt.app.agent_master.callbacks.outdir import _list_dirs

        with pytest.raises(FileNotFoundError):
            _list_dirs(str(tmp_path / "nope"))

    def test_permission_error_returns_empty(self, monkeypatch, tmp_path):
        from pycsamt.app.agent_master.callbacks.outdir import _list_dirs

        def boom(path):
            raise PermissionError("denied")

        monkeypatch.setattr("os.listdir", boom)
        assert _list_dirs(str(tmp_path)) == []


class TestDirListing:
    def test_empty_dir_shows_placeholder(self, tmp_path):
        from pycsamt.app.agent_master.callbacks.outdir import _dir_listing

        children = _dir_listing(str(tmp_path))
        assert "No subdirectories" in str(children[0])

    def test_populated_dir_lists_buttons(self, tmp_path):
        from pycsamt.app.agent_master.callbacks.outdir import _dir_listing

        (tmp_path / "sub1").mkdir()
        (tmp_path / "sub2").mkdir()
        children = _dir_listing(str(tmp_path))
        assert len(children) == 2
        text = str(children[0]) + str(children[1])
        assert "sub1" in text
        assert "sub2" in text


class TestRegisterOutdir:
    def test_register_outdir_is_callable(self):
        from pycsamt.app.agent_master.callbacks.outdir import (
            register_outdir,
        )

        assert callable(register_outdir)

    def test_expected_outputs_wired(self):
        from pycsamt.app.agent_master._ids import IDs
        from pycsamt.app.agent_master.app import create_app

        app = create_app()
        cb_outputs = str(app.callback_map)
        assert IDs.MODAL_OUTPUT_BROWSE in cb_outputs
        assert IDs.OUTPUT_DIR in cb_outputs


def _unwrap(entry):
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def _find(agent_app, input_id, output_hint):
    matches = [
        k
        for k, entry in agent_app.callback_map.items()
        if entry["inputs"] and entry["inputs"][0]["id"] == input_id and output_hint in k
    ]
    assert len(matches) == 1, (input_id, output_hint, matches)
    return _unwrap(agent_app.callback_map[matches[0]])


def _set_pattern_trigger(component_id_dict, prop="n_clicks"):
    import json

    import dash._callback_context as cc
    from dash._utils import AttributeDict

    prop_id = json.dumps(component_id_dict, sort_keys=True) + f".{prop}"
    cc.context_value.set(
        AttributeDict(triggered_inputs=[{"prop_id": prop_id, "value": 1}])
    )


def _clear_trigger():
    import dash._callback_context as cc

    cc.context_value.set({})


class TestOpenOutputModal:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(agent_app, IDs.BTN_OUTPUT_BROWSE, IDs.MODAL_OUTPUT_BROWSE)

    def test_prevent_update_without_click(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(None, None)

    def test_opens_at_current_value_when_valid_dir(self, agent_app, tmp_path):
        fn = self._fn(agent_app)
        is_open, store, path_text, listing = fn(1, str(tmp_path))
        assert is_open is True
        assert store == {"path": str(tmp_path)}
        assert path_text == str(tmp_path)

    def test_falls_back_to_home_for_invalid_path(self, agent_app):
        import os

        fn = self._fn(agent_app)
        is_open, store, path_text, _listing = fn(1, "/definitely/not/a/dir")
        assert is_open is True
        assert store == {"path": os.path.expanduser("~")}


class TestEnterSubdir:
    def teardown_method(self, _method):
        _clear_trigger()

    def _fn(self, agent_app):
        return _find(
            agent_app,
            '{"index":["ALL"],"name":["ALL"],"type":"am-outdir-entry"}',
            "",
        )

    def test_prevent_update_without_clicks(self, agent_app, tmp_path):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn([None], {"path": str(tmp_path)})

    def test_prevent_update_for_nonexistent_target(self, agent_app, tmp_path):
        from dash.exceptions import PreventUpdate

        _set_pattern_trigger({"index": 0, "name": "ghost", "type": "am-outdir-entry"})
        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn([1], {"path": str(tmp_path)})

    def test_navigates_into_subdir(self, agent_app, tmp_path):
        (tmp_path / "child").mkdir()
        _set_pattern_trigger({"index": 0, "name": "child", "type": "am-outdir-entry"})
        fn = self._fn(agent_app)
        store, path_text, _listing, status = fn([1], {"path": str(tmp_path)})
        assert store == {"path": str(tmp_path / "child")}
        assert path_text == str(tmp_path / "child")
        assert status == ""


class TestGoUp:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(agent_app, IDs.BTN_OUTPUT_UP, IDs.OUTPUT_BROWSE_STORE)

    def test_prevent_update_without_click(self, agent_app, tmp_path):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(None, {"path": str(tmp_path)})

    def test_prevent_update_at_filesystem_root(self, agent_app, tmp_path):
        from dash.exceptions import PreventUpdate

        root = tmp_path.anchor  # e.g. "C:\\" on Windows, "/" on POSIX
        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(1, {"path": root})

    def test_navigates_to_parent(self, agent_app, tmp_path):
        child = tmp_path / "child"
        child.mkdir()
        fn = self._fn(agent_app)
        store, path_text, _listing = fn(1, {"path": str(child)})
        assert store == {"path": str(tmp_path)}
        assert path_text == str(tmp_path)


class TestMkdir:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(agent_app, IDs.BTN_OUTPUT_MKDIR, IDs.OUTPUT_BROWSE_STORE)

    def test_prevent_update_without_click(self, agent_app, tmp_path):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(None, "newdir", {"path": str(tmp_path)})

    def test_prevent_update_without_name(self, agent_app, tmp_path):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(1, "   ", {"path": str(tmp_path)})

    def test_creates_directory(self, agent_app, tmp_path):
        fn = self._fn(agent_app)
        store, path_text, _listing, status, cleared = fn(
            1, "newdir", {"path": str(tmp_path)}
        )
        assert store == {"path": str(tmp_path / "newdir")}
        assert path_text == str(tmp_path / "newdir")
        assert (tmp_path / "newdir").is_dir()
        assert "Created: newdir" in str(status)
        assert cleared == ""

    def test_mkdir_failure_reports_error(self, agent_app, tmp_path, monkeypatch):
        from pathlib import Path

        from dash import no_update

        def boom(self, parents=True, exist_ok=True):
            raise OSError("disk full")

        monkeypatch.setattr(Path, "mkdir", boom)
        fn = self._fn(agent_app)
        store, path_text, listing, status, name_back = fn(
            1, "newdir", {"path": str(tmp_path)}
        )
        assert store is no_update
        assert path_text is no_update
        assert listing is no_update
        assert "disk full" in str(status)
        assert name_back == "newdir"


class TestConfirmOutputDir:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(agent_app, IDs.BTN_OUTPUT_CONFIRM, IDs.OUTPUT_DIR)

    def test_prevent_update_without_click(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(None, {"path": "/tmp"})

    def test_confirms_and_closes(self, agent_app):
        fn = self._fn(agent_app)
        path, is_open = fn(1, {"path": "/some/dir"})
        assert path == "/some/dir"
        assert is_open is False
