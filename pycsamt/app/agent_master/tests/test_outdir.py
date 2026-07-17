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
