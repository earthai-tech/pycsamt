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
