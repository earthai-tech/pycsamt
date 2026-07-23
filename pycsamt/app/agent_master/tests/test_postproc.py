# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.agent_master.callbacks.postproc pure helpers."""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("dash", reason="dash required")


def _dup_edi(tmp_path: Path, src: Path, stem: str) -> Path:
    dst = tmp_path / f"{stem}.edi"
    dst.write_text(src.read_text(encoding="utf-8"), encoding="utf-8")
    return dst


def _sites(tmp_path: Path, simulated_edi: Path):
    from pycsamt.seg.edi import EDIFile
    from pycsamt.site.base import Site, Sites

    p1 = _dup_edi(tmp_path, simulated_edi, "T01")
    p2 = _dup_edi(tmp_path, simulated_edi, "T02")
    s1 = Site(EDIFile(p1))
    s2 = Site(EDIFile(p2))
    return Sites([s1.edi, s2.edi])


class TestExport:
    def test_writes_edi_files(self, tmp_path, simulated_edi):
        from pycsamt.app.agent_master.callbacks.postproc import _export

        sites = _sites(tmp_path, simulated_edi)
        dest = tmp_path / "out"
        paths = _export(sites, dest)
        assert len(paths) == 2
        assert all(Path(p).suffix == ".edi" for p in paths)


class TestCountValidImpedance:
    def test_counts_stations_with_z(self, tmp_path, simulated_edi):
        from pycsamt.app.agent_master.callbacks.postproc import (
            _count_valid_impedance,
        )

        sites = _sites(tmp_path, simulated_edi)
        assert _count_valid_impedance(sites) == 2


class TestValidateExport:
    def test_valid_export_reports_all_stations(self, tmp_path, simulated_edi):
        from pycsamt.app.agent_master.callbacks.postproc import (
            _export,
            _validate_export,
        )

        sites = _sites(tmp_path, simulated_edi)
        dest = tmp_path / "out2"
        _export(sites, dest)
        n_valid, n_total, detail = _validate_export(dest)
        assert n_valid == 2
        assert n_total == 2
        assert detail == ""

    def test_empty_folder_reports_zero_valid(self, tmp_path):
        from pycsamt.app.agent_master.callbacks.postproc import (
            _validate_export,
        )

        dest = tmp_path / "empty_out"
        dest.mkdir()
        n_valid, n_total, detail = _validate_export(dest)
        assert n_valid == 0
        assert n_total == 0

    def test_exception_reports_validation_error(self, tmp_path, monkeypatch):
        from pycsamt.app.agent_master.callbacks import postproc as pp_mod

        def boom(*a, **k):
            raise RuntimeError("boom")

        # _validate_export resolves ensure_sites via a local import, so
        # patch it on the actual target module.
        import pycsamt.emtools._core as core_mod

        monkeypatch.setattr(core_mod, "ensure_sites", boom)
        n_valid, n_total, detail = pp_mod._validate_export(tmp_path)
        assert n_valid == -1
        assert n_total == 0
        assert "validation error" in detail


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


class TestOpenModal:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(agent_app, IDs.STORE_POSTPROC, IDs.MODAL_POSTPROC)

    def test_prevent_update_without_jid(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(None)
        with pytest.raises(PreventUpdate):
            fn({})

    def test_opens_with_workflow_message(self, agent_app):
        fn = self._fn(agent_app)
        is_open, msg = fn({"jid": "j1", "workflow": "static-shift"})
        assert is_open is True
        assert "static-shift" in str(msg)


class TestShowExportPath:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(
            agent_app, IDs.BTN_POSTPROC_EXPORT, IDs.POSTPROC_COLLAPSE
        )

    def test_none_clicks_collapsed(self, agent_app):
        fn = self._fn(agent_app)
        assert fn(None) is False

    def test_odd_clicks_open(self, agent_app):
        fn = self._fn(agent_app)
        assert fn(1) is True

    def test_even_clicks_closed(self, agent_app):
        fn = self._fn(agent_app)
        assert fn(2) is False


class TestApplyToSession:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(agent_app, IDs.BTN_POSTPROC_APPLY, IDs.STORE_EDI)

    def test_prevent_update_without_click(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(None, {}, {}, {})

    def test_sites_not_in_memory_returns_error(self, agent_app, monkeypatch):
        import pycsamt.app.agent_master.callbacks.chat as chat_mod
        from dash import no_update

        monkeypatch.setattr(chat_mod, "_CORR_CACHE", {})
        fn = self._fn(agent_app)
        edi, status, is_open = fn(1, {"jid": "missing"}, {}, {"path": "/x"})
        assert edi is no_update
        assert is_open is True
        assert "not in memory" in str(status)

    def test_export_exception_returns_error(
        self, agent_app, monkeypatch, tmp_path
    ):
        import pycsamt.app.agent_master.callbacks.chat as chat_mod
        import pycsamt.app.agent_master.callbacks.postproc as pp_mod
        from dash import no_update

        monkeypatch.setattr(chat_mod, "_CORR_CACHE", {"j1": object()})
        monkeypatch.setattr(
            pp_mod, "_count_valid_impedance", lambda sites: 2
        )

        def boom(sites, dest):
            raise RuntimeError("disk full")

        monkeypatch.setattr(pp_mod, "_export", boom)
        fn = self._fn(agent_app)
        edi, status, is_open = fn(
            1, {"jid": "j1", "output_dir": str(tmp_path)}, {}, {}
        )
        assert edi is no_update
        assert is_open is True
        assert "disk full" in str(status)

    def test_invalid_export_leaves_session_untouched(
        self, agent_app, monkeypatch, tmp_path
    ):
        import pycsamt.app.agent_master.callbacks.chat as chat_mod
        import pycsamt.app.agent_master.callbacks.postproc as pp_mod
        from dash import no_update

        monkeypatch.setattr(chat_mod, "_CORR_CACHE", {"j1": object()})
        monkeypatch.setattr(
            pp_mod, "_count_valid_impedance", lambda sites: 0
        )
        monkeypatch.setattr(
            pp_mod, "_export", lambda sites, dest: ["a.edi"]
        )
        monkeypatch.setattr(
            pp_mod, "_validate_export", lambda dest: (0, 0, "no Z blocks")
        )
        fn = self._fn(agent_app)
        edi, status, is_open = fn(
            1, {"jid": "j1", "output_dir": str(tmp_path)}, {}, {"path": "old"}
        )
        assert edi is no_update
        assert is_open is True
        assert "no Z blocks" in str(status)

    def test_successful_apply_updates_session(
        self, agent_app, monkeypatch, tmp_path
    ):
        import pycsamt.app.agent_master.callbacks.chat as chat_mod
        import pycsamt.app.agent_master.callbacks.postproc as pp_mod

        cache = {"j1": object()}
        monkeypatch.setattr(chat_mod, "_CORR_CACHE", cache)
        monkeypatch.setattr(
            pp_mod, "_count_valid_impedance", lambda sites: 2
        )
        monkeypatch.setattr(
            pp_mod, "_export", lambda sites, dest: ["a.edi", "b.edi"]
        )
        monkeypatch.setattr(
            pp_mod, "_validate_export", lambda dest: (2, 2, "")
        )
        fn = self._fn(agent_app)
        edi, status, is_open = fn(
            1,
            {"jid": "j1", "output_dir": str(tmp_path), "workflow": "shift"},
            {},
            {"path": "old", "sites": "x", "groups": "y"},
        )
        assert "j1" not in cache
        assert edi["path"] != "old"
        assert "sites" not in edi
        assert "groups" not in edi
        assert is_open is False
        assert "2/2" in str(status)


class TestConfirmExport:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(agent_app, IDs.BTN_POSTPROC_OK, IDs.POSTPROC_STATUS)

    def test_prevent_update_without_click(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(None, "/dest", {})

    def test_empty_path_shows_warning(self, agent_app):
        fn = self._fn(agent_app)
        status, is_open, collapse_open = fn(1, "   ", {})
        assert "Enter a folder path" in str(status)
        assert is_open is True
        assert collapse_open is True

    def test_sites_not_in_memory_returns_error(
        self, agent_app, monkeypatch, tmp_path
    ):
        import pycsamt.app.agent_master.callbacks.chat as chat_mod

        monkeypatch.setattr(chat_mod, "_CORR_CACHE", {})
        fn = self._fn(agent_app)
        status, is_open, collapse_open = fn(
            1, str(tmp_path), {"jid": "missing"}
        )
        assert "not in memory" in str(status)
        assert is_open is True
        assert collapse_open is True

    def test_successful_export(self, agent_app, monkeypatch, tmp_path):
        import pycsamt.app.agent_master.callbacks.chat as chat_mod
        import pycsamt.app.agent_master.callbacks.postproc as pp_mod

        cache = {"j1": object()}
        monkeypatch.setattr(chat_mod, "_CORR_CACHE", cache)
        monkeypatch.setattr(
            pp_mod, "_count_valid_impedance", lambda sites: 2
        )
        monkeypatch.setattr(
            pp_mod, "_export", lambda sites, dest: ["a.edi", "b.edi"]
        )
        monkeypatch.setattr(
            pp_mod, "_validate_export", lambda dest: (2, 2, "")
        )
        fn = self._fn(agent_app)
        status, is_open, collapse_open = fn(
            1, str(tmp_path), {"jid": "j1"}
        )
        assert "j1" not in cache
        assert "Exported 2/2" in str(status)
        assert is_open is False
        assert collapse_open is False
