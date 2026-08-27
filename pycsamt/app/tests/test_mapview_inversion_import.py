# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.mapview.callbacks.inversion_import."""

from __future__ import annotations

import base64
import os
import shutil

import numpy as np
import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")


def _b64_payload(text: str) -> str:
    encoded = base64.b64encode(text.encode()).decode()
    return f"data:text/plain;base64,{encoded}"


def _pcsf_file_payload(tmp_path, suffix, *, kind="grid2d", name="demo"):
    """Build a small real PCSF/PCSM file and return its (filename, b64
    data-URI content) pair, ready to feed a staging/confirm callback."""
    from pycsamt.format import write_pcsf, write_pcsm
    from pycsamt.format.schema import (
        Grid2DGeometry,
        Grid3DGeometry,
        PCSFModel,
        StationTable,
        UnstructuredMeshGeometry,
    )

    if kind == "grid3d":
        nxg, nyg, nzg = 3, 3, 3
        model = PCSFModel(
            geometry=Grid3DGeometry(
                x=np.linspace(0, 300, nxg),
                y=np.linspace(0, 300, nyg),
                z=np.linspace(0, 200, nzg),
                x_nodes=np.linspace(-50, 350, nxg + 1),
                y_nodes=np.linspace(-50, 350, nyg + 1),
            ),
            resistivity=np.full((nzg, nyg, nxg), 100.0),
            source_backend="modem3d",
            stations=StationTable(
                name=["23-01-001", "23-01-002", "23-01-003"],
                x=np.linspace(-150, 150, 3),
                y=np.zeros(3),
                z=np.array([12.0, 11.5, 13.0]),
            ),
        )
    elif kind == "mesh_unstructured":
        model = PCSFModel(
            geometry=UnstructuredMeshGeometry(
                nodes=np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
                connectivity=np.array([[0, 1, 2]], dtype=np.int64),
                region_ids=np.array([0], dtype=np.int32),
            ),
            resistivity=np.array([100.0]),
            source_backend="mare2dem",
            stations=StationTable(
                # Interior of the (0,0)-(1,0)-(0,1) triangle, not on an
                # edge/vertex, so point-location is unambiguous.
                name=["S0"], x=np.array([0.3]), y=np.array([0.0]), z=np.array([0.0])
            ),
        )
    else:
        nx, nz = 4, 3
        x = np.linspace(0, 300, nx)
        z = np.linspace(0, 200, nz)
        model = PCSFModel(
            geometry=Grid2DGeometry(x=x, z=z),
            resistivity=np.full((nz, nx), 100.0),
            source_backend="generic",
            stations=StationTable(
                name=[f"S{i}" for i in range(nx)],
                x=x,
                y=np.zeros(nx),
                z=np.array([12.0, 11.5, 13.0, 10.0]),
            ),
        )

    filename = f"{name}{suffix}"
    path = tmp_path / filename
    if suffix == ".pcsf":
        write_pcsf(model, path)
    else:
        write_pcsm(model, path)
    content = (
        f"data:application/octet-stream;base64,"
        f"{base64.b64encode(path.read_bytes()).decode()}"
    )
    return filename, content


class TestIsPcsfName:
    def test_recognises_pcsf_pcsm_and_gz(self):
        from pycsamt.app.mapview.callbacks.inversion_import import _is_pcsf_name

        assert _is_pcsf_name("model.pcsf")
        assert _is_pcsf_name("MODEL.PCSF")
        assert _is_pcsf_name("model.pcsm")
        assert _is_pcsf_name("model.pcsm.gz")

    def test_rejects_other_extensions(self):
        from pycsamt.app.mapview.callbacks.inversion_import import _is_pcsf_name

        assert not _is_pcsf_name("model.rho")
        assert not _is_pcsf_name("model.dat")
        assert not _is_pcsf_name("model.gz")


class TestDecodeToTempfile:
    def test_writes_single_file_preserving_extension(self):
        from pycsamt.app.mapview.callbacks.inversion_import import (
            _decode_to_tempfile,
        )

        path = _decode_to_tempfile("result.pcsm.gz", _b64_payload("payload"))
        try:
            assert path.endswith("result.pcsm.gz")
            with open(path) as fh:
                assert fh.read() == "payload"
        finally:
            shutil.rmtree(os.path.dirname(path), ignore_errors=True)

    def test_sanitizes_nested_path_to_basename(self):
        from pycsamt.app.mapview.callbacks.inversion_import import (
            _decode_to_tempfile,
        )

        path = _decode_to_tempfile("sub/dir/model.pcsf", _b64_payload("x"))
        try:
            assert os.path.basename(path) == "model.pcsf"
            assert not os.path.isdir(os.path.join(os.path.dirname(path), "sub"))
        finally:
            shutil.rmtree(os.path.dirname(path), ignore_errors=True)


class TestPeekKindSafe:
    def test_grid2d_pcsf_reports_kind(self, tmp_path):
        from pycsamt.app.mapview.callbacks.inversion_import import _peek_kind_safe

        _, content = _pcsf_file_payload(tmp_path, ".pcsf", kind="grid2d")
        assert _peek_kind_safe("demo.pcsf", content) == "grid2d"

    def test_grid3d_pcsf_reports_kind(self, tmp_path):
        from pycsamt.app.mapview.callbacks.inversion_import import _peek_kind_safe

        _, content = _pcsf_file_payload(tmp_path, ".pcsf", kind="grid3d")
        assert _peek_kind_safe("demo.pcsf", content) == "grid3d"

    def test_bad_content_returns_none(self):
        from pycsamt.app.mapview.callbacks.inversion_import import _peek_kind_safe

        assert _peek_kind_safe("bad.pcsf", "not-a-valid-data-uri") is None


class TestRegisterInversionImport:
    def test_register_is_callable(self):
        from pycsamt.app.mapview.callbacks.inversion_import import (
            register_inversion_import,
        )

        assert callable(register_inversion_import)

    def test_expected_outputs_wired(self):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.app import create_app

        app = create_app()
        cb_outputs = str(app.callback_map)
        assert IDs.BTN_INV_CONFIRM in cb_outputs
        assert IDs.INV_STATUS in cb_outputs
        assert IDs.INV_PCSF_UPLOAD in cb_outputs
        assert IDs.INV_CANDIDATES_STORE in cb_outputs
        assert IDs.INV_CANDIDATE_PICKER in cb_outputs
        assert IDs.INV_RESOLVED_STORE in cb_outputs


def _find_callback(app, *, input_id):
    for spec in app.callback_map.values():
        if input_id in str(spec.get("inputs")):
            cb = spec["callback"]
            return getattr(cb, "__wrapped__", cb)
    raise AssertionError(f"no callback found with input {input_id!r}")


class TestCaptureCandidates:
    def _capture(self, app):
        from pycsamt.app.mapview._ids import IDs

        return _find_callback(app, input_id=IDs.INV_PCSF_UPLOAD)

    def _set_trigger(self, monkeypatch, triggered_id):
        """``ctx.triggered_id`` needs a real Dash callback dispatch to
        be populated; unit-testing the unwrapped function needs the
        module's bound ``ctx`` name replaced directly instead (see
        memory: monkeypatch the consuming module's import, not Dash's
        own ``dash.ctx``, which ``has_context`` guards outside a real
        callback)."""
        import types

        import pycsamt.app.mapview.callbacks.inversion_import as inv_mod

        monkeypatch.setattr(
            inv_mod, "ctx", types.SimpleNamespace(triggered_id=triggered_id)
        )

    def test_drop_files_filters_by_extension(self, monkeypatch):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.app import create_app

        self._set_trigger(monkeypatch, IDs.INV_PCSF_UPLOAD)
        app = create_app()
        capture = self._capture(app)
        result = capture(
            [_b64_payload("a"), _b64_payload("b")],
            {},
            ["line01.pcsf", "notes.txt"],
        )
        assert result == {"filenames": ["line01.pcsf"], "contents": [_b64_payload("a")]}

    def test_single_drop_content_is_a_bare_string(self, monkeypatch):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.app import create_app

        self._set_trigger(monkeypatch, IDs.INV_PCSF_UPLOAD)
        app = create_app()
        capture = self._capture(app)
        result = capture(_b64_payload("a"), {}, "line01.pcsm")
        assert result == {"filenames": ["line01.pcsm"], "contents": [_b64_payload("a")]}

    def test_no_upload_returns_empty(self, monkeypatch):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.app import create_app

        self._set_trigger(monkeypatch, IDs.INV_PCSF_UPLOAD)
        app = create_app()
        capture = self._capture(app)
        assert capture(None, {}, None) == {"filenames": [], "contents": []}

    def test_folder_scan_source_reads_from_folder_store(self, monkeypatch):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.app import create_app

        self._set_trigger(monkeypatch, IDs.INV_FOLDER_STORE)
        app = create_app()
        capture = self._capture(app)
        folder_data = {
            "filenames": ["line01.pcsf", "line02.pcsm.gz"],
            "contents": [_b64_payload("a"), _b64_payload("b")],
        }
        result = capture(None, folder_data, None)
        assert result == folder_data


class TestClassifyCandidates:
    def _classify(self, app):
        from pycsamt.app.mapview._ids import IDs

        return _find_callback(app, input_id=IDs.INV_CANDIDATES_STORE)

    def test_empty_candidates_hides_picker(self):
        from pycsamt.app.mapview.app import create_app

        app = create_app()
        classify = self._classify(app)
        options, value, style, count, status = classify({})
        assert options == [] and value is None
        assert style == {"display": "none"}

    def test_single_grid2d_candidate_auto_selected(self, tmp_path):
        from pycsamt.app.mapview.app import create_app

        filename, content = _pcsf_file_payload(tmp_path, ".pcsf", kind="grid2d")
        app = create_app()
        classify = self._classify(app)
        options, value, style, count, status = classify(
            {"filenames": [filename], "contents": [content]}
        )
        assert len(options) == 1
        assert options[0]["value"] == value == "0"
        assert not options[0].get("disabled")
        assert style == {"display": "block"}
        assert "Ready to import" in status

    def test_grid3d_candidate_is_selectable_and_auto_selected(self, tmp_path):
        """A native ModEM grid3d volume is a real, supported curtain
        source, same as grid2d/multiline."""
        from pycsamt.app.mapview.app import create_app

        filename, content = _pcsf_file_payload(tmp_path, ".pcsf", kind="grid3d")
        app = create_app()
        classify = self._classify(app)
        options, value, style, count, status = classify(
            {"filenames": [filename], "contents": [content]}
        )
        assert len(options) == 1
        assert options[0]["value"] == value == "0"
        assert not options[0].get("disabled")
        assert "grid3d" in options[0]["label"]
        assert "Ready to import" in status

    def test_mesh_unstructured_candidate_is_selectable_and_auto_selected(
        self, tmp_path
    ):
        """A MARE2DEM mesh is also a real, supported curtain source
        (point-location sampling on its triangulation)."""
        from pycsamt.app.mapview.app import create_app

        filename, content = _pcsf_file_payload(
            tmp_path, ".pcsf", kind="mesh_unstructured"
        )
        app = create_app()
        classify = self._classify(app)
        options, value, style, count, status = classify(
            {"filenames": [filename], "contents": [content]}
        )
        assert len(options) == 1
        assert options[0]["value"] == value == "0"
        assert not options[0].get("disabled")
        assert "mesh_unstructured" in options[0]["label"]
        assert "Ready to import" in status

    def test_unreadable_candidate_is_disabled_and_not_auto_selected(self, tmp_path):
        from pycsamt.app.mapview.app import create_app

        app = create_app()
        classify = self._classify(app)
        options, value, style, count, status = classify(
            {"filenames": ["corrupt.pcsf"], "contents": ["not-a-valid-data-uri"]}
        )
        assert options[0]["disabled"] is True
        assert "unreadable" in options[0]["label"]
        assert value is None
        assert "None of these are importable" in status

    def test_two_candidates_one_unreadable_only_good_one_selectable(self, tmp_path):
        from pycsamt.app.mapview.app import create_app

        f2, c2 = _pcsf_file_payload(tmp_path, ".pcsf", kind="grid2d", name="line01")
        app = create_app()
        classify = self._classify(app)
        options, value, style, count, status = classify(
            {
                "filenames": [f2, "corrupt.pcsf"],
                "contents": [c2, "not-a-valid-data-uri"],
            }
        )
        assert count == "2 files found"
        # exactly one importable -> auto-selected
        enabled = [o for o in options if not o.get("disabled")]
        assert len(enabled) == 1
        assert value == enabled[0]["value"]

    def test_two_selectable_candidates_requires_manual_pick(self, tmp_path):
        """grid2d + grid3d are both real, importable curtain sources --
        with more than one selectable candidate, nothing is
        auto-picked; the user must choose."""
        from pycsamt.app.mapview.app import create_app

        f2, c2 = _pcsf_file_payload(tmp_path, ".pcsf", kind="grid2d", name="line01")
        f3, c3 = _pcsf_file_payload(tmp_path, ".pcsf", kind="grid3d", name="volume")
        app = create_app()
        classify = self._classify(app)
        options, value, style, count, status = classify(
            {"filenames": [f2, f3], "contents": [c2, c3]}
        )
        enabled = [o for o in options if not o.get("disabled")]
        assert len(enabled) == 2
        assert value is None
        assert "2 importable" in status


class TestResolvePick:
    def _resolve(self, app):
        from pycsamt.app.mapview._ids import IDs

        return _find_callback(app, input_id=IDs.INV_CANDIDATE_PICKER)

    def test_resolves_selected_index_to_single_file_store(self):
        from pycsamt.app.mapview.app import create_app

        app = create_app()
        resolve = self._resolve(app)
        candidates = {
            "filenames": ["a.pcsf", "b.pcsf"],
            "contents": [_b64_payload("A"), _b64_payload("B")],
        }
        result = resolve("1", candidates)
        assert result == {"filenames": ["b.pcsf"], "contents": [_b64_payload("B")]}

    def test_none_value_resolves_empty(self):
        from pycsamt.app.mapview.app import create_app

        app = create_app()
        resolve = self._resolve(app)
        assert resolve(None, {"filenames": ["a.pcsf"], "contents": ["x"]}) == {}


class TestConfirm:
    def _confirm(self, app):
        from pycsamt.app.mapview._ids import IDs

        return _find_callback(app, input_id=IDs.BTN_INV_CONFIRM)

    def test_grid2d_pcsf_file_imports_successfully(self, tmp_path):
        from pycsamt.app.mapview.app import create_app
        from pycsamt.app.mapview.cache import get_view

        filename, content = _pcsf_file_payload(tmp_path, ".pcsf", kind="grid2d")
        app = create_app()
        confirm = self._confirm(app)
        staged = {"filenames": [filename], "contents": [content]}
        store, feedback, is_open, badge, badge_cls = confirm(
            1, staged, False, "sess-pcsf", "replace", "light"
        )
        assert store["n_stations"] == 4
        assert "Imported 4 station(s)" in feedback
        assert filename in feedback
        assert is_open is False
        assert get_view("sess-pcsf") is not None

    def test_pcsm_gz_file_imports_successfully(self, tmp_path):
        from pycsamt.app.mapview.app import create_app

        filename, content = _pcsf_file_payload(tmp_path, ".pcsm.gz", kind="grid2d")
        app = create_app()
        confirm = self._confirm(app)
        staged = {"filenames": [filename], "contents": [content]}
        store, feedback, *_ = confirm(1, staged, False, "sess-pcsm-gz", "replace", "light")
        assert store["n_stations"] == 4
        assert "Imported 4 station(s)" in feedback

    def test_grid3d_pcsf_file_imports_successfully(self, tmp_path):
        """A native ModEM grid3d volume is a real, supported source --
        this exercises the same registration + curtain-slicing path as
        pycsamt/map/tests/test_pcsf_import.py::TestLoadPcsfLinesGrid3D,
        through the actual UI confirm() callback."""
        from pycsamt.app.mapview.app import create_app

        filename, content = _pcsf_file_payload(tmp_path, ".pcsf", kind="grid3d")
        app = create_app()
        confirm = self._confirm(app)
        staged = {"filenames": [filename], "contents": [content]}
        store, feedback, is_open, badge, badge_cls = confirm(
            1, staged, False, "sess-grid3d", "replace", "light"
        )
        assert store["n_stations"] == 3
        assert "Imported 3 station(s)" in feedback
        assert is_open is False

    def test_mesh_unstructured_pcsf_file_imports_successfully(self, tmp_path):
        """A MARE2DEM mesh is a real, supported source too -- this
        exercises the same point-location curtain-slicing path as
        pycsamt/map/tests/test_pcsf_import.py::TestLoadPcsfLinesMeshUnstructured,
        through the actual UI confirm() callback."""
        from pycsamt.app.mapview.app import create_app

        filename, content = _pcsf_file_payload(
            tmp_path, ".pcsf", kind="mesh_unstructured"
        )
        app = create_app()
        confirm = self._confirm(app)
        staged = {"filenames": [filename], "contents": [content]}
        store, feedback, is_open, badge, badge_cls = confirm(
            1, staged, False, "sess-mesh", "replace", "light"
        )
        assert store["n_stations"] == 1
        assert "Imported 1 station(s)" in feedback
        assert is_open is False

    def test_nothing_staged_shows_hint(self):
        from pycsamt.app.mapview.app import create_app

        app = create_app()
        confirm = self._confirm(app)
        store, feedback, *_ = confirm(1, {}, False, "sess-empty", "replace", "light")
        from dash import no_update

        assert store is no_update
        assert "Browse to a folder" in feedback

    def test_no_session_id_shows_hint(self, tmp_path):
        from pycsamt.app.mapview.app import create_app

        filename, content = _pcsf_file_payload(tmp_path, ".pcsf", kind="grid2d")
        app = create_app()
        confirm = self._confirm(app)
        staged = {"filenames": [filename], "contents": [content]}
        store, feedback, *_ = confirm(1, staged, False, None, "replace", "light")
        from dash import no_update

        assert store is no_update
        assert "Session not initialised" in feedback
