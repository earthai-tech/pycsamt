# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.agents.occam2d_agent`.

Mocks ``OccamData``/``OccamConfig``/``OccamMesh``/``OccamModel``/
``OccamStartup`` (imported locally inside ``execute()``) with fast,
deterministic fakes so the whole four-file pipeline runs without needing
real EDI geometry: the ``models.occam2d`` import guard, ``ensure_sites``
exception, the ``OccamData.from_edi`` exception, the mesh/model/startup
success and exception branches, the "file not found after write" warning,
and the LLM interpretation path.
"""

from __future__ import annotations

import sys
import types
from pathlib import Path

import pytest

from pycsamt.agents.occam2d_agent import Occam2DAgent

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


class _FakeWritable:
    def __init__(self, *, write_file=True, **attrs):
        self.__dict__.update(attrs)
        self._write_file = write_file

    def write(self, path):
        if self._write_file:
            Path(path).write_text("fake")


class _FakeOccamModel(_FakeWritable):
    def __init__(self, *, description=None, n_params=50, **attrs):
        super().__init__(n_params=n_params, description=description, **attrs)

    @classmethod
    def from_mesh(cls, mesh, config=None):
        return cls(n_params=50)


def _mock_occam2d(
    monkeypatch,
    *,
    data_raises=False,
    mesh_raises=False,
    model_raises=False,
    startup_raises=False,
    mesh_writes_file=True,
):
    import pycsamt.models.occam2d as occam_pkg
    import pycsamt.models.occam2d.config as config_mod
    import pycsamt.models.occam2d.mesh as mesh_mod
    import pycsamt.models.occam2d.model as model_mod
    import pycsamt.models.occam2d.startup as startup_mod

    monkeypatch.setattr(config_mod, "OccamConfig", lambda **k: k)

    class _FakeOccamData(_FakeWritable):
        @classmethod
        def from_edi(cls, sites, *, modes, config, title):
            if data_raises:
                raise RuntimeError("from_edi boom")
            return cls(n_sites=3, n_frequencies=10, n_data=120)

    monkeypatch.setattr(occam_pkg, "OccamData", _FakeOccamData)

    class _FakeOccamMesh(_FakeWritable):
        @classmethod
        def from_data(cls, occ_data):
            if mesh_raises:
                raise RuntimeError("mesh boom")
            return cls(
                n_xcells=20, n_zcells=15, write_file=mesh_writes_file
            )

    monkeypatch.setattr(mesh_mod, "OccamMesh", _FakeOccamMesh)

    class _MaybeBoomOccamModel(_FakeOccamModel):
        @classmethod
        def from_mesh(cls, mesh, config=None):
            if model_raises:
                raise RuntimeError("model boom")
            return cls(n_params=50)

    monkeypatch.setattr(model_mod, "OccamModel", _MaybeBoomOccamModel)

    class _FakeOccamStartup(_FakeWritable):
        @classmethod
        def from_model(cls, occ_model, config=None):
            if startup_raises:
                raise RuntimeError("startup boom")
            return cls()

    monkeypatch.setattr(startup_mod, "OccamStartup", _FakeOccamStartup)


def _passthrough_ensure_sites(monkeypatch, sites=object()):
    import pycsamt.emtools._core as core

    monkeypatch.setattr(core, "ensure_sites", lambda x, **k: sites)


def test_ensure_sites_exception_fails(monkeypatch):
    import pycsamt.emtools._core as core

    monkeypatch.setattr(
        core,
        "ensure_sites",
        lambda *a, **k: (_ for _ in ()).throw(ValueError("bad sites")),
    )
    agent = Occam2DAgent()
    result = agent.execute({"path": "/nope"})
    assert result.status == "failed"
    assert "bad sites" in result.error


def test_models_occam2d_import_error_fails(monkeypatch, tmp_output):
    _passthrough_ensure_sites(monkeypatch)
    fake_mod = types.ModuleType("pycsamt.models.occam2d")
    monkeypatch.setitem(sys.modules, "pycsamt.models.occam2d", fake_mod)
    agent = Occam2DAgent()
    result = agent.execute(
        {"sites": object(), "output_dir": str(tmp_output)}
    )
    assert result.status == "failed"
    assert "pycsamt.models.occam2d not available" in result.error


def test_occam_data_exception_fails(monkeypatch, tmp_output):
    _passthrough_ensure_sites(monkeypatch)
    _mock_occam2d(monkeypatch, data_raises=True)
    agent = Occam2DAgent()
    result = agent.execute(
        {"sites": object(), "output_dir": str(tmp_output)}
    )
    assert result.status == "failed"
    assert "OccamData.from_edi failed" in result.error


def test_mesh_exception_is_recorded_as_warning(monkeypatch, tmp_output):
    _passthrough_ensure_sites(monkeypatch)
    _mock_occam2d(monkeypatch, mesh_raises=True)
    agent = Occam2DAgent()
    result = agent.execute(
        {"sites": object(), "output_dir": str(tmp_output)}
    )
    assert result.status == "success"
    assert any("OccamMesh.from_data failed" in w for w in result.warnings)
    assert result["mesh_path"] is None


def test_model_exception_falls_back_to_default_description(
    monkeypatch, tmp_output
):
    _passthrough_ensure_sites(monkeypatch)
    _mock_occam2d(monkeypatch, model_raises=True)
    agent = Occam2DAgent()
    result = agent.execute(
        {"sites": object(), "output_dir": str(tmp_output)}
    )
    assert any("OccamModel.from_mesh failed" in w for w in result.warnings)
    # startup is still built from the default-description fallback model
    assert result["startup_path"] is not None


def test_startup_exception_is_recorded_as_warning(monkeypatch, tmp_output):
    _passthrough_ensure_sites(monkeypatch)
    _mock_occam2d(monkeypatch, startup_raises=True)
    agent = Occam2DAgent()
    result = agent.execute(
        {"sites": object(), "output_dir": str(tmp_output)}
    )
    assert any("OccamStartup.from_model failed" in w for w in result.warnings)
    assert result["startup_path"] is None


def test_file_not_found_after_write_warns(monkeypatch, tmp_output):
    _passthrough_ensure_sites(monkeypatch)
    _mock_occam2d(monkeypatch, mesh_writes_file=False)
    agent = Occam2DAgent()
    result = agent.execute(
        {"sites": object(), "output_dir": str(tmp_output)}
    )
    assert any("mesh file not found after write" in w for w in result.warnings)


def test_happy_path_all_four_files_and_llm(monkeypatch, tmp_output):
    _passthrough_ensure_sites(monkeypatch)
    _mock_occam2d(monkeypatch)
    agent = Occam2DAgent(api_key="fake-key")
    agent.query_llm = lambda *a, **k: "mocked interpretation"
    result = agent.execute(
        {"sites": object(), "output_dir": str(tmp_output)}
    )
    assert result.status == "success"
    assert result["n_stations"] == 3
    assert result["n_periods"] == 10
    assert result["n_data"] == 120
    for key in ("data_path", "mesh_path", "model_path", "startup_path"):
        assert result[key] is not None
        assert Path(result[key]).exists()
    assert result.llm_interpretation == "mocked interpretation"
