# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.agents.modem_agent`.

Mocks ``ModEmData``/``ModEmConfig``/``InputBuilder`` (imported locally
inside ``execute()``) with fast, deterministic fakes so every branch of
the ModEM3D data-file pipeline can be exercised without real EDI
geometry: the ``ensure_sites`` exception, the ``models.modem`` import
guard, the list-valued ``component_types`` branch, the ``period_range``
warning branch, the ``ModEmData.from_edi``/``write`` exception, the
"data file not found after write" warning, the ``InputBuilder`` failure
warning, and the LLM interpretation path.
"""

from __future__ import annotations

import sys
import types
from pathlib import Path

import pytest

from pycsamt.agents.modem_agent import ModEmAgent

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


def _passthrough_ensure_sites(monkeypatch, sites=object()):
    import pycsamt.emtools._core as core

    monkeypatch.setattr(core, "ensure_sites", lambda x, **k: sites)


def _mock_modem(
    monkeypatch,
    *,
    data_raises=False,
    write_file=True,
    builder_raises=False,
):
    import pycsamt.models.modem as modem_pkg
    import pycsamt.models.modem.builder as builder_mod
    import pycsamt.models.modem.config as config_mod

    monkeypatch.setattr(config_mod, "ModEmConfig", lambda **k: k)

    class _FakeModEmData:
        def __init__(self, n_sites=3, n_periods=10):
            self.n_sites = n_sites
            self.n_periods = n_periods

        @classmethod
        def from_edi(cls, sites, config=None, **kw):
            if data_raises:
                raise RuntimeError("from_edi boom")
            return cls()

        def write(self, path):
            if write_file:
                Path(path).write_text("fake")

    monkeypatch.setattr(modem_pkg, "ModEmData", _FakeModEmData)

    class _FakeInputBuilder:
        def __init__(self, config=None, **kw):
            pass

        def build_from_data(self, data, workdir, **kw):
            if builder_raises:
                raise RuntimeError("builder boom")
            model_p = Path(workdir) / "ModEM_Model.rho"
            cov_p = Path(workdir) / "ModEM.cov"
            ctrl_p = Path(workdir) / "ModEM.inv"
            for p in (model_p, cov_p, ctrl_p):
                p.write_text("fake")
            return {"model": model_p, "covariance": cov_p, "control": ctrl_p}

    monkeypatch.setattr(builder_mod, "InputBuilder", _FakeInputBuilder)


def test_ensure_sites_exception_fails(monkeypatch):
    import pycsamt.emtools._core as core

    monkeypatch.setattr(
        core,
        "ensure_sites",
        lambda *a, **k: (_ for _ in ()).throw(ValueError("bad sites")),
    )
    agent = ModEmAgent()
    result = agent.execute({"path": "/nope"})
    assert result.status == "failed"
    assert "bad sites" in result.error


def test_models_modem_import_error_fails(monkeypatch, tmp_output):
    _passthrough_ensure_sites(monkeypatch)
    fake_mod = types.ModuleType("pycsamt.models.modem")
    monkeypatch.setitem(sys.modules, "pycsamt.models.modem", fake_mod)
    agent = ModEmAgent()
    result = agent.execute(
        {"sites": object(), "output_dir": str(tmp_output)}
    )
    assert result.status == "failed"
    assert "pycsamt.models.modem not available" in result.error


def test_modem_data_exception_fails(monkeypatch, tmp_output):
    _passthrough_ensure_sites(monkeypatch)
    _mock_modem(monkeypatch, data_raises=True)
    agent = ModEmAgent()
    result = agent.execute(
        {"sites": object(), "output_dir": str(tmp_output)}
    )
    assert result.status == "failed"
    assert "ModEmData.from_edi failed" in result.error


def test_data_file_not_found_after_write_warns(monkeypatch, tmp_output):
    _passthrough_ensure_sites(monkeypatch)
    _mock_modem(monkeypatch, write_file=False)
    agent = ModEmAgent()
    result = agent.execute(
        {"sites": object(), "output_dir": str(tmp_output)}
    )
    assert any(
        "ModEM data file not found after write" in w for w in result.warnings
    )
    assert result.status == "needs_review"


def test_builder_exception_is_recorded_as_warning(monkeypatch, tmp_output):
    _passthrough_ensure_sites(monkeypatch)
    _mock_modem(monkeypatch, builder_raises=True)
    agent = ModEmAgent()
    result = agent.execute(
        {"sites": object(), "output_dir": str(tmp_output)}
    )
    assert result.status == "success"
    assert any(
        "Starting model / covariance / control not written" in w
        for w in result.warnings
    )
    assert result["model_path"] is None
    assert result["cov_path"] is None
    assert result["ctrl_path"] is None


def test_list_component_types_and_period_range_warning(
    monkeypatch, tmp_output
):
    _passthrough_ensure_sites(monkeypatch)
    _mock_modem(monkeypatch)
    agent = ModEmAgent()
    result = agent.execute(
        {
            "sites": object(),
            "output_dir": str(tmp_output),
            "component_types": ["Full_Impedance", "Off_Diagonal_Impedance"],
            "period_range": [0.01, 100.0],
        }
    )
    assert result.status == "success"
    assert any(
        "Period-range filtering should be applied" in w
        for w in result.warnings
    )


def test_happy_path_all_files_and_llm(monkeypatch, tmp_output):
    _passthrough_ensure_sites(monkeypatch)
    _mock_modem(monkeypatch)
    agent = ModEmAgent(api_key="fake-key")
    agent.query_llm = lambda *a, **k: "mocked interpretation"
    result = agent.execute(
        {"sites": object(), "output_dir": str(tmp_output)}
    )
    assert result.status == "success"
    assert result["n_stations"] == 3
    assert result["n_periods"] == 10
    for key in ("data_path", "model_path", "cov_path", "ctrl_path"):
        assert result[key] is not None
        assert Path(result[key]).exists()
    assert result.llm_interpretation == "mocked interpretation"
