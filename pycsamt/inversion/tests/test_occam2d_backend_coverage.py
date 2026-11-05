# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage tests for pycsamt.inversion.backends.occam2d.Occam2DBackend.

test_inversion_api.py already covers the "no data source" dry-run path
(command assembly, missing-runner-files reporting) and the small
``_missing_runner_files``/``_command_string`` helpers directly. This file
fills the remaining gap: the data-source build path through the real,
well-tested ``pycsamt.models.occam2d.InputBuilder`` (mirroring the synthetic
EDI-style site fixture used in
``pycsamt/models/occam2d/tests/test_occam2d_builder.py``), the
``run_external=True`` runner-execution branches, the ``OccamResult``
result-loading branches, and the ``_occam_config``/``_build_options``/
``_builder_files`` helpers.
"""

from __future__ import annotations

import sys

import numpy as np
import pytest

import pycsamt.inversion.backends.occam2d as occam2d_mod
import pycsamt.models.occam2d as occam2d_models
from pycsamt.inversion.backends.occam2d import Occam2DBackend
from pycsamt.inversion.config import InversionConfig
from pycsamt.inversion.results import InversionResult
from pycsamt.models.occam2d import InputBuilder, OccamConfig, OccamRunner


class _ZContainer:
    freq = np.array([10.0, 1.0])
    resistivity = np.ones((2, 2, 2), dtype=float) * 100.0
    resistivity_err = np.ones((2, 2, 2), dtype=float) * 5.0
    phase = np.ones((2, 2, 2), dtype=float) * 45.0
    phase_err = np.ones((2, 2, 2), dtype=float)


class _EDIStyleSite:
    def __init__(self, name, lon):
        self.station = name
        self.coords = (0.0, lon, 0.0)
        self.Z = _ZContainer()


def _make_sites(n=2):
    return [_EDIStyleSite(f"S{i:03d}", i * 0.01) for i in range(n)]


# ---------------------------------------------------------------------------
# Data-source build path (real InputBuilder)
# ---------------------------------------------------------------------------


def test_run_with_data_source_builds_and_reports_ready(tmp_path):
    cfg = InversionConfig(
        method="mt",
        dimension="2d",
        backend="occam2d",
        workdir=str(tmp_path),
        data=_make_sites(2),
        backend_options={
            "modes": ["TE", "TM"],
            "freq_min": 1.0,
            "freq_max": 10.0,
        },
    )
    result = Occam2DBackend(cfg).run()

    assert isinstance(result, InversionResult)
    assert result.backend == "occam2d"
    assert result.status in {"ready", "loaded"}
    assert result.files["data"].endswith("OccamDataFile.dat")
    assert result.files["mesh"].endswith("Occam2DMesh")
    assert result.files["model"].endswith("Occam2DModel")
    assert result.files["startup"].endswith("Startup")
    assert result.metadata["command"] == "Occam2D Startup"
    assert (tmp_path / "OccamDataFile.dat").exists()


def test_run_data_source_build_failure_sets_needs_review(tmp_path):
    cfg = InversionConfig(
        method="mt",
        dimension="2d",
        backend="occam2d",
        workdir=str(tmp_path),
        data=[],  # empty source -> OccamData.from_edi raises
    )
    result = Occam2DBackend(cfg).run()
    assert result.status == "needs_review"
    assert any("Occam2D preparation failed" in w for w in result.warnings)


def test_run_explicit_data_argument_overrides_config_data(tmp_path):
    cfg = InversionConfig(
        method="mt",
        dimension="2d",
        backend="occam2d",
        workdir=str(tmp_path),
        data=None,
    )
    result = Occam2DBackend(cfg).run(data=_make_sites(2))
    assert result.status in {"ready", "loaded"}
    assert (tmp_path / "OccamDataFile.dat").exists()


# ---------------------------------------------------------------------------
# ImportError branch
# ---------------------------------------------------------------------------


def test_run_raises_import_error_when_occam2d_package_unavailable(
    monkeypatch, tmp_path
):
    monkeypatch.setitem(sys.modules, "pycsamt.models.occam2d", None)
    cfg = InversionConfig(
        method="mt",
        dimension="2d",
        backend="occam2d",
        workdir=str(tmp_path),
    )
    with pytest.raises(ImportError, match="Occam2D backend requires"):
        Occam2DBackend(cfg).run()


# ---------------------------------------------------------------------------
# run_external=True branches
# ---------------------------------------------------------------------------


def test_run_external_success_sets_executed_status(monkeypatch, tmp_path):
    for name in (
        "OccamDataFile.dat",
        "Occam2DMesh",
        "Occam2DModel",
        "Startup",
    ):
        (tmp_path / name).write_text("", encoding="utf-8")

    monkeypatch.setattr(OccamRunner, "run", lambda self, **kw: 0)
    cfg = InversionConfig(
        method="mt",
        dimension="2d",
        backend="occam2d",
        workdir=str(tmp_path),
        run_external=True,
        backend_options={"binary_name": "Occam2D"},
    )
    result = Occam2DBackend(cfg).run()
    assert result.metadata["executed"] is True
    assert result.metadata["exit_code"] == 0
    assert result.status in {"executed", "loaded"}


@pytest.mark.parametrize(
    ("runner_opts", "expected"),
    [({}, 25), ({"runner": {"max_iter": 7}}, 7)],
)
def test_run_external_passes_iteration_budget(
    monkeypatch, tmp_path, runner_opts, expected
):
    # InversionConfig.max_iter is the caller's budget (what the agent sets);
    # it used to be dropped, leaving the startup file's default of 100.
    for name in ("OccamDataFile.dat", "Occam2DMesh", "Occam2DModel", "Startup"):
        (tmp_path / name).write_text("", encoding="utf-8")
    seen = {}

    def _run(self, **kw):
        seen.update(kw)
        return 0

    monkeypatch.setattr(OccamRunner, "run", _run)
    cfg = InversionConfig(
        method="mt",
        dimension="2d",
        backend="occam2d",
        workdir=str(tmp_path),
        run_external=True,
        max_iter=25,
        backend_options={"binary_name": "Occam2D", **runner_opts},
    )
    Occam2DBackend(cfg).run()
    assert seen["max_iter"] == expected


def test_run_external_nonzero_exit_sets_needs_review(monkeypatch, tmp_path):
    for name in (
        "OccamDataFile.dat",
        "Occam2DMesh",
        "Occam2DModel",
        "Startup",
    ):
        (tmp_path / name).write_text("", encoding="utf-8")

    monkeypatch.setattr(OccamRunner, "run", lambda self, **kw: 2)
    cfg = InversionConfig(
        method="mt",
        dimension="2d",
        backend="occam2d",
        workdir=str(tmp_path),
        run_external=True,
    )
    result = Occam2DBackend(cfg).run()
    assert result.metadata["executed"] is True
    assert result.metadata["exit_code"] == 2
    assert result.status == "needs_review"
    assert any("exited with code 2" in w for w in result.warnings)


def test_run_external_runner_exception_sets_needs_review(monkeypatch, tmp_path):
    for name in (
        "OccamDataFile.dat",
        "Occam2DMesh",
        "Occam2DModel",
        "Startup",
    ):
        (tmp_path / name).write_text("", encoding="utf-8")

    def _boom(self, **kw):
        raise RuntimeError("binary crashed")

    monkeypatch.setattr(OccamRunner, "run", _boom)
    cfg = InversionConfig(
        method="mt",
        dimension="2d",
        backend="occam2d",
        workdir=str(tmp_path),
        run_external=True,
    )
    result = Occam2DBackend(cfg).run()
    assert result.metadata["executed"] is False
    assert result.status == "needs_review"
    assert any("Occam2D runner failed" in w for w in result.warnings)


# ---------------------------------------------------------------------------
# OccamResult loading branches
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(("reached", "expected"), [(True, True), (False, False)])
def test_loaded_result_reports_iteration_and_native_convergence(
    monkeypatch, tmp_path, reached, expected
):
    # A run that stalls above the target is loaded, but not converged.
    from types import SimpleNamespace

    class FakeLoadedResult:
        def __init__(self, workdir, config=None, **kw):
            self.rho_2d = np.ones((2, 2))
            self.final_rms = 4.3
            self.best_iter = SimpleNamespace(
                iteration=12, misfit_reached=reached, target_misfit=1.0
            )

    monkeypatch.setattr(occam2d_models, "InversionResult", FakeLoadedResult)
    cfg = InversionConfig(
        method="mt", dimension="2d", backend="occam2d",
        workdir=str(tmp_path), data=_make_sites(2),
    )
    result = Occam2DBackend(cfg).run()
    assert result.n_iter == 12
    assert result.metadata["converged"] is expected
    assert result.metadata["target_misfit"] == 1.0


def test_run_result_loaded_when_rho_2d_present(monkeypatch, tmp_path):
    class FakeLoadedResult:
        def __init__(self, workdir, config=None, **kw):
            self.rho_2d = np.ones((2, 2))
            self.final_rms = 1.234

    monkeypatch.setattr(occam2d_models, "InversionResult", FakeLoadedResult)
    cfg = InversionConfig(
        method="mt",
        dimension="2d",
        backend="occam2d",
        workdir=str(tmp_path),
        data=_make_sites(2),
    )
    result = Occam2DBackend(cfg).run()
    assert result.status == "loaded"
    assert result.rms == pytest.approx(1.234)


def test_run_result_loading_exception_appends_warning(monkeypatch, tmp_path):
    class FakeBrokenResult:
        def __init__(self, workdir, config=None, **kw):
            raise RuntimeError("cannot parse iteration files")

    monkeypatch.setattr(occam2d_models, "InversionResult", FakeBrokenResult)
    cfg = InversionConfig(
        method="mt",
        dimension="2d",
        backend="occam2d",
        workdir=str(tmp_path),
    )
    result = Occam2DBackend(cfg).run()
    assert any("Occam2D result loading skipped" in w for w in result.warnings)
    assert np.isnan(result.rms)


# ---------------------------------------------------------------------------
# _occam_config
# ---------------------------------------------------------------------------


def test_occam_config_from_backend_options_fields(tmp_path):
    cfg = InversionConfig(
        method="mt",
        dimension="2d",
        backend="occam2d",
        workdir=str(tmp_path),
        backend_options={"binary_name": "MyOccam2D", "n_layers": 40},
    )
    occam_cfg = occam2d_mod._occam_config(OccamConfig, cfg)
    assert occam_cfg.binary_name == "MyOccam2D"
    assert occam_cfg.n_layers == 40


def test_occam_config_dict_branch():
    cfg = InversionConfig(
        method="mt",
        dimension="2d",
        backend="occam2d",
        backend_options={"config": {"binary_name": "FromDict"}},
    )
    occam_cfg = occam2d_mod._occam_config(OccamConfig, cfg)
    assert isinstance(occam_cfg, OccamConfig)
    assert occam_cfg.binary_name == "FromDict"


def test_occam_config_instance_branch():
    provided = OccamConfig(binary_name="AlreadyBuilt")
    cfg = InversionConfig(
        method="mt",
        dimension="2d",
        backend="occam2d",
        backend_options={"config": provided},
    )
    occam_cfg = occam2d_mod._occam_config(OccamConfig, cfg)
    assert occam_cfg is provided


def test_occam_config_raw_passthrough_branch():
    sentinel = object()
    cfg = InversionConfig(
        method="mt",
        dimension="2d",
        backend="occam2d",
        backend_options={"config": sentinel},
    )
    assert occam2d_mod._occam_config(OccamConfig, cfg) is sentinel


# ---------------------------------------------------------------------------
# _build_options
# ---------------------------------------------------------------------------


def test_build_options_excludes_reserved_and_config_field_keys():
    options = {
        "config": {"a": 1},
        "runner": {"b": 2},
        "files": {"c": 3},
        "binary_path": "x",
        "data_file": "y",
        "mesh_file": "z",
        "model_file": "m",
        "startup_file": "s",
        "binary_name": "n",
        "modes": ["TE"],
        "n_layers": 30,
    }
    out = occam2d_mod._build_options(options)
    assert out == {"modes": ["TE"], "n_layers": 30}


# ---------------------------------------------------------------------------
# _builder_files
# ---------------------------------------------------------------------------


def test_builder_files_uses_object_path_when_present(tmp_path):
    builder = InputBuilder(_make_sites(2), workdir=tmp_path).build()
    files = occam2d_mod._builder_files(builder)
    assert files["data"].endswith("OccamDataFile.dat")
    assert files["startup"].endswith("Startup")


def test_builder_files_falls_back_to_workdir_defaults_when_no_path(tmp_path):
    class _NoPath:
        pass

    builder = InputBuilder(_make_sites(2), workdir=tmp_path)
    builder.config = OccamConfig()
    builder.data = _NoPath()
    builder.mesh = _NoPath()
    builder.model = _NoPath()
    builder.startup = _NoPath()
    files = occam2d_mod._builder_files(builder)
    assert files["data"] == str(tmp_path / "OccamDataFile.dat")


# ---------------------------------------------------------------------------
# _missing_runner_files -- absent-key branch
# ---------------------------------------------------------------------------


def test_missing_runner_files_key_absent_entirely():
    missing = occam2d_mod._missing_runner_files({"data": None})
    assert "mesh" in missing
    assert "model" in missing
    assert "startup" in missing
    assert "data" in missing
