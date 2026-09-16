# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage tests for pycsamt.inversion.backends.modem.ModEMBackend.

test_backend_compatibility.py / test_inversion_api.py already exercise the
registry-listing and generic no-data-source path shared by every backend.
This file fills the remaining gap: the real ``InputBuilder`` data-source
path (mirroring the site fixture used in
``pycsamt/models/modem/tests/test_modem_phase7.py``), the configured-files
fallback when no data source is supplied, the ``run_external=True``
runner-execution branches, the ``InversionResult`` result-loading branches,
and the small ``_modem_config``/``_build_options``/``_configured_files``/
``_required_runner_keys``/``_missing_runner_files``/``_relative_to_workdir``/
``_resolve_workdir_path``/``_has_loaded_result`` helpers.
"""

from __future__ import annotations

import sys

import numpy as np
import pytest

import pycsamt.inversion.backends.modem as modem_backend_mod
import pycsamt.models.modem as modem_models
from pycsamt.inversion.backends.modem import ModEMBackend
from pycsamt.inversion.config import InversionConfig
from pycsamt.inversion.results import InversionResult
from pycsamt.models.modem.config import ModEmConfig


def _make_site(name, x_offset, y_offset=0.0, n_freq=5):
    freqs = np.logspace(2, -1, n_freq)
    rho = 100.0
    omega = 2 * np.pi * freqs
    mu0 = 4 * np.pi * 1e-7
    z_mag = np.sqrt(omega * mu0 * rho)
    z_val = z_mag * (1.0 + 1.0j) / np.sqrt(2)
    z_arr = np.zeros((n_freq, 2, 2), dtype=complex)
    z_arr[:, 0, 1] = z_val
    z_arr[:, 1, 0] = -z_val

    class _Site:
        pass

    s = _Site()
    s.name = name
    s.coords = (x_offset, y_offset, 0.0)
    s.freq = freqs
    s.z = z_arr
    s.z_err = np.abs(z_arr) * 0.05
    return s


def _make_sites(n=3):
    return [_make_site(f"S{i:02d}", i * 1000.0) for i in range(n)]


# ---------------------------------------------------------------------------
# Data-source build path (real InputBuilder)
# ---------------------------------------------------------------------------


def test_run_with_data_source_builds_and_reports_ready(tmp_path):
    cfg = InversionConfig(
        method="mt",
        dimension="3d",
        backend="modem",
        workdir=str(tmp_path),
        data=_make_sites(),
        backend_options={"nz": 6, "n_airlayers": 2},
    )
    result = ModEMBackend(cfg).run()

    assert isinstance(result, InversionResult)
    assert result.backend == "modem"
    assert result.status in {"ready", "loaded"}
    assert result.files["data"].endswith("data.dat")
    assert result.files["model"].endswith("m0.ws")
    assert result.files["control"].endswith("control.inv")
    assert result.files["covariance"].endswith("covariance.cov")
    assert result.metadata["command"]
    assert result.metadata["mode"] == "3d"
    assert (tmp_path / "data.dat").exists()


def test_run_data_source_build_failure_sets_needs_review(tmp_path):
    cfg = InversionConfig(
        method="mt",
        dimension="3d",
        backend="modem",
        workdir=str(tmp_path),
        data=[],  # empty source -> InputBuilder.build raises
    )
    result = ModEMBackend(cfg).run()
    assert result.status == "needs_review"
    assert any("ModEM preparation failed" in w for w in result.warnings)


def test_run_explicit_data_argument_overrides_config_data(tmp_path):
    cfg = InversionConfig(
        method="mt",
        dimension="3d",
        backend="modem",
        workdir=str(tmp_path),
        data=None,
        backend_options={"nz": 6, "n_airlayers": 2},
    )
    result = ModEMBackend(cfg).run(data=_make_sites())
    assert result.status in {"ready", "loaded"}
    assert (tmp_path / "data.dat").exists()


def test_run_2d_data_source_has_no_covariance_key(tmp_path):
    cfg = InversionConfig(
        method="mt",
        dimension="2d",
        backend="modem",
        workdir=str(tmp_path),
        data=_make_sites(),
    )
    result = ModEMBackend(cfg).run()
    assert "covariance" not in result.files
    assert result.metadata["mode"] == "2d"


# ---------------------------------------------------------------------------
# No-data-source fallback (_configured_files)
# ---------------------------------------------------------------------------


def test_run_without_data_source_configures_files_only(tmp_path):
    cfg = InversionConfig(
        method="mt",
        dimension="3d",
        backend="modem",
        workdir=str(tmp_path),
    )
    result = ModEMBackend(cfg).run()
    assert any(
        "InputBuilder is not available" in w for w in result.warnings
    )
    assert result.files["data"].endswith("ModEMData.dat")
    assert result.files["covariance"].endswith("ModEM.cov")
    assert any("missing files" in w for w in result.warnings)
    assert result.status == "prepared"


def test_run_without_data_source_ready_when_files_present(tmp_path):
    for name in (
        "ModEMData.dat",
        "ModEM_Model.rho",
        "ModEM.inv",
        "ModEM.cov",
    ):
        (tmp_path / name).write_text("", encoding="utf-8")
    cfg = InversionConfig(
        method="mt",
        dimension="3d",
        backend="modem",
        workdir=str(tmp_path),
    )
    result = ModEMBackend(cfg).run()
    assert result.status in {"ready", "loaded"}
    assert result.metadata["command"]


def test_configured_files_uses_alias_options(tmp_path):
    cfg = InversionConfig(
        method="mt",
        dimension="2d",
        backend="modem",
        workdir=str(tmp_path),
        backend_options={
            "data_file": "obs.dat",
            "model_file": "start.ws",
            "control_file": "run.inv",
        },
    )
    result = ModEMBackend(cfg).run()
    assert result.files["data"].endswith("obs.dat")
    assert result.files["model"].endswith("start.ws")
    assert result.files["control"].endswith("run.inv")


# ---------------------------------------------------------------------------
# ImportError branch
# ---------------------------------------------------------------------------


def test_run_raises_import_error_when_modem_package_unavailable(
    monkeypatch, tmp_path
):
    # ``from ...models import modem`` is a single-level fromlist import, so
    # the ``None``-in-sys.modules sentinel alone is not enough: Python's
    # import machinery first checks ``hasattr(pycsamt.models, "modem")``
    # and, since the real submodule is already imported elsewhere in the
    # test session, short-circuits there without consulting sys.modules.
    # The parent-package attribute must be removed too.
    import pycsamt.models as models_pkg

    monkeypatch.delattr(models_pkg, "modem", raising=False)
    monkeypatch.setitem(sys.modules, "pycsamt.models.modem", None)
    cfg = InversionConfig(
        method="mt",
        dimension="3d",
        backend="modem",
        workdir=str(tmp_path),
    )
    with pytest.raises(ImportError, match="ModEM backend requires"):
        ModEMBackend(cfg).run()


# ---------------------------------------------------------------------------
# ModEmRunner not available branch
# ---------------------------------------------------------------------------


def test_run_runner_unavailable_sets_warning(monkeypatch, tmp_path):
    monkeypatch.setattr(modem_models, "ModEmRunner", None)
    cfg = InversionConfig(
        method="mt",
        dimension="3d",
        backend="modem",
        workdir=str(tmp_path),
        data=_make_sites(),
        backend_options={"nz": 6, "n_airlayers": 2},
    )
    result = ModEMBackend(cfg).run()
    assert any(
        "ModEmRunner is not available" in w for w in result.warnings
    )
    assert result.metadata["command"] is None


# ---------------------------------------------------------------------------
# run_external=True branches
# ---------------------------------------------------------------------------


def test_run_external_success_loads_result(monkeypatch, tmp_path):
    for name in (
        "ModEMData.dat",
        "ModEM_Model.rho",
        "ModEM.inv",
        "ModEM.cov",
    ):
        (tmp_path / name).write_text("", encoding="utf-8")

    sentinel = object()
    monkeypatch.setattr(
        modem_models.ModEmRunner, "run", lambda self, *a, **kw: sentinel
    )
    cfg = InversionConfig(
        method="mt",
        dimension="3d",
        backend="modem",
        workdir=str(tmp_path),
        run_external=True,
    )
    result = ModEMBackend(cfg).run()
    assert result.metadata["executed"] is True
    assert result.status == "loaded"
    assert result.native is sentinel


def test_run_external_no_loaded_result_sets_executed(monkeypatch, tmp_path):
    for name in (
        "ModEMData.dat",
        "ModEM_Model.rho",
        "ModEM.inv",
        "ModEM.cov",
    ):
        (tmp_path / name).write_text("", encoding="utf-8")

    monkeypatch.setattr(
        modem_models.ModEmRunner, "run", lambda self, *a, **kw: None
    )
    monkeypatch.setattr(modem_models, "InversionResult", None)
    cfg = InversionConfig(
        method="mt",
        dimension="3d",
        backend="modem",
        workdir=str(tmp_path),
        run_external=True,
    )
    result = ModEMBackend(cfg).run()
    assert result.metadata["executed"] is True
    assert result.status == "executed"


def test_run_external_runner_exception_sets_needs_review(monkeypatch, tmp_path):
    for name in (
        "ModEMData.dat",
        "ModEM_Model.rho",
        "ModEM.inv",
        "ModEM.cov",
    ):
        (tmp_path / name).write_text("", encoding="utf-8")

    def _boom(self, *a, **kw):
        raise RuntimeError("binary crashed")

    monkeypatch.setattr(modem_models.ModEmRunner, "run", _boom)
    cfg = InversionConfig(
        method="mt",
        dimension="3d",
        backend="modem",
        workdir=str(tmp_path),
        run_external=True,
    )
    result = ModEMBackend(cfg).run()
    assert result.status == "needs_review"
    assert any("ModEM runner failed" in w for w in result.warnings)


def test_run_external_passes_runner_options(monkeypatch, tmp_path):
    for name in (
        "ModEMData.dat",
        "ModEM_Model.rho",
        "ModEM.inv",
        "ModEM.cov",
    ):
        (tmp_path / name).write_text("", encoding="utf-8")

    captured = {}

    def _fake_run(self, model, data, control, covariance=None, **kw):
        captured.update(kw)
        return None

    monkeypatch.setattr(modem_models.ModEmRunner, "run", _fake_run)
    monkeypatch.setattr(modem_models, "InversionResult", None)
    cfg = InversionConfig(
        method="mt",
        dimension="3d",
        backend="modem",
        workdir=str(tmp_path),
        run_external=True,
        backend_options={
            "runner": {
                "use_mpi": True,
                "n_procs": 8,
                "extra_args": ["--verbose"],
                "timeout": 120,
                "load_result": False,
            }
        },
    )
    result = ModEMBackend(cfg).run()
    assert captured["use_mpi"] is True
    assert captured["n_procs"] == 8
    assert captured["extra_args"] == ["--verbose"]
    assert captured["timeout"] == 120
    assert captured["load_result"] is False
    # metadata reflects the static ModEmConfig, not the per-call runner
    # override, so it stays at the dataclass default here.
    assert result.metadata["use_mpi"] is False
    assert result.metadata["n_procs"] == 4


# ---------------------------------------------------------------------------
# InversionResult loading branches
# ---------------------------------------------------------------------------


def test_run_result_loaded_when_models_present(monkeypatch, tmp_path):
    class FakeLoadedResult:
        def __init__(self, workdir, config=None, **kw):
            self.models = [object()]
            self.final_rms = 1.234

    monkeypatch.setattr(modem_models, "InversionResult", FakeLoadedResult)
    cfg = InversionConfig(
        method="mt",
        dimension="3d",
        backend="modem",
        workdir=str(tmp_path),
        data=_make_sites(),
        backend_options={"nz": 6, "n_airlayers": 2},
    )
    result = ModEMBackend(cfg).run()
    assert result.status == "loaded"
    assert result.rms == pytest.approx(1.234)


def test_run_result_loading_exception_appends_warning(monkeypatch, tmp_path):
    class FakeBrokenResult:
        def __init__(self, workdir, config=None, **kw):
            raise RuntimeError("cannot parse log")

    monkeypatch.setattr(modem_models, "InversionResult", FakeBrokenResult)
    cfg = InversionConfig(
        method="mt",
        dimension="3d",
        backend="modem",
        workdir=str(tmp_path),
    )
    result = ModEMBackend(cfg).run()
    assert any(
        "ModEM result loading skipped" in w for w in result.warnings
    )
    assert np.isnan(result.rms)


def test_run_result_not_skipped_when_status_needs_review(monkeypatch, tmp_path):
    calls = []

    class TrackingResult:
        def __init__(self, workdir, config=None, **kw):
            calls.append(1)
            raise RuntimeError("should not be called")

    monkeypatch.setattr(modem_models, "InversionResult", TrackingResult)
    cfg = InversionConfig(
        method="mt",
        dimension="3d",
        backend="modem",
        workdir=str(tmp_path),
        data=[],  # forces needs_review before reaching result-loading block
    )
    result = ModEMBackend(cfg).run()
    assert result.status == "needs_review"
    assert calls == []


# ---------------------------------------------------------------------------
# _modem_config
# ---------------------------------------------------------------------------


def test_modem_config_from_backend_options_fields(tmp_path):
    cfg = InversionConfig(
        method="mt",
        dimension="3d",
        backend="modem",
        workdir=str(tmp_path),
        backend_options={"nz": 12, "target_rms": 1.2},
    )
    modem_cfg = modem_backend_mod._modem_config(modem_models, cfg)
    assert modem_cfg.nz == 12
    assert modem_cfg.target_rms == 1.2
    assert modem_cfg.mode == "3d"


def test_modem_config_dict_branch():
    cfg = InversionConfig(
        method="mt",
        dimension="2d",
        backend="modem",
        backend_options={"config": {"binary_2d": "MyMod2DMT"}},
    )
    modem_cfg = modem_backend_mod._modem_config(modem_models, cfg)
    assert isinstance(modem_cfg, ModEmConfig)
    assert modem_cfg.binary_2d == "MyMod2DMT"
    assert modem_cfg.mode == "2d"


def test_modem_config_instance_branch():
    provided = ModEmConfig(mode="2d", binary_3d="AlreadyBuilt")
    cfg = InversionConfig(
        method="mt",
        dimension="3d",
        backend="modem",
        backend_options={"config": provided},
    )
    modem_cfg = modem_backend_mod._modem_config(modem_models, cfg)
    assert modem_cfg is provided
    assert modem_cfg.mode == "3d"  # overwritten to match cfg.dimension


def test_modem_config_raw_passthrough_with_mode_attr():
    class _Fake:
        mode = "old"

    fake = _Fake()
    cfg = InversionConfig(
        method="mt",
        dimension="3d",
        backend="modem",
        backend_options={"config": fake},
    )
    out = modem_backend_mod._modem_config(modem_models, cfg)
    assert out is fake
    assert out.mode == "3d"


def test_modem_config_raw_passthrough_without_mode_attr():
    sentinel = object()
    cfg = InversionConfig(
        method="mt",
        dimension="3d",
        backend="modem",
        backend_options={"config": sentinel},
    )
    assert modem_backend_mod._modem_config(modem_models, cfg) is sentinel


# ---------------------------------------------------------------------------
# _build_options
# ---------------------------------------------------------------------------


def test_build_options_excludes_reserved_and_modem_field_keys():
    options = {
        "config": {"a": 1},
        "runner": {"b": 2},
        "files": {"c": 3},
        "mode": "3d",
        "binary_2d": "x",
        "binary_3d": "y",
        "use_mpi": True,
        "n_procs": 4,
        "mpi_command": "mpirun",
        "data_filename": "obs.dat",
        "nz": 30,
    }
    out = modem_backend_mod._build_options(options)
    assert out == {"data_filename": "obs.dat", "nz": 30}


# ---------------------------------------------------------------------------
# _configured_files
# ---------------------------------------------------------------------------


def test_configured_files_defaults_2d_no_covariance(tmp_path):
    modem_cfg = ModEmConfig(mode="2d")
    files = modem_backend_mod._configured_files(tmp_path, modem_cfg, {})
    assert "covariance" not in files
    assert files["data"] == str(tmp_path / modem_cfg.data_file)


def test_configured_files_defaults_3d_has_covariance(tmp_path):
    modem_cfg = ModEmConfig(mode="3d")
    files = modem_backend_mod._configured_files(tmp_path, modem_cfg, {})
    assert files["covariance"] == str(tmp_path / modem_cfg.covariance_file)


def test_configured_files_explicit_files_mapping_wins(tmp_path):
    modem_cfg = ModEmConfig(mode="3d")
    files = modem_backend_mod._configured_files(
        tmp_path, modem_cfg, {"files": {"data": "custom.dat"}}
    )
    assert files["data"] == str(tmp_path / "custom.dat")


# ---------------------------------------------------------------------------
# _required_runner_keys / _missing_runner_files
# ---------------------------------------------------------------------------


def test_required_runner_keys_2d_vs_3d():
    assert modem_backend_mod._required_runner_keys(
        ModEmConfig(mode="2d")
    ) == ("model", "data", "control")
    assert modem_backend_mod._required_runner_keys(
        ModEmConfig(mode="3d")
    ) == ("model", "data", "control", "covariance")


def test_missing_runner_files_key_absent_entirely():
    missing = modem_backend_mod._missing_runner_files(
        {"data": None}, required=("data", "model", "control")
    )
    assert "data" in missing
    assert "model" in missing
    assert "control" in missing


def test_missing_runner_files_nonexistent_path(tmp_path):
    missing = modem_backend_mod._missing_runner_files(
        {"data": str(tmp_path / "nope.dat")}, required=("data",)
    )
    assert len(missing) == 1
    assert "data=" in missing[0]


def test_missing_runner_files_none_when_all_present(tmp_path):
    p = tmp_path / "data.dat"
    p.write_text("", encoding="utf-8")
    missing = modem_backend_mod._missing_runner_files(
        {"data": str(p)}, required=("data",)
    )
    assert missing == []


# ---------------------------------------------------------------------------
# _relative_to_workdir / _resolve_workdir_path
# ---------------------------------------------------------------------------


def test_relative_to_workdir_inside(tmp_path):
    p = tmp_path / "sub" / "data.dat"
    rel = modem_backend_mod._relative_to_workdir(str(p), tmp_path)
    assert rel == str(p.relative_to(tmp_path))


def test_relative_to_workdir_outside_returns_original(tmp_path, tmp_path_factory):
    other = tmp_path_factory.mktemp("elsewhere") / "data.dat"
    rel = modem_backend_mod._relative_to_workdir(str(other), tmp_path)
    assert rel == str(other)


def test_resolve_workdir_path_relative_and_absolute(tmp_path):
    rel_result = modem_backend_mod._resolve_workdir_path("x.dat", tmp_path)
    assert rel_result == tmp_path / "x.dat"

    abs_path = tmp_path / "abs.dat"
    abs_result = modem_backend_mod._resolve_workdir_path(abs_path, tmp_path)
    assert abs_result == abs_path


# ---------------------------------------------------------------------------
# _has_loaded_result
# ---------------------------------------------------------------------------


def test_has_loaded_result_variants():
    class _Empty:
        pass

    assert modem_backend_mod._has_loaded_result(_Empty()) is False

    class _WithModels:
        models = [1]

    assert modem_backend_mod._has_loaded_result(_WithModels()) is True

    class _WithLog:
        log = "something"

    assert modem_backend_mod._has_loaded_result(_WithLog()) is True

    class _WithData:
        data = {"a": 1}

    assert modem_backend_mod._has_loaded_result(_WithData()) is True

    class _WithRho2D:
        rho_2d = np.zeros((2, 2))

    assert modem_backend_mod._has_loaded_result(_WithRho2D()) is True
