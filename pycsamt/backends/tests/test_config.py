# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.backends._config.BackendConfig — the resolution
chain (session override -> env var -> config file -> auto-detect),
the singleton contract, and config-file persistence.

All filesystem access is redirected to ``tmp_path`` via a patched
``Path.home`` so these tests never touch the real
``~/.pycsamt/config.json``, and the process-wide singleton is reset
around every test so state does not leak into other test modules.
"""

from __future__ import annotations

import json
from unittest import mock

import pytest

from pycsamt.backends import _config
from pycsamt.backends._config import BackendConfig


@pytest.fixture(autouse=True)
def reset_singleton():
    saved = BackendConfig._instance
    BackendConfig._instance = None
    yield
    BackendConfig._instance = saved


@pytest.fixture(autouse=True)
def clean_env(monkeypatch):
    monkeypatch.delenv("PYCSAMT_AI_BACKEND", raising=False)


# ─── singleton ──────────────────────────────────────────────────────────────


def test_singleton_identity():
    assert BackendConfig() is BackendConfig()


def test_fresh_instance_has_no_resolved_backend():
    cfg = BackendConfig()
    assert cfg._backend_name is None


# ─── set / backend_name ─────────────────────────────────────────────────────


def test_set_valid_name():
    cfg = BackendConfig()
    cfg.set("torch")
    assert cfg.backend_name == "torch"


def test_set_is_case_insensitive_and_strips():
    cfg = BackendConfig()
    cfg.set("  TensorFlow  ")
    assert cfg.backend_name == "tensorflow"


def test_set_invalid_raises():
    cfg = BackendConfig()
    with pytest.raises(ValueError, match="Unknown backend"):
        cfg.set("keras")


def test_set_auto_triggers_resolve():
    cfg = BackendConfig()
    with mock.patch.object(cfg, "_resolve", return_value="torch"):
        cfg.set("auto")
        assert cfg.backend_name == "torch"


def test_reset_clears_cache():
    cfg = BackendConfig()
    cfg.set("torch")
    cfg.reset()
    with mock.patch.object(cfg, "_resolve", return_value="tensorflow"):
        assert cfg.backend_name == "tensorflow"


def test_backend_name_resolves_lazily_and_caches():
    cfg = BackendConfig()
    with mock.patch.object(cfg, "_resolve", return_value="torch") as m:
        assert cfg.backend_name == "torch"
        assert cfg.backend_name == "torch"
        m.assert_called_once()


def test_backend_name_re_resolves_when_cached_value_is_auto():
    """A cached value of 'auto' must never be returned as-is — it is
    re-resolved on every access until a concrete name is cached."""
    cfg = BackendConfig()
    cfg._backend_name = "auto"
    with mock.patch.object(cfg, "_resolve", return_value="torch") as m:
        assert cfg.backend_name == "torch"
        m.assert_called_once()


# ─── _resolve chain ──────────────────────────────────────────────────────────


def test_resolve_prefers_env_var(monkeypatch):
    monkeypatch.setenv("PYCSAMT_AI_BACKEND", "torch")
    cfg = BackendConfig()
    assert cfg._resolve() == "torch"


def test_resolve_env_var_case_insensitive(monkeypatch):
    monkeypatch.setenv("PYCSAMT_AI_BACKEND", "  TensorFlow  ")
    cfg = BackendConfig()
    assert cfg._resolve() == "tensorflow"


def test_resolve_env_var_auto_falls_through_to_config_file(monkeypatch):
    monkeypatch.setenv("PYCSAMT_AI_BACKEND", "auto")
    cfg = BackendConfig()
    with mock.patch.object(cfg, "_read_config_file", return_value="tensorflow"):
        assert cfg._resolve() == "tensorflow"


def test_resolve_env_var_invalid_falls_through(monkeypatch):
    monkeypatch.setenv("PYCSAMT_AI_BACKEND", "bogus")
    cfg = BackendConfig()
    with mock.patch.object(cfg, "_read_config_file", return_value="tensorflow"):
        assert cfg._resolve() == "tensorflow"


def test_resolve_falls_back_to_config_file():
    cfg = BackendConfig()
    with mock.patch.object(cfg, "_read_config_file", return_value="tensorflow"):
        assert cfg._resolve() == "tensorflow"


def test_resolve_config_file_auto_falls_through_to_autodetect():
    cfg = BackendConfig()
    with (
        mock.patch.object(cfg, "_read_config_file", return_value="auto"),
        mock.patch.object(cfg, "_auto_detect", return_value="torch"),
    ):
        assert cfg._resolve() == "torch"


def test_resolve_config_file_invalid_falls_through_to_autodetect():
    cfg = BackendConfig()
    with (
        mock.patch.object(cfg, "_read_config_file", return_value="bogus"),
        mock.patch.object(cfg, "_auto_detect", return_value="torch"),
    ):
        assert cfg._resolve() == "torch"


def test_resolve_falls_back_to_auto_detect():
    cfg = BackendConfig()
    with (
        mock.patch.object(cfg, "_read_config_file", return_value=None),
        mock.patch.object(cfg, "_auto_detect", return_value="none"),
    ):
        assert cfg._resolve() == "none"


# ─── _read_config_file ────────────────────────────────────────────────────────


def test_read_config_file_missing_home(tmp_path):
    with mock.patch.object(_config.Path, "home", return_value=tmp_path / "nohome"):
        assert BackendConfig._read_config_file() is None


def test_read_config_file_valid(tmp_path):
    cfg_dir = tmp_path / ".pycsamt"
    cfg_dir.mkdir()
    (cfg_dir / "config.json").write_text(json.dumps({"ai_backend": "torch"}))
    with mock.patch.object(_config.Path, "home", return_value=tmp_path):
        assert BackendConfig._read_config_file() == "torch"


def test_read_config_file_bad_json(tmp_path):
    cfg_dir = tmp_path / ".pycsamt"
    cfg_dir.mkdir()
    (cfg_dir / "config.json").write_text("{not valid json")
    with mock.patch.object(_config.Path, "home", return_value=tmp_path):
        assert BackendConfig._read_config_file() is None


def test_read_config_file_missing_key_returns_none(tmp_path):
    cfg_dir = tmp_path / ".pycsamt"
    cfg_dir.mkdir()
    (cfg_dir / "config.json").write_text(json.dumps({"other": 1}))
    with mock.patch.object(_config.Path, "home", return_value=tmp_path):
        assert BackendConfig._read_config_file() is None


# ─── _auto_detect ──────────────────────────────────────────────────────────────


def test_auto_detect_returns_first_available():
    with mock.patch(
        "pycsamt.backends._detect.detect_available_backends",
        return_value=["torch", "tensorflow"],
    ):
        assert BackendConfig._auto_detect() == "torch"


def test_auto_detect_returns_none_when_empty():
    with mock.patch(
        "pycsamt.backends._detect.detect_available_backends",
        return_value=[],
    ):
        assert BackendConfig._auto_detect() == "none"


# ─── write_config_file ─────────────────────────────────────────────────────────


def test_write_config_file_creates_new(tmp_path):
    with mock.patch.object(_config.Path, "home", return_value=tmp_path):
        BackendConfig().write_config_file("torch")
    data = json.loads((tmp_path / ".pycsamt" / "config.json").read_text())
    assert data == {"ai_backend": "torch"}


def test_write_config_file_merges_existing(tmp_path):
    cfg_dir = tmp_path / ".pycsamt"
    cfg_dir.mkdir()
    (cfg_dir / "config.json").write_text(json.dumps({"other_key": 1}))
    with mock.patch.object(_config.Path, "home", return_value=tmp_path):
        BackendConfig().write_config_file("TensorFlow")
    data = json.loads((cfg_dir / "config.json").read_text())
    assert data == {"other_key": 1, "ai_backend": "tensorflow"}


def test_write_config_file_bad_existing_json_is_overwritten(tmp_path):
    cfg_dir = tmp_path / ".pycsamt"
    cfg_dir.mkdir()
    (cfg_dir / "config.json").write_text("{not valid json")
    with mock.patch.object(_config.Path, "home", return_value=tmp_path):
        BackendConfig().write_config_file("torch")
    data = json.loads((cfg_dir / "config.json").read_text())
    assert data == {"ai_backend": "torch"}


def test_write_config_file_invalid_name_raises():
    with pytest.raises(ValueError, match="Cannot persist unknown backend"):
        BackendConfig().write_config_file("keras")
