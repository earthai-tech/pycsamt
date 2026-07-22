# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for the pycsamt.backends public API (__init__.py): get_backend,
set_backend, auto_detect, list_backends, get_backend_instance.

The module-level ``_CFG`` singleton is swapped for a mock in every test
so these tests exercise only the __init__.py dispatch logic, are
independent of whichever ML frameworks happen to be installed, and never
touch the real ``~/.pycsamt/config.json``.
"""

from __future__ import annotations

from unittest import mock

import pytest

import pycsamt.backends as backends_pkg
from pycsamt.backends import (
    auto_detect,
    get_backend,
    get_backend_instance,
    list_backends,
    set_backend,
)


@pytest.fixture
def fake_cfg():
    original = backends_pkg._CFG
    fake = mock.MagicMock()
    backends_pkg._CFG = fake
    yield fake
    backends_pkg._CFG = original


# ─── get_backend ────────────────────────────────────────────────────────────


def test_get_backend_returns_cfg_backend_name(fake_cfg):
    fake_cfg.backend_name = "torch"
    assert get_backend() == "torch"


# ─── set_backend ────────────────────────────────────────────────────────────


def test_set_backend_delegates_to_cfg(fake_cfg):
    set_backend("torch")
    fake_cfg.set.assert_called_once_with("torch")
    fake_cfg.write_config_file.assert_not_called()


def test_set_backend_persist_writes_config(fake_cfg):
    fake_cfg.backend_name = "tensorflow"
    set_backend("tensorflow", persist=True)
    fake_cfg.set.assert_called_once_with("tensorflow")
    fake_cfg.write_config_file.assert_called_once_with("tensorflow")


def test_set_backend_no_persist_by_default(fake_cfg):
    fake_cfg.backend_name = "torch"
    set_backend("torch", persist=False)
    fake_cfg.write_config_file.assert_not_called()


# ─── auto_detect ────────────────────────────────────────────────────────────


def test_auto_detect_resets_then_returns_name(fake_cfg):
    fake_cfg.backend_name = "torch"
    assert auto_detect() == "torch"
    fake_cfg.reset.assert_called_once()


# ─── list_backends ──────────────────────────────────────────────────────────


def test_list_backends_uses_probe_backend():
    with mock.patch.object(
        backends_pkg, "probe_backend", side_effect=lambda n: n == "torch"
    ):
        assert list_backends() == {"torch": True, "tensorflow": False}


def test_list_backends_both_available():
    with mock.patch.object(backends_pkg, "probe_backend", return_value=True):
        assert list_backends() == {"torch": True, "tensorflow": True}


def test_list_backends_none_available():
    with mock.patch.object(backends_pkg, "probe_backend", return_value=False):
        assert list_backends() == {"torch": False, "tensorflow": False}


# ─── get_backend_instance ───────────────────────────────────────────────────


def test_get_backend_instance_torch(fake_cfg):
    fake_cfg.backend_name = "torch"
    from pycsamt.backends._torch import TorchBackend

    assert isinstance(get_backend_instance(), TorchBackend)


def test_get_backend_instance_tensorflow(fake_cfg):
    """Constructing a TensorFlowBackend does not itself import tensorflow
    (all framework imports are deferred to method bodies), so this must
    work even when TensorFlow is not installed / not importable."""
    fake_cfg.backend_name = "tensorflow"
    from pycsamt.backends._tensorflow import TensorFlowBackend

    assert isinstance(get_backend_instance(), TensorFlowBackend)


def test_get_backend_instance_none(fake_cfg):
    fake_cfg.backend_name = "none"
    assert get_backend_instance() is None


def test_get_backend_instance_unknown_raises(fake_cfg):
    fake_cfg.backend_name = "bogus"
    with pytest.raises(RuntimeError, match="Unknown backend"):
        get_backend_instance()
