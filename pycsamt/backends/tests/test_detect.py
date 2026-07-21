# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.backends._detect — framework-independent, fully mocked
so results do not depend on which ML frameworks happen to be installed in
the environment running the suite.
"""

from __future__ import annotations

import sys
import types
from unittest import mock

import pytest

from pycsamt.backends import _detect


# ─── probe_backend ──────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "alias,resolved_pkg",
    [
        ("torch", "torch"),
        ("pytorch", "torch"),
        ("tensorflow", "tensorflow"),
        ("tf", "tensorflow"),
        ("keras", "keras"),
        ("TORCH", "torch"),  # case-insensitive
        ("TensorFlow", "tensorflow"),
    ],
)
def test_probe_backend_resolves_alias(alias, resolved_pkg):
    with mock.patch.object(
        _detect.importlib.util, "find_spec"
    ) as m_find_spec:
        m_find_spec.return_value = object()  # any non-None spec
        result = _detect.probe_backend(alias)
        m_find_spec.assert_called_once_with(resolved_pkg)
        assert result is True


def test_probe_backend_unknown_name_passthrough():
    """Names without an alias entry are looked up as-is."""
    with mock.patch.object(
        _detect.importlib.util, "find_spec"
    ) as m_find_spec:
        m_find_spec.return_value = None
        result = _detect.probe_backend("some_unrelated_package")
        m_find_spec.assert_called_once_with("some_unrelated_package")
        assert result is False


def test_probe_backend_not_found_returns_false():
    assert _detect.probe_backend("this_package_definitely_does_not_exist_xyz") is False


def test_probe_backend_real_stdlib_module_found():
    # os is always importable — exercises the "found" path with the real
    # importlib.util.find_spec implementation (no mocking).
    assert _detect.probe_backend("os") is True


# ─── detect_available_backends ─────────────────────────────────────────────


@pytest.mark.parametrize(
    "torch_avail,tf_avail,expected",
    [
        (True, True, ["torch", "tensorflow"]),
        (True, False, ["torch"]),
        (False, True, ["tensorflow"]),
        (False, False, []),
    ],
)
def test_detect_available_backends_combinations(torch_avail, tf_avail, expected):
    def fake_probe(name):
        return {"torch": torch_avail, "tensorflow": tf_avail}[name]

    with mock.patch.object(_detect, "probe_backend", side_effect=fake_probe):
        assert _detect.detect_available_backends() == expected


# ─── get_backend_versions ───────────────────────────────────────────────────


def _fake_module(version: str) -> types.ModuleType:
    mod = types.ModuleType("fake")
    mod.__version__ = version
    return mod


def test_get_backend_versions_both_installed():
    def fake_probe(name):
        return True

    with (
        mock.patch.object(_detect, "probe_backend", side_effect=fake_probe),
        mock.patch.dict(
            sys.modules,
            {
                "torch": _fake_module("2.1.0"),
                "tensorflow": _fake_module("2.13.1"),
            },
        ),
    ):
        versions = _detect.get_backend_versions()
        assert versions == {"torch": "2.1.0", "tensorflow": "2.13.1"}


def test_get_backend_versions_neither_installed():
    with mock.patch.object(_detect, "probe_backend", return_value=False):
        versions = _detect.get_backend_versions()
        assert versions == {"torch": None, "tensorflow": None}


def test_get_backend_versions_probe_true_but_import_fails():
    """probe_backend can be True (spec found) while the actual import still
    raises — e.g. a broken native extension. get_backend_versions must
    swallow the exception and leave the version as None."""

    def fake_probe(name):
        return True

    with (
        mock.patch.object(_detect, "probe_backend", side_effect=fake_probe),
        mock.patch.dict(sys.modules, {"torch": None, "tensorflow": None}),
    ):
        # sys.modules[name] = None forces `import name` to raise ImportError
        versions = _detect.get_backend_versions()
        assert versions == {"torch": None, "tensorflow": None}


def test_get_backend_versions_torch_only():
    def fake_probe(name):
        return name == "torch"

    with (
        mock.patch.object(_detect, "probe_backend", side_effect=fake_probe),
        mock.patch.dict(sys.modules, {"torch": _fake_module("1.13.0")}),
    ):
        versions = _detect.get_backend_versions()
        assert versions == {"torch": "1.13.0", "tensorflow": None}
