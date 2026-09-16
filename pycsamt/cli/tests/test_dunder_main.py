# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for ``python -m pycsamt.cli`` (pycsamt/cli/__main__.py)."""

from __future__ import annotations

import runpy
import sys

import pytest


def test_dunder_main_help_exits_zero(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["pycsamt", "--help"])
    with pytest.raises(SystemExit) as exc:
        runpy.run_module("pycsamt.cli.__main__", run_name="__main__")
    assert exc.value.code == 0
    out = capsys.readouterr().out
    assert "pycsamt" in out.lower() or "usage" in out.lower()


def test_dunder_main_no_args_prints_help(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["pycsamt"])
    with pytest.raises(SystemExit) as exc:
        runpy.run_module("pycsamt.cli.__main__", run_name="__main__")
    assert exc.value.code == 0
    out = capsys.readouterr().out
    assert len(out) > 0


def test_dunder_main_not_executed_on_plain_import():
    import pycsamt.cli.__main__ as mod

    assert hasattr(mod, "main")
