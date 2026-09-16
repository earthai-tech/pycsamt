# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for ``python -m pycsamt.assistant.rag`` (__main__.py)."""

from __future__ import annotations

import runpy
import sys

import pytest


def test_dunder_main_help_exits_zero(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["pycsamt.assistant.rag", "--help"])
    with pytest.raises(SystemExit) as exc:
        runpy.run_module(
            "pycsamt.assistant.rag.__main__", run_name="__main__"
        )
    assert exc.value.code == 0
    out = capsys.readouterr().out
    assert "usage" in out.lower()


def test_dunder_main_no_args_exits_nonzero(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["pycsamt.assistant.rag"])
    with pytest.raises(SystemExit) as exc:
        runpy.run_module(
            "pycsamt.assistant.rag.__main__", run_name="__main__"
        )
    assert exc.value.code != 0


def test_dunder_main_stats_command_runs(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(
        sys,
        "argv",
        ["pycsamt.assistant.rag", "--out", str(tmp_path), "stats"],
    )
    with pytest.raises(SystemExit) as exc:
        runpy.run_module(
            "pycsamt.assistant.rag.__main__", run_name="__main__"
        )
    assert exc.value.code == 1
    out = capsys.readouterr().out
    assert "No persisted index" in out


def test_dunder_main_not_executed_on_plain_import():
    import pycsamt.assistant.rag.__main__ as mod

    assert hasattr(mod, "main")
