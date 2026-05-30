# -*- coding: utf-8 -*-
"""Tests for OccamRunner."""

import pytest
from pathlib import Path


def test_runner_imports():
    from pycsamt.models.occam2d.runner import OccamRunner
    assert OccamRunner is not None


def test_runner_no_binary_raises(tmp_path):
    """discover_binary should raise when binary is absent and compile fails."""
    from pycsamt.models.occam2d.runner import OccamRunner
    import shutil
    runner = OccamRunner(workdir=tmp_path)
    # If gfortran is not on PATH, compile() will raise RuntimeError
    # If it is, it should still fail because _source/ is bundled — just
    # verify that discover_binary() raises something meaningful.
    with pytest.raises((FileNotFoundError, RuntimeError)):
        runner.discover_binary(auto_compile=False)


def test_runner_is_not_running_by_default(tmp_path):
    from pycsamt.models.occam2d.runner import OccamRunner
    runner = OccamRunner(workdir=tmp_path)
    assert not runner.is_running
