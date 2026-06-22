# -*- coding: utf-8 -*-
"""Tests for InputBuilder."""

import pytest
from pathlib import Path


def test_builder_imports():
    from pycsamt.models.occam2d.builder import InputBuilder
    assert InputBuilder is not None


def test_builder_not_ready_before_build(tmp_path):
    from pycsamt.models.occam2d.builder import InputBuilder
    builder = InputBuilder(source=None, workdir=tmp_path)
    assert not builder.is_ready


def test_builder_summary_not_built(tmp_path):
    from pycsamt.models.occam2d.builder import InputBuilder
    builder = InputBuilder(source=None, workdir=tmp_path)
    assert "not yet built" in builder.summary()
