# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Targeted coverage for :mod:`pycsamt.pipeline._config`.

``test_pipeline.py::TestPipelineConfigFiles`` already exercises the
normal load_yaml/load_json/load_py/pipeline_to_yaml round trips with
PyYAML installed. This file adds:

* the ``ImportError`` branch of :func:`load_yaml` (PyYAML unavailable),
* the ``output_dir`` branch of :func:`pipeline_to_yaml`,
* its own PyYAML-unavailable fallback (hand-crafted YAML text).
"""

from __future__ import annotations

import sys

import pytest

from pycsamt.pipeline._config import load_yaml, pipeline_to_yaml
from pycsamt.pipeline._steps import Step


@pytest.fixture()
def no_pyyaml(monkeypatch):
    """Make ``import yaml`` raise ImportError inside the module under test.

    Setting ``sys.modules['yaml'] = None`` is the standard trick: the
    import system special-cases a ``None`` entry and raises ImportError
    instead of a fresh import attempt, without needing PyYAML to be
    actually uninstalled.
    """
    monkeypatch.setitem(sys.modules, "yaml", None)
    yield


def test_load_yaml_raises_import_error_without_pyyaml(tmp_path, no_pyyaml):
    p = tmp_path / "cfg.yaml"
    p.write_text("name: x\n", encoding="utf-8")
    with pytest.raises(ImportError, match="PyYAML is required"):
        load_yaml(p)


def test_pipeline_to_yaml_output_dir_included():
    steps = [("notch", Step("NR001"))]
    text = pipeline_to_yaml(steps, name="wf", output_dir="results/")
    assert "output_dir" in text
    assert "results/" in text


def test_pipeline_to_yaml_fallback_without_pyyaml(no_pyyaml):
    steps = [
        ("notch", Step("NR001", mains_hz=50)),
        ("align", Step("FREQ004")),
    ]
    text = pipeline_to_yaml(steps, name="fallback_wf", output_dir="out/")
    assert "name: 'fallback_wf'" in text
    assert "output_dir: 'out/'" in text
    assert "NR001" in text
    assert "FREQ004" in text
    assert "mains_hz: 50" in text


def test_pipeline_to_yaml_fallback_no_output_dir_no_params(no_pyyaml):
    steps = [("align", Step("FREQ004"))]
    text = pipeline_to_yaml(steps, name="wf2")
    assert "output_dir" not in text
    assert "params" not in text
    assert "FREQ004" in text
