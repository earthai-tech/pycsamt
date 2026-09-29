# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Targeted coverage for :mod:`pycsamt.pipeline._output` error paths.

The happy paths (setup, step_plot_dir, save_pipeline_config, save_text,
repr) are already covered in ``test_pipeline.py::TestOutputDir`` and the
full-pipeline run in ``test_pipeline_integration.py``. This file adds the
failure branches (each wrapped in a ``try/except -> warnings.warn``) and
the config-fallback / api-override paths that neither of those exercise.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from pycsamt.pipeline._output import OutputDir


class _RaisingFigure:
    def savefig(self, *a, **k):
        raise RuntimeError("disk full")


class _OkFigure:
    def __init__(self):
        self.calls = []

    def savefig(self, *a, **k):
        self.calls.append((a, k))


def test_save_figure_success_returns_path(tmp_path):
    od = OutputDir(tmp_path / "out")
    od.setup()
    fig = _OkFigure()
    p = od.save_figure(fig, "myplot", 1, "notch")
    assert p is not None
    assert p.name == "myplot.png"
    assert fig.calls


def test_save_figure_failure_warns_and_returns_none(tmp_path):
    od = OutputDir(tmp_path / "out")
    od.setup()
    with pytest.warns(UserWarning, match="Failed to save figure"):
        result = od.save_figure(_RaisingFigure(), "bad", 1, "notch")
    assert result is None


def test_save_figure_accepts_api_override(tmp_path):
    from pycsamt.api.pipe import PipelineAPIConfig

    od = OutputDir(tmp_path / "out")
    od.setup()
    fig = _OkFigure()
    cfg = PipelineAPIConfig(plot_fmt="svg", plot_dpi=72)
    p = od.save_figure(fig, "override", 2, "select_band", api=cfg)
    assert p is not None
    assert p.suffix == ".svg"


def test_write_edis_failure_warns_and_returns_empty(tmp_path, monkeypatch):
    import pycsamt.site.export as site_export

    def _boom(*a, **k):
        raise RuntimeError("export failed")

    monkeypatch.setattr(site_export, "write_sites", _boom)
    od = OutputDir(tmp_path / "out")
    od.setup()
    with pytest.warns(UserWarning, match="EDI export failed"):
        out = od.write_edis(sites=[])
    assert out == []


def test_save_pipeline_config_failure_warns_and_returns_none(
    tmp_path, monkeypatch
):
    od = OutputDir(tmp_path / "out")
    od.setup()

    def _boom(self, *a, **k):
        raise OSError("no space left")

    monkeypatch.setattr(Path, "write_text", _boom)
    with pytest.warns(UserWarning, match="Could not save pipeline.yaml"):
        result = od.save_pipeline_config("name: x\n")
    assert result is None


def test_save_text_failure_warns_and_returns_none(tmp_path, monkeypatch):
    od = OutputDir(tmp_path / "out")
    od.setup()

    def _boom(self, *a, **k):
        raise OSError("no space left")

    monkeypatch.setattr(Path, "write_text", _boom)
    with pytest.warns(UserWarning, match="Could not save summary.txt"):
        result = od.save_text("hello", "summary.txt")
    assert result is None


def test_root_none_falls_back_to_global_singleton_output_root():
    from pycsamt.api.pipe import PYCSAMT_PIPE

    od = OutputDir(None)
    assert od.root == Path(PYCSAMT_PIPE.output_root)


def test_cfg_uses_injected_api():
    from pycsamt.api.pipe import PipelineAPIConfig

    cfg = PipelineAPIConfig(output_root="somewhere")
    od = OutputDir(None, api=cfg)
    assert od._cfg() is cfg
    assert od.root == Path("somewhere")
