# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Supplementary coverage tests for :mod:`pycsamt.agents.metrics`.

Complements ``test_metrics_agent.py`` (which already exercises every
metric kind against real 3-EDI data) with: the ``ensure_sites`` exception
guard, the per-kind compute-exception guard, the ``_m_summary``
sub-metric exception guard, ``_fmt_deg``, and the empty-text branch of
``looks_like_metric_query``.
"""

from __future__ import annotations

import pytest

from pycsamt.agents.metrics import MetricsAgent, _fmt_deg, looks_like_metric_query

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


def test_looks_like_metric_query_empty_text():
    assert looks_like_metric_query("") is False


def test_fmt_deg():
    assert _fmt_deg(42.4) == "42°"


def test_ensure_sites_exception_fails(monkeypatch):
    import pycsamt.emtools._core as core

    monkeypatch.setattr(
        core,
        "ensure_sites",
        lambda *a, **k: (_ for _ in ()).throw(ValueError("bad sites")),
    )
    agent = MetricsAgent()
    result = agent.execute({"path": "/nope", "kinds": ["strike"]})
    assert result.status == "failed"
    assert "Could not load survey data" in result.error


def test_per_kind_compute_exception_is_recorded(monkeypatch, edi_dir):
    agent = MetricsAgent()
    monkeypatch.setattr(
        agent,
        "_m_strike",
        lambda sites, warnings: (_ for _ in ()).throw(
            RuntimeError("strike boom")
        ),
    )
    result = agent.execute({"path": str(edi_dir), "kinds": ["strike"]})
    assert result.status == "success"
    assert result.data["values"]["strike"] == "could not be computed"
    assert any("strike boom" in w for w in result.warnings)


def test_summary_sub_metric_exception_is_recorded(monkeypatch, edi_dir):
    agent = MetricsAgent()
    monkeypatch.setattr(
        agent,
        "_m_quality",
        lambda sites, warnings: (_ for _ in ()).throw(
            RuntimeError("quality boom")
        ),
    )
    result = agent.execute({"path": str(edi_dir), "kinds": ["summary"]})
    assert result.status == "success"
    assert any("quality boom" in w for w in result.warnings)
