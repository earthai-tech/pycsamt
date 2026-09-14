# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.agents.report`.

Targets: the figure-copy loop's isinstance/success branches, the
markdown-write exception guard, the HTML ImportError and generic
exception guards, every ``_section_*`` builder's "results present"
branch (static_shift / phase_analysis / forward were previously
untested), the LLM-backed section paths, and ``_build_markdown``'s
no-results fallback, figure-embedding, and warnings-appendix blocks.
"""

from __future__ import annotations

import sys

import pytest

from pycsamt.agents._base import AgentResult
from pycsamt.agents.report import ReportAgent, _build_markdown

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


def _full_results(tmp_path):
    fig_qc = tmp_path / "confidence_section.png"
    fig_qc.write_bytes(b"fake png")
    fig_ss = tmp_path / "ss_summary.png"
    fig_ss.write_bytes(b"fake png")

    load = AgentResult(
        status="success",
        summary="loaded",
        data={
            "n_stations": 5,
            "summary_stats": {
                "global_t_min_s": 1e-3,
                "global_t_max_s": 100.0,
                "mean_qc_score": 87.0,
            },
        },
    )
    qc = AgentResult(
        status="success",
        summary="qc",
        data={
            "n_flagged": 1,
            "flagged_stations": ["s3"],
            "figure_paths": {"confidence_section": str(fig_qc)},
        },
        warnings=["QC dropped a bad frequency"],
    )
    static_shift = AgentResult(
        status="success",
        summary="ss",
        data={
            "delta_stats": {"mean": 0.12, "n_shifted": 3},
            "figure_paths": {"ss_summary": str(fig_ss)},
        },
    )
    phase_analysis = AgentResult(
        status="success",
        summary="pt",
        data={
            "n_1d": 10,
            "n_2d": 20,
            "n_3d": 5,
            "strike_consensus": 42.5,
            "strike_iqr": 8.1,
        },
    )
    forward = AgentResult(
        status="success",
        summary="fwd",
        data={"rms": 0.15},
    )
    return {
        "load": load,
        "qc": qc,
        "static_shift": static_shift,
        "phase_analysis": phase_analysis,
        "forward": forward,
    }


def test_full_report_all_sections_and_figures(tmp_path):
    results = _full_results(tmp_path)
    outdir = tmp_path / "report_out"
    agent = ReportAgent()
    result = agent.execute(
        {"results": results, "output_dir": str(outdir)}
    )
    assert result.status == "success"
    md = result["report_md"]
    assert "5 stations were loaded" in md
    assert "Static-shift correction was applied" in md
    assert "Phase tensor analysis found predominantly" in md
    assert "1-D MT forward model was computed" in md
    assert "![qc_confidence_section]" in md
    assert "![static_shift_ss_summary]" in md
    assert "Appendix: Processing Warnings" in md
    assert "QC dropped a bad frequency" in md
    assert result["figure_refs"]
    assert (outdir / "qc_confidence_section.png").exists()


def test_non_agentresult_entries_are_skipped(tmp_path):
    # Plain dicts still satisfy the section builders' `.get(...)` calls,
    # but are not AgentResult instances, so the figure-copy loop's
    # isinstance guard skips them.
    results = {"load": {}, "qc": {}}
    agent = ReportAgent()
    result = agent.execute(
        {"results": results, "output_dir": str(tmp_path / "out")}
    )
    assert result.status == "success"
    assert result["figure_refs"] == {}


def test_no_results_sections_fall_back_to_default_text(tmp_path):
    agent = ReportAgent()
    result = agent.execute(
        {"results": {}, "output_dir": str(tmp_path / "out")}
    )
    assert result.status == "success"
    assert "No data loading results available" in result["sections"]["loading"]
    assert "No QC analysis results available" in result["sections"]["qc"]
    assert (
        "No static-shift correction results available"
        in result["sections"]["static_shift"]
    )
    assert (
        "No phase tensor analysis results available"
        in result["sections"]["phase_analysis"]
    )
    assert (
        "No forward modelling results available"
        in result["sections"]["forward"]
    )


def test_build_markdown_empty_section_uses_placeholder():
    md = _build_markdown(
        "Title",
        {"loading": "", "qc": "text"},
        {},
        {},
    )
    assert "*No results available for this section.*" in md


def test_markdown_write_exception_is_recorded(tmp_path):
    outdir = tmp_path / "out"
    outdir.mkdir()
    # Pre-create the target path as a directory so Path.write_text raises.
    (outdir / "survey_report.md").mkdir()
    agent = ReportAgent()
    result = agent.execute({"results": {}, "output_dir": str(outdir)})
    assert result["report_path_md"] is None
    assert any(
        "Could not write markdown report" in w for w in result.warnings
    )


def test_html_import_error_is_recorded(tmp_path, monkeypatch):
    monkeypatch.setitem(sys.modules, "markdown", None)
    agent = ReportAgent()
    result = agent.execute(
        {"results": {}, "output_dir": str(tmp_path / "out")}
    )
    assert result["report_html"] is None
    assert any(
        "markdown package not installed" in w for w in result.warnings
    )


def test_html_generic_exception_is_recorded(tmp_path, monkeypatch):
    import markdown as md_pkg

    monkeypatch.setattr(
        md_pkg,
        "markdown",
        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("md boom")),
    )
    agent = ReportAgent()
    result = agent.execute(
        {"results": {}, "output_dir": str(tmp_path / "out")}
    )
    assert result["report_html"] is None
    assert any("HTML report failed" in w for w in result.warnings)


def test_llm_backed_sections_use_query_llm(tmp_path):
    results = _full_results(tmp_path)
    agent = ReportAgent(api_key="fake-key")
    agent.query_llm = lambda *a, **k: "LLM narrative text."
    result = agent.execute(
        {"results": results, "output_dir": str(tmp_path / "out")}
    )
    assert result.status == "success"
    for key in (
        "loading",
        "qc",
        "static_shift",
        "phase_analysis",
        "forward",
        "recommendations",
    ):
        assert result["sections"][key] == "LLM narrative text."


def test_llm_empty_response_falls_back_to_base_text(tmp_path):
    results = _full_results(tmp_path)
    agent = ReportAgent(api_key="fake-key")
    agent.query_llm = lambda *a, **k: ""
    result = agent.execute(
        {"results": results, "output_dir": str(tmp_path / "out")}
    )
    assert "Static-shift correction was applied" in result["sections"][
        "static_shift"
    ]
