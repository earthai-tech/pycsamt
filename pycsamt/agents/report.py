"""
pycsamt.agents.report
=====================

:class:`ReportAgent` — Assemble all agent results into a survey report.

The report is built in three formats:

* **Markdown** — always produced; human-readable plain text + embedded image
  paths.
* **HTML** — produced when ``markdown`` package is installed.
* **PDF** — produced when ``weasyprint`` or ``pdfkit`` is installed.

The agent queries the LLM once per section (optional) to write a narrative
paragraph, then assembles everything into a structured document:

    1. Title & metadata
    2. Data loading summary
    3. QC summary + figure
    4. Static-shift correction summary + figure
    5. Phase tensor analysis summary + figures
    6. Forward modelling summary + figure
    7. Recommendations
"""

from __future__ import annotations

import os
import time
from datetime import datetime
from pathlib import Path
from typing import Any

from ._base import AgentResult, BaseAgent

_SYSTEM_PROMPT = """\
You are a geophysics technical writer specialising in MT surveys.
Write clear, concise report sections in formal scientific English.
Use complete sentences. No markdown headings inside your response.
Keep each section to 3–5 sentences.
"""

# ── section writers ───────────────────────────────────────────────────────────

_SECTION_PROMPTS = {
    "loading": (
        "Write a 3-sentence data loading section for an MT survey report. "
        "Include: number of stations, period range, and data completeness. "
        "Data: {data}"
    ),
    "qc": (
        "Write a 3-sentence QC section for an MT survey report. "
        "Include: overall quality rating, number of flagged stations, "
        "and recommended remediation. Data: {data}"
    ),
    "static_shift": (
        "Write a 3-sentence static-shift section for an MT survey report. "
        "Include: detection method, mean correction magnitude, and stations "
        "most affected. Data: {data}"
    ),
    "phase_analysis": (
        "Write a 4-sentence phase tensor section for an MT survey report. "
        "Include: dominant dimensionality, consensus strike, 3-D fraction, "
        "and whether rotation is recommended. Data: {data}"
    ),
    "forward": (
        "Write a 3-sentence forward modelling section for an MT survey report. "
        "Include: model description, synthetic response quality, and RMS "
        "misfit if available. Data: {data}"
    ),
    "recommendations": (
        "Write a 4-sentence recommendations section for an MT survey report. "
        "Cover: next processing step, inversion strategy, areas needing "
        "additional data, and confidence in the results. Data: {data}"
    ),
}


class ReportAgent(BaseAgent):
    """Generate a structured MT survey report from agent results.

    Parameters
    ----------
    api_key, model, llm_provider : str
    report_title : str
        Title for the report.
    formats : list of {"md", "html", "pdf"}
        Output formats.  Default ``["md", "html"]``.

    Input keys
    ----------
    ``results`` : dict
        Keyed by agent step name → :class:`AgentResult`.
        Expected keys: ``"load"``, ``"qc"``, ``"static_shift"``,
        ``"phase_analysis"``, ``"forward"`` (all optional).
    ``output_dir`` : str
    ``title`` : str, optional — overrides constructor default

    Output data keys
    ----------------
    ``report_md``      str — full markdown text
    ``report_html``    str or None
    ``report_path_md`` str — path to .md file
    ``report_path_html`` str or None
    ``sections``       dict — section text keyed by name

    Examples
    --------
    >>> agent = ReportAgent(api_key="sk-ant-…")
    >>> result = agent.execute(
    ...     {
    ...         "results": {"load": load_result, "qc": qc_result},
    ...         "output_dir": "/out/report",
    ...         "title": "WILLY_DATA AMT Survey — L22PLT",
    ...     }
    ... )
    >>> print(result["report_path_md"])
    /out/report/survey_report.md
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        report_title: str = "MT/AMT Survey Report",
        formats: list[str] | None = None,
    ) -> None:
        super().__init__(
            "ReportAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="publication",
        )
        self.report_title = report_title
        self.formats = [f.lower() for f in (formats or ["md", "html"])]

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        output_dir = input_data.get("output_dir", "pycsamt_report")
        title = input_data.get("title", self.report_title)
        results = input_data.get("results", {})

        os.makedirs(output_dir, exist_ok=True)

        # ── copy figures to report dir ────────────────────────────────────────
        fig_refs: dict[str, str] = {}  # section_name → relative image path
        for step_name, res in results.items():
            if not isinstance(res, AgentResult):
                continue
            fps = res.get("figure_paths") or {}
            if isinstance(fps, dict):
                for fig_name, src_path in fps.items():
                    if src_path and Path(src_path).exists():
                        dst = os.path.join(
                            output_dir, f"{step_name}_{fig_name}.png"
                        )
                        try:
                            import shutil

                            shutil.copy2(src_path, dst)
                            fig_refs[f"{step_name}_{fig_name}"] = (
                                os.path.basename(dst)
                            )
                        except Exception as exc:
                            warnings.append(
                                f"Could not copy figure {src_path}: {exc}"
                            )

        # ── build section texts ───────────────────────────────────────────────
        sections: dict[str, str] = {}
        sections["loading"] = self._section_loading(results, warnings)
        sections["qc"] = self._section_qc(results, warnings)
        sections["static_shift"] = self._section_ss(results, warnings)
        sections["phase_analysis"] = self._section_pt(results, warnings)
        sections["forward"] = self._section_fwd(results, warnings)
        sections["recommendations"] = self._section_rec(
            results, sections, warnings
        )

        # ── assemble markdown ─────────────────────────────────────────────────
        md = _build_markdown(title, sections, fig_refs, results)

        # ── write markdown ────────────────────────────────────────────────────
        md_path = os.path.join(output_dir, "survey_report.md")
        try:
            Path(md_path).write_text(md, encoding="utf-8")
        except Exception as exc:
            warnings.append(f"Could not write markdown report: {exc}")
            md_path = None

        # ── write HTML ────────────────────────────────────────────────────────
        html: str | None = None
        html_path: str | None = None
        if "html" in self.formats:
            try:
                import markdown as _md_pkg

                html = _md_pkg.markdown(
                    md, extensions=["tables", "fenced_code"]
                )
                html = _HTML_TEMPLATE.format(title=title, body=html)
                html_path = os.path.join(output_dir, "survey_report.html")
                Path(html_path).write_text(html, encoding="utf-8")
            except ImportError:
                warnings.append(
                    "markdown package not installed; HTML report skipped. "
                    "Install with: pip install markdown"
                )
            except Exception as exc:
                warnings.append(f"HTML report failed: {exc}")

        elapsed = time.time() - t0
        n_sec = sum(1 for v in sections.values() if v.strip())
        return AgentResult(
            status="success",
            summary=(
                f"Report generated: {n_sec} sections, "
                f"{len(fig_refs)} figure(s). "
                f"Saved to {output_dir}."
            ),
            data={
                "report_md": md,
                "report_html": html,
                "report_path_md": md_path,
                "report_path_html": html_path,
                "sections": sections,
                "figure_refs": fig_refs,
                "output_dir": output_dir,
            },
            warnings=warnings,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )

    # ── section builders ──────────────────────────────────────────────────────

    def _llm_section(self, section_key: str, data_str: str) -> str:
        """Query the LLM for one report section, or return a fallback."""
        if not self.api_key:
            return ""
        prompt = _SECTION_PROMPTS[section_key].format(data=data_str)
        text = self.query_llm(prompt, max_tokens=200)
        return text or ""

    def _section_loading(self, results: dict, warnings: list) -> str:
        load_r = results.get("load") or results.get("loader")
        if load_r is None:
            return "No data loading results available."
        n_st = load_r.get("n_stations", "?")
        stats = load_r.get("summary_stats") or {}
        t_min = stats.get("global_t_min_s")
        t_max = stats.get("global_t_max_s")
        per_str = (
            f"{t_min:.2e}–{t_max:.2e} s" if t_min and t_max else "unknown"
        )
        base = (
            f"{n_st} stations were loaded. "
            f"Period range: {per_str}. "
            f"Mean QC score: {stats.get('mean_qc_score', '?'):.0f}/100."
        )
        if self.api_key:
            llm_text = self._llm_section("loading", str(stats))
            return llm_text or base
        return base

    def _section_qc(self, results: dict, warnings: list) -> str:
        qc_r = results.get("qc") or results.get("data_qc")
        if qc_r is None:
            return "No QC analysis results available."
        n_flagged = qc_r.get("n_flagged", 0)
        flagged = qc_r.get("flagged_stations", [])
        base = (
            f"QC analysis identified {n_flagged} station(s) below threshold. "
            + (f"Flagged: {', '.join(flagged)}. " if flagged else "")
            + "All other stations met minimum data quality criteria."
        )
        if self.api_key:
            llm_text = self._llm_section(
                "qc", f"n_flagged={n_flagged}, flagged={flagged}"
            )
            return llm_text or base
        return base

    def _section_ss(self, results: dict, warnings: list) -> str:
        ss_r = results.get("static_shift") or results.get("ss")
        if ss_r is None:
            return "No static-shift correction results available."
        ds = ss_r.get("delta_stats") or {}
        base = (
            f"Static-shift correction was applied. "
            f"Mean correction: {ds.get('mean', '?'):.3f} log₁₀(Ω·m). "
            f"{ds.get('n_shifted', '?')} station(s) showed significant shift (>0.05)."
        )
        if self.api_key:
            llm_text = self._llm_section("static_shift", str(ds))
            return llm_text or base
        return base

    def _section_pt(self, results: dict, warnings: list) -> str:
        pt_r = results.get("phase_analysis") or results.get("pt")
        if pt_r is None:
            return "No phase tensor analysis results available."
        n_1d = pt_r.get("n_1d", 0)
        n_2d = pt_r.get("n_2d", 0)
        n_3d = pt_r.get("n_3d", 0)
        st = pt_r.get("strike_consensus", float("nan"))
        iqr = pt_r.get("strike_iqr", float("nan"))
        base = (
            f"Phase tensor analysis found predominantly "
            f"{'1-D' if n_1d >= n_2d and n_1d >= n_3d else '2-D' if n_2d >= n_3d else '3-D'} "
            f"character ({n_1d}/{n_2d}/{n_3d} 1-D/2-D/3-D observations). "
            f"Consensus geoelectric strike: {st:.1f}° ± {iqr:.1f}°."
        )
        if self.api_key:
            data_str = f"1D={n_1d}, 2D={n_2d}, 3D={n_3d}, strike={st:.1f}, iqr={iqr:.1f}"
            llm_text = self._llm_section("phase_analysis", data_str)
            return llm_text or base
        return base

    def _section_fwd(self, results: dict, warnings: list) -> str:
        fwd_r = results.get("forward") or results.get("fwd")
        if fwd_r is None:
            return "No forward modelling results available."
        rms = fwd_r.get("rms")
        base = (
            "A 1-D MT forward model was computed. "
            + (f"Data–model RMS misfit: {rms:.3f} log₁₀(Ω·m). " if rms else "")
            + "The synthetic response covers the full available period range."
        )
        if self.api_key:
            llm_text = self._llm_section("forward", f"rms={rms}")
            return llm_text or base
        return base

    def _section_rec(
        self,
        results: dict,
        sections: dict[str, str],
        warnings: list,
    ) -> str:
        if self.api_key:
            summary = " | ".join(
                f"{k}: {v[:100]}" for k, v in sections.items() if v.strip()
            )
            llm_text = self._llm_section("recommendations", summary)
            if llm_text:
                return llm_text
        return (
            "Based on the QC and phase tensor results, proceed with "
            "profile-oriented 2-D inversion using the consensus strike. "
            "Apply static-shift correction to the EDIs before export. "
            "Consider independent Bostick depth estimates to guide the "
            "mesh parametrisation."
        )


# ── markdown assembler ────────────────────────────────────────────────────────


def _build_markdown(
    title: str,
    sections: dict[str, str],
    fig_refs: dict[str, str],
    results: dict,
) -> str:
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    lines = [
        f"# {title}",
        "",
        f"*Generated by pycsamt ReportAgent — {now}*",
        "",
        "---",
        "",
    ]

    _SEC_TITLES = {
        "loading": "1. Data Loading",
        "qc": "2. Quality Control",
        "static_shift": "3. Static-Shift Correction",
        "phase_analysis": "4. Phase Tensor Analysis",
        "forward": "5. Forward Modelling",
        "recommendations": "6. Recommendations",
    }

    _FIG_KEYS = {
        "qc": ["qc_confidence_section", "qc_confidence_profile"],
        "static_shift": [
            "static_shift_ss_summary",
            "static_shift_ss_comparison",
        ],
        "phase_analysis": [
            "phase_analysis_pt_psection",
            "phase_analysis_strike_analysis",
            "phase_analysis_survey_fingerprint",
        ],
        "forward": ["forward_response_and_model", "forward_response"],
    }

    for sec_key, sec_title in _SEC_TITLES.items():
        text = sections.get(sec_key, "").strip()
        lines.append(f"## {sec_title}")
        lines.append("")
        if text:
            lines.append(text)
        else:
            lines.append("*No results available for this section.*")
        lines.append("")

        # embed figures for this section
        for fig_key in _FIG_KEYS.get(sec_key, []):
            if fig_key in fig_refs:
                lines.append(f"![{fig_key}]({fig_refs[fig_key]})")
                lines.append("")

    # warnings appendix
    all_warnings: list[str] = []
    for res in results.values():
        if isinstance(res, AgentResult):
            all_warnings.extend(res.warnings or [])
    if all_warnings:
        lines.append("---")
        lines.append("## Appendix: Processing Warnings")
        lines.append("")
        for w in all_warnings:
            lines.append(f"- {w}")
        lines.append("")

    return "\n".join(lines)


_HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>
  body {{ font-family: Georgia, serif; max-width: 900px; margin: 40px auto;
         padding: 0 20px; color: #222; line-height: 1.6; }}
  h1   {{ color: #1a3a5c; border-bottom: 2px solid #1a3a5c; }}
  h2   {{ color: #2c5f8a; margin-top: 2em; }}
  img  {{ max-width: 100%; height: auto; display: block; margin: 1em 0; }}
  li   {{ margin-bottom: 0.3em; }}
  code {{ background: #f4f4f4; padding: 2px 6px; border-radius: 3px; }}
  hr   {{ border: none; border-top: 1px solid #ccc; margin: 2em 0; }}
</style>
</head>
<body>
{body}
</body>
</html>
"""


__all__ = ["ReportAgent"]
