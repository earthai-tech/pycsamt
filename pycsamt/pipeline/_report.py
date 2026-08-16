"""Pipeline report generators.

:func:`make_text_report` produces a compact plain-text summary.
:func:`make_html_report` produces a self-contained HTML report with per-step
cards, timing, plot thumbnails (linked, not embedded), and a processing log.
"""

from __future__ import annotations

import datetime
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._steps import StepResult


# ---------------------------------------------------------------------------
# Plain-text report
# ---------------------------------------------------------------------------


def make_text_report(
    pipeline_name: str,
    step_results: list[StepResult],
    elapsed_sec: float,
    outdir: Path | None,
    n_sites_in: int,
    n_sites_out: int,
) -> str:
    """Return a compact plain-text pipeline summary."""
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    width = 72
    sep = "─" * width

    lines = [
        sep,
        "  pyCSAMT Pipeline Report",
        f"  Pipeline : {pipeline_name}",
        f"  Run at   : {now}",
        f"  Sites    : {n_sites_in} in → {n_sites_out} out",
        f"  Total    : {elapsed_sec:.2f}s",
        sep,
        "",
        "  Step results",
        "  " + "─" * (width - 2),
        f"  {'#':>3}  {'Name':<22} {'Code':<8} {'Status':<5} "
        f"{'Time(s)':>7}  {'Sites':>10}  {'Plots':>5}  {'Cached':>6}",
        "  " + "─" * (width - 2),
    ]
    for sr in step_results:
        status = "OK " if sr.ok else "ERR"
        sites_str = f"{sr.n_sites_in}→{sr.n_sites_out}"
        cached_str = "yes" if sr.cached else "-"
        lines.append(
            f"  {sr.step_idx:>3}  {sr.step_name:<22} "
            f"{sr.step_code:<8} {status:<5} "
            f"{sr.elapsed_sec:>7.2f}  {sites_str:>10}  "
            f"{len(sr.plots):>5}  {cached_str:>6}"
        )
    lines.append("  " + "─" * (width - 2))

    # Errors section
    errors = [sr for sr in step_results if not sr.ok]
    if errors:
        lines += ["", "  Errors"]
        for sr in errors:
            lines.append(f"  [{sr.step_code}] {sr.step_name}: {sr.error}")

    # Output directory
    if outdir is not None:
        lines += [
            "",
            f"  Output → {outdir}",
            f"    EDIs   : {outdir / 'processed'}",
            f"    Plots  : {outdir / 'plots'}",
        ]

    lines += ["", sep, ""]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------

# pyCSAMT brand tokens (docs/source/_static/css/custom.css) — the fast tier
# gets a palette swap only, no new markup: real ink/surface/status colors in
# place of the old generic blue/green/red, kept cheap to render and small.
# Status colors (good/warning/critical) are the dataviz skill's validated,
# brand-independent set — see pycsamt/pipeline/_dashboard.py for the full
# provenance and the richer "dashboard" report tier built on the same tokens.
_CSS = """
<style>
  body { font-family: system-ui, sans-serif; margin: 2em; color: #24324b; }
  h1 { border-bottom: 2px solid #3e65b0; padding-bottom: 0.4em; }
  h2 { color: #3e65b0; margin-top: 1.5em; }
  table { border-collapse: collapse; width: 100%; }
  th, td { border: 1px solid #dde3ee; padding: 0.4em 0.7em; text-align: left; }
  th { background: #3e65b0; color: white; }
  tr:nth-child(even) { background: #f6f8fc; }
  .ok   { color: #0ca30c; font-weight: bold; }
  .err  { color: #d03b3b; font-weight: bold; }
  .step-card { border: 1px solid #dde3ee; border-radius: 6px;
               padding: 1em; margin: 1em 0; }
  .step-card h3 { margin: 0 0 0.5em; }
  .step-meta { font-size: 0.85em; color: #5c677d; }
  .plots { display: flex; flex-wrap: wrap; gap: 0.5em; margin-top: 0.5em; }
  .plots a img { max-height: 120px; border: 1px solid #dde3ee;
                 border-radius: 4px; }
  .error-box { background: #fdecea; border-left: 4px solid #d03b3b;
               padding: 0.5em 1em; font-size: 0.9em; }
  footer { margin-top: 3em; font-size: 0.8em; color: #5c677d;
           border-top: 1px solid #dde3ee; padding-top: 1em; }
</style>
"""


def _plot_thumbs(plots: list[Path], outdir: Path | None) -> str:
    """Return HTML <a><img></a> tags for a list of plot paths."""
    if not plots or outdir is None:
        return ""
    tags = []
    for p in plots:
        try:
            rel = p.relative_to(outdir)
        except ValueError:
            rel = p
        tags.append(f'<a href="{rel}" target="_blank"><img src="{rel}" /></a>')
    return f'<div class="plots">{"".join(tags)}</div>'


def make_html_report(
    pipeline_name: str,
    step_results: list[StepResult],
    elapsed_sec: float,
    outdir: Path | None,
    n_sites_in: int,
    n_sites_out: int,
    pipeline_yaml: str = "",
) -> str:
    """Return a self-contained HTML pipeline report."""
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    n_ok = sum(1 for sr in step_results if sr.ok)
    n_err = len(step_results) - n_ok
    status_cls = "ok" if n_err == 0 else "err"
    status_txt = "SUCCESS" if n_err == 0 else f"FAILED ({n_err} error(s))"

    # ── Summary table ─────────────────────────────────────────────────
    summary_rows = (
        f"<tr><th>Pipeline</th><td>{pipeline_name}</td></tr>"
        f"<tr><th>Run at</th><td>{now}</td></tr>"
        f"<tr><th>Status</th>"
        f"<td><span class='{status_cls}'>{status_txt}</span></td></tr>"
        f"<tr><th>Sites</th><td>{n_sites_in} → {n_sites_out}</td></tr>"
        f"<tr><th>Elapsed</th><td>{elapsed_sec:.2f} s</td></tr>"
        f"<tr><th>Steps</th><td>{len(step_results)} ({n_ok} ok, {n_err} err)</td></tr>"
    )

    # ── Step cards ────────────────────────────────────────────────────
    step_cards_html = []
    for sr in step_results:
        status_icon = "&#10003;" if sr.ok else "&#10007;"
        status_c = "ok" if sr.ok else "err"
        params_str = (
            ", ".join(f"{k}={v!r}" for k, v in sr.params.items())
            if sr.params
            else "—"
        )
        thumbs = _plot_thumbs(sr.plots, outdir)
        cached_badge = " &nbsp;|&nbsp; <b>Cached:</b> yes" if sr.cached else ""
        err_html = ""
        if not sr.ok:
            err_html = f'<div class="error-box"><b>Error:</b> {sr.error}</div>'
        card = (
            f'<div class="step-card">'
            f'<h3><span class="{status_c}">{status_icon}</span> '
            f"({sr.step_idx}) {sr.step_name} "
            f"<small>[{sr.step_code}]</small></h3>"
            f'<div class="step-meta">'
            f"<b>Label:</b> {sr.step_label} &nbsp;|&nbsp; "
            f"<b>Sites:</b> {sr.n_sites_in} → {sr.n_sites_out} &nbsp;|&nbsp; "
            f"<b>Time:</b> {sr.elapsed_sec:.2f} s &nbsp;|&nbsp; "
            f"<b>Plots:</b> {len(sr.plots)}{cached_badge}</div>"
            f"<p><b>Params:</b> <code>{params_str}</code></p>"
            f"{err_html}"
            f"{thumbs}"
            f"</div>"
        )
        step_cards_html.append(card)

    # ── Config YAML block ─────────────────────────────────────────────
    yaml_block = ""
    if pipeline_yaml:
        yaml_block = (
            "<h2>Pipeline configuration</h2>"
            f"<pre><code>{pipeline_yaml}</code></pre>"
        )

    # ── Assemble ──────────────────────────────────────────────────────
    html = (
        "<!DOCTYPE html>"
        '<html lang="en">'
        "<head>"
        '<meta charset="utf-8">'
        f"<title>pyCSAMT Pipeline — {pipeline_name}</title>"
        f"{_CSS}"
        "</head>"
        "<body>"
        f"<h1>pyCSAMT Pipeline Report</h1>"
        f"<table>{summary_rows}</table>"
        "<h2>Step results</h2>"
        + "".join(step_cards_html)
        + yaml_block
        + "<footer>Generated by pyCSAMT — "
        + now
        + "</footer>"
        "</body></html>"
    )
    return html


# ---------------------------------------------------------------------------
# Notebook (_repr_html_) — compact, scoped, no plot thumbnails/config block
# ---------------------------------------------------------------------------

# Same palette as _CSS above, but every selector is scoped under
# .pycsamt-result: this HTML is injected directly into a Jupyter output
# cell's DOM (not a standalone document), so unscoped selectors like the
# ones in _CSS ("table { width: 100% }", "th { background: ... }") would
# leak out and restyle every other table/output in the notebook.
_NOTEBOOK_CSS = """
<style>
  .pycsamt-result { font-family: system-ui, sans-serif; color: #222; }
  .pycsamt-result table { border-collapse: collapse; width: auto; }
  .pycsamt-result th, .pycsamt-result td {
    border: 1px solid #ddd; padding: 0.3em 0.6em; text-align: left;
    font-size: 0.9em;
  }
  .pycsamt-result th { background: #4a6fa5; color: white; }
  .pycsamt-result tr:nth-child(even) { background: #f5f7fb; }
  .pycsamt-result .ok  { color: #2d7a3e; font-weight: bold; }
  .pycsamt-result .err { color: #c0392b; font-weight: bold; }
  .pycsamt-result .cached { color: #8a6d00; }
  .pycsamt-result .title { font-weight: bold; margin-bottom: 0.4em; }
</style>
"""


def make_notebook_html(
    pipeline_name: str,
    step_results: list[StepResult],
    elapsed_sec: float,
    outdir: Path | None,
    n_sites_in: int,
    n_sites_out: int,
) -> str:
    """Return a compact, scoped HTML fragment for ``PipelineResult._repr_html_``.

    A quick-glance step table, not a substitute for :func:`make_html_report`
    (no plot thumbnails, no embedded config YAML).
    """
    n_ok = sum(1 for sr in step_results if sr.ok)
    n_err = len(step_results) - n_ok
    status_cls = "ok" if n_err == 0 else "err"
    status_txt = "OK" if n_err == 0 else f"{n_err} error(s)"

    rows = []
    for sr in step_results:
        status_c = "ok" if sr.ok else "err"
        status_t = "OK" if sr.ok else "ERR"
        cached = '<span class="cached">yes</span>' if sr.cached else "&ndash;"
        rows.append(
            "<tr>"
            f"<td>{sr.step_idx}</td>"
            f"<td>{sr.step_name}</td>"
            f"<td>{sr.step_code}</td>"
            f'<td><span class="{status_c}">{status_t}</span></td>'
            f"<td>{sr.elapsed_sec:.2f}s</td>"
            f"<td>{sr.n_sites_in}&rarr;{sr.n_sites_out}</td>"
            f"<td>{cached}</td>"
            "</tr>"
        )

    outdir_line = f"<br><b>Output:</b> {outdir}" if outdir is not None else ""

    return (
        f"{_NOTEBOOK_CSS}"
        '<div class="pycsamt-result">'
        f'<div class="title">PipelineResult &#39;{pipeline_name}&#39; '
        f'&mdash; <span class="{status_cls}">{status_txt}</span></div>'
        f"<b>Sites:</b> {n_sites_in}&rarr;{n_sites_out} &nbsp;|&nbsp; "
        f"<b>Steps:</b> {len(step_results)} ({n_ok} ok, {n_err} err) "
        f"&nbsp;|&nbsp; <b>Time:</b> {elapsed_sec:.2f}s"
        f"{outdir_line}"
        "<table><tr><th>#</th><th>Step</th><th>Code</th><th>Status</th>"
        "<th>Time</th><th>Sites</th><th>Cached</th></tr>"
        + "".join(rows)
        + "</table></div>"
    )


__all__ = ["make_text_report", "make_html_report", "make_notebook_html"]
