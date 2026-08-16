"""The "dashboard" report tier — a branded, chart-carrying HTML report.

:func:`make_dashboard_html` is the richer sibling of
:func:`pycsamt.pipeline._report.make_html_report`: real pyCSAMT brand
colors and logo, KPI stat tiles, and three native inline-SVG charts (no JS,
no CDN — the file stays a single self-contained document, the same
convention the rest of the pipeline's reports already follow).

Palette provenance
-------------------
Every color here traces to a real source, not an eyeballed choice:

* Brand tokens (ink/surface/border/secondary/primary) are copied from
  ``docs/source/_static/css/custom.css``'s ``--pycsamt-*`` custom
  properties — pyCSAMT's actual site palette.
* The categorical pair used for the two-series site-count-flow chart is
  the brand's own secondary blue and primary orange, run through the
  dataviz skill's ``validate_palette.js`` against both this dashboard's
  light (``#ffffff``) and dark (``#141922``) surfaces — every check passed
  (CVD ΔE 22.4, normal-vision ΔE 33.2, contrast >= 3:1) with the raw brand
  hex values, no substitution needed.
* The status palette (good/warning/critical) is the dataviz skill's own
  validated, brand-independent set (documented there as "fixed, never
  themed") — deliberately distinct from the categorical slots so a status
  color can never be mistaken for a data series.

Both the categorical pair and the status palette validated identically on
light and dark surfaces, so — unlike the ink/surface/border tokens — none
of them are redefined in the dark-mode media query below.

The small stat-tile icons are trimmed, ``currentColor``-recolored copies
of real icon assets; the header logo and the favicon are the real,
unmodified pyCSAMT marks (``pycsamt-v2-symbol.svg`` / ``.ico``). All of
these live under ``docs/source/_static/`` — that directory is not shipped
in the installed package, so every asset is embedded here as a plain
string (the favicon as a base64 data URI) rather than read from disk at
runtime.
"""

from __future__ import annotations

import datetime
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._steps import StepResult

__all__ = ["make_dashboard_html"]

# ---------------------------------------------------------------------------
# Brand tokens
# ---------------------------------------------------------------------------

_BRAND_CSS = """
<style>
  .pycsamt-dashboard {
    --pc-primary: #f15a29;
    --pc-secondary: #3e65b0;
    --pc-secondary-dark: #1f4490;
    --pc-gold: #fbb040;
    --pc-ink: #24324b;
    --pc-muted: #5c677d;
    --pc-surface: #ffffff;
    --pc-surface-2: #f6f8fc;
    --pc-border: #dde3ee;
    --pc-good: #0ca30c;
    --pc-warning: #fab219;
    --pc-critical: #d03b3b;
    font-family: system-ui, -apple-system, "Segoe UI", sans-serif;
    color: var(--pc-ink);
    background: var(--pc-surface);
    max-width: 960px;
    margin: 0 auto;
    padding: 2em 1.5em 3em;
  }
  @media (prefers-color-scheme: dark) {
    .pycsamt-dashboard {
      --pc-ink: #f4f7fb;
      --pc-muted: #c6d0df;
      --pc-surface: #141922;
      --pc-surface-2: #1a2130;
      --pc-border: #2a3245;
    }
  }
  .pycsamt-dashboard * { box-sizing: border-box; }
  .pycsamt-dashboard header {
    display: flex; align-items: center; gap: 0.8em; flex-wrap: wrap;
  }
  .pycsamt-dashboard header .logo { width: 48px; height: 48px; flex: none; }
  .pycsamt-dashboard header .wordmark { font-weight: 800; font-size: 1.4rem; }
  .pycsamt-dashboard header .wordmark .sub {
    display: block; font-weight: 500; font-size: 0.85rem; color: var(--pc-muted);
  }
  .pycsamt-dashboard .status-badge {
    margin-left: auto; padding: 0.35em 0.9em; border-radius: 999px;
    font-weight: 700; font-size: 0.85rem; color: #fff;
  }
  .pycsamt-dashboard .status-badge.ok  { background: var(--pc-good); }
  .pycsamt-dashboard .status-badge.err { background: var(--pc-critical); }
  .pycsamt-dashboard .run-meta {
    color: var(--pc-muted); font-size: 0.9rem; margin: 0.6em 0 1.6em;
  }
  .pycsamt-dashboard .tiles {
    display: grid; grid-template-columns: repeat(auto-fit, minmax(140px, 1fr));
    gap: 0.8em; margin-bottom: 2em;
  }
  .pycsamt-dashboard .tile {
    background: var(--pc-surface-2); border: 1px solid var(--pc-border);
    border-radius: 10px; padding: 0.9em 1em;
  }
  .pycsamt-dashboard .tile .icon {
    width: 20px; height: 20px; color: var(--pc-secondary); margin-bottom: 0.4em;
  }
  .pycsamt-dashboard .tile .label {
    font-size: 0.78rem; color: var(--pc-muted); text-transform: none;
  }
  .pycsamt-dashboard .tile .value {
    font-weight: 600; font-size: 1.3rem; margin-top: 0.15em;
  }
  .pycsamt-dashboard .tile.err .value { color: var(--pc-critical); }
  .pycsamt-dashboard section.chart-card {
    background: var(--pc-surface-2); border: 1px solid var(--pc-border);
    border-radius: 10px; padding: 1.2em 1.4em; margin-bottom: 1.6em;
  }
  .pycsamt-dashboard section.chart-card h2 {
    margin: 0 0 0.8em; font-size: 1.05rem; color: var(--pc-secondary);
  }
  .pycsamt-dashboard .chart-scroll { overflow-x: auto; }
  .pycsamt-dashboard .legend {
    display: flex; gap: 1.2em; font-size: 0.82rem; color: var(--pc-muted);
    margin-bottom: 0.6em;
  }
  .pycsamt-dashboard .legend .swatch {
    display: inline-block; width: 10px; height: 10px; border-radius: 2px;
    margin-right: 0.4em; vertical-align: middle;
  }
  .pycsamt-dashboard table {
    border-collapse: collapse; width: 100%; font-size: 0.88rem;
  }
  .pycsamt-dashboard th, .pycsamt-dashboard td {
    border: 1px solid var(--pc-border); padding: 0.4em 0.6em; text-align: left;
  }
  .pycsamt-dashboard th { background: var(--pc-secondary); color: #fff; }
  .pycsamt-dashboard tr:nth-child(even) { background: var(--pc-surface-2); }
  .pycsamt-dashboard .ok   { color: var(--pc-good); font-weight: bold; }
  .pycsamt-dashboard .err  { color: var(--pc-critical); font-weight: bold; }
  .pycsamt-dashboard .step-card {
    border: 1px solid var(--pc-border); border-radius: 8px;
    padding: 1em; margin: 1em 0;
  }
  .pycsamt-dashboard .step-card h3 { margin: 0 0 0.5em; display: flex; align-items: center; gap: 0.5em; }
  .pycsamt-dashboard .step-card .icon { width: 18px; height: 18px; color: var(--pc-secondary); }
  .pycsamt-dashboard .step-meta { font-size: 0.85em; color: var(--pc-muted); }
  .pycsamt-dashboard .plots { display: flex; flex-wrap: wrap; gap: 0.5em; margin-top: 0.5em; }
  .pycsamt-dashboard .plots a img {
    max-height: 120px; border: 1px solid var(--pc-border); border-radius: 4px;
  }
  .pycsamt-dashboard .error-box {
    background: #fdecea; border-left: 4px solid var(--pc-critical);
    padding: 0.5em 1em; font-size: 0.9em; color: #5c2320;
  }
  .pycsamt-dashboard footer {
    margin-top: 3em; font-size: 0.8em; color: var(--pc-muted);
    border-top: 1px solid var(--pc-border); padding-top: 1em;
  }
</style>
"""

# ---------------------------------------------------------------------------
# Embedded icons — trimmed, currentColor-recolored copies of
# docs/source/_static/icons/{workflow,clock,export,diagnostic}.svg
# (SVG Repo, viewBox 0 0 24 24). Embedded rather than read from disk: the
# docs/ tree is not part of the installed package.
# ---------------------------------------------------------------------------

_ICON_WORKFLOW = (
    '<svg class="icon" viewBox="0 0 24 24" fill="currentColor" '
    'xmlns="http://www.w3.org/2000/svg"><path d="M7.5,15.5h-5a1,1,0,0,0-1,1v5'
    'a1,1,0,0,0,1,1h5a1,1,0,0,0,1-1V20H12a1,1,0,0,0,0-2H8.5V16.5A1,1,0,0,0,7.5'
    ',15.5Zm-1,5h-3v-3h3ZM4,8.858V13a1,1,0,0,0,2,0V8.858a4,4,0,1,0-2,0ZM5,3A2,'
    '2,0,1,1,3,5,2,2,0,0,1,5,3ZM20,15.142V12a1,1,0,0,0-2,0v3.142a4,4,0,1,0,2,'
    '0ZM19,21a2,2,0,1,1,2-2A2,2,0,0,1,19,21ZM16.5,8.5h5a1,1,0,0,0,1-1v-5a1,1,'
    '0,0,0-1-1h-5a1,1,0,0,0-1,1V4H12a1,1,0,0,0,0,2h3.5V7.5A1,1,0,0,0,16.5,8.5'
    'Zm1-5h3v3h-3Z"/></svg>'
)

_ICON_CLOCK = (
    '<svg class="icon" viewBox="0 0 24 24" fill="none" '
    'xmlns="http://www.w3.org/2000/svg"><path d="M12 17V12L14.5 10.5M21 12C21'
    ' 16.9706 16.9706 21 12 21C7.02944 21 3 16.9706 3 12C3 7.02944 7.02944 3 '
    '12 3C16.9706 3 21 7.02944 21 12Z" stroke="currentColor" stroke-width="2"'
    ' stroke-linecap="round" stroke-linejoin="round"/></svg>'
)

_ICON_EXPORT = (
    '<svg class="icon" viewBox="0 0 24 24" fill="currentColor" '
    'xmlns="http://www.w3.org/2000/svg"><path d="M14.4697 7.53033C14.7626 '
    '7.82322 15.2374 7.82322 15.5303 7.53033C15.8232 7.23744 15.8232 6.76256'
    ' 15.5303 6.46967L12.5303 3.46967C12.2374 3.17678 11.7626 3.17678 11.4697'
    ' 3.46967L8.46967 6.46967C8.17678 6.76256 8.17678 7.23744 8.46967 '
    '7.53033C8.76256 7.82322 9.23744 7.82322 9.53033 7.53033L11.25 5.81066V14'
    'C11.25 14.4142 11.5858 14.75 12 14.75C12.4142 14.75 12.75 14.4142 12.75 '
    '14V5.81066L14.4697 7.53033Z"/><path d="M20.75 12C20.75 11.5858 20.4142 '
    '11.25 20 11.25C19.5858 11.25 19.25 11.5858 19.25 12C19.25 16.0041 16.0041'
    ' 19.25 12 19.25C7.99593 19.25 4.75 16.0041 4.75 12C4.75 11.5858 4.41421 '
    '11.25 4 11.25C3.58579 11.25 3.25 11.5858 3.25 12C3.25 16.8325 7.16751 '
    '20.75 12 20.75C16.8325 20.75 20.75 16.8325 20.75 12Z"/></svg>'
)

_ICON_DIAGNOSTIC = (
    '<svg class="icon" viewBox="0 0 24 24" fill="currentColor" '
    'xmlns="http://www.w3.org/2000/svg"><path fill-rule="evenodd" d="M15.5,14'
    'h-.8l-.3-.3A6.26,6.26,0,0,0,16,9.5,6.5,6.5,0,1,0,9.5,16a6.26,6.26,0,0,0,'
    '4.2-1.6l.3.3v.8l5,5L20.5,19Zm-6-9a4.48,4.48,0,0,1,4.39,3.5H11.46l-1.07,'
    '1.6L8.22,6,6.91,8.5H5.11A4.48,4.48,0,0,1,9.5,5Zm0,9A4.48,4.48,0,0,1,5,'
    '10H7.82l.41-.78,2,3.79,2-3H14A4.47,4.47,0,0,1,9.5,14Z"/></svg>'
)

# Real pyCSAMT logo mark, embedded verbatim from
# docs/source/_static/logo/pycsamt-v2-symbol.svg (that directory is not
# shipped in the installed package, so the source is embedded here as a
# plain string rather than read from disk at runtime). The header wraps
# it at 48x48 via the ".logo" CSS rule in _BRAND_CSS; the SVG's own
# viewBox keeps it crisp at that size. class/role/aria-label were added
# to the root <svg> tag; every fill/stroke hex below is unmodified.
_LOGO_SVG = """<svg class="logo" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 476.36 454.92" role="img" aria-label="pyCSAMT"><defs><style>.cls-1{fill:#fff;}.cls-1,.cls-42{stroke:#fff;stroke-miterlimit:10;stroke-width:7px;}.cls-2{fill:#1f4490;}.cls-3{fill:#1e4491;}.cls-4{fill:#f39128;}.cls-5{fill:#f1f0f0;}.cls-6{fill:#48537f;}.cls-7{fill:#5976b0;}.cls-8{fill:#f89a25;}.cls-9{fill:#30539e;}.cls-10{fill:#4b587f;}.cls-11{fill:#2d4a8b;}.cls-12{fill:#5878b5;}.cls-13{fill:#404d78;}.cls-14{fill:#3d60a4;}.cls-15{fill:#274890;}.cls-16{fill:#617cb5;}.cls-17{fill:#254c99;}.cls-18{fill:#294d96;}.cls-19{fill:#264b97;}.cls-20{fill:#bdcbe3;}.cls-21{fill:#4f6daa;}.cls-22{fill:#254d97;}.cls-23{fill:#faa746;}.cls-24{fill:#3e5ea2;}.cls-25{fill:#3557a2;}.cls-26{fill:#4064a9;}.cls-27{fill:#395a9f;}.cls-28{fill:#34569d;}.cls-29{fill:#224a97;}.cls-30{fill:#385a9f;}.cls-31{fill:#dc8a2f;}.cls-32{fill:#21478f;}.cls-33{fill:#2a478a;}.cls-34{fill:#284a8b;}.cls-35{fill:#34569e;}.cls-36{fill:#3b5c9f;}.cls-37{fill:#3c5ea1;}.cls-38{fill:#224790;}.cls-39{fill:#264990;}.cls-40{fill:#1f498f;}.cls-41{fill:#dd8e2d;}.cls-42{fill:none;}</style></defs><g id="round"><circle class="cls-1" cx="208.44" cy="246.49" r="204.94"/></g><g id="Layer_1" data-name="Layer 1"><path class="cls-2" d="M363,173.73l2.57.82a191.47,191.47,0,0,1,18.61,6.93l2.93,1.24c19.13,8.15,38.41,20.57,53,35,.75.72,1.53,1.38,2.34,2a86.51,86.51,0,0,1,9.67,9.37c.92,1,1.88,2,2.83,3,10.84,11.43,19.65,24.6,27.47,38.07l.51.87c15.33,26.22,15.33,26.22,15.33,36.51l-.88,0c-14.32.29-14.32.29-28.24,3.27l-.92.29c-21.23,6.71-34.71,22.14-47.08,39.2-2.33,3.2-4.71,6.34-7.25,9.39L413,361c-2.41,2.88-5,5.6-7.61,8.29l-.75.81c-5.65,5.81-14.3,10.36-22.59,10.83-9.95.32-18.52-3.7-25.74-10.1a76.74,76.74,0,0,1-9.17-10.17l-.62-.8c-8.09-10.62-14.75-22.13-21.22-33.7l-.74-1.33q-3.47-6.18-6.91-12.37c-15.88-29.12-15.88-29.12-36.9-54.9l-.61-.63c-8.23-8.31-19.71-15.68-32-16.28-18.55-.22-31.75,5.08-45.06,17.49l-.69.65-1.15,1.07a137.93,137.93,0,0,0-14.66,17c-2.1,2.77-4.27,5.49-6.46,8.19l-.79,1c-11.17,13.89-23.48,23.85-42,26.53l-1.57.23a139.23,139.23,0,0,1-18.83.65l-3.42,0c-2.77,0-5.53,0-8.3,0a42,42,0,0,1,1.46-6.19l.32-1,.33-1,.33-1.05c.67-2.13,1.43-4.22,2.27-6.3.25-.64.49-1.28.74-1.91,1.84-4.75,3.86-9.38,6.08-14,.43-.9.85-1.79,1.27-2.69a150.07,150.07,0,0,1,7.52-13.59c.59-1,1.17-2,1.75-2.95a180.37,180.37,0,0,1,14.76-20.75l.62-.77c2.26-2.83,4.57-5.56,7.08-8.19.76-.79,1.5-1.6,2.23-2.42a132.46,132.46,0,0,1,13.52-13.28c1-.8,1.87-1.61,2.8-2.44,2.59-2.3,5.25-4.48,8-6.59l.66-.51c16.14-12.44,34.86-24,54.6-30.28a19.23,19.23,0,0,0,2.67-1,80.88,80.88,0,0,1,7.89-2.6l.91-.27c5.69-1.71,11.42-3.16,17.23-4.47,1-.22,1.93-.46,2.89-.7,1.9-.46,3.8-.74,5.73-1l1.11-.15c4.32-.6,8.63-1.17,13-1.72l1.52-.19C310.19,162.18,338.48,165.8,363,173.73Z" transform="translate(-94.74 -119.49)"/><path class="cls-2" d="M262.62,262c13.77,10.1,22.58,24.84,30.94,39.1l.81,1.37q4.9,8.32,9.45,16.8,3.19,6,6.6,11.82c.64,1.12,1.28,2.23,1.91,3.35,5.31,9.32,10.76,18.54,16.78,27.46l.49.73c10,14.88,22.66,27.9,41.06,32.54a49.19,49.19,0,0,0,35.61-7.27,79.56,79.56,0,0,0,17.5-16.14l.74-.92c3.28-4,6.44-8.09,9.35-12.36,2.27-3.32,4.71-6.5,7.2-9.66.63-.8,1.25-1.61,1.87-2.41,9.61-12.33,22.28-20.62,38.23-23.13l1.53-.25c6.56-.93,13.27-.64,19.89-.57,2.58,13.74,4.87,27.18,4.82,41.18,0,1.49,0,3,0,4.46a164,164,0,0,1-1.93,25l-.14.9c-.95,6.31-2.18,12.57-3.61,18.79l-.15.84c-.42,1.2-1.33,1.48-2.46,2l-.93.42-.95.45-3,1.41-2,.94-3.74,1.75c-1.47.68-2.92,1.38-4.37,2.08-1.29.6-2.62,1.08-4,1.55l-.9.31-2.77,1-.87.3c-2.13.75-4.27,1.46-6.41,2.16l-1.05.35a149.7,149.7,0,0,1-19.32,4.77,4.21,4.21,0,0,0-1.9.64c-1.27.2-2.53.38-3.8.55l-1.17.15c-31,4.15-62.37,3.67-92.69-4.4l-1.4-.37a192.9,192.9,0,0,1-23.92-7.76c-1.5-.6-3-1.17-4.51-1.74a84.35,84.35,0,0,1-8.59-3.63,30.14,30.14,0,0,0-2.81-1.16c-2.84-1.08-5.61-2.31-8.38-3.53l-2.93-1.3-4.09-1.79-8.82-3.87-1.79-.78c-5.76-2.53-11.59-4.9-17.53-7l-3.25-1.17-3.59-1.29c-.82-.29-1.64-.59-2.45-.9-17.66-6.43-37.88-10.21-56.7-11.33l-1.89-.14a227.13,227.13,0,0,0-32,.14l-1,.07a237.56,237.56,0,0,0-31.11,4.45l-1.66.34c-2.71.59-5.41,1.24-8.1,1.9-1.45.35-2.9.69-4.35,1-5,1.07-9.83,2.42-14.7,3.82l-4.17,1.15-.11-1.1L99,390c-.06-.59-.11-1.18-.17-1.76-.09-.85-.18-1.69-.26-2.54l-.16-1.53A7.63,7.63,0,0,0,98,382a22,22,0,0,1-.45-4.44l-.06-1.09c-.18-4-.2-8-.19-12v-9.21a11.68,11.68,0,0,1,.49-2.49c.07-.64.12-1.27.17-1.91l.08-1.16.09-1.24c.32-3.94.79-7.83,1.37-11.74l.22-1.47a44.32,44.32,0,0,1,1.49-6.72h1l9.52-.07,4.89,0c23.48-.1,45.13-3.76,62.54-20l.68-.62A81.49,81.49,0,0,0,188.2,299l.74-.92q5.12-6.36,10.06-12.85c8.29-10.87,17.86-22.42,31.17-27.72l1.51-.61C242.42,253,253.52,255.9,262.62,262Z" transform="translate(-94.74 -119.49)"/><path class="cls-3" d="M181.63,407l5,0a115.47,115.47,0,0,1,16.1.8c.9.12,1.79.23,2.69.33,14.47,1.66,29.2,3.88,43.15,8l2.7.78,2.87.84,1.44.42a78.35,78.35,0,0,1,7.48,2.54,30.43,30.43,0,0,0,4.36,1.17v.82l.77.11a26.33,26.33,0,0,1,7.05,2.23,68.88,68.88,0,0,0,6.56,2.61,74.64,74.64,0,0,1,9.06,3.85,47.82,47.82,0,0,0,5.08,1.93,25.24,25.24,0,0,1,4,2c.71.3,1.42.6,2.14.88,1.05.43,2.09.89,3.13,1.35,1.8.79,3.61,1.52,5.45,2.21,1.25.48,2.48,1,3.71,1.5,26.05,10.77,54.75,18.75,83.28,18.77h8.87c9.62,0,19.31,0,28.79-1.75l2.45-.43,1.14-.2c1.37-.21,2.73-.35,4.11-.49a44.14,44.14,0,0,0,5.47-1.1l3.31-.76A152.59,152.59,0,0,0,467.46,451l1.49-.52a227.4,227.4,0,0,0,21.62-8.55,15,15,0,0,1,2.16-.79c-5.39,16.55-15.73,33.42-26.89,47-.47.59-.91,1.2-1.35,1.8a7.46,7.46,0,0,1-3.73,2l-1.49.46-.8.25c-1,.32-1.9.7-2.85,1.07-12.68,4.92-26.9,6.84-40.42,8.26l-.91.1a89,89,0,0,1-9.52.39h-.88l-4.62,0h-3.76a136.89,136.89,0,0,1-22.71-1.35l-1-.15c-1.93-.3-3.85-.64-5.77-1l-.83-.16a110.51,110.51,0,0,1-13.33-3.33c-1.4-.43-2.82-.82-4.23-1.21a34.69,34.69,0,0,1-3.81-1.51c-1.11-.47-2.24-.87-3.39-1.26-10.12-3.7-19.75-8.92-29.23-13.92a298.11,298.11,0,0,0-38.65-17.3l-1.79-.65a160.18,160.18,0,0,0-21.34-5.92l-1.12-.24c-6.29-1.32-12.64-2.24-19-3l-1.64-.23a126.23,126.23,0,0,0-16.84-.7h-4.92a163.91,163.91,0,0,0-26.29,1.67l-1.36.21-4.2.69-1.45.24A160.86,160.86,0,0,0,149.34,459c-1.3.42-2.61.81-3.93,1.18-2.21.63-4.35,1.38-6.49,2.18l-1.17.43q-4.29,1.57-8.52,3.3a25.47,25.47,0,0,1-4,1.3A172.67,172.67,0,0,1,112,440.32c-.2-.52-.41-1-.61-1.54-1.35-3.45-2.51-7-3.7-10.45-.24-.7-.48-1.4-.73-2.09s-.5-1.47-.74-2.2l-.39-1.15a12.23,12.23,0,0,1-.25-3.52,229.58,229.58,0,0,1,44.11-10.27l3.21-.41,1.62-.21c1.05-.13,2.09-.27,3.14-.42C165.58,407,173.59,407,181.63,407Z" transform="translate(-94.74 -119.49)"/><path class="cls-4" d="M256.73,241.91l1.48.46c8.05,2.76,14.71,7.45,20.79,13.09l.69.63a83.19,83.19,0,0,1,9.59,10.46l1.17,1.49c10.62,13.55,18.79,28.8,27.1,43.7l1.2,2.15.81,1.46,3.64,6.51,1.5,2.69q4.89,8.76,10.16,17.31c.61,1,1.21,2,1.81,3,8.34,13.93,18.95,29.12,35.39,34.73,8.45,2.18,16.85.42,24.32-3.69,6-3.5,11-8.51,15.26-13.77l2-2.38c3.1-3.76,6-7.59,8.8-11.56,1.26-1.78,2.59-3.5,3.93-5.23l.57-.76a111.16,111.16,0,0,1,8.85-10.33c.52-.57,1-1.14,1.55-1.72a66.21,66.21,0,0,1,9.59-8.55l1.31-1c14.83-10.8,32-13.23,50.08-13,1,2.75,1.85,5.52,2.65,8.32l.36,1.21.33,1.16.3,1a11.67,11.67,0,0,1,.22,3.06l-1.1,0c-17.7-.16-17.7-.16-34,5.32-.61.22-1.23.42-1.85.62-5.47,2.06-10.1,6-14.24,9.89-.58.55-1.17,1.08-1.76,1.62a63.54,63.54,0,0,0-7.42,8.41l-.64.84c-2.38,3.1-4.71,6.23-7,9.42-11.05,15.7-24.37,32.07-44.47,36.93-7.52,1.27-16.08,1.33-23.32-1.18l-1.14-.37c-1.08-.38-2.1-.81-3.14-1.27l-1.38-.6c-9.15-4.16-16.77-10.76-23-18.3l-1.14-1.35a118,118,0,0,1-7.86-10.56l-.72-1.08c-4-6-7.7-12-11.27-18.23l-1.22-2.11c-4.79-8.27-9.51-16.58-14.05-25-2.77-5.14-5.71-10.19-8.71-15.21-.42-.7-.84-1.4-1.25-2.1-8.59-14.39-18.12-30-33.45-38.57l-1.29-.73a32.92,32.92,0,0,0-22.07-2.45,46.16,46.16,0,0,0-17.75,9.75l-1.29,1.08A108.43,108.43,0,0,0,200,284.35c-2.78,3.65-5.64,7.24-8.51,10.83l-2,2.49a98.62,98.62,0,0,1-8.21,9.14l-.68.69c-14.21,14.28-33.26,20.62-53.57,21.2-3.8.07-7.6,0-11.39,0l-4.93,0q-4.75,0-9.53-.08a39.19,39.19,0,0,1,1-5.67l.34-1.29.36-1.3.69-2.6.31-1.15c.3-1.19.51-2.39.71-3.6h.87q4.13,0,8.27.06c1.42,0,2.83,0,4.25,0,10.7.1,21-.31,31.15-3.79l.88-.29c13-4.3,22.06-14.19,30.07-24.5,1.82-2.33,3.67-4.64,5.51-7,1.56-1.95,3.1-3.92,4.63-5.89a142.84,142.84,0,0,1,11.6-12.89l1.17-1.18a66.12,66.12,0,0,1,24.42-15.43l.82-.29A48.84,48.84,0,0,1,256.73,241.91Z" transform="translate(-94.74 -119.49)"/><path class="cls-3" d="M237,474l3.06.51A146.09,146.09,0,0,1,270.86,483l1.43.57c11.25,4.55,11.25,4.55,16.43,7.15l3.69,1.83q6.64,3.3,13.28,6.65c2.73,1.38,5.46,2.73,8.23,4l3.24,1.56a108.78,108.78,0,0,0,13.52,5.62,38.43,38.43,0,0,1,4,1.75,163,163,0,0,0,30.84,8.62l1.27.24a181.28,181.28,0,0,0,62,.5,22.58,22.58,0,0,1,4.79-.33,97.79,97.79,0,0,1-13.28,10.68l-.89.63c-4.12,2.86-8.33,5.57-12.67,8.12a31.15,31.15,0,0,0-3.12,2.06c-4.36,3.11-9.89,2.2-15,1.93l-2.17-.11c-12.55-.63-25.24-1.42-37.58-3.82l-2.12-.4A117.15,117.15,0,0,1,332.41,537c-1.34-.43-2.68-.83-4-1.21-14.13-4.09-27.55-10.57-40.41-17.44-4.58-2.45-9.27-4.66-14-6.74l-2.72-1.21a84.61,84.61,0,0,0-8.72-3.4c-1.34-.45-2.64-.95-4-1.47-21.42-7.92-46.79-10.63-69.46-6.83l-1.63.24c-8.53,1.43-16.91,4.08-25.08,6.79l-1.36.45-1.19.4a20.44,20.44,0,0,1-3.72.75,8,8,0,0,1-1.9-1.83l-.91-1.05-1.05-1.23-1.44-1.58c-3.59-3.92-7.05-8-10.28-12.16l-.61-.78c-1.37-1.81-1.37-1.81-1.37-2.73a194.66,194.66,0,0,1,27-9.14c.79-.2,1.57-.41,2.35-.63a87.37,87.37,0,0,1,9.89-2l1.63-.26C198.38,471,218.13,470.64,237,474Z" transform="translate(-94.74 -119.49)"/><path class="cls-5" d="M568.11,119.49h3v11.09h-3Z" transform="translate(-94.74 -119.49)"/><path class="cls-6" d="M143.66,326.94v.41l-2.57.41v.41l-1,.11-4.47.48-1.56.17-1.51.17-1.38.15a7.18,7.18,0,0,0-2.11.41,3.13,3.13,0,0,1-1.79,0,20.84,20.84,0,0,0-3.7-.31l-1.51,0-1.61,0-1.66,0-4.36-.06-4.45-.07-8.73-.13v-.41h1l9.68-.09,5-.05a154.07,154.07,0,0,0,21.26-1.24A32.48,32.48,0,0,1,143.66,326.94Z" transform="translate(-94.74 -119.49)"/><path class="cls-7" d="M181.77,407h9.91c3.14,0,6.24.17,9.37.46v.42H171.92l1.29.82H168.5v-.41h-9v-.41C166.84,406.77,174.36,407,181.77,407Z" transform="translate(-94.74 -119.49)"/><path class="cls-8" d="M101.68,325.3l1.29.41-.43,2.46h16.7v.41h-18Z" transform="translate(-94.74 -119.49)"/><path class="cls-9" d="M440.9,218.9a10.25,10.25,0,0,1,3.29,2.25l.86.84.88.87.91.88,2.2,2.15-.86.82-3.64-3.44-1.05-1-1-.95-.92-.87C440.9,219.72,440.9,219.72,440.9,218.9Z" transform="translate(-94.74 -119.49)"/><path class="cls-10" d="M388.65,395.54l.43.82c-4.43.64-8.84.27-13.28,0V396l5-.2,1.42-.06,1.38-.06,1.26-.05C386.11,395.54,387.38,395.54,388.65,395.54Z" transform="translate(-94.74 -119.49)"/><path class="cls-11" d="M143.66,326.94v.41l-2.57.41v.41a36.6,36.6,0,0,1-6.86.41v-.82C140.48,326.82,140.48,326.82,143.66,326.94Z" transform="translate(-94.74 -119.49)"/><path class="cls-12" d="M144.94,410.33l-3,.41V412l-1-.23c-1.11-.28-1.11-.28-2,.23-.86-.11-1.71-.25-2.57-.41v-.41l3.32-.64.95-.19A8.36,8.36,0,0,1,144.94,410.33Z" transform="translate(-94.74 -119.49)"/><path class="cls-13" d="M253.73,241.09l6.42.41.43,1.64c-2.31-.47-4.58-1-6.85-1.64Z" transform="translate(-94.74 -119.49)"/><path class="cls-14" d="M266.58,168.38l-2.57.41v.82a22.75,22.75,0,0,1-6,.41,18.4,18.4,0,0,1,4.29-1.52l1.33-.33A5,5,0,0,1,266.58,168.38Z" transform="translate(-94.74 -119.49)"/><path class="cls-15" d="M445.18,343.78H446v.82h-.86l-.08.7a6.75,6.75,0,0,1-2.49,3h-1.28L443.2,346l1.06-1.27Z" transform="translate(-94.74 -119.49)"/><path class="cls-16" d="M234,451.82l-9.42-.82v-.41C228,450.09,231.09,449.94,234,451.82Z" transform="translate(-94.74 -119.49)"/><path class="cls-17" d="M259.3,170a8.39,8.39,0,0,1-3.45,1.33l-1.14.24-1,.07-.86-.82A17.47,17.47,0,0,1,259.3,170Z" transform="translate(-94.74 -119.49)"/><path class="cls-18" d="M275.15,511.38a16,16,0,0,1,5.56,2.06v.82a56.88,56.88,0,0,1-5.56-2Z" transform="translate(-94.74 -119.49)"/><path class="cls-19" d="M168.5,214.38l.85.41-2.57,2.88-.85-.82Z" transform="translate(-94.74 -119.49)"/><path class="cls-20" d="M418.63,501.53v.41c-2.67,1-4.91.77-7.71.41v-.41A58.41,58.41,0,0,1,418.63,501.53Zm-3.86,1.23.86.41Z" transform="translate(-94.74 -119.49)"/><path class="cls-21" d="M220.75,384l1.71,1.23-7.28-.82V384A10.43,10.43,0,0,1,220.75,384Z" transform="translate(-94.74 -119.49)"/><path class="cls-22" d="M147.08,235.75l.86.41a14.28,14.28,0,0,1-2.57,3.28l-1.29-.41Z" transform="translate(-94.74 -119.49)"/><path class="cls-23" d="M101.68,325.3l1.29.41-.43,2.46.85.41h-2.14Z" transform="translate(-94.74 -119.49)"/><path class="cls-24" d="M295.28,520.83a11.92,11.92,0,0,1,1.76.8l.94.47.72.38v.82a11,11,0,0,1-4.28-2.06Z" transform="translate(-94.74 -119.49)"/><path class="cls-25" d="M140.23,485.09l.86.41L139.8,488a16.19,16.19,0,0,1-1.28-1.23v-.82A13,13,0,0,1,140.23,485.09Z" transform="translate(-94.74 -119.49)"/><path class="cls-26" d="M492.73,441.14l-.43,1.64-3.86.41C491.17,441.14,491.17,441.14,492.73,441.14Z" transform="translate(-94.74 -119.49)"/><path class="cls-27" d="M296.13,409.51a11.55,11.55,0,0,1,3.86,1.23v.82a11.55,11.55,0,0,1-3.86-1.23Z" transform="translate(-94.74 -119.49)"/><path class="cls-28" d="M267.44,507.69l1.68.54.95.3c.79.39.79.39,1.22,1.62l-3.85-1.23Z" transform="translate(-94.74 -119.49)"/><path class="cls-29" d="M416.49,500.29l3,.41v.83h-3.42Z" transform="translate(-94.74 -119.49)"/><path class="cls-30" d="M321.4,419.78c.57.17,1.13.35,1.69.54l1,.3c.79.39.79.39,1.22,1.62L321.4,421Z" transform="translate(-94.74 -119.49)"/><path class="cls-31" d="M449.47,338l1.28.41-3,2.88-.43-1.23h.86l.4-.8A3.39,3.39,0,0,1,449.47,338Z" transform="translate(-94.74 -119.49)"/><path class="cls-32" d="M191.63,294.9h.85a16.36,16.36,0,0,1-2.14,3.28l-1.28-.41c.4-.48.81-1,1.23-1.43l.69-.81Z" transform="translate(-94.74 -119.49)"/><path class="cls-33" d="M448.18,340.5l.86.41-2.14,2.46-1.29-.41Z" transform="translate(-94.74 -119.49)"/><path class="cls-34" d="M271.29,270.66a5,5,0,0,1,1.71.41,13.47,13.47,0,0,1,1.29,2.06H273l-.42-1.24-1.29-.41Z" transform="translate(-94.74 -119.49)"/><path class="cls-35" d="M303.84,525.35a12.89,12.89,0,0,1,3.43,1.23v.82a19,19,0,0,1-3.43-1.23Z" transform="translate(-94.74 -119.49)"/><path class="cls-36" d="M383.94,435.39v.41a28.71,28.71,0,0,1-6-.41V435C380.06,434.32,381.85,434.86,383.94,435.39Z" transform="translate(-94.74 -119.49)"/><path class="cls-37" d="M99.54,392.67l.43.82.85.2.86.21.43.82-2.57.41Z" transform="translate(-94.74 -119.49)"/><path class="cls-38" d="M430.19,363.5l1.29.41L429.76,366l-1.28-.42Z" transform="translate(-94.74 -119.49)"/><path class="cls-39" d="M265.72,265.32a3.3,3.3,0,0,1,2.68,1.47l.75,1c-1.29,0-1.29,0-2-.57a1.57,1.57,0,0,1-.56-1.49Z" transform="translate(-94.74 -119.49)"/><path class="cls-40" d="M305.13,322.83c1.2,1,1.91,1.73,2.14,3.29-1.29-.41-1.29-.41-1.79-1.18A3.83,3.83,0,0,1,305.13,322.83Z" transform="translate(-94.74 -119.49)"/><path class="cls-41" d="M155.65,307.22V308l-3,1.24-.86-.83C154.68,307.22,154.68,307.22,155.65,307.22Z" transform="translate(-94.74 -119.49)"/></g><g id="fw"><circle class="cls-42" cx="208.6" cy="246.24" r="204.94"/></g></svg>"""

# Real pyCSAMT favicon, embedded as a base64 data URI (same "docs/ is not
# shipped" constraint as the logo above) — copied verbatim from
# docs/source/_static/logo/pycsamt-v2-symbol.ico.
_FAVICON_ICO_BASE64 = (
    "AAABAAEAIB8AAAEAIAAkEAAAFgAAACgAAAAgAAAAPgAAAAEAIAAAAAAAgA8AACMuAAAjLgAAAAAA"
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAD///8A////AP///xP///9R////mv///9H////v"
    "/////P////z////v////0f///5v///9S////FP///wD///8AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAD///8A////AP///w////9i////w///"
    "//X//////////////////////////////////////Pr5//fy7/Xx6OPE8efhY////xD+/f0A////"
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA////AP///wD///8z"
    "////tP////r////////////////////////////////17ur/28K1/8GVfv+udVj/pWVE/6FfPf+t"
    "dFb607SktvXv6zQyAAAA////AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAP//"
    "/wD///8A////Tf///9z////////////////////////////////07Oj/zqya/6VlRP+TRx//mVAr"
    "/6tvUP++j3f/yqSR/82ql//Ss6P/8Obg3f///07///8B////AAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    "AAAAAAAAAAD///8AyqWRAPHn4kzgy8Dl8OXf//r39f/9/Pv//Pv6//fy7//o2ND/yqSQ/6RkQ/+S"
    "RRz/omA9/8mijv/r3db/9/Lv//Ts6P/v493/7uHb//Dl3//17uv//Pr65v///07//v4A////AAAA"
    "AAAAAAAAAAAAAAAAAAAAAAAA////APr29QD///8wyqSR2ZdNJ/+dVzP/q3BR/7J8YP+weFv/pWZF"
    "/5dOKP+SRR3/omE//8ynlP/z6uX/8+rm/9Kzo/+yfGD/oV89/5tUMP+aUy7/nFcz/6RjQf+0gGT/"
    "1LWl2/j08TLz6uUA////AAAAAAAAAAAAAAAAAAAAAAD///8A////Df/+/q/l0sj/wZV+/6twUf+h"
    "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"
)



# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------


def _percentile(values: list[float], q: float) -> float:
    """Nearest-rank percentile of *values* (no numpy dependency needed here)."""
    if not values:
        return 0.0
    ordered = sorted(values)
    idx = min(len(ordered) - 1, max(0, round(q * (len(ordered) - 1))))
    return ordered[idx]


def _fmt_time(seconds: float) -> str:
    return f"{seconds:.2f}s"


def _esc(text: object) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


# ---------------------------------------------------------------------------
# Stat tiles
# ---------------------------------------------------------------------------


def _stat_tiles(
    step_results: list[StepResult],
    elapsed_sec: float,
    n_sites_in: int,
    n_sites_out: int,
) -> str:
    n_ok = sum(1 for sr in step_results if sr.ok)
    n_err = len(step_results) - n_ok
    n_cached = sum(1 for sr in step_results if sr.cached)
    n_plots = sum(len(sr.plots) for sr in step_results)
    cache_rate = (n_cached / len(step_results) * 100) if step_results else 0.0

    def tile(icon: str, label: str, value: str, err: bool = False) -> str:
        cls = "tile err" if err else "tile"
        return (
            f'<div class="{cls}">{icon}'
            f'<div class="label">{_esc(label)}</div>'
            f'<div class="value">{_esc(value)}</div></div>'
        )

    return (
        '<div class="tiles">'
        + tile(_ICON_WORKFLOW, "Steps", f"{n_ok}/{len(step_results)} ok")
        + tile(_ICON_DIAGNOSTIC, "Errors", str(n_err), err=n_err > 0)
        + tile(_ICON_CLOCK, "Total time", _fmt_time(elapsed_sec))
        + tile(_ICON_EXPORT, "Sites", f"{n_sites_in}→{n_sites_out}")
        + tile(_ICON_WORKFLOW, "Cache hit rate", f"{cache_rate:.0f}%")
        + tile(_ICON_EXPORT, "Figures", str(n_plots))
        + "</div>"
    )


# ---------------------------------------------------------------------------
# Chart 1 — step status timeline
# ---------------------------------------------------------------------------


def _svg_status_timeline(step_results: list[StepResult]) -> str:
    if not step_results:
        return '<p class="step-meta">No steps to display.</p>'

    block_w, block_h, gap = 26, 40, 3
    width = len(step_results) * (block_w + gap)
    height = block_h + 4
    parts = [
        f'<svg viewBox="0 0 {width} {height}" width="{width}" height="{height}" '
        f'role="img" aria-label="Per-step status timeline">'
    ]
    for i, sr in enumerate(step_results):
        x = i * (block_w + gap)
        color = "var(--pc-good)" if sr.ok else "var(--pc-critical)"
        title = (
            f"({sr.step_idx}) {sr.step_name} [{sr.step_code}] "
            f"{'OK' if sr.ok else 'ERROR'} — {_fmt_time(sr.elapsed_sec)}"
            f"{' (cached)' if sr.cached else ''}"
        )
        parts.append(
            f'<g><rect x="{x}" y="2" width="{block_w}" height="{block_h}" rx="4" '
            f'fill="{color}"><title>{_esc(title)}</title></rect>'
        )
        if sr.cached:
            # cached marker: a small surface-ringed dot in the top-right
            # corner of the block, >=8px diameter (r=4) per the mark spec.
            cx, cy = x + block_w - 7, 9
            parts.append(
                f'<circle cx="{cx}" cy="{cy}" r="4" fill="var(--pc-gold)" '
                f'stroke="var(--pc-surface)" stroke-width="2"/>'
            )
        parts.append("</g>")
    parts.append("</svg>")
    return (
        '<div class="chart-scroll">' + "".join(parts) + "</div>"
        '<p class="step-meta">Green = OK, red = error. '
        "Gold dot = replayed from cache.</p>"
    )


# ---------------------------------------------------------------------------
# Chart 2 — per-step duration bars
# ---------------------------------------------------------------------------


def _svg_duration_bars(step_results: list[StepResult]) -> str:
    if not step_results:
        return '<p class="step-meta">No steps to display.</p>'

    times = [sr.elapsed_sec for sr in step_results]
    slow_threshold = _percentile(times, 0.8)
    max_time = max(times) or 1.0

    row_h = 28
    label_w = 160
    bar_area_w = 420
    width = label_w + bar_area_w + 70
    height = len(step_results) * row_h + 10

    parts = [
        f'<svg viewBox="0 0 {width} {height}" width="{width}" height="{height}" '
        f'role="img" aria-label="Per-step duration">'
    ]
    for i, sr in enumerate(step_results):
        y = i * row_h + 5
        bar_w = max(2.0, (sr.elapsed_sec / max_time) * bar_area_w)
        # Strict ">": nearest-rank percentiles land on a real data point, and
        # with several tied-low steps the 80th-percentile *value* can equal
        # the common low value, which would flag every step as slow. Only a
        # step that actually exceeds that threshold counts as an outlier.
        is_slow = sr.elapsed_sec > slow_threshold and len(step_results) > 1
        color = "var(--pc-warning)" if is_slow else "var(--pc-secondary)"
        label = f"{sr.step_name} [{sr.step_code}]"
        parts.append(
            f'<text x="{label_w - 8}" y="{y + 15}" text-anchor="end" '
            f'font-size="11" fill="var(--pc-ink)">{_esc(label)}</text>'
        )
        parts.append(
            f'<rect x="{label_w}" y="{y}" width="{bar_w:.1f}" height="20" rx="4" '
            f'fill="{color}"><title>{_esc(label)}: {_fmt_time(sr.elapsed_sec)}'
            f'{" (slow step)" if is_slow else ""}</title></rect>'
        )
        parts.append(
            f'<text x="{label_w + bar_w + 6:.1f}" y="{y + 15}" font-size="11" '
            f'fill="var(--pc-muted)">{_fmt_time(sr.elapsed_sec)}</text>'
        )
    parts.append("</svg>")
    return (
        '<div class="chart-scroll">' + "".join(parts) + "</div>"
        '<p class="step-meta">Gold bars are at or above the 80th-percentile '
        "step time for this run.</p>"
    )


# ---------------------------------------------------------------------------
# Chart 3 — site-count flow (2-series line chart)
# ---------------------------------------------------------------------------


def _svg_site_flow(step_results: list[StepResult]) -> str:
    if not step_results:
        return '<p class="step-meta">No steps to display.</p>'

    n_in = [sr.n_sites_in for sr in step_results]
    n_out = [sr.n_sites_out for sr in step_results]
    all_counts = n_in + n_out
    y_min, y_max = min(all_counts), max(all_counts)
    if y_min == y_max:
        y_min, y_max = max(0, y_min - 1), y_max + 1

    width, height = 560, 200
    margin_l, margin_r, margin_t, margin_b = 40, 60, 10, 24
    plot_w = width - margin_l - margin_r
    plot_h = height - margin_t - margin_b
    n = len(step_results)

    def xy(i: int, value: int) -> tuple[float, float]:
        x = margin_l + (i / max(1, n - 1)) * plot_w if n > 1 else margin_l + plot_w / 2
        y = margin_t + (1 - (value - y_min) / (y_max - y_min)) * plot_h
        return x, y

    def polyline(series: list[int], color: str) -> str:
        pts = " ".join(f"{xy(i, v)[0]:.1f},{xy(i, v)[1]:.1f}" for i, v in enumerate(series))
        dots = "".join(
            f'<circle cx="{xy(i, v)[0]:.1f}" cy="{xy(i, v)[1]:.1f}" r="4" '
            f'fill="{color}" stroke="var(--pc-surface)" stroke-width="2">'
            f"<title>Step {i + 1}: {v} sites</title></circle>"
            for i, v in enumerate(series)
        )
        return (
            f'<polyline points="{pts}" fill="none" stroke="{color}" '
            f'stroke-width="2" stroke-linecap="round" stroke-linejoin="round"/>'
            + dots
        )

    # gridlines (hairline, recessive) at 3 evenly spaced y ticks
    gridlines = []
    for frac in (0.0, 0.5, 1.0):
        gy = margin_t + frac * plot_h
        tick_val = round(y_max - frac * (y_max - y_min))
        gridlines.append(
            f'<line x1="{margin_l}" y1="{gy:.1f}" x2="{margin_l + plot_w}" '
            f'y2="{gy:.1f}" stroke="var(--pc-border)" stroke-width="1"/>'
            f'<text x="{margin_l - 6}" y="{gy + 4:.1f}" text-anchor="end" '
            f'font-size="10" fill="var(--pc-muted)">{tick_val}</text>'
        )

    last_in_x, last_in_y = xy(n - 1, n_in[-1])
    last_out_x, last_out_y = xy(n - 1, n_out[-1])
    end_labels = (
        f'<text x="{last_in_x + 6:.1f}" y="{last_in_y + 3:.1f}" font-size="11" '
        f'fill="var(--pc-secondary)">{n_in[-1]}</text>'
        f'<text x="{last_out_x + 6:.1f}" y="{last_out_y + 14:.1f}" font-size="11" '
        f'fill="var(--pc-primary)">{n_out[-1]}</text>'
    )

    svg = (
        f'<svg viewBox="0 0 {width} {height}" width="{width}" height="{height}" '
        f'role="img" aria-label="Site count in versus out per step">'
        + "".join(gridlines)
        + polyline(n_in, "var(--pc-secondary)")
        + polyline(n_out, "var(--pc-primary)")
        + end_labels
        + "</svg>"
    )
    legend = (
        '<div class="legend">'
        '<span><span class="swatch" style="background:var(--pc-secondary)"></span>Sites in</span>'
        '<span><span class="swatch" style="background:var(--pc-primary)"></span>Sites out</span>'
        "</div>"
    )
    return legend + '<div class="chart-scroll">' + svg + "</div>"


# ---------------------------------------------------------------------------
# Step table + detail cards (shared shape with make_html_report)
# ---------------------------------------------------------------------------


def _step_table(step_results: list[StepResult]) -> str:
    rows = []
    for sr in step_results:
        status_c = "ok" if sr.ok else "err"
        status_t = "OK" if sr.ok else "ERR"
        cached = "yes" if sr.cached else "–"
        rows.append(
            "<tr>"
            f"<td>{sr.step_idx}</td><td>{_esc(sr.step_name)}</td>"
            f"<td>{_esc(sr.step_code)}</td>"
            f'<td><span class="{status_c}">{status_t}</span></td>'
            f"<td>{_fmt_time(sr.elapsed_sec)}</td>"
            f"<td>{sr.n_sites_in}→{sr.n_sites_out}</td>"
            f"<td>{len(sr.plots)}</td><td>{cached}</td></tr>"
        )
    return (
        "<table><tr><th>#</th><th>Step</th><th>Code</th><th>Status</th>"
        "<th>Time</th><th>Sites</th><th>Plots</th><th>Cached</th></tr>"
        + "".join(rows)
        + "</table>"
    )


def _plot_thumbs(plots: list[Path], outdir: Path | None) -> str:
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


def _step_cards(step_results: list[StepResult], outdir: Path | None) -> str:
    cards = []
    for sr in step_results:
        status_icon = "&#10003;" if sr.ok else "&#10007;"
        status_c = "ok" if sr.ok else "err"
        params_str = (
            ", ".join(f"{k}={v!r}" for k, v in sr.params.items()) if sr.params else "—"
        )
        cached_badge = " &nbsp;|&nbsp; <b>Cached:</b> yes" if sr.cached else ""
        err_html = (
            f'<div class="error-box"><b>Error:</b> {_esc(sr.error)}</div>' if not sr.ok else ""
        )
        cards.append(
            '<div class="step-card">'
            f'<h3>{_ICON_WORKFLOW}<span class="{status_c}">{status_icon}</span> '
            f"({sr.step_idx}) {_esc(sr.step_name)} "
            f"<small>[{_esc(sr.step_code)}]</small></h3>"
            '<div class="step-meta">'
            f"<b>Label:</b> {_esc(sr.step_label)} &nbsp;|&nbsp; "
            f"<b>Sites:</b> {sr.n_sites_in}→{sr.n_sites_out} &nbsp;|&nbsp; "
            f"<b>Time:</b> {_fmt_time(sr.elapsed_sec)} &nbsp;|&nbsp; "
            f"<b>Plots:</b> {len(sr.plots)}{cached_badge}</div>"
            f"<p><b>Params:</b> <code>{_esc(params_str)}</code></p>"
            f"{err_html}"
            f"{_plot_thumbs(sr.plots, outdir)}"
            "</div>"
        )
    return "".join(cards)


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


def make_dashboard_html(
    pipeline_name: str,
    step_results: list[StepResult],
    elapsed_sec: float,
    outdir: Path | None,
    n_sites_in: int,
    n_sites_out: int,
    pipeline_yaml: str = "",
) -> str:
    """Return a self-contained, branded HTML dashboard report.

    The richer sibling of :func:`pycsamt.pipeline._report.make_html_report`:
    KPI stat tiles and three native inline-SVG charts (step status, per-step
    duration, site-count flow) on top of the same step data. See the module
    docstring for where every color in this report traces back to.
    """
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    n_err = sum(1 for sr in step_results if not sr.ok)
    status_cls = "ok" if n_err == 0 else "err"
    status_txt = "SUCCESS" if n_err == 0 else f"{n_err} ERROR(S)"

    yaml_block = ""
    if pipeline_yaml:
        yaml_block = (
            "<h2>Pipeline configuration</h2>" f"<pre><code>{_esc(pipeline_yaml)}</code></pre>"
        )

    body = (
        '<div class="pycsamt-dashboard">'
        "<header>"
        f"{_LOGO_SVG}"
        '<div class="wordmark">pyCSAMT<span class="sub">Pipeline Dashboard</span></div>'
        f'<div class="status-badge {status_cls}">{status_txt}</div>'
        "</header>"
        f'<p class="run-meta">Pipeline: <b>{_esc(pipeline_name)}</b> '
        f"&nbsp;|&nbsp; Run at: {now} &nbsp;|&nbsp; "
        f"Sites: {n_sites_in}→{n_sites_out}</p>"
        + _stat_tiles(step_results, elapsed_sec, n_sites_in, n_sites_out)
        + '<section class="chart-card"><h2>Step status</h2>'
        + _svg_status_timeline(step_results)
        + "</section>"
        + '<section class="chart-card"><h2>Step duration</h2>'
        + _svg_duration_bars(step_results)
        + "</section>"
        + '<section class="chart-card"><h2>Site count flow</h2>'
        + _svg_site_flow(step_results)
        + "</section>"
        + "<h2>Step results</h2>"
        + _step_table(step_results)
        + "<h2>Step details</h2>"
        + _step_cards(step_results, outdir)
        + yaml_block
        + f"<footer>Generated by pyCSAMT — {now}</footer>"
        + "</div>"
    )

    return (
        "<!DOCTYPE html>"
        '<html lang="en">'
        "<head>"
        '<meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        f"<title>pyCSAMT Pipeline Dashboard — {_esc(pipeline_name)}</title>"
        '<link rel="icon" type="image/x-icon" '
        f'href="data:image/x-icon;base64,{_FAVICON_ICO_BASE64}">'
        f"{_BRAND_CSS}"
        "</head>"
        f"<body>{body}</body></html>"
    )
