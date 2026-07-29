"""
pycsamt.agents.loader
=====================

:class:`MTLoaderAgent` — Load MT/AMT/CSAMT data into a pycsamt
:class:`~pycsamt.site.Sites` object.

Accepts any input supported by :func:`~pycsamt.emtools._core.ensure_sites`:

* A single EDI file path.
* A directory of EDI / AVG / J files.
* A list of file paths.
* An existing :class:`~pycsamt.site.Sites` or :class:`~pycsamt.seg.EDICollection`.

After loading the agent runs a per-station quality scan and returns:

* ``data["sites"]``         — the :class:`Sites` object.
* ``data["station_names"]`` — ordered list of station names.
* ``data["n_stations"]``    — total station count.
* ``data["quality_table"]`` — ``pandas.DataFrame`` with per-station scores.
* ``data["summary_stats"]`` — dict of survey-level statistics.

The quality table columns:

================  =====================================================
Column            Meaning
================  =====================================================
station           Station name
has_z             Whether the Z impedance tensor is present
has_tipper        Whether the Tipper block is present
has_coords        Whether lat / lon coordinates are available
n_freq            Number of frequencies with finite Z
t_min_s           Shortest available period (s)
t_max_s           Longest available period (s)
snr_proxy         Median |Zxy| / std(|Zxy|) across frequencies (proxy)
qc_score          Integer 0–100 composite quality score
================  =====================================================
"""

from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Any

import numpy as np

from ..emtools._core import (
    _get_t_block,
    _get_z_block,
    _iter_items,
    _name,
    ensure_sites,
)
from ._base import AgentResult, BaseAgent

logger = logging.getLogger(__name__)


# ── system prompt for LLM interpretation ─────────────────────────────────────

_SYSTEM_PROMPT = """\
You are an expert MT/AMT/CSAMT data quality analyst.
Given a per-station data quality summary, write 2–3 concise sentences that:
1. State the overall data quality.
2. Flag any stations or frequency ranges that need attention.
3. Recommend the next processing step.
Reply in plain English — no bullet points, no markdown.
"""


class MTLoaderAgent(BaseAgent):
    """Load MT data from any pycsamt-supported format and assess quality.

    Parameters
    ----------
    api_key : str or None
    model, llm_provider : str
    recursive : bool
        When loading a directory, recurse into sub-directories.
    on_dup : str
        Duplicate-station handling: ``"replace"`` (default) or ``"skip"``.

    Examples
    --------
    >>> agent = MTLoaderAgent()
    >>> result = agent.execute({"path": "/data/AMT/WILLY_DATA/L22PLT"})
    >>> result.status
    'success'
    >>> result["n_stations"]
    25
    >>> result["quality_table"].head()
         station  has_z  ...  qc_score
    0  22-22BF    True  ...        88
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        recursive: bool = True,
        on_dup: str = "replace",
    ) -> None:
        super().__init__(
            "MTLoaderAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="pseudosection",
        )
        self.recursive = recursive
        self.on_dup = on_dup

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()

        # ── resolve input path ────────────────────────────────────────────────
        path = (
            input_data.get("path")
            or input_data.get("data_path")
            or input_data.get("sites")  # already a Sites object
        )
        if path is None:
            return AgentResult.failed(
                "No data path supplied. Pass {'path': '/path/to/EDIs'}.",
                hint="Provide 'path' key pointing to EDI files or a directory.",
                elapsed=time.time() - t0,
            )

        # ── load ──────────────────────────────────────────────────────────────
        try:
            sites = ensure_sites(
                path,
                recursive=self.recursive,
                on_dup=self.on_dup,
                verbose=0,
            )
        except Exception as exc:
            return AgentResult.failed(
                f"Failed to load data from {path!r}: {exc}",
                hint=(
                    "Check that the path contains valid EDI / AVG / J files. "
                    "Use recursive=True to search sub-directories."
                ),
                elapsed=time.time() - t0,
            )

        # ── quality scan ──────────────────────────────────────────────────────
        rows, warnings = _quality_scan(sites)

        if not rows:
            return AgentResult.failed(
                "No stations with valid impedance data found.",
                hint="Ensure the files contain ZXXR/ZXXI … blocks or RES/PHS data.",
                elapsed=time.time() - t0,
            )

        try:
            import pandas as pd

            qt = pd.DataFrame(rows)
        except ImportError:
            qt = rows  # fall back to list of dicts if pandas not available

        station_names = [r["station"] for r in rows]
        n_st = len(station_names)

        # ── summary statistics ────────────────────────────────────────────────
        scores = np.array([r["qc_score"] for r in rows], float)
        has_z = sum(1 for r in rows if r["has_z"])
        has_coords = sum(1 for r in rows if r["has_coords"])
        n_freq_all = [r["n_freq"] for r in rows if r["n_freq"] > 0]

        summary_stats = {
            "n_stations": n_st,
            "n_with_z": has_z,
            "n_with_coords": has_coords,
            "n_with_tipper": sum(1 for r in rows if r["has_tipper"]),
            "mean_qc_score": float(np.mean(scores)) if scores.size else 0.0,
            "min_qc_score": float(np.min(scores)) if scores.size else 0.0,
            "median_n_freq": float(np.median(n_freq_all)) if n_freq_all else 0.0,
            "global_t_min_s": min(
                (r["t_min_s"] for r in rows if r["t_min_s"] is not None),
                default=None,
            ),
            "global_t_max_s": max(
                (r["t_max_s"] for r in rows if r["t_max_s"] is not None),
                default=None,
            ),
        }

        # ── LLM interpretation ────────────────────────────────────────────────
        interpretation: str | None = None
        if self.api_key:
            prompt = _build_interpretation_prompt(summary_stats, rows, warnings)
            interpretation = self.query_llm(prompt, max_tokens=200)

        elapsed = time.time() - t0
        mean_score = summary_stats["mean_qc_score"]
        status = (
            "success"
            if mean_score >= 60
            else "needs_review"
            if mean_score >= 30
            else "failed"
        )
        summary = (
            (
                f"Loaded {n_st} station(s) from {_path_label(path)}. "
                f"Mean QC score {mean_score:.0f}/100. "
                f"Period range: "
                f"{summary_stats['global_t_min_s']:.2e}–"
                f"{summary_stats['global_t_max_s']:.2e} s."
            )
            if summary_stats["global_t_min_s"]
            else (f"Loaded {n_st} station(s); no finite periods found.")
        )

        return AgentResult(
            status=status,
            summary=summary,
            data={
                "sites": sites,
                "station_names": station_names,
                "n_stations": n_st,
                "quality_table": qt,
                "summary_stats": summary_stats,
                "path": str(path),
            },
            warnings=warnings,
            llm_interpretation=interpretation,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


# ── quality scan ──────────────────────────────────────────────────────────────


def _quality_scan(
    sites: Any,
) -> tuple[list[dict], list[str]]:
    """Iterate over *sites* and compute a per-station quality record.

    Returns
    -------
    rows : list of dict
    warnings : list of str
    """
    rows: list[dict] = []
    warnings: list[str] = []

    for i, ed in enumerate(_iter_items(sites)):
        nm = _name(ed, i)

        # ── Z tensor ─────────────────────────────────────────────────────────
        Z_obj, z, fr = _get_z_block(ed)
        has_z = z is not None and fr is not None

        # ── tipper ───────────────────────────────────────────────────────────
        T_obj, t, fr_t = _get_t_block(ed)
        has_tipper = t is not None

        # ── coordinates ──────────────────────────────────────────────────────
        coords = getattr(ed, "coords", None)
        if coords is None:
            edi_inner = getattr(ed, "edi", None)
            if edi_inner is not None:
                coords = getattr(edi_inner, "coords", None)
        has_coords = (
            coords is not None
            and len(coords) >= 2
            and coords[0] is not None
            and coords[1] is not None
        )

        # ── frequency / period stats ──────────────────────────────────────────
        n_freq = 0
        t_min_s = None
        t_max_s = None
        snr_proxy = np.nan

        if has_z and fr is not None:
            per = 1.0 / np.where(fr == 0, np.nan, fr)
            finite_mask = np.isfinite(per)
            finite_per = per[finite_mask]
            n_freq = int(finite_mask.sum())
            if n_freq > 0:
                t_min_s = float(np.nanmin(finite_per))
                t_max_s = float(np.nanmax(finite_per))
                # SNR proxy: |Zxy| / std(|Zxy|) across freqs — higher = cleaner
                zxy = np.abs(z[finite_mask, 0, 1])
                if zxy.size > 2 and np.std(zxy) > 0:
                    snr_proxy = float(np.mean(zxy) / np.std(zxy))

        # ── QC score 0-100 ────────────────────────────────────────────────────
        score = _compute_qc_score(
            has_z=has_z,
            has_coords=has_coords,
            has_tipper=has_tipper,
            n_freq=n_freq,
            snr_proxy=snr_proxy,
        )

        # ── per-station warnings ──────────────────────────────────────────────
        if not has_z:
            warnings.append(f"{nm}: no impedance (Z) data found.")
        if not has_coords:
            warnings.append(f"{nm}: no geographic coordinates.")
        if n_freq < 5:
            warnings.append(f"{nm}: only {n_freq} finite frequency band(s).")

        rows.append(
            {
                "station": nm,
                "has_z": has_z,
                "has_tipper": has_tipper,
                "has_coords": has_coords,
                "n_freq": n_freq,
                "t_min_s": t_min_s,
                "t_max_s": t_max_s,
                "snr_proxy": round(float(snr_proxy), 2)
                if np.isfinite(snr_proxy)
                else None,
                "qc_score": score,
            }
        )

    return rows, warnings


def _compute_qc_score(
    *,
    has_z: bool,
    has_coords: bool,
    has_tipper: bool,
    n_freq: int,
    snr_proxy: float,
) -> int:
    """Return an integer QC score 0–100.

    Weights:
      - Z tensor present             : 40 pts
      - Coordinates present          : 20 pts
      - ≥ 10 finite frequencies      : 20 pts  (linear up to 10)
      - SNR proxy ≥ 3                : 10 pts  (capped at 10)
      - Tipper present               :  5 pts
      - No extreme SNR warning       :  5 pts
    """
    score = 0
    if has_z:
        score += 40
    if has_coords:
        score += 20
    # frequency count contribution: linear 0→10 pts for n_freq 0→10+
    score += min(20, n_freq * 2)
    if np.isfinite(snr_proxy) and snr_proxy > 0:
        score += min(10, int(snr_proxy * 3))
    if has_tipper:
        score += 5
    # bonus: clean SNR (>5 = no aliasing suspected)
    if np.isfinite(snr_proxy) and snr_proxy >= 5:
        score += 5
    return min(100, score)


# ── helpers ───────────────────────────────────────────────────────────────────


def _path_label(path: Any) -> str:
    """Return a short display label for *path*."""
    try:
        p = Path(str(path))
        return p.name or str(p)
    except Exception:
        return str(path)[:60]


def _build_interpretation_prompt(
    stats: dict,
    rows: list[dict],
    warnings: list[str],
) -> str:
    low_score = [r["station"] for r in rows if r["qc_score"] < 50]
    no_z = [r["station"] for r in rows if not r["has_z"]]
    return (
        f"Survey quality summary:\n"
        f"  Stations: {stats['n_stations']}, "
        f"  with Z: {stats['n_with_z']}, "
        f"  with coords: {stats['n_with_coords']}\n"
        f"  Mean QC score: {stats['mean_qc_score']:.0f}/100\n"
        f"  Period range: {stats['global_t_min_s']:.2e} – "
        f"{stats['global_t_max_s']:.2e} s\n"
        f"  Low-score stations (<50): {low_score or 'none'}\n"
        f"  Stations without Z: {no_z or 'none'}\n"
        f"  Warnings ({len(warnings)}): "
        f"{warnings[:3] if warnings else 'none'}\n\n"
        f"Briefly assess overall quality and recommend the next step."
    )


__all__ = ["MTLoaderAgent"]
