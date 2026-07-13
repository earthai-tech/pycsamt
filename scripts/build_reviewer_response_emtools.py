"""Build reviewer-response diagnostics for bundled WILLY AMT lines.

The script intentionally uses only public pycsamt.emtools functions so it
can be handed to a student and adapted to their own EDI folders.  It writes
CSV tables and PNG figures under results/reviewer_response_emtools/.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pycsamt.emtools import (
    confidence_gated_emap_filter,
    correct_ss_ama,
    edit_frequencies_by_confidence,
    emi_mitigation_report,
    ensure_sites,
    estimate_ss_ama,
    frequency_confidence_table,
    groom_bailey_table,
    plot_confidence_band_summary,
    plot_confidence_profile,
    plot_frequency_confidence_psection,
    plot_skew_traffic_psection,
    plot_ss_delta_profile,
    plot_strike_profile,
    pre2d_inversion_assessment,
    qc_flags,
    station_confidence_table,
)
from pycsamt.emtools.remove_noise import snr_table


DEFAULT_LINES = {
    "L18PLT": ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT",
    "L22PLT": ROOT / "data" / "AMT" / "WILLY_DATA" / "L22PLT",
}


def _write_table(table: Any, path: Path) -> pd.DataFrame:
    """Write a DataFrame-like object and return a plain DataFrame."""
    df = getattr(table, "df", table)
    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)
    df.to_csv(path, index=False)
    return df


def _save_plot(result: Any, path: Path) -> None:
    """Save a figure returned either as a Figure or an Axes."""
    fig = result
    if hasattr(result, "get_figure"):
        fig = result.get_figure()
    if fig is None:
        fig = plt.gcf()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _safe_plot(label: str, func: Any, path: Path, *args: Any, **kwargs: Any) -> None:
    error_path = path.with_suffix(".error.txt")
    if error_path.exists():
        error_path.unlink()
    try:
        _save_plot(func(*args, **kwargs), path)
    except Exception as exc:  # pragma: no cover - report-building guard
        error_path.write_text(
            f"{label} failed: {exc}\n",
            encoding="utf-8",
        )


def _summarise_line(line: str, edi_dir: Path, out_dir: Path) -> None:
    line_dir = out_dir / line
    fig_dir = line_dir / "figures"
    line_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    sites = ensure_sites(str(edi_dir), recursive=False, verbose=0)

    conf_freq = _write_table(
        frequency_confidence_table(sites, method="composite"),
        line_dir / "frequency_confidence.csv",
    )
    _write_table(
        station_confidence_table(sites, method="composite"),
        line_dir / "station_confidence.csv",
    )
    _write_table(qc_flags(sites), line_dir / "qc_flags.csv")
    _write_table(snr_table(sites), line_dir / "snr_table.csv")
    _write_table(
        emi_mitigation_report(
            sites,
            remote_reference_attempted=False,
            remote_reference_reason=(
                "No remote-reference time series are included with the "
                "bundled WILLY EDI example; mitigation is documented at "
                "the transfer-function level."
            ),
            mains_hz=50.0,
        ),
        line_dir / "emi_mitigation_report.csv",
    )

    edited = edit_frequencies_by_confidence(
        sites,
        mode="recover",
        before_sites=sites,
        ci_hi=0.95,
        ci_lo=0.85,
        reject="mask",
        interpolation="linear",
    )
    _write_table(edited.report, line_dir / "frequency_edit_report.csv")
    _write_table(edited.decisions, line_dir / "frequency_edit_decisions.csv")

    emap = confidence_gated_emap_filter(
        sites,
        before_sites=sites,
        method="flma",
        component="xy",
        ci_hi=0.95,
        ci_lo=0.85,
        window=5,
    )
    _write_table(emap.report, line_dir / "confidence_gated_emap_report.csv")
    _write_table(
        emap.decisions,
        line_dir / "confidence_gated_emap_decisions.csv",
    )

    rows = []
    for half_window in (1, 2, 3, 4):
        factors = estimate_ss_ama(
            sites,
            sort_by="lon",
            half_window=half_window,
            weights="tri",
            max_skew=6.0,
        )
        factors_df = _write_table(
            factors,
            line_dir / f"ss_ama_factors_half_window_{half_window}.csv",
        )
        for _, row in factors_df.iterrows():
            rows.append(
                {
                    "line": line,
                    "half_window": half_window,
                    "station": row.get("station"),
                    "delta_log10_rho": row.get("delta_log10_rho"),
                    "fac_rho": row.get("fac_rho"),
                    "fac_z": row.get("fac_z"),
                    "n_used": row.get("n_used"),
                }
            )
    pd.DataFrame(rows).to_csv(
        line_dir / "ss_ama_window_sensitivity.csv",
        index=False,
    )

    gb = _write_table(
        groom_bailey_table(sites, min_freq=4, robust=True),
        line_dir / "groom_bailey_table.csv",
    )
    gb_ok = bool((gb.get("status", pd.Series(dtype=str)) == "ok").any())
    _write_table(
        pre2d_inversion_assessment(
            sites,
            rotation_applied=False,
            groom_bailey_attempted=True,
            groom_bailey_applied=False,
            groom_bailey_reason=(
                "Groom-Bailey table was computed as a diagnostic. "
                "Correction was not applied in this screening script; "
                f"valid station fits present: {gb_ok}."
            ),
        ),
        line_dir / "pre2d_inversion_assessment.csv",
    )

    _safe_plot(
        "confidence profile",
        plot_confidence_profile,
        fig_dir / "confidence_profile.png",
        sites,
        method="composite",
    )
    _safe_plot(
        "frequency confidence pseudosection",
        plot_frequency_confidence_psection,
        fig_dir / "frequency_confidence_psection.png",
        sites,
        method="composite",
    )
    _safe_plot(
        "confidence band summary",
        plot_confidence_band_summary,
        fig_dir / "confidence_band_summary.png",
        sites,
        method="composite",
    )
    ss_corrected = correct_ss_ama(
        sites,
        sort_by="lon",
        half_window=3,
        weights="tri",
        max_skew=6.0,
        inplace=False,
    )
    _safe_plot(
        "static-shift delta profile",
        plot_ss_delta_profile,
        fig_dir / "ss_delta_profile.png",
        sites,
        ss_corrected,
    )
    _safe_plot(
        "strike profile",
        plot_strike_profile,
        fig_dir / "strike_profile.png",
        sites,
    )
    _safe_plot(
        "skew traffic pseudosection",
        plot_skew_traffic_psection,
        fig_dir / "skew_traffic_psection.png",
        sites,
    )

    summary = {
        "line": line,
        "edi_dir": str(edi_dir),
        "n_frequency_confidence_rows": len(conf_freq),
        "n_safe": int((conf_freq["confidence"] >= 0.95).sum())
        if "confidence" in conf_freq
        else 0,
        "n_recoverable": int(
            (
                (conf_freq["confidence"] >= 0.85)
                & (conf_freq["confidence"] < 0.95)
            ).sum()
        )
        if "confidence" in conf_freq
        else 0,
        "n_reject": int((conf_freq["confidence"] < 0.85).sum())
        if "confidence" in conf_freq
        else 0,
        "emap_preserved": emap.n_preserved,
        "emap_blended": emap.n_blended,
        "emap_filtered": emap.n_filtered,
    }
    pd.DataFrame([summary]).to_csv(line_dir / "line_summary.csv", index=False)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build reviewer-response tables and figures.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=ROOT / "results" / "reviewer_response_emtools",
        help="Output directory.",
    )
    parser.add_argument(
        "--line",
        action="append",
        default=[],
        help=(
            "Line specification NAME=EDI_DIR. May be repeated. "
            "Defaults to bundled L18PLT and L22PLT."
        ),
    )
    args = parser.parse_args()

    lines = dict(DEFAULT_LINES)
    for spec in args.line:
        if "=" not in spec:
            raise SystemExit("--line must look like NAME=EDI_DIR")
        name, value = spec.split("=", 1)
        lines[name.strip()] = Path(value).expanduser().resolve()

    args.out.mkdir(parents=True, exist_ok=True)
    summaries = []
    for line, edi_dir in lines.items():
        if not edi_dir.exists():
            raise SystemExit(f"EDI directory not found for {line}: {edi_dir}")
        _summarise_line(line, edi_dir, args.out)
        summaries.append(pd.read_csv(args.out / line / "line_summary.csv"))
    pd.concat(summaries, ignore_index=True).to_csv(
        args.out / "all_lines_summary.csv",
        index=False,
    )
    print(f"Wrote reviewer-response artifacts to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
