"""Generate figures and sample outputs for the inspect-and-QC tutorial."""

from __future__ import annotations

import contextlib
import io
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
IMAGE_DIR = (
    ROOT
    / "docs"
    / "source"
    / "images"
    / "tutorials"
    / "inspect_and_qc_survey"
)


def _import_pycsamt():
    sys.path.insert(0, str(ROOT))
    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        from pycsamt.api import read_edis
        from pycsamt.emtools.qc import (
            build_qc_table,
            frequency_confidence_table,
            plot_confidence_band_summary,
            plot_confidence_profile,
            plot_frequency_confidence_psection,
            qc_flags,
            station_confidence_table,
        )

    return {
        "read_edis": read_edis,
        "build_qc_table": build_qc_table,
        "frequency_confidence_table": frequency_confidence_table,
        "plot_confidence_band_summary": plot_confidence_band_summary,
        "plot_confidence_profile": plot_confidence_profile,
        "plot_frequency_confidence_psection": plot_frequency_confidence_psection,
        "qc_flags": qc_flags,
        "station_confidence_table": station_confidence_table,
    }


def _style_axis(ax: plt.Axes) -> None:
    ax.set_facecolor("#fbfbf7")
    ax.grid(True, color="#d7d7cf", linewidth=0.75, alpha=0.65)
    for spine in ax.spines.values():
        spine.set_color("#3d4752")
        spine.set_linewidth(0.8)


def _save(fig: plt.Figure, name: str) -> None:
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(IMAGE_DIR / name, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _make_inventory_plot(qc_df, station_ci_df) -> None:
    plot_df = qc_df.merge(
        station_ci_df[["station", "confidence"]],
        on="station",
        how="left",
    )
    x = np.arange(len(plot_df))
    fig, ax1 = plt.subplots(figsize=(11.5, 4.2))
    ax2 = ax1.twinx()

    ax1.bar(
        x,
        plot_df["snr_med"],
        color="#5aa9a4",
        edgecolor="#285c59",
        linewidth=0.5,
        alpha=0.82,
        label="median SNR proxy",
    )
    ax2.plot(
        x,
        plot_df["confidence"],
        color="#c44536",
        marker="D",
        markersize=4.5,
        linewidth=1.4,
        label="confidence",
    )
    ax2.axhspan(0.0, 0.60, color="#f3b3a6", alpha=0.18, zorder=0)
    ax2.axhline(0.60, color="#c44536", linestyle="--", linewidth=1.0)
    ax2.axhline(0.80, color="#234b6d", linestyle=":", linewidth=1.0)

    ax1.set_xticks(x[::2])
    ax1.set_xticklabels(plot_df["station"].iloc[::2], rotation=45, ha="right")
    ax1.set_ylabel("Median SNR proxy")
    ax2.set_ylabel("Station confidence")
    ax1.set_xlabel("Station")
    ax1.set_title("L18PLT station inventory: coverage is complete, review focuses on confidence")
    ax2.set_ylim(0.35, 1.02)
    _style_axis(ax1)
    for spine in ax2.spines.values():
        spine.set_color("#3d4752")
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper right", fontsize=8)
    _save(fig, "station_inventory_qc.png")


def _make_confidence_grid(functions, sites) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(11.2, 8.4), constrained_layout=True)
    ax = functions["plot_confidence_profile"](
        sites,
        method="composite",
        ci_hi=0.95,
        ci_lo=0.50,
        station_labels=True,
        station_label_step=2,
        spacing_m=200.0,
        recursive=False,
        ax=axes[0],
    )
    ax.set_title("Station confidence along the profile")
    for line in ax.lines:
        line.set_marker("o")
        line.set_markersize(4.5)
    _style_axis(ax)

    ax = functions["plot_confidence_band_summary"](
        sites,
        method="composite",
        ci_hi=0.95,
        ci_lo=0.50,
        recursive=False,
        ax=axes[1],
    )
    ax.set_title("Line-wide confidence by period")
    _style_axis(ax)
    _save(fig, "confidence_review_grid.png")


def _make_frequency_psection(functions, sites) -> None:
    fig, ax = plt.subplots(figsize=(11.2, 5.0))
    ax = functions["plot_frequency_confidence_psection"](
        sites,
        method="composite",
        ci_hi=0.95,
        ci_lo=0.50,
        cmap="viridis",
        station_label_step=2,
        recursive=False,
        ax=ax,
    )
    ax.set_title("Frequency-level confidence pseudosection")
    _style_axis(ax)
    _save(fig, "frequency_confidence_psection.png")


def main() -> int:
    functions = _import_pycsamt()
    survey = functions["read_edis"](
        DATA_DIR,
        recursive=False,
        strict=False,
        progress=False,
    )
    sites = survey.collection
    inventory = survey.df.to_pandas(copy=True)

    qc = functions["build_qc_table"](
        sites,
        include_skew=True,
        recursive=False,
        api=True,
    )
    qc_df = qc.to_pandas(copy=True)
    flagged = functions["qc_flags"](
        sites,
        min_frac_ok=0.75,
        min_snr_med=3.0,
        max_skew_med=6.0,
        recursive=False,
    )
    station_ci = functions["station_confidence_table"](
        sites,
        method="composite",
        relerr_threshold=0.20,
        offdiag_tolerance_log10=0.35,
        diagonal_leakage_max=0.35,
        phase_jump_tolerance_deg=90.0,
        spatial_tolerance_log10=0.60,
        spacing_m=200.0,
        recursive=False,
        api=True,
    )
    station_ci_df = station_ci.to_pandas(copy=True)
    freq_ci = functions["frequency_confidence_table"](
        sites,
        method="composite",
        ci_hi=0.95,
        ci_lo=0.50,
        recursive=False,
        api=True,
    )
    freq_ci_df = freq_ci.to_pandas(copy=True)

    _make_inventory_plot(qc_df, station_ci_df)
    _make_confidence_grid(functions, sites)
    _make_frequency_psection(functions, sites)

    print(f"survey_summary: {len(inventory)} stations x {inventory['n_freq'].iloc[0]} frequencies")
    print("inventory_head:")
    print(inventory.head(3).to_string(index=False))
    print("qc_head:")
    print(
        qc_df[
            ["station", "n_freq", "frac_ok", "snr_med", "skew_med", "pmin", "pmax"]
        ]
        .head(5)
        .to_string(index=False)
    )
    print("flag_counts:")
    print(flagged["flags"].value_counts().head().to_string())
    low_ci = station_ci_df[station_ci_df["confidence"] < 0.60]
    high_ci = station_ci_df[station_ci_df["confidence"] >= 0.80]
    print(f"confidence_groups: review={len(low_ci)}, first_trial={len(high_ci)}")
    print("low_confidence:")
    print(
        low_ci[["station", "confidence", "coverage", "phase", "spatial"]]
        .to_string(index=False)
    )
    weak_freq = freq_ci_df[freq_ci_df["confidence"] < 0.50]
    print(f"weak_frequency_rows: {len(weak_freq)}")
    print(
        weak_freq[["station", "frequency_hz", "period_s", "confidence", "flags"]]
        .head(8)
        .to_string(index=False)
    )
    print(f"images: {IMAGE_DIR.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
