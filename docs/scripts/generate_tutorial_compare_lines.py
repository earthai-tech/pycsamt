"""Generate figures and sample outputs for the survey-line comparison tutorial."""

from __future__ import annotations

import contextlib
import io
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = ROOT / "data" / "AMT" / "WILLY_DATA"
LINES = {
    "L18PLT": DATA_ROOT / "L18PLT",
    "L22PLT": DATA_ROOT / "L22PLT",
}
IMAGE_DIR = (
    ROOT
    / "docs"
    / "source"
    / "images"
    / "tutorials"
    / "compare_survey_lines_for_qc"
)


def _import_pycsamt():
    sys.path.insert(0, str(ROOT))
    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        pass

    return locals()


def _style_axis(ax: plt.Axes) -> None:
    ax.set_facecolor("#fbfbf7")
    ax.grid(True, color="#d8dad1", linewidth=0.7, alpha=0.65)
    for spine in ax.spines.values():
        spine.set_color("#39434d")
        spine.set_linewidth(0.8)


def _save(fig: plt.Figure, name: str) -> None:
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(IMAGE_DIR / name, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _load_line(functions, name: str, path: Path) -> dict:
    survey = functions["read_edis"](
        path,
        recursive=False,
        strict=False,
        on_dup="replace",
        progress=False,
    )
    sites = survey.collection
    inventory = survey.summary().to_pandas(copy=True)
    qc = functions["build_qc_table"](
        sites,
        include_skew=True,
        recursive=False,
        api=True,
    ).to_pandas(copy=True)
    confidence = functions["station_confidence_table"](
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
    ).to_pandas(copy=True)
    inventory["line"] = name
    qc["line"] = name
    confidence["line"] = name
    return {
        "survey": survey,
        "sites": sites,
        "inventory": inventory,
        "qc": qc,
        "confidence": confidence,
    }


def _line_summary(records: dict[str, dict]) -> pd.DataFrame:
    rows = []
    for line, rec in records.items():
        inv = rec["inventory"]
        qc = rec["qc"]
        ci = rec["confidence"]
        rows.append(
            {
                "line": line,
                "stations": len(inv),
                "median_n_freq": int(inv["n_freq"].median()),
                "min_freq_hz": 1.0 / float(qc["pmax"].max()),
                "max_freq_hz": 1.0 / float(qc["pmin"].min()),
                "median_frac_ok": float(qc["frac_ok"].median()),
                "median_confidence": float(ci["confidence"].median()),
                "tipper_files": int(inv["tipper"].sum()),
                "spectra_files": int(inv["spectra"].sum()),
            }
        )
    return pd.DataFrame(rows)


def _inventory_plot(summary: pd.DataFrame) -> None:
    labels = summary["line"].tolist()
    x = np.arange(len(labels))
    width = 0.24
    fig, ax = plt.subplots(figsize=(8.6, 4.4))
    metrics = [
        ("stations", "#2f6f8f"),
        ("tipper_files", "#d5962c"),
        ("spectra_files", "#3f8f61"),
    ]
    for idx, (column, color) in enumerate(metrics):
        ax.bar(
            x + (idx - 1) * width,
            summary[column],
            width=width,
            label=column.replace("_", " "),
            color=color,
            edgecolor="#27323a",
            linewidth=0.45,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Count")
    ax.set_title("Line inventory comparison")
    ax.legend(fontsize=8, ncol=3, loc="upper right")
    _style_axis(ax)
    _save(fig, "line_inventory_comparison.png")


def _frequency_overlap_plot(summary: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(9.4, 3.8))
    colors = {"L18PLT": "#2f6f8f", "L22PLT": "#c85745"}
    for idx, row in summary.iterrows():
        ax.hlines(
            y=idx,
            xmin=row["min_freq_hz"],
            xmax=row["max_freq_hz"],
            color=colors.get(row["line"], "#7c4d79"),
            linewidth=9,
            alpha=0.82,
        )
        ax.scatter(
            [row["min_freq_hz"], row["max_freq_hz"]],
            [idx, idx],
            color="#27323a",
            s=28,
            zorder=3,
        )
    shared_min = summary["min_freq_hz"].max()
    shared_max = summary["max_freq_hz"].min()
    ax.axvspan(
        shared_min,
        shared_max,
        color="#d5962c",
        alpha=0.18,
        label="shared band",
    )
    ax.set_xscale("log")
    ax.set_yticks(np.arange(len(summary)))
    ax.set_yticklabels(summary["line"])
    ax.set_xlabel("Frequency (Hz)")
    ax.set_title("Usable frequency-band overlap")
    ax.legend(fontsize=8, loc="lower right")
    _style_axis(ax)
    _save(fig, "frequency_overlap_l18_l22.png")


def _confidence_plot(records: dict[str, dict]) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 4.4))
    values = [
        records[line]["confidence"]["confidence"].to_numpy(dtype=float)
        for line in records
    ]
    labels = list(records)
    parts = ax.violinplot(values, showmeans=True, showextrema=True)
    colors = ["#2f6f8f", "#c85745"]
    for body, color in zip(parts["bodies"], colors):
        body.set_facecolor(color)
        body.set_edgecolor("#27323a")
        body.set_alpha(0.72)
    for key in ("cmeans", "cmins", "cmaxes", "cbars"):
        parts[key].set_color("#27323a")
        parts[key].set_linewidth(0.9)
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_ylabel("Station confidence")
    ax.set_ylim(0.0, 1.05)
    ax.set_title("Station-confidence distribution by line")
    _style_axis(ax)
    _save(fig, "qc_confidence_comparison.png")


def _decision_plot(summary: pd.DataFrame) -> None:
    shared_min = summary["min_freq_hz"].max()
    shared_max = summary["max_freq_hz"].min()
    decisions = pd.DataFrame(
        [
            ["Station count close enough for shared review?", "yes"],
            [
                "Common frequency band available?",
                f"yes, {shared_min:.2g}-{shared_max:.2g} Hz",
            ],
            ["Same first-pass notch settings?", "yes, start with 50 Hz"],
            ["Use identical static-shift factors?", "no, estimate per line"],
            ["Reuse one pipeline config directly?", "yes, after band check"],
        ],
        columns=["question", "decision"],
    )

    fig, ax = plt.subplots(figsize=(10.4, 3.2))
    ax.axis("off")
    table = ax.table(
        cellText=decisions.values,
        colLabels=decisions.columns,
        cellLoc="left",
        colLoc="left",
        loc="center",
        colWidths=[0.68, 0.32],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9.2)
    table.scale(1, 1.55)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#39434d")
        cell.set_linewidth(0.5)
        if row == 0:
            cell.set_facecolor("#2f6f8f")
            cell.get_text().set_color("white")
            cell.get_text().set_weight("bold")
        elif col == 1:
            cell.set_facecolor("#eef4ed")
        else:
            cell.set_facecolor("#fbfbf7")
    ax.set_title("Processing decision table", pad=12)
    _save(fig, "processing_decision_table.png")


def main() -> int:
    functions = _import_pycsamt()
    records = {
        line: _load_line(functions, line, path) for line, path in LINES.items()
    }
    summary = _line_summary(records)

    _inventory_plot(summary)
    _frequency_overlap_plot(summary)
    _confidence_plot(records)
    _decision_plot(summary)

    combined_inventory = pd.concat(
        [rec["inventory"] for rec in records.values()],
        ignore_index=True,
    )
    combined_qc = pd.concat(
        [rec["qc"] for rec in records.values()],
        ignore_index=True,
    )
    combined_ci = pd.concat(
        [rec["confidence"] for rec in records.values()],
        ignore_index=True,
    )

    print("line_summary:")
    print(
        summary[
            [
                "line",
                "stations",
                "median_n_freq",
                "min_freq_hz",
                "max_freq_hz",
                "median_frac_ok",
                "median_confidence",
            ]
        ].to_string(
            index=False,
            formatters={
                "min_freq_hz": "{:.3g}".format,
                "max_freq_hz": "{:.3g}".format,
                "median_frac_ok": "{:.3f}".format,
                "median_confidence": "{:.3f}".format,
            },
        )
    )
    print("inventory_head:")
    print(
        combined_inventory[["line", "station", "n_freq", "tipper", "spectra"]]
        .head(6)
        .to_string(index=False)
    )
    print("qc_by_line:")
    qc_by_line = combined_qc.groupby("line").agg(
        stations=("station", "count"),
        frac_ok_median=("frac_ok", "median"),
        snr_med_median=("snr_med", "median"),
        skew_med_median=("skew_med", "median"),
    )
    print(qc_by_line.to_string(float_format=lambda value: f"{value:.3f}"))
    print("confidence_by_line:")
    ci_by_line = combined_ci.groupby("line").agg(
        stations=("station", "count"),
        confidence_min=("confidence", "min"),
        confidence_median=("confidence", "median"),
        confidence_max=("confidence", "max"),
    )
    print(ci_by_line.to_string(float_format=lambda value: f"{value:.3f}"))
    shared_min = summary["min_freq_hz"].max()
    shared_max = summary["max_freq_hz"].min()
    print("decision:")
    print(f"shared_band_hz={shared_min:.3g}-{shared_max:.3g}")
    print("reuse_config=yes_for_first_qc")
    print("static_shift=estimate_per_line")
    print(f"images: {IMAGE_DIR.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
