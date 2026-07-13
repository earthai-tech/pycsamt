"""Generate figures and sample outputs for the pipeline-config tutorial."""

from __future__ import annotations

import contextlib
import io
import shutil
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
    / "run_pipeline_from_config"
)
BUILD_DIR = ROOT / "docs" / "_generated" / "tutorial_pipeline"
CONFIG_PATH = BUILD_DIR / "config" / "l18_first_qc.yaml"


CONFIG_TEXT = """name: l18_first_qc
output_dir: docs/_generated/tutorial_pipeline/results/l18_first_qc

steps:
  - name: notch
    code: NR001
    params:
      mains_hz: 50.0
      n_harm: 30
      tol_hz: 0.08

  - name: drop_duplicates
    code: FREQ002

  - name: select_band
    code: FREQ001
    params:
      band_hz: [1.0, 10000.0]

  - name: align_grid
    code: FREQ004

  - name: qc_snapshot
    code: QC001
"""


def _import_pycsamt():
    sys.path.insert(0, str(ROOT))
    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        from pycsamt.api import read_edis
        from pycsamt.pipeline import Pipeline, configure_pipe

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


def _step_status_plot(result) -> None:
    labels = [sr.step_name for sr in result.step_results]
    elapsed = np.asarray([sr.elapsed_sec for sr in result.step_results], dtype=float)
    plots = np.asarray([len(sr.plots) for sr in result.step_results], dtype=float)
    colors = ["#2f6f8f" if sr.ok else "#c85745" for sr in result.step_results]
    y = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    ax.barh(y, elapsed, color=colors, edgecolor="#27323a", linewidth=0.45)
    for yi, seconds, count in zip(y, elapsed, plots):
        ax.text(
            seconds + max(elapsed) * 0.025,
            yi,
            f"{seconds:.2f}s, {int(count)} plots",
            va="center",
            fontsize=9,
            color="#27323a",
        )
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Elapsed time (s)")
    ax.set_title("Pipeline step status for L18PLT")
    _style_axis(ax)
    _save(fig, "pipeline_step_status.png")


def _artifact_plot(result) -> None:
    counts = {
        "processed EDIs": len(result.processed_paths),
        "QC plots": len(result.plots),
        "reports": int((result.outdir / "summary.txt").exists())
        + int((result.outdir / "report.html").exists()),
        "saved configs": int((result.outdir / "pipeline.yaml").exists()),
    }
    labels = list(counts)
    values = list(counts.values())
    colors = ["#3f8f61", "#2f6f8f", "#d5962c", "#7c4d79"]

    fig, ax = plt.subplots(figsize=(8.6, 4.0))
    ax.bar(labels, values, color=colors, edgecolor="#27323a", linewidth=0.5)
    ax.set_ylabel("File count")
    ax.set_title("Artifacts written by the config run")
    ax.tick_params(axis="x", rotation=15)
    for idx, value in enumerate(values):
        ax.text(idx, value + 0.4, str(value), ha="center", fontsize=9)
    _style_axis(ax)
    _save(fig, "pipeline_output_artifacts.png")


def _workflow_plot(pipe) -> None:
    table = pipe.describe()
    labels = table["label"].tolist()
    codes = table["code"].tolist()
    categories = table["category"].tolist()
    palette = {
        "noise_removal": "#5f8fb4",
        "frequency": "#d5962c",
        "qc": "#3f8f61",
    }

    fig, ax = plt.subplots(figsize=(10.8, 2.8))
    ax.set_axis_off()
    for idx, (label, code, category) in enumerate(zip(labels, codes, categories)):
        x = idx / max(len(labels) - 1, 1)
        color = palette.get(category, "#7c4d79")
        ax.scatter(x, 0.56, s=900, color=color, edgecolor="#27323a", linewidth=1.0)
        ax.text(x, 0.56, code, color="white", ha="center", va="center", fontsize=9, weight="bold")
        ax.text(x, 0.18, label, ha="center", va="center", fontsize=9, color="#27323a")
        if idx < len(labels) - 1:
            ax.annotate(
                "",
                xy=((idx + 0.75) / (len(labels) - 1), 0.56),
                xytext=((idx + 0.25) / (len(labels) - 1), 0.56),
                arrowprops={"arrowstyle": "->", "color": "#27323a", "lw": 1.3},
            )
    ax.set_xlim(-0.08, 1.08)
    ax.set_ylim(0.0, 1.0)
    ax.set_title("Configured processing chain", color="#27323a", pad=10)
    _save(fig, "pipeline_configured_chain.png")


def _compact_summary(text: str) -> str:
    return text.replace("→", "->")


def main() -> int:
    functions = _import_pycsamt()
    if BUILD_DIR.exists():
        shutil.rmtree(BUILD_DIR)
    CONFIG_PATH.parent.mkdir(parents=True, exist_ok=True)
    CONFIG_PATH.write_text(CONFIG_TEXT, encoding="utf-8")

    functions["configure_pipe"](show_progress=False)
    pipe = functions["Pipeline"].from_yaml(CONFIG_PATH)
    survey = functions["read_edis"](
        DATA_DIR,
        recursive=False,
        strict=False,
        on_dup="replace",
        progress=False,
    )
    result = pipe.run(survey.collection)

    _workflow_plot(pipe)
    _step_status_plot(result)
    _artifact_plot(result)

    table = pipe.describe()[["label", "code", "category", "params"]]
    compact_tree = []
    for path in sorted(result.outdir.rglob("*")):
        rel = path.relative_to(result.outdir)
        if len(rel.parts) == 1 or rel.parts[0] in {"plots", "processed"}:
            compact_tree.append(rel)
        if len(compact_tree) >= 24:
            break

    print("config_path:")
    print(CONFIG_PATH.relative_to(ROOT))
    print("describe_table:")
    print(table.to_string())
    print("survey_summary:")
    print(survey.summary())
    print("result_summary:")
    print(_compact_summary(result.summary()))
    print("step_summaries:")
    for sr in result.step_results:
        print(_compact_summary(sr.summary_line()))
    print("result_fields:")
    print(f"ok={result.ok}")
    print(f"n_errors={result.n_errors}")
    print(f"outdir={result.outdir}")
    print(f"processed_count={len(result.processed_paths)}")
    print(f"plot_count={len(result.plots)}")
    print("output_tree:")
    for rel in compact_tree:
        print(rel)
    print(f"images: {IMAGE_DIR.relative_to(ROOT)}")

    shutil.rmtree(BUILD_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
