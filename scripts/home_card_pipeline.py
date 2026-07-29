"""Home-page card image: "Pipelines, agents & apps".

Left panel: a real :class:`pycsamt.pipeline.Pipeline` run (the same 5-step
chain as the pipeline gallery) on the bundled demo sites, drawn as a
per-step timing bar chart from the returned
:class:`~pycsamt.pipeline.PipelineResult`.

Right panel: a real :class:`pycsamt.agents.AgentCoordinator` **dry-run**
preview (fully offline) of a context -> load -> QC agent chain, drawn as a
flow of steps exactly as the agents gallery does.

Output: docs/source/_static/images/home/card-pipeline.png

Usage (any cwd):
    python scripts/home_card_pipeline.py
"""

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

os.environ.setdefault("PYCSAMT_DOCS_BUILD", "1")
os.environ.setdefault("PYCSAMT_DOCS_REPO_ROOT", str(ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "docs" / "examples" / "pipeline"))

from _pipe_data import demo_sites, quiet_logs, scratch_dir  # noqa: E402

from pycsamt.agents import (  # noqa: E402
    AgentCoordinator,
    ContextInputAgent,
    DataQCAgent,
    MTLoaderAgent,
)
from pycsamt.api.agents import AGENT_CONFIG  # noqa: E402
from pycsamt.pipeline import Pipeline, Step, configure_pipe  # noqa: E402

BLUE, ORANGE, GOLD, SLATE = "#3e65b0", "#f15a29", "#d99114", "#5c677d"
INK = "#1c2941"
GREEN = "#2e9e5b"

STYLE = {
    "font.size": 9,
    "axes.linewidth": 0.8,
    "axes.edgecolor": "#c8d0dc",
    "xtick.color": "#5c677d",
    "ytick.color": "#5c677d",
    "axes.labelcolor": INK,
    "axes.titlesize": 10.5,
    "axes.titleweight": "bold",
    "axes.titlecolor": INK,
}


def run_pipeline():
    configure_pipe(show_progress=False, plot_dpi=72, plot_fmt="png")
    pipe = Pipeline(
        [
            ("drop_dup", Step("FREQ002")),
            ("select_band", Step("FREQ001", band_hz=(0.01, 10_000))),
            ("align", Step("FREQ004")),
            ("notch", Step("NR001", mains_hz=50)),
            ("qc_snap", Step("QC001")),
        ]
    )
    with quiet_logs():
        result = pipe.run(demo_sites(), outdir=str(scratch_dir()))
    return result


def agent_preview():
    with AGENT_CONFIG.offline():
        coordinator = AgentCoordinator("home_preview", verbose=False)
        coordinator.add_step(
            "context",
            ContextInputAgent(),
            description="parse the request",
        )
        coordinator.add_step(
            "load",
            MTLoaderAgent(),
            input_fn=lambda r: {
                "path": (r["context"].get("config") or {}).get("data_path", "")
            },
            description="load the survey",
        )
        coordinator.add_step(
            "qc",
            DataQCAgent(),
            input_fn=lambda r: {"sites": r["load"].get("sites")},
            description="flag noisy stations",
        )
        preview = coordinator.execute(
            {"request": "QC the EDI line and report"}, dry_run=True
        )
    return list(preview["steps"])


def main():
    plt.rcParams.update(STYLE)
    fig, (ax_l, ax_r) = plt.subplots(
        1,
        2,
        figsize=(6.4, 3.0),
        dpi=200,
        gridspec_kw={"width_ratios": [1.3, 1.0]},
    )
    fig.patch.set_facecolor("white")

    # -- pipeline: per-step timings ------------------------------------------
    result = run_pipeline()
    steps = list(result.step_results)
    labels = [sr.step_name for sr in steps]
    times = np.array([sr.elapsed_sec for sr in steps])
    colors = [BLUE, ORANGE, GOLD, SLATE, "#7a6fb3"][: len(steps)]
    y = np.arange(len(steps))[::-1]

    ax_l.barh(y, times, height=0.62, color=colors, zorder=3)
    for yi, sr, t in zip(y, steps, times):
        ax_l.text(
            t + times.max() * 0.03,
            yi,
            f"✓ {sr.n_sites_in}→{sr.n_sites_out} sites",
            va="center",
            fontsize=6.8,
            color=GREEN if sr.ok else ORANGE,
            fontweight="bold",
        )
    ax_l.set_yticks(y)
    ax_l.set_yticklabels(labels, fontsize=7.5)
    ax_l.set_xlim(0, times.max() * 1.45)
    ax_l.set_xlabel("elapsed (s)", fontsize=8.5)
    ax_l.set_title(
        f"Pipeline.run — {len(steps)} steps · {result.elapsed_sec:.1f}s",
        pad=6,
    )
    ax_l.grid(True, axis="x", ls=":", lw=0.4, color="#c8d0dc", alpha=0.7)
    ax_l.tick_params(labelsize=7.5, length=2.5)

    # -- agents: offline dry-run chain ---------------------------------------
    chain = agent_preview()
    ax_r.axis("off")
    ax_r.set_xlim(0, 1)
    ax_r.set_ylim(0, 1)
    ax_r.set_title("agents — dry-run plan", pad=6)

    n = len(chain)
    box_h = 0.16
    gap = (0.9 - n * box_h) / max(n - 1, 1)
    for i, step in enumerate(chain):
        y0 = 0.92 - box_h - i * (box_h + gap)
        offline = str(step.get("llm", "")).lower().startswith("no")
        edge = BLUE if offline else ORANGE
        face = "#eef3fb" if offline else "#fdf1ea"
        ax_r.add_patch(
            mpatches.FancyBboxPatch(
                (0.12, y0),
                0.76,
                box_h,
                boxstyle="round,pad=0.012,rounding_size=0.025",
                facecolor=face,
                edgecolor=edge,
                linewidth=1.4,
                zorder=3,
            )
        )
        ax_r.text(
            0.5,
            y0 + box_h * 0.66,
            f"{step['step']}. {step['name']}",
            ha="center",
            va="center",
            fontsize=8.2,
            fontweight="bold",
            color=edge,
        )
        ax_r.text(
            0.5,
            y0 + box_h * 0.28,
            f"{step['agent']} · {step['llm']}",
            ha="center",
            va="center",
            fontsize=6.4,
            color="#5c677d",
        )
        if i < n - 1:
            ax_r.annotate(
                "",
                xy=(0.5, y0 - gap + 0.012),
                xytext=(0.5, y0 - 0.012),
                arrowprops=dict(arrowstyle="-|>", color="#8a94a6", lw=1.4),
            )

    fig.tight_layout(pad=0.6, w_pad=1.0)
    out = ROOT / "docs/source/_static/images/home/card-pipeline.png"
    fig.savefig(out, dpi=200, facecolor="white", bbox_inches="tight")
    print(
        f"saved: {out} ({out.stat().st_size / 1024:.0f} KB)  "
        f"pipeline ok={result.ok}  agents={len(chain)} steps"
    )


if __name__ == "__main__":
    main()
