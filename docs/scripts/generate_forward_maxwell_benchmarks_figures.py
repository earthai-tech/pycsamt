"""Generate executed figures for user_guide/forward/maxwell_benchmarks.rst."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.forward.maxwell import MaxwellMesh, ReceiverSet
from pycsamt.forward.maxwell.mt2d import MT2DAdapter
from pycsamt.forward.maxwell.benchmarks import (
    BenchmarkThresholds,
    half_space_benchmark,
    layered_earth_benchmark,
    run_benchmarks,
)

IMAGE_DIR = ROOT / "docs/source/images/user_guide/forward"


def _save(fig: plt.Figure, name: str) -> None:
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(IMAGE_DIR / name, dpi=190, bbox_inches="tight")
    plt.close(fig)


def _utilization(outcome, thresholds: BenchmarkThresholds) -> tuple[float, float, float]:
    metrics = outcome.metrics
    return (
        metrics.normalized_rms / thresholds.maximum_normalized_rms,
        metrics.maximum_amplitude_relative_error
        / thresholds.maximum_amplitude_relative_error,
        metrics.maximum_phase_error_deg / thresholds.maximum_phase_error_deg,
    )


def make_benchmark_scorecard() -> tuple[bool, float]:
    dz = 20.0 * 1.25 ** np.arange(30)
    z_edges = np.concatenate([[0.0], np.cumsum(dz)])
    x_edges = np.linspace(0, 180_000, 19)
    mesh = MaxwellMesh(x_edges, z_edges)
    stations = ReceiverSet([[60_000.0, 0.0], [120_000.0, 0.0]], ["S00", "S01"])
    adapter = MT2DAdapter(verbose=False)

    default_thresholds = BenchmarkThresholds()
    half = half_space_benchmark(mesh, stations, [10.0, 1.0], resistivity_ohm_m=150.0)
    interface1, interface2 = z_edges[10], z_edges[18]
    layered = layered_earth_benchmark(
        mesh, stations, np.geomspace(10, 0.5, 5),
        [150.0, 20.0, 800.0], [interface1, interface2 - interface1],
    )
    report = run_benchmarks(adapter, [half, layered])

    strict_thresholds = BenchmarkThresholds(maximum_phase_error_deg=0.1)
    strict_half = half_space_benchmark(
        mesh, stations, [10.0, 1.0], resistivity_ohm_m=150.0,
        thresholds=strict_thresholds,
    )
    strict_outcome = strict_half.evaluate(adapter.solve(strict_half.problem))

    cases = [
        ("half-space\n(default limits)", report.outcomes[0], default_thresholds),
        ("layered-earth\n(default limits)", report.outcomes[1], default_thresholds),
        ("half-space\n(phase limit 0.1 deg)", strict_outcome, strict_thresholds),
    ]
    metric_labels = ["NRMS", "amplitude", "phase"]
    colours = ["#457b9d", "#2a9d8f", "#e76f51"]

    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    width = 0.24
    x = np.arange(len(cases))
    for mi, (label, colour) in enumerate(zip(metric_labels, colours)):
        heights = [_utilization(outcome, thresholds)[mi] for _, outcome, thresholds in cases]
        ax.bar(x + (mi - 1) * width, heights, width, label=label, color=colour)
    ax.axhline(1.0, color="black", ls="--", lw=1.3)
    ax.text(len(cases) - 0.5, 1.12, "pass / fail boundary", ha="right", fontsize=9)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([label for label, _, _ in cases])
    ax.set_ylabel("achieved metric / threshold (log scale)")
    ax.set_title("Benchmark scorecard: threshold utilization, not just pass/fail")
    ax.legend(title="metric", fontsize=9)
    _save(fig, "maxwell_benchmarks_scorecard.png")
    return report.passed, float(_utilization(strict_outcome, strict_thresholds)[2])


def main() -> int:
    passed, strict_phase_utilization = make_benchmark_scorecard()
    print("default-threshold report passed:", passed)
    print("strict phase utilization:", round(strict_phase_utilization, 3))
    print("figures generated: 1")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
