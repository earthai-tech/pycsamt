"""Generate executed figures for user_guide/forward/maxwell_caching_and_batch.rst."""

from __future__ import annotations

import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.forward.maxwell import MaxwellMesh, MaxwellProblem, ReceiverSet
from pycsamt.forward.maxwell.mt2d import MT2DAdapter
from pycsamt.forward.maxwell.modem3d import ModEm3DAdapter
from pycsamt.forward.maxwell.batch import solve_batch, BatchPolicy

IMAGE_DIR = ROOT / "docs/source/images/user_guide/forward"


def _save(fig: plt.Figure, name: str) -> None:
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(IMAGE_DIR / name, dpi=190, bbox_inches="tight")
    plt.close(fig)


def _timed_batch(problems, backend, max_workers: int) -> float:
    start = time.perf_counter()
    report = solve_batch(problems, backend, policy=BatchPolicy(max_workers=max_workers))
    elapsed = time.perf_counter() - start
    assert len(report.solved) == len(problems)
    return elapsed


def make_concurrency_comparison() -> dict[str, float]:
    mt2d_mesh = MaxwellMesh(np.linspace(0, 10_000, 41), np.linspace(0, 5_000, 31))
    mt2d_adapter = MT2DAdapter(verbose=False)
    mt2d_problems = [
        MaxwellProblem(
            mt2d_mesh, np.full(mt2d_mesh.shape, 0.01), np.geomspace(1000, 1, 8),
            ReceiverSet([[x, 0.0]], [f"S{i:02d}"]), ("zxy", "zyx"),
        )
        for i, x in enumerate(np.linspace(1000, 9000, 8))
    ]

    n, h = 8, 300.0
    modem_mesh = MaxwellMesh(np.arange(n + 1) * h, np.arange(11) * h, np.arange(n + 1) * h)
    modem_adapter = ModEm3DAdapter()
    modem_problems = [
        MaxwellProblem(
            modem_mesh, np.full(modem_mesh.shape, 1.0 / 100.0), [1.0, 0.5],
            ReceiverSet([[n * h / 2, n * h / 2, 0.0]], [f"S{i:02d}"]), ("zxy", "zyx"),
            metadata={"case": i},
        )
        for i in range(6)
    ]

    timings = {
        "mt2d_seq": _timed_batch(mt2d_problems, mt2d_adapter, 1),
        "mt2d_par": _timed_batch(mt2d_problems, mt2d_adapter, 4),
        "modem3d_seq": _timed_batch(modem_problems, modem_adapter, 1),
        "modem3d_par": _timed_batch(modem_problems, modem_adapter, 4),
    }

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    groups = ["MT2DAdapter\n(8 in-process solves)", "ModEm3DAdapter\n(6 external-process solves)"]
    sequential = [timings["mt2d_seq"], timings["modem3d_seq"]]
    parallel = [timings["mt2d_par"], timings["modem3d_par"]]
    x = np.arange(len(groups))
    width = 0.35
    ax.bar(x - width / 2, sequential, width, label="max_workers=1", color="#457b9d")
    ax.bar(x + width / 2, parallel, width, label="max_workers=4", color="#e76f51")
    for i, (s, p) in enumerate(zip(sequential, parallel)):
        ax.annotate(f"{s / p:.2f}x", (x[i], max(s, p) * 1.03), ha="center", fontsize=10)
    ax.set_xticks(x)
    ax.set_xticklabels(groups)
    ax.set_ylabel("wall-clock time (s)")
    ax.set_title("solve_batch concurrency: in-process GIL vs. external subprocesses")
    ax.legend()
    _save(fig, "maxwell_caching_and_batch_concurrency.png")
    return timings


def main() -> int:
    timings = make_concurrency_comparison()
    for key, value in timings.items():
        print(key, round(value, 3))
    print("figures generated: 1")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
