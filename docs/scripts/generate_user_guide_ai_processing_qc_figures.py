"""Generate figures for the EMQCScorer user guide.

Run from the repository root::

    python docs/scripts/generate_user_guide_ai_processing_qc_figures.py \
        --save-dir docs/source/images/user_guide/ai_processing

Exercises :class:`~pycsamt.ai.processing.qc.EMQCScorer` on pyCSAMT's
bundled AMT line ``data/AMT/WILLY_DATA/L18PLT`` (28 stations, real
impedance-error tensors -- the same line used by
:mod:`pycsamt.emtools.qc`'s user guide).
"""

from __future__ import annotations

import argparse
from pathlib import Path

from pycsamt.ai.processing import (
    EMQCScorer,
    plot_qc_feature_heatmap,
    plot_qc_heatmap,
    plot_qc_summary,
)
from pycsamt.emtools import ensure_sites


def _repository_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _save(fig, save_dir: Path, name: str) -> None:
    save_dir.mkdir(parents=True, exist_ok=True)
    path = save_dir / name
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"  saved {path.relative_to(_repository_root())}")


def load_l18(verbose: int = 1):
    data_dir = _repository_root() / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
    sites = ensure_sites(
        data_dir, recursive=True, on_dup="replace", strict=True,
        verbose=verbose,
    )
    return sites


def run_qc(sites, save_dir: Path) -> None:
    print("\n=== EMQCScorer ===")
    scorer = EMQCScorer(random_state=0)
    scorer.fit(sites)
    table = scorer.score_table(sites)
    print(f"  rows: {len(table)}")
    print(f"  stations: {table['station'].nunique()}")
    print(f"  flagged bad: {(table['flag'] == 0).sum()} "
          f"/ {len(table)}")
    print(table[["station", "freq", "score", "flag"]].head(5)
          .to_string(index=False))

    fig = plot_qc_summary(
        table, score_threshold=0.5,
        suptitle="EMQCScorer -- L18PLT quality scores",
    )
    _save(fig, save_dir, "qc_summary.png")

    ax = plot_qc_heatmap(table, title="L18PLT -- QC score heat-map")
    _save(ax.figure, save_dir, "qc_heatmap.png")

    fig = plot_qc_feature_heatmap(
        table, title="L18PLT -- QC feature diagnostics"
    )
    _save(fig, save_dir, "qc_feature_heatmap.png")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--save-dir",
        type=Path,
        default=_repository_root()
        / "docs"
        / "source"
        / "images"
        / "user_guide"
        / "ai_processing",
    )
    args = parser.parse_args()

    import matplotlib

    matplotlib.use("Agg")

    sites = load_l18()
    run_qc(sites, args.save_dir)


if __name__ == "__main__":
    main()
