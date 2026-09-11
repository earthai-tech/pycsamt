"""Generate figures for the DimensionalityClassifier user guide.

Run from the repository root::

    python docs/scripts/generate_user_guide_ai_processing_classify_figures.py \
        --save-dir docs/source/images/user_guide/ai_processing

Exercises :class:`~pycsamt.ai.processing.classify.DimensionalityClassifier`
on pyCSAMT's bundled AMT line ``data/AMT/WILLY_DATA/L18PLT`` (28
stations, the same line used throughout :mod:`pycsamt.emtools.qc`'s
user guide).
"""

from __future__ import annotations

import argparse
from pathlib import Path

from pycsamt.ai.processing import (
    DimensionalityClassifier,
    plot_dimensionality_summary,
    plot_training_history,
)
from pycsamt.emtools import ensure_sites
from pycsamt.emtools.dimensionality import phase_features_table
from pycsamt.emtools.strike import estimate_strike_phase_tensor


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


def run_classify(sites, save_dir: Path) -> None:
    print("\n=== DimensionalityClassifier ===")
    feats = phase_features_table(sites)
    print(f"  feature rows: {len(feats)}")
    print(f"  stations: {feats['station'].nunique()}")

    # Classical, per-station phase-tensor strike (emtools.strike) supplies
    # the strike-head training target.  Without it, from_features_table's
    # default self-training path (rule-based dimensionality labels only)
    # never supervises the strike head at all -- see the note below.
    consensus = estimate_strike_phase_tensor(sites)
    feats = feats.merge(
        consensus[["station", "ang"]], on="station", how="left"
    ).rename(columns={"ang": "strike_target"})
    print("  classical (phase-tensor) strike per station -- head/tail:")
    print(
        consensus[["station", "ang", "iqr"]]
        .head(3)
        .to_string(index=False)
    )

    clf = DimensionalityClassifier.from_features_table(
        feats, strike_col="strike_target", epochs=120, seed=0,
        verbose=False,
    )
    print(f"  {clf!r}")
    hist = clf.history_
    if hist.get("val_loss"):
        print(f"  final val_loss: {hist['val_loss'][-1]:.5f}")

    table = clf.predict_table(sites)
    counts = table["dim_label"].value_counts()
    print("  predicted class counts:")
    print(counts.to_string())
    n2d = int((table["dim"] == 1).sum())
    if n2d:
        strike_mean = table.loc[table["dim"] == 1, "strike"].mean()
        strike_std = table.loc[table["dim"] == 1, "strike"].std()
        print(f"  predicted strike over {n2d} 2-D samples: "
              f"{strike_mean:.1f} +/- {strike_std:.1f} deg")
        print(f"  classical consensus (all stations): "
              f"{consensus['ang'].mean():.1f} deg")

    fig = plot_dimensionality_summary(
        table, suptitle="DimensionalityClassifier -- L18PLT"
    )
    _save(fig, save_dir, "dimensionality_summary.png")

    ax = plot_training_history(
        clf, title="DimensionalityClassifier training"
    )
    _save(ax.figure, save_dir, "dimensionality_training.png")


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
    run_classify(sites, args.save_dir)


if __name__ == "__main__":
    main()
