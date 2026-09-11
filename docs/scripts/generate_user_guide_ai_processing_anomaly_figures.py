"""Generate figures for the AnomalyDetector user guide.

Run from the repository root::

    python docs/scripts/generate_user_guide_ai_processing_anomaly_figures.py \
        --save-dir docs/source/images/user_guide/ai_processing

Exercises :class:`~pycsamt.ai.processing.anomaly.AnomalyDetector` on
pyCSAMT's bundled AMT line ``data/AMT/WILLY_DATA/L18PLT`` (28 stations,
the same line used throughout :mod:`pycsamt.emtools.qc`'s user guide).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from pycsamt.ai.processing import (
    AnomalyDetector,
    plot_anomaly_summary,
    plot_training_history,
    prepare_z_features,
)
from pycsamt.emtools import ensure_sites


def _repository_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _save(fig, save_dir: Path, name: str) -> None:
    save_dir.mkdir(parents=True, exist_ok=True)
    path = save_dir / name
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"  saved {path.relative_to(_repository_root())}")


def station_names(sites) -> list[str]:
    """Real station IDs in the same iteration order as prepare_z_features."""
    from pycsamt.emtools._core import _iter_items, _name

    return [_name(ed, i) for i, ed in enumerate(_iter_items(sites))]


def load_l18(verbose: int = 1):
    data_dir = _repository_root() / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
    sites = ensure_sites(
        data_dir, recursive=True, on_dup="replace", strict=True,
        verbose=verbose,
    )
    return sites


def run_anomaly(sites, save_dir: Path) -> None:
    print("\n=== AnomalyDetector ===")
    X = prepare_z_features(sites, n_components=4)
    labels = station_names(sites)
    n_sites, n_comp, n_freq = X.shape
    X_flat = X.reshape(n_sites, n_comp * n_freq)
    X_flat = np.nan_to_num(X_flat, nan=0.0)
    print(f"  feature matrix: {X_flat.shape}")

    det = AnomalyDetector(latent_dim=8, channels=(32, 16))
    det.fit(X_flat, epochs=150, seed=0, verbose=False)
    scores = det.transform(X_flat)
    flags = det.flag_anomalies(X_flat)
    print(f"  {det!r}")
    print(f"  threshold_ (95th pct): {det.threshold_:.4f}")
    print(f"  flagged anomalous: {int(flags.sum())} / {n_sites}")
    if flags.any():
        flagged = [labels[i] for i in np.where(flags)[0]]
        print(f"  stations: {flagged}")

    fig = plot_anomaly_summary(
        scores, station_labels=labels, threshold=det.threshold_,
        suptitle="AnomalyDetector -- L18PLT profile-level scores",
    )
    _save(fig, save_dir, "anomaly_summary.png")

    ax = plot_training_history(det, title="AnomalyDetector training")
    _save(ax.figure, save_dir, "anomaly_training.png")


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
    run_anomaly(sites, args.save_dir)


if __name__ == "__main__":
    main()
