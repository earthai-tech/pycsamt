"""Generate figures for the DistortionTypeClassifier user guide.

Run from the repository root::

    python docs/scripts/generate_user_guide_ai_processing_distortion_figures.py \
        --save-dir docs/source/images/user_guide/ai_processing

Exercises :class:`~pycsamt.ai.processing.distortion.DistortionTypeClassifier`
on pyCSAMT's bundled AMT line ``data/AMT/WILLY_DATA/L18PLT`` (28
stations, the same line used throughout :mod:`pycsamt.emtools.qc`'s
user guide), then routes the predicted regimes to the real
correction tools -- :func:`~pycsamt.emtools.ss.correct_ss_ama` and
:func:`~pycsamt.emtools.gb.apply_groom_bailey` -- to demonstrate the
triage value proposition end to end.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from pycsamt.ai.processing import (
    DistortionTypeClassifier,
    build_distortion_features_table,
    plot_distortion_feature_space,
    plot_distortion_summary,
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


def run_distortion(sites, save_dir: Path) -> None:
    print("\n=== DistortionTypeClassifier ===")
    feats = build_distortion_features_table(sites)
    print(f"  station rows: {len(feats)}")
    print(f"  delta_log10_rho: min={feats['delta_log10_rho'].min():.3f} "
          f"max={feats['delta_log10_rho'].max():.3f} "
          f"std={feats['delta_log10_rho'].std():.3f}")
    print(f"  twist_deg: min={feats['twist_deg'].min():.1f} "
          f"max={feats['twist_deg'].max():.1f} "
          f"std={feats['twist_deg'].std():.1f}")
    print(f"  shear: min={feats['shear'].min():.3f} "
          f"max={feats['shear'].max():.3f} "
          f"std={feats['shear'].std():.3f}")

    from pycsamt.ai.processing.distortion import (
        _DISTORTION_LABELS,
        _rule_labels,
    )

    y_rule = _rule_labels(
        feats["delta_log10_rho"].to_numpy(),
        feats["twist_deg"].to_numpy(),
        feats["shear"].to_numpy(),
    )
    print("  rule-based label counts (deterministic self-training "
          "target):")
    import pandas as pd
    print(pd.Series([_DISTORTION_LABELS[i] for i in y_rule])
          .value_counts().to_string())

    clf = DistortionTypeClassifier.from_features_table(
        feats, epochs=200, seed=0, verbose=False,
    )
    print(f"  {clf!r}")
    hist = clf.history_
    if hist.get("val_loss"):
        best_ep = hist["val_loss"].index(min(hist["val_loss"])) + 1
        print(f"  final val_loss: {hist['val_loss'][-1]:.5f}, "
              f"best epoch: {best_ep}/{len(hist['val_loss'])}")

    table = clf.predict_table(sites)
    print(table[["station", "regime_label", "confidence"]]
          .head(5).to_string(index=False))
    print("  regime counts:")
    print(table["regime_label"].value_counts().to_string())

    fig = plot_distortion_summary(
        table, history=clf,
        suptitle="DistortionTypeClassifier -- L18PLT",
    )
    _save(fig, save_dir, "distortion_summary.png")

    ax = plot_distortion_feature_space(
        table, title="Feature space -- delta_log10_rho vs. twist_deg",
    )
    _save(ax.figure, save_dir, "distortion_feature_space.png")

    # Triage in action: route each regime to the real correction tool
    # the classifier is meant to point at, not a new one of its own.
    from pycsamt.emtools.gb import apply_groom_bailey
    from pycsamt.emtools.ss import correct_ss_ama

    ss_only = table.loc[
        table["regime_label"] == "static_shift_only", "station"
    ].tolist()
    distorted = table.loc[
        table["regime_label"] == "distorted", "station"
    ].tolist()
    print(f"  routed to correct_ss_ama (static_shift_only): {len(ss_only)}")
    print(f"  routed to apply_groom_bailey (distorted): {len(distorted)}")

    if ss_only:
        ss_sites = sites.select(names=ss_only)
        ss_corrected = correct_ss_ama(ss_sites, inplace=False)
        print(f"  correct_ss_ama: {len(list(ss_sites))} -> "
              f"{len(list(ss_corrected))} stations corrected")

    if distorted:
        gb_sites = sites.select(names=distorted)
        gb_corrected = apply_groom_bailey(gb_sites, inplace=False)
        print(f"  apply_groom_bailey: {len(list(gb_sites))} -> "
              f"{len(list(gb_corrected))} stations corrected")


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
    run_distortion(sites, args.save_dir)


if __name__ == "__main__":
    main()
