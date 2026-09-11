"""Generate figures for the UncertaintyCalibrator user guide.

Run from the repository root::

    python docs/scripts/generate_user_guide_ai_processing_uncertainty_figures.py \
        --save-dir docs/source/images/user_guide/ai_processing

Exercises :class:`~pycsamt.ai.processing.uncertainty.UncertaintyCalibrator`
on pyCSAMT's bundled AMT line ``data/AMT/WILLY_DATA/L18PLT`` (28
stations) -- its recalibration target is the field ``z_err`` already
present on every station there, with no dataset-specific requirement
the way :class:`~pycsamt.ai.processing.imputer.EMImputer` needs real
gaps.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from pycsamt.ai.processing import (
    UncertaintyCalibrator,
    build_uncertainty_features_table,
    plot_uncertainty_summary,
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


def run_uncertainty(sites, save_dir: Path) -> None:
    print("\n=== UncertaintyCalibrator ===")
    feats = build_uncertainty_features_table(sites)
    print(f"  feature rows: {len(feats)}")
    print(f"  z_err_frac stats: min={feats['z_err_frac'].min():.4f} "
          f"max={feats['z_err_frac'].max():.4f} "
          f"median={feats['z_err_frac'].median():.4f}")

    # Held-out cells (not stations) for a quantitative RMSE / R^2 check,
    # the same row-level split UncertaintyCalibrator.fit() itself uses
    # internally for its own validation loss.
    rng = np.random.default_rng(0)
    n = len(feats)
    idx = rng.permutation(n)
    n_val = max(1, int(n * 0.2))
    val_idx, train_idx = idx[:n_val], idx[n_val:]
    train_df = feats.iloc[train_idx].reset_index(drop=True)
    val_df = feats.iloc[val_idx].reset_index(drop=True)

    cal_val = UncertaintyCalibrator(hidden=(32, 16))
    cal_val.fit(train_df, epochs=150, seed=0, verbose=False)
    y_true = val_df["z_err_frac"].to_numpy()
    y_pred = cal_val.transform(val_df)
    rmse = float(np.sqrt(np.mean((y_true - y_pred) ** 2)))
    ss_res = float(np.sum((y_true - y_pred) ** 2))
    ss_tot = float(np.sum((y_true - y_true.mean()) ** 2)) + 1e-24
    r2 = 1.0 - ss_res / ss_tot
    print(f"  held-out ({len(val_df)} cells): RMSE={rmse:.4f} R^2={r2:.3f}")

    baseline = np.full_like(y_true, train_df["z_err_frac"].median())
    rmse_baseline = float(np.sqrt(np.mean((y_true - baseline) ** 2)))
    print(f"  RMSE (median baseline): {rmse_baseline:.4f}")

    # Deployed calibrator: fit on the full feature table.
    cal = UncertaintyCalibrator(hidden=(32, 16))
    cal.fit(feats, epochs=150, seed=0, verbose=False)
    print(f"  {cal!r}")
    hist = cal.history_
    if hist.get("val_loss"):
        print(f"  final val_loss: {hist['val_loss'][-1]:.5f}")

    table = cal.predict_table(sites)
    print(
        table[
            ["station", "freq", "z_err_frac", "z_err_frac_calibrated"]
        ].head(5).to_string(index=False)
    )

    agg = table.groupby("station")[
        ["z_err_frac", "z_err_frac_calibrated"]
    ].median()
    print("  median fractional error per station (head):")
    print(agg.head(5).to_string())

    fig = plot_uncertainty_summary(
        table, y_true, y_pred, history=cal,
        suptitle="UncertaintyCalibrator -- L18PLT",
    )
    _save(fig, save_dir, "uncertainty_summary.png")

    # apply(): Z.z untouched, Z.z_err rescaled.
    calibrated_sites = cal.apply(sites, inplace=False)
    from pycsamt.emtools._core import _get_z_block, _iter_items

    ed_in = next(_iter_items(sites))
    ed_out = next(_iter_items(calibrated_sites))
    _, z_in, _fr_in, ze_in = _get_z_block(ed_in, with_errors=True)
    _, z_out, _fr_out, ze_out = _get_z_block(ed_out, with_errors=True)
    same_z = bool(np.allclose(z_in, z_out, equal_nan=True))
    changed_err = not bool(np.allclose(ze_in, ze_out, equal_nan=True))
    print(f"  apply(): Z.z unchanged={same_z}, Z.z_err changed={changed_err}")
    print(f"  z_err before[:3, 0, 1]: {ze_in[:3, 0, 1]}")
    print(f"  z_err after  [:3, 0, 1]: {ze_out[:3, 0, 1]}")


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
    run_uncertainty(sites, args.save_dir)


if __name__ == "__main__":
    main()
