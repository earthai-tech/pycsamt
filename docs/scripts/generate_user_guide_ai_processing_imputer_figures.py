"""Generate figures for the EMImputer user guide.

Run from the repository root::

    python docs/scripts/generate_user_guide_ai_processing_imputer_figures.py \
        --save-dir docs/source/images/user_guide/ai_processing

Exercises :class:`~pycsamt.ai.processing.imputer.EMImputer` on
``data/MT/broken-hill`` (21 real "ultra-wide-band" MT soundings,
AlQahtani et al. 2026, CC-BY -- see that folder's ``README.md`` for the
required citation), because several of its stations have *genuinely*
missing frequency rows -- a real gap-filling problem, unlike pyCSAMT's
other bundled AMT line (``data/AMT/WILLY_DATA/L18PLT``), which has
none. A separate synthetic-mask validation (hiding a known fraction of
the *observed* Broken Hill cells) still supplies the quantitative
RMSE / R^2 check that real gaps, having no known ground truth, cannot.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from pycsamt.ai.processing import (
    EMImputer,
    plot_imputer_reconstruction,
    plot_imputer_summary,
    plot_imputer_validation,
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


def load_broken_hill(verbose: int = 1):
    data_dir = _repository_root() / "data" / "MT" / "broken-hill" / "edis"
    sites = ensure_sites(
        data_dir, recursive=True, on_dup="replace", strict=True,
        verbose=verbose,
    )
    return sites


def run_imputer(sites, save_dir: Path) -> None:
    print("\n=== EMImputer ===")
    labels = station_names(sites)
    X = prepare_z_features(sites, n_components=4)
    n_sites, _n_comp, n_freq = X.shape
    print(f"  X shape: {X.shape}")

    from pycsamt.emtools._core import _get_z_block, _iter_items

    freq = None
    for ed in _iter_items(sites):
        result = _get_z_block(ed, with_errors=False)
        freq = result[2] if len(result) >= 3 else None
        if freq is not None:
            break

    missing_mask = ~np.isfinite(X)
    row_missing = missing_mask.any(axis=1)  # (n_sites, n_freq)
    per_site = row_missing.sum(axis=1)
    gapped = [
        (labels[i], int(per_site[i]))
        for i in range(n_sites)
        if per_site[i] > 0
    ]
    print(f"  genuinely missing rows: {int(row_missing.sum())} "
          f"across {len(gapped)} / {n_sites} stations")
    for name, n in gapped:
        print(f"    {name}: {n} missing rows (of {n_freq})")

    gap_table = pd.DataFrame(
        {
            "station": np.repeat(labels, n_freq),
            "freq": np.tile(freq, n_sites),
            "missing": row_missing.astype(int).ravel(),
        }
    )

    # Quantitative validation: real gaps have no known ground truth, so
    # hold out a random 12% of the *genuinely observed* cells instead
    # (seeded), fit with them additionally hidden, then compare the
    # reconstruction against the real, known values at those cells.
    rng = np.random.default_rng(0)
    observed = np.isfinite(X)
    holdout = (rng.random(X.shape) < 0.12) & observed
    X_holdout = X.copy()
    X_holdout[holdout] = np.nan
    print(f"  held-out cells for validation: {int(holdout.sum())}")

    imp_val = EMImputer(channels=(32, 64, 32))
    imp_val.fit(
        X_holdout, mask_frac=0.15, epochs=150, seed=0, verbose=False
    )
    X_recon_holdout = imp_val.transform(X_holdout)
    y_true = X[holdout]
    y_pred = X_recon_holdout[holdout]
    rmse_model = float(np.sqrt(np.mean((y_true - y_pred) ** 2)))

    chan_mean = np.nanmean(X_holdout, axis=(0, 2), keepdims=True)
    y_pred_naive = np.broadcast_to(chan_mean, X.shape)[holdout]
    rmse_naive = float(np.sqrt(np.mean((y_true - y_pred_naive) ** 2)))
    print(f"  RMSE (EMImputer), pooled:     {rmse_model:.4f}")
    print(f"  RMSE (channel-mean baseline): {rmse_naive:.4f}")
    print("  (pooled RMSE mixes log-amplitude and degree-valued phase "
          "channels -- see per-channel breakdown below)")

    comp_ids = np.broadcast_to(
        np.arange(X.shape[1])[None, :, None], X.shape
    )[holdout]
    for ci, name in enumerate(
        ["log|Zxy|", "phase_xy", "log|Zyx|", "phase_yx"]
    ):
        cm = comp_ids == ci
        rm = float(np.sqrt(np.mean((y_true[cm] - y_pred[cm]) ** 2)))
        rb = float(
            np.sqrt(np.mean((y_true[cm] - y_pred_naive[cm]) ** 2))
        )
        print(f"    {name:>9s} (n={int(cm.sum())}): "
              f"RMSE model={rm:.3f}  baseline={rb:.3f}")

    # The deployed imputer: fit on the full array (only the genuine
    # gaps are NaN -- none of the synthetic holdout above), then fill
    # those real gaps via apply().
    imp = EMImputer(channels=(32, 64, 32))
    imp.fit(X, mask_frac=0.15, epochs=150, seed=0, verbose=False)
    print(f"  {imp!r}")
    hist = imp.history_
    if hist.get("val_loss"):
        print(f"  final val_loss: {hist['val_loss'][-1]:.5f}")

    X_filled = imp.transform(X)
    filled_sites = imp.apply(sites, inplace=False)

    # Verify: observed cells byte-identical, gapped station now filled.
    gap_name = gapped[0][0] if gapped else labels[0]
    gi = labels.index(gap_name)
    for i, (ed_in, ed_out) in enumerate(
        zip(_iter_items(sites), _iter_items(filled_sites))
    ):
        if i != gi:
            continue
        _, z_in, _fr_in = _get_z_block(ed_in, with_errors=False)[:3]
        _, z_out, _fr_out = _get_z_block(ed_out, with_errors=False)[:3]
        obs = ~np.isnan(z_in).any(axis=(1, 2))
        same_observed = bool(
            np.allclose(z_in[obs], z_out[obs], equal_nan=True)
        )
        was_missing = np.isnan(z_in).any(axis=(1, 2))
        offdiag_filled = int(
            (
                was_missing
                & ~np.isnan(z_out[:, 0, 1])
                & ~np.isnan(z_out[:, 1, 0])
            ).sum()
        )
        diag_still_nan = int(
            (
                was_missing
                & (np.isnan(z_out[:, 0, 0]) | np.isnan(z_out[:, 1, 1]))
            ).sum()
        )
        print(f"  {gap_name}: observed rows unchanged = {same_observed}; "
              f"of {int(was_missing.sum())} missing rows, "
              f"{offdiag_filled} now have Zxy/Zyx filled, "
              f"{diag_still_nan} still have Zxx/Zyy = NaN "
              f"(n_components=4 -> diagonal not modelled)")
        break

    fig = plot_imputer_summary(
        gap_table, freq, X, X_filled,
        history=imp, station_labels=labels,
        n_show=2,
        suptitle="EMImputer -- Broken Hill MT survey",
    )
    _save(fig, save_dir, "imputer_summary.png")

    show_idx = [labels.index(n) for n, _ in gapped[:4]] or [0]
    fig = plot_imputer_reconstruction(
        freq, X, X_filled, station_labels=labels, sites=show_idx,
        suptitle="Reconstructed cells -- Broken Hill stations with real gaps",
    )
    _save(fig, save_dir, "imputer_reconstruction.png")

    fig = plot_imputer_validation(
        y_true, y_pred, component_ids=comp_ids,
        component_labels=["log|Zxy|", "phase_xy", "log|Zyx|", "phase_yx"],
        suptitle="Held-out cell reconstruction, by component",
    )
    _save(fig, save_dir, "imputer_validation.png")


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

    sites = load_broken_hill()
    run_imputer(sites, args.save_dir)


if __name__ == "__main__":
    main()
