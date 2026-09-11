"""Generate figures for the EMDenoiser user guide.

Run from the repository root::

    python docs/scripts/generate_user_guide_ai_processing_denoise_figures.py \
        --save-dir docs/source/images/user_guide/ai_processing

Exercises :class:`~pycsamt.ai.processing.denoise.EMDenoiser` on
pyCSAMT's bundled AMT line ``data/AMT/WILLY_DATA/L18PLT`` (28 stations,
the same line used throughout :mod:`pycsamt.emtools.qc`'s user guide).

The denoiser has no "known-clean" field twin to validate against, so
this script follows the standard denoising-autoencoder protocol: the
28 real station spectra are treated as the clean reference, a fixed
``numpy.random.default_rng(0)`` Gaussian-noise corruption is added on
top, and the trained network's reconstruction is compared against both
the corrupted input and the real clean reference.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from pycsamt.ai.processing import (
    EMDenoiser,
    plot_denoise_spectra,
    plot_denoise_summary,
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


def run_denoise(sites, save_dir: Path) -> None:
    print("\n=== EMDenoiser ===")
    freq = None
    X_clean = prepare_z_features(sites, n_components=4)
    labels = station_names(sites)
    print(f"  X_clean shape: {X_clean.shape}")

    rng = np.random.default_rng(0)
    noise = rng.normal(
        scale=0.35 * X_clean.std(axis=(0, 2), keepdims=True),
        size=X_clean.shape,
    ).astype("float32")
    X_noisy = (X_clean + noise).astype("float32")

    den = EMDenoiser(channels=(32, 64, 32))
    den.fit(X_clean, noise_level=0.15, epochs=80, seed=0, verbose=False)
    X_den = den.transform(X_noisy)
    print(f"  {den!r}")
    hist = den.history_
    if hist.get("val_loss"):
        print(f"  final val_loss: {hist['val_loss'][-1]:.5f}")
        print(f"  best val_loss:  {min(hist['val_loss']):.5f}")

    rmse_noisy = float(np.sqrt(np.mean((X_noisy - X_clean) ** 2)))
    rmse_den = float(np.sqrt(np.mean((X_den - X_clean) ** 2)))
    print(f"  RMSE vs. clean reference -- noisy:    {rmse_noisy:.4f}")
    print(f"  RMSE vs. clean reference -- denoised: {rmse_den:.4f}")
    for name, ci in [
        ("log|Zxy|", 0), ("phase_xy", 1), ("log|Zyx|", 2), ("phase_yx", 3),
    ]:
        rn = float(np.sqrt(np.mean((X_noisy[:, ci] - X_clean[:, ci]) ** 2)))
        rd = float(np.sqrt(np.mean((X_den[:, ci] - X_clean[:, ci]) ** 2)))
        print(f"    {name:>9s}: RMSE noisy={rn:.3f} denoised={rd:.3f}")

    from pycsamt.ai.processing.plot import _roughness
    r_noisy = _roughness(X_noisy)
    r_den = _roughness(X_den)
    pct = 100.0 * (r_noisy - r_den) / (r_noisy + 1e-24)
    pct_site = np.nanmedian(pct, axis=1)
    print(f"  roughness reduction per station: "
          f"min={pct_site.min():.1f}% max={pct_site.max():.1f}%")

    # frequency grid: reuse the site's own grid (first site) for the x-axis
    from pycsamt.emtools._core import _get_z_block, _iter_items

    for ed in _iter_items(sites):
        result = _get_z_block(ed, with_errors=False)
        freq = result[2] if len(result) >= 3 else None
        if freq is not None:
            break

    fig = plot_denoise_summary(
        freq, X_noisy, X_den, history=den, station_labels=labels,
        n_show=3, suptitle="EMDenoiser -- L18PLT (synthetic noise +35%)",
    )
    _save(fig, save_dir, "denoise_summary.png")

    fig = plot_denoise_spectra(
        freq, X_noisy, X_den, station_labels=labels, n_show=4,
        suptitle="Before / after spectra -- 4 L18PLT stations",
    )
    _save(fig, save_dir, "denoise_spectra.png")


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
    run_denoise(sites, args.save_dir)


if __name__ == "__main__":
    main()
