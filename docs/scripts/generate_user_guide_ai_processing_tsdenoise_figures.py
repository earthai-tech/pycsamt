"""Generate figures for the TimeSeriesDenoiser (MMF-SVM-K-SVD) user guide.

Run from the repository root::

    python docs/scripts/generate_user_guide_ai_processing_tsdenoise_figures.py \
        --save-dir docs/source/images/user_guide/ai_processing

Two complementary demonstrations, mirroring the structure of
:mod:`pycsamt.ai.processing.denoise`'s own user-guide script:

1. A fully synthetic record with a *known* clean reference (a
   band-limited colored-noise background with two injected transients
   of known shape, amplitude, and location) -- the only way to report
   a genuine before/after SNR / NCC number, since no real field
   recording carries an independent clean twin.
2. pyCSAMT's bundled real long-period MT time series
   ``data/MT/TS/kap103as.ts`` (station kap103, 5 channels, 5 s
   sampling, ~27 days -- the same file used throughout
   :doc:`../user_guide/transformers` and the IoT edge-QC examples),
   applied to a real ~580-unit interference burst in the EY channel
   (station background std ~1) to show the pipeline generalises beyond
   the synthetic case.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import matplotlib.pyplot as plt

from pycsamt.ai.processing import (
    KSVDDenoiser,
    SignalQualityClassifier,
    TimeSeriesDenoiser,
    compute_entropy_features,
    generate_synthetic_library,
    mmf_split,
    plot_ts_denoise_mmf_split,
    plot_ts_denoise_segments,
    plot_ts_denoise_summary,
    snr_db,
    time_domain_ncc,
)
from pycsamt.ts import TSData, read_ts


def _repository_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _save(fig, save_dir: Path, name: str) -> None:
    save_dir.mkdir(parents=True, exist_ok=True)
    path = save_dir / name
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"  saved {path.relative_to(_repository_root())}")


# ------------------------------------------------------------- synthetic


def make_synthetic_record(seed: int = 42):
    """A 30-minute, 1 Hz synthetic record with a known clean reference.

    The background is band-limited colored noise (weak, strongly
    random -- the qualitative description of a high-quality MT signal
    used throughout Gui et al., 2024), built by shaping white noise's
    spectrum as ``1/f^0.6`` and rescaling to a background standard
    deviation of 2.0 (arbitrary units). Two transients are injected on
    top: a charge/discharge-like exponential decay (amplitude 40,
    starting at t=600 s, decay constant 25 s) and a square-wave block
    (amplitude 30, t=1215-1245 s, deliberately offset from a
    classification-window boundary so the transition edges fall
    *inside* a window rather than exactly on its border).
    """
    rng = np.random.default_rng(seed)
    n = 1800
    dt = 1.0
    t = np.arange(n) * dt

    white = rng.normal(size=n)
    freqs = np.fft.rfftfreq(n)
    freqs[0] = freqs[1]
    spec = np.fft.rfft(white) / (freqs**0.6)
    clean = np.fft.irfft(spec, n=n)
    clean = clean / clean.std() * 2.0

    noise = np.zeros(n)
    c, tau = 600.0, 25.0
    env = 40.0 * np.exp(-np.maximum(t - c, 0.0) / tau)
    env[t < c] = 0.0
    noise += env
    noise[1215:1245] -= 30.0

    noisy = clean + noise
    return t, dt, clean, noisy


def run_synthetic_demo(save_dir: Path) -> None:
    print("\n=== TimeSeriesDenoiser -- synthetic record with known ground truth ===")
    t, dt, clean, noisy = make_synthetic_record()

    ts = TSData(data={"EX": noisy}, dt=dt, station="synthetic_demo")
    den = TimeSeriesDenoiser(mmf_size=121, win_seconds=60.0, random_state=0)
    out = den.apply(ts)
    denoised = out.get("EX")
    diag = den.diagnostics_

    n_noisy = int((diag["label"] == "noisy").sum())
    print(f"  {n_noisy}/{len(diag)} windows flagged noisy")
    print(diag.to_string(index=False))

    low, high = mmf_split(noisy, size=121)
    ax = plot_ts_denoise_segments(
        t, high, diag,
        title="MMF high-frequency residual -- SVM classification only "
              "(before K-SVD)",
    )
    _save(ax.figure, save_dir, "tsdenoise_segments_synthetic.png")

    print(f"  SNR  noisy    -> clean ref: {snr_db(clean, noisy):.2f} dB")
    print(f"  SNR  denoised -> clean ref: {snr_db(clean, denoised):.2f} dB")
    print(f"  NCC  noisy    -> clean ref: {time_domain_ncc(clean, noisy):.4f}")
    print(f"  NCC  denoised -> clean ref: {time_domain_ncc(clean, denoised):.4f}")

    # false-positive cost: any flagged window with no injected transient
    injected = np.zeros(len(t), dtype=bool)
    injected[550:750] = True  # generous margin around the decaying burst
    injected[1200:1260] = True
    for _, row in diag[diag["label"] == "noisy"].iterrows():
        lo, hi = int(row["win_start"] / dt), int(row["win_stop"] / dt)
        if not injected[lo:hi].any():
            sub_snr_before = snr_db(clean[lo:hi], noisy[lo:hi])
            sub_snr_after = snr_db(clean[lo:hi], denoised[lo:hi])
            print(
                f"  false-positive window [{row['win_start']:.0f},"
                f" {row['win_stop']:.0f}) s: SNR {sub_snr_before:.1f} ->"
                f" {sub_snr_after:.1f} dB (genuinely clean signal, "
                f"needlessly re-coded)"
            )

    # amplitude at which the square block stops being flagged (see
    # "Parameters and limitations")
    for amp in (30.0, 35.0):
        noise = np.zeros(len(t))
        env = 40.0 * np.exp(-np.maximum(t - 600.0, 0.0) / 25.0)
        env[t < 600.0] = 0.0
        noise += env
        noise[1215:1245] -= amp
        ts_amp = TSData(data={"EX": clean + noise}, dt=dt, station="amp_sweep")
        den_amp = TimeSeriesDenoiser(mmf_size=121, win_seconds=60.0, random_state=0)
        den_amp.apply(ts_amp)
        hit = den_amp.diagnostics_
        hit = hit[(hit["win_start"] <= 1230) & (1230 < hit["win_stop"])]
        print(
            f"  square amplitude {amp:.0f} (background std ~2): window "
            f"label={hit['label'].values[0]!r} score={hit['score'].values[0]:.3f}"
        )

    fig = plot_ts_denoise_summary(
        t, noisy, denoised, diag, low=low, high=high,
        suptitle="Synthetic record: charge/discharge burst + square block",
    )
    _save(fig, save_dir, "tsdenoise_synthetic_summary.png")


def run_entropy_features_demo(save_dir: Path) -> None:
    print("\n=== compute_entropy_features -- clean vs. interference windows ===")
    t, dt, clean, noisy = make_synthetic_record()
    _, high = mmf_split(noisy, size=121)

    windows = {
        "clean": high[0:60],
        "charge/discharge burst": high[600:660],
        "square block": high[1200:1260],
    }
    feats = {name: compute_entropy_features(w[None, :])[0]
             for name, w in windows.items()}
    for name, F in feats.items():
        print(f"  {name:>25s}: SE={F[0]:.3f} FE={F[1]:.3f} "
              f"AE={F[2]:.3f} BD={F[3]:.3f}")

    labels = ["Sample\nentropy", "Fuzzy\nentropy", "Approximate\nentropy",
              "Box-counting\ndimension"]
    colors = {"clean": "#2ca02c", "charge/discharge burst": "#d62728",
              "square block": "#9467bd"}
    x = np.arange(len(labels))
    width = 0.25
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    for i, (name, F) in enumerate(feats.items()):
        ax.bar(x + (i - 1) * width, F, width, label=name, color=colors[name])
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Feature value")
    ax.set_title("compute_entropy_features -- one example window per class")
    ax.legend(fontsize=8, frameon=False)
    ax.grid(True, axis="y", ls=":", lw=0.4, color="gray", alpha=0.5)
    ax.set_axisbelow(True)
    fig.tight_layout()
    _save(fig, save_dir, "tsdenoise_entropy_features.png")


def run_ksvd_direct_demo(ts_full: TSData, save_dir: Path) -> None:
    print("\n=== KSVDDenoiser -- direct example on a real segment ===")
    dt = ts_full.dt
    s, e = 339500, 342500  # ~4.2-hour real EY segment containing the burst
    ey = ts_full.get("EY")[s:e]
    t = np.arange(ey.size) * dt
    peak = int(np.argmax(np.abs(ey)))

    low, high = mmf_split(ey, size=241)
    print(f"  segment length {high.size} samples -> "
          f"~{(high.size - 32) // 16 + 1} overlapping 32-sample patches")

    for n_atoms in (64, 32, 16, 8, 4):
        den = KSVDDenoiser(n_atoms=n_atoms, random_state=0).fit(high)
        K = den._D.shape[1]
        denoised = den.transform(high)
        quiet_raw = high[700:1000]
        quiet_den = denoised[700:1000]
        print(f"  n_atoms={n_atoms:>2d} (K={K:>2d}): "
              f"denoised.std={denoised.std():.4f}  "
              f"quiet-region std raw={quiet_raw.std():.4f} "
              f"-> denoised={quiet_den.std():.4f}")

    den = KSVDDenoiser(random_state=0).fit(high)
    noise_profile = den.noise_profile(high)
    denoised = den.transform(high)
    print(f"  default n_atoms=64 -> K={den._D.shape[1]} atoms from "
          f"~{(high.size - 32) // 16 + 1} patches")
    print(f"  at the peak: raw={high[peak]:.2f} "
          f"noise_profile={noise_profile[peak]:.2f} "
          f"denoised={denoised[peak]:.2f}")

    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(9.0, 4.6), sharex=True,
                                    layout="constrained")
    ax0.plot(t, high, color="#7f7f7f", lw=0.8, label="High-frequency residual")
    ax0.plot(t, noise_profile, color="#d62728", lw=1.0, ls="--",
              label="K-SVD noise_profile")
    ax0.set_ylabel("Amplitude")
    ax0.set_title("(a) Residual vs. its K-SVD sparse reconstruction "
                   "(default n_atoms=64)")
    ax0.legend(fontsize=8, frameon=False, loc="upper right")

    ax1.plot(t, denoised, color="#1f77b4", lw=0.8)
    ax1.axhline(0, color="0.3", lw=0.6)
    ax1.set_ylabel("Amplitude")
    ax1.set_xlabel("Time (s)")
    ax1.set_title("(b) denoised = residual - noise_profile")

    _save(fig, save_dir, "tsdenoise_ksvd_direct.png")


def run_classifier_accuracy_demo() -> None:
    print("\n=== SignalQualityClassifier -- synthetic-library holdout accuracy ===")
    segs, y = generate_synthetic_library(win_len=60, n_per_class=150, seed=11)
    rng = np.random.default_rng(12)
    idx = rng.permutation(len(segs))
    n_train = int(0.7 * len(segs))
    tr, te = idx[:n_train], idx[n_train:]

    clf = SignalQualityClassifier(random_state=0).fit(segs[tr], y[tr])
    pred = clf.predict(segs[te])
    acc = float((pred == y[te]).mean())
    print(f"  holdout accuracy ({len(te)} windows, 4 noise types): {acc:.4f}")


# ------------------------------------------------------------------ real


def load_kap103(verbose: int = 1) -> TSData:
    path = (
        _repository_root() / "data" / "MT" / "TS" / "kap103as.ts" / "kap103as.ts"
    )
    ts = read_ts(str(path))
    if verbose:
        print(f"  loaded {ts!r}")
        print(f"  channels={ts.channels()} dt={ts.dt}s "
              f"duration={ts.duration / 86400:.1f} days")
    return ts


def run_real_data_demo(ts_full: TSData, save_dir: Path) -> None:
    print("\n=== TimeSeriesDenoiser -- real field recording (kap103as.ts) ===")
    dt = ts_full.dt
    s, e = 341100, 341900  # ~66-minute window around a real EY burst
    ey_raw = ts_full.get("EY")[s:e]
    t = np.arange(ey_raw.size) * dt

    peak = int(np.argmax(np.abs(ey_raw)))
    robust_std = 1.4826 * np.median(np.abs(ey_raw - np.median(ey_raw)))
    print(f"  EY peak amplitude in this window: {ey_raw[peak]:.1f} "
          f"(robust background std ~{robust_std:.1f} -- the plain std, "
          f"{np.std(ey_raw):.1f}, is itself inflated by this burst)")

    for size in (31, 61, 121, 241, 361):
        low, _ = mmf_split(ey_raw, size=size)
        print(f"  mmf_size={size:>3d} samples ({size * dt:>4.0f} s): "
              f"low range [{np.nanmin(low):.2f}, {np.nanmax(low):.2f}]")

    low_tuned, high_tuned = mmf_split(ey_raw, size=241)
    fig = plot_ts_denoise_mmf_split(
        t, ey_raw, low_tuned, high_tuned,
        suptitle="kap103as.ts -- EY, mmf_size=241 (1205 s)",
    )
    _save(fig, save_dir, "tsdenoise_mmf_split_real.png")

    ts_sub = TSData(dt=dt, station=ts_full.station)
    ts_sub.add_channel("EY", ey_raw)

    den = TimeSeriesDenoiser(mmf_size=241, win_seconds=300.0, random_state=0)
    out = den.apply(ts_sub)
    denoised = out.get("EY")
    diag = den.diagnostics_
    print(diag.to_string(index=False))
    print(f"  peak amplitude: raw={ey_raw[peak]:.2f} "
          f"denoised={denoised[peak]:.2f}")

    low, high = mmf_split(ey_raw, size=241)
    fig = plot_ts_denoise_summary(
        t, ey_raw, denoised, diag, low=low, high=high,
        suptitle="kap103as.ts -- EY, real interference burst",
    )
    _save(fig, save_dir, "tsdenoise_real_summary.png")


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

    run_synthetic_demo(args.save_dir)
    run_entropy_features_demo(args.save_dir)
    run_classifier_accuracy_demo()

    ts_full = load_kap103()
    run_ksvd_direct_demo(ts_full, args.save_dir)
    run_real_data_demo(ts_full, args.save_dir)


if __name__ == "__main__":
    main()
