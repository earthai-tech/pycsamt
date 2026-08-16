"""Generate figures and sample outputs for the static-shift tutorial."""

from __future__ import annotations

import contextlib
import io
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
CSAMT_DIR = ROOT / "data" / "CSAMT"
MT_DIR = ROOT / "data" / "MT" / "kap03lmt_edis"
IMAGE_DIR = (
    ROOT / "docs" / "source" / "images" / "tutorials" / "correct_static_shift"
)


def _import_pycsamt():
    sys.path.insert(0, str(ROOT))
    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        from pycsamt.api import read_edis
        from pycsamt.emtools.ss import (
            _get_z_block,
            _iter_items,
            _name,
            _rho_det_from_z,
            apply_ss_factors,
            correct_ss_ama,
            detect_near_surface,
            estimate_ss_ama,
            estimate_ss_bilateral,
            estimate_ss_loess,
            estimate_ss_refmedian,
            plot_ss_radar,
            plot_ss_station_curves,
        )

    return locals()


def _style_axis(ax: plt.Axes) -> None:
    ax.set_facecolor("#fbfbf7")
    ax.grid(True, color="#d6d8d1", linewidth=0.7, alpha=0.65)
    for spine in ax.spines.values():
        spine.set_color("#39434d")
        spine.set_linewidth(0.8)


def _save(fig: plt.Figure, name: str) -> None:
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(IMAGE_DIR / name, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _factor_plot(
    factor_df,
    filename: str = "ama_factor_profile.png",
    title: str = "AMA static-shift factors, exploratory full-band estimate",
    xtick_step: int = 2,
) -> None:
    x = np.arange(len(factor_df))
    colors = np.where(
        factor_df["delta_log10_rho"].to_numpy(dtype=float) >= 0,
        "#2f6f8f",
        "#c85745",
    )
    fig, ax = plt.subplots(figsize=(11.2, 4.2))
    ax.bar(
        x,
        factor_df["delta_log10_rho"],
        color=colors,
        edgecolor="#27323a",
        linewidth=0.45,
        alpha=0.86,
    )
    ax.axhline(0.0, color="#27323a", linewidth=1.0)
    ax.axhline(0.30, color="#c85745", linestyle="--", linewidth=1.0)
    ax.axhline(-0.30, color="#c85745", linestyle="--", linewidth=1.0)
    ax.set_xticks(x[::xtick_step])
    ax.set_xticklabels(
        factor_df["station"].iloc[::xtick_step], rotation=45, ha="right"
    )
    ax.set_ylabel(r"$\Delta \log_{10}\rho_a$")
    ax.set_xlabel("Station")
    ax.set_title(title)
    _style_axis(ax)
    _save(fig, filename)


def _method_comparison_plot(tables) -> None:
    fig, ax = plt.subplots(figsize=(11.2, 4.3))
    styles = {
        "AMA": ("#2f6f8f", "o"),
        "LOESS": ("#7c4d79", "s"),
        "reference median": ("#d5962c", "^"),
        "bilateral": ("#3f8f61", "D"),
    }
    for label, df in tables.items():
        color, marker = styles[label]
        ax.plot(
            df["station"],
            df["fac_z"],
            color=color,
            marker=marker,
            markersize=3.8,
            linewidth=1.25,
            label=label,
        )
    ax.axhline(1.0, color="#27323a", linewidth=1.0)
    ax.axhspan(0.5, 2.0, color="#89b58a", alpha=0.14)
    ax.set_xticks(np.arange(0, len(tables["AMA"]), 2))
    ax.set_xticklabels(
        tables["AMA"]["station"].iloc[::2], rotation=45, ha="right"
    )
    ax.set_ylabel("Impedance scale factor")
    ax.set_xlabel("Station")
    ax.set_title("Estimator comparison for the same period band")
    ax.legend(fontsize=8, ncol=2)
    _style_axis(ax)
    _save(fig, "estimator_factor_comparison.png")


def _logrho_matrix(functions, sites):
    rows = []
    labels = []
    freqs_ref = None
    for i, ed in enumerate(functions["_iter_items"](sites)):
        _, z, fr = functions["_get_z_block"](ed)
        if z is None:
            continue
        rho = functions["_rho_det_from_z"](z, fr)
        labels.append(functions["_name"](ed, i))
        if freqs_ref is None:
            freqs_ref = np.asarray(fr, dtype=float)
        rows.append(np.log10(np.maximum(rho, 1e-24)))
    return (
        np.asarray(rows, dtype=float),
        np.asarray(freqs_ref, dtype=float),
        labels,
    )


def _compare_plots(functions, sites, corrected, factors) -> None:
    before_rho, freqs, labels = _logrho_matrix(functions, sites)
    after_rho, _, _ = _logrho_matrix(functions, corrected)
    delta = after_rho - before_rho
    periods = 1.0 / freqs
    order = np.argsort(np.log10(periods))

    fig, axes = plt.subplots(
        2, 1, figsize=(11.2, 8.4), constrained_layout=True
    )
    ax = axes[0]
    x = np.arange(len(labels))
    station_delta = np.nanmedian(delta, axis=1)
    colors = np.where(station_delta >= 0.0, "#c85745", "#2f6f8f")
    ax.bar(x, station_delta, color=colors, edgecolor="#27323a", linewidth=0.45)
    ax.axhline(0.0, color="#27323a", linewidth=1.0)
    ax.set_xticks(x[::2])
    ax.set_xticklabels(labels[::2], rotation=45, ha="right")
    ax.set_ylabel(r"$\Delta \log_{10}\rho_a$")
    ax.set_xlabel("Station")
    ax.set_title("Applied determinant-resistivity correction by station")
    _style_axis(ax)

    ax = axes[1]
    im = ax.imshow(
        delta[:, order].T,
        aspect="auto",
        origin="lower",
        cmap="RdBu_r",
        vmin=-1.5,
        vmax=1.5,
    )
    yticks = np.linspace(0, len(order) - 1, 6).astype(int)
    ax.set_yticks(yticks)
    ax.set_yticklabels([f"{np.log10(periods[order[i]]):.2g}" for i in yticks])
    ax.set_xticks(x[::2])
    ax.set_xticklabels(labels[::2], rotation=45, ha="right")
    ax.set_ylabel(r"$\log_{10}T$ (s)")
    ax.set_xlabel("Station")
    ax.set_title("The applied factor is period independent")
    cb = fig.colorbar(im, ax=ax, pad=0.015)
    cb.set_label(r"$\Delta \log_{10}\rho_a$ after - before")
    _style_axis(ax)
    _save(fig, "before_after_delta_grid.png")

    station = str(
        factors.reindex(
            factors["delta_log10_rho"].abs().sort_values(ascending=False).index
        )["station"].iloc[0]
    )
    ax = functions["plot_ss_station_curves"](
        sites,
        corrected,
        station=station,
        pband=None,
    )
    ax.set_title(f"{station}: apparent-resistivity curve before and after")
    _style_axis(ax)
    _save(ax.figure, "station_curve_before_after.png")
    return station


def _radar_before_after(functions, sites, corrected, station: str, name: str) -> None:
    """Two-panel raw-vs-AMA-corrected static-shift radar for one station.

    Both panels share one radial (log10 rho_a) scale so the correction's
    actual amplitude shift is visible -- matplotlib polar axes otherwise
    autoscale each panel independently, which would hide a uniform shift.
    """
    fig, axes = plt.subplots(
        1, 2, figsize=(9.6, 4.8), subplot_kw={"polar": True}
    )
    functions["plot_ss_radar"](
        sites, station=station, rotate="none", ax=axes[0],
    )
    axes[0].set_title(f"{station} (raw)", pad=14)
    functions["plot_ss_radar"](
        corrected, station=station, rotate="none", ax=axes[1],
    )
    axes[1].set_title(f"{station} (AMA-corrected)", pad=14)

    r_all = np.concatenate(
        [
            np.concatenate([line.get_ydata() for line in ax.lines])
            for ax in axes
        ]
    )
    r_all = r_all[np.isfinite(r_all)]
    pad = 0.08 * (r_all.max() - r_all.min() + 1e-9)
    rmin, rmax = r_all.min() - pad, r_all.max() + pad
    for ax in axes:
        ax.set_ylim(rmin, rmax)

    fig.subplots_adjust(wspace=0.55)
    _save(fig, name)


def _csamt_static_shift(functions):
    """Real, freshly-computed Tongkeng CSAMT static-shift estimate."""
    sites = functions["read_edis"](
        CSAMT_DIR, recursive=False, strict=False, progress=False,
    ).collection
    ss = functions["estimate_ss_ama"](
        sites, recursive=False, api=True,
    ).to_pandas(copy=True)
    _factor_plot(
        ss,
        filename="csamt_tongkeng_factor_profile.png",
        title="AMA static-shift factors, Tongkeng CSAMT (default settings)",
        xtick_step=1,
    )
    ax = functions["plot_ss_radar"](
        sites, station="csa000", rotate="none",
    )
    ax.set_title("csa000 (Tongkeng CSAMT, raw)", pad=14)
    _save(ax.figure, "radar_csamt_tongkeng_csa000_raw.png")
    return ss


def _mt_static_shift(functions):
    """Real, freshly-computed KAP03 MT static-shift trial."""
    sites = functions["read_edis"](
        MT_DIR, recursive=False, strict=True, progress=False,
    ).collection
    ss = functions["estimate_ss_ama"](
        sites,
        sort_by="name",
        half_window=3,
        max_skew=None,
        recursive=False,
        api=True,
    ).to_pandas(copy=True)
    _factor_plot(
        ss,
        filename="mt_kap03_factor_profile.png",
        title="AMA static-shift factors, KAP03 MT (sort_by='name')",
        xtick_step=2,
    )
    ax = functions["plot_ss_radar"](
        sites, station="kap112", rotate="none",
    )
    ax.set_title("kap112 (KAP03 MT, raw)", pad=14)
    _save(ax.figure, "radar_mt_kap03_kap112_raw.png")
    return ss


def main() -> int:
    functions = _import_pycsamt()
    survey = functions["read_edis"](
        DATA_DIR,
        recursive=False,
        strict=False,
        progress=False,
    )
    sites = survey.collection

    conservative = functions["estimate_ss_ama"](
        sites,
        half_window=3,
        weights="tri",
        pband=None,
        max_skew=6.0,
        recursive=False,
        api=True,
    ).to_pandas(copy=True)

    factors = functions["estimate_ss_ama"](
        sites,
        sort_by="name",
        half_window=3,
        weights="tri",
        pband=None,
        max_skew=None,
        recursive=False,
        api=True,
    ).to_pandas(copy=True)
    corrected = functions["apply_ss_factors"](
        sites,
        factors,
        key="fac_z",
        inplace=False,
        recursive=False,
        verbose=0,
    )
    one_step = functions["correct_ss_ama"](
        sites,
        sort_by="name",
        half_window=3,
        weights="tri",
        pband=None,
        max_skew=None,
        inplace=False,
        recursive=False,
        verbose=0,
    )

    pband = (0.01, 10.0)
    loess = functions["estimate_ss_loess"](
        sites,
        half_window=3,
        poly=1,
        it=2,
        pband=pband,
        max_skew=6.0,
        api=True,
    ).to_pandas(copy=True)
    bilateral = functions["estimate_ss_bilateral"](
        sites,
        half_window=4,
        pband=pband,
        max_skew=6.0,
        summary="median",
        api=True,
    ).to_pandas(copy=True)
    refmedian = functions["estimate_ss_refmedian"](
        sites,
        pband=pband,
        max_skew=6.0,
        api=True,
    ).to_pandas(copy=True)

    _factor_plot(factors)
    _method_comparison_plot(
        {
            "AMA": factors,
            "LOESS": loess,
            "reference median": refmedian,
            "bilateral": bilateral,
        }
    )
    max_delta_station = _compare_plots(functions, sites, corrected, factors)
    _radar_before_after(
        functions,
        sites,
        corrected,
        max_delta_station,
        "radar_amt_station_raw_vs_corrected.png",
    )

    csamt_ss = _csamt_static_shift(functions)
    mt_ss = _mt_static_shift(functions)

    before_rho, freqs, labels = _logrho_matrix(functions, sites)
    after_rho, _, _ = _logrho_matrix(functions, corrected)
    delta = np.nanmedian(after_rho - before_rho, axis=1)
    print(f"survey_summary: {len(labels)} stations x {len(freqs)} frequencies")
    print("conservative_ama_head:")
    print(conservative.head(5).to_string(index=False))
    print("conservative_n_used:")
    print(conservative["n_used"].describe().to_string())
    print("exploratory_ama_head:")
    print(factors.head(6).to_string(index=False))
    print("exploratory_fac_z_range:")
    print(factors[["fac_rho", "fac_z", "n_used"]].describe().to_string())
    large = factors[
        (factors["fac_rho"] < 0.5)
        | (factors["fac_rho"] > 2.0)
        | (factors["n_used"] < 5)
    ]
    print(f"large_review_count: {len(large)}")
    print(
        large[["station", "delta_log10_rho", "fac_rho", "fac_z", "n_used"]]
        .head(8)
        .to_string(index=False)
    )
    print("applied_delta_head:")
    for station, value in zip(labels[:6], delta[:6]):
        print(f"{station:8s} {value:+.6f}")
    print(f"one_step_sites_type: {type(one_step).__name__}")
    print("method_head_fac_z:")
    joined = factors[["station", "fac_z"]].rename(columns={"fac_z": "ama"})
    joined = joined.merge(
        loess[["station", "fac_z"]].rename(columns={"fac_z": "loess"})
    )
    joined = joined.merge(
        refmedian[["station", "fac_z"]].rename(columns={"fac_z": "refmedian"})
    )
    joined = joined.merge(
        bilateral[["station", "fac_z"]].rename(columns={"fac_z": "bilateral"})
    )
    print(joined.head(6).to_string(index=False))
    print(f"max_delta_station: {max_delta_station}")

    print("csamt_tongkeng_ama:")
    print(
        csamt_ss[["station", "delta_log10_rho", "fac_rho", "fac_z", "n_used"]]
        .to_string(index=False)
    )
    print(f"csamt_delta_log10_rho_nunique: {csamt_ss['delta_log10_rho'].nunique()}")
    print("csamt_fac_rho_describe:")
    print(csamt_ss["fac_rho"].describe().to_string())
    csamt_ns = functions["detect_near_surface"](
        functions["read_edis"](
            CSAMT_DIR, recursive=False, strict=False, progress=False,
        ).collection,
        f_split=32.0,
        recursive=False,
    )
    print("csamt_distortion_type_counts:")
    print(csamt_ns["distortion_type"].value_counts().to_string())

    print("mt_kap03_ama_head:")
    print(
        mt_ss[["station", "delta_log10_rho", "fac_rho", "fac_z", "n_used"]]
        .head(8)
        .to_string(index=False)
    )
    print("mt_kap03_fac_z_describe:")
    print(mt_ss["fac_z"].describe().to_string())
    kap112_fac_z = float(mt_ss.loc[mt_ss["station"] == "kap112", "fac_z"].iloc[0])
    print(f"mt_kap112_fac_z: {kap112_fac_z:.6f}")

    print(f"images: {IMAGE_DIR.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
