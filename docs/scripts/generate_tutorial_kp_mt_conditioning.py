"""Generate figures and sample outputs for the KP MT conditioning tutorial."""

from __future__ import annotations

import contextlib
import copy
import io
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Ellipse

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data" / "MT" / "kap03lmt_edis"
IMAGE_DIR = (
    ROOT
    / "docs"
    / "source"
    / "images"
    / "tutorials"
    / "condition_mt_line_with_tipper_and_rotation"
)
STATIONS = ["kap103", "kap112", "kap136", "kap169"]
MU0_FACTOR = 0.2


def _import_pycsamt():
    sys.path.insert(0, str(ROOT))
    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        pass

    return locals()


def _style_axis(ax: plt.Axes) -> None:
    ax.set_facecolor("#fbfbf7")
    ax.grid(True, color="#d8dad1", linewidth=0.7, alpha=0.65)
    for spine in ax.spines.values():
        spine.set_color("#39434d")
        spine.set_linewidth(0.8)


def _save(fig: plt.Figure, name: str) -> None:
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(IMAGE_DIR / name, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _site_map(sites) -> dict[str, object]:
    out = {}
    for i, site in enumerate(sites):
        names = [
            getattr(site, "station", None),
            getattr(site, "id", None),
            getattr(getattr(site, "edi", None), "station", None),
        ]
        fallback = getattr(site, "name", None)
        if fallback is not None:
            names.append(fallback)
        names.append(str(i))
        for name in names:
            if name is not None:
                out[str(name)] = site
    return out


def _get_site(sites, station: str):
    by_name = _site_map(sites)
    if station in by_name:
        return by_name[station]
    for name, site in by_name.items():
        if name.lower() == station.lower() or station.lower() in name.lower():
            return site
    raise KeyError(station)


def _edi(site):
    return getattr(site, "edi", site)


def _rho_phase(
    site, comp: tuple[int, int]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    z_obj = _edi(site).Z
    freq = np.asarray(z_obj.freq, dtype=float)
    z = np.asarray(z_obj.z, dtype=complex)[:, comp[0], comp[1]]
    rho = MU0_FACTOR * np.abs(z) ** 2 / np.maximum(freq, 1e-30)
    phase = np.angle(z, deg=True)
    return freq, rho, phase


def _tipper_arrays(site) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    tip = _edi(site).Tip
    freq = np.asarray(tip.freq, dtype=float)
    tx = np.asarray(tip.tipper[:, 0, 0], dtype=complex)
    ty = np.asarray(tip.tipper[:, 0, 1], dtype=complex)
    return freq, tx, ty


def _plot_raw_tensor(sites) -> None:
    comps = {
        "xx": (0, 0, "#7c4d79"),
        "xy": (0, 1, "#2f6f8f"),
        "yx": (1, 0, "#c85745"),
        "yy": (1, 1, "#3f8f61"),
    }
    fig, axes = plt.subplots(
        2, 2, figsize=(12.0, 8.0), constrained_layout=True
    )
    for ax, station in zip(axes.flat, STATIONS):
        site = _get_site(sites, station)
        for label, (i, j, color) in comps.items():
            freq, rho, phase = _rho_phase(site, (i, j))
            period = 1.0 / freq
            ax.plot(
                period,
                rho,
                color=color,
                marker="o",
                markersize=3,
                linewidth=1.0,
                label=f"rho {label}",
            )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.invert_xaxis()
        ax.set_title(station)
        ax.set_xlabel("Period (s)")
        ax.set_ylabel("Apparent resistivity (ohm m)")
        _style_axis(ax)
    axes[0, 0].legend(fontsize=8, ncol=2)
    _save(fig, "kp_raw_tensor_components.png")

    fig, axes = plt.subplots(
        2, 2, figsize=(12.0, 8.0), constrained_layout=True
    )
    for ax, station in zip(axes.flat, STATIONS):
        site = _get_site(sites, station)
        for label, (i, j, color) in comps.items():
            freq, rho, phase = _rho_phase(site, (i, j))
            period = 1.0 / freq
            ax.plot(
                period,
                phase,
                color=color,
                marker="o",
                markersize=3,
                linewidth=1.0,
                label=f"phi {label}",
            )
        ax.set_xscale("log")
        ax.invert_xaxis()
        ax.set_title(station)
        ax.set_xlabel("Period (s)")
        ax.set_ylabel("Phase (deg)")
        _style_axis(ax)
    axes[0, 0].legend(fontsize=8, ncol=2)
    _save(fig, "kp_raw_phase_components.png")


def _plot_tipper(sites) -> None:
    fig, axes = plt.subplots(
        2, 2, figsize=(12.0, 7.6), constrained_layout=True
    )
    for ax, station in zip(axes.flat, STATIONS):
        site = _get_site(sites, station)
        freq, tx, ty = _tipper_arrays(site)
        period = 1.0 / freq
        ax.plot(
            period,
            np.abs(tx),
            color="#2f6f8f",
            marker="o",
            markersize=3,
            label="|Tx|",
        )
        ax.plot(
            period,
            np.abs(ty),
            color="#c85745",
            marker="s",
            markersize=3,
            label="|Ty|",
        )
        ax.plot(
            period,
            np.abs(tx + 1j * ty),
            color="#3f8f61",
            linewidth=1.1,
            label="combined",
        )
        ax.set_xscale("log")
        ax.invert_xaxis()
        ax.set_title(station)
        ax.set_xlabel("Period (s)")
        ax.set_ylabel("Tipper amplitude")
        _style_axis(ax)
    axes[0, 0].legend(fontsize=8)
    _save(fig, "kp_raw_tipper_components.png")


def _plot_qc_frequency(functions, sites) -> pd.DataFrame:
    freq_ci = functions["frequency_confidence_table"](
        sites,
        method="composite",
        ci_hi=0.9,
        ci_lo=0.5,
        recursive=False,
        api=True,
    ).to_pandas(copy=True)
    pivot = freq_ci.pivot_table(
        index="period_s",
        columns="station",
        values="confidence",
        aggfunc="median",
    ).sort_index()
    stations = list(pivot.columns)
    fig, ax = plt.subplots(figsize=(11.4, 5.8))
    im = ax.imshow(
        pivot.to_numpy(),
        aspect="auto",
        origin="lower",
        cmap="viridis",
        vmin=0,
        vmax=1,
    )
    ax.set_xticks(np.arange(0, len(stations), 2))
    ax.set_xticklabels(stations[::2], rotation=45, ha="right")
    y_idx = np.linspace(
        0, len(pivot.index) - 1, min(8, len(pivot.index))
    ).astype(int)
    ax.set_yticks(y_idx)
    ax.set_yticklabels([f"{pivot.index[i]:.3g}" for i in y_idx])
    ax.set_ylabel("Period (s)")
    ax.set_title("Frequency confidence by station")
    cb = fig.colorbar(im, ax=ax, pad=0.015)
    cb.set_label("Confidence")
    _save(fig, "kp_qc_frequency_confidence.png")
    return freq_ci


def _plot_bad_frequency_mask(freq_ci: pd.DataFrame) -> None:
    weak = freq_ci.assign(weak=freq_ci["confidence"] < 0.5)
    summary = (
        weak.groupby("period_s")
        .agg(
            weak_count=("weak", "sum"),
            total=("weak", "count"),
            median_confidence=("confidence", "median"),
        )
        .reset_index()
    )
    summary["weak_fraction"] = summary["weak_count"] / summary["total"]
    fig, ax1 = plt.subplots(figsize=(9.8, 4.6))
    ax1.bar(
        summary["period_s"],
        summary["weak_fraction"],
        width=summary["period_s"] * 0.12,
        color="#c85745",
        edgecolor="#27323a",
        linewidth=0.35,
        alpha=0.75,
        label="weak station fraction",
    )
    ax1.set_xscale("log")
    ax1.invert_xaxis()
    ax1.set_xlabel("Period (s)")
    ax1.set_ylabel("Weak fraction")
    _style_axis(ax1)
    ax2 = ax1.twinx()
    ax2.plot(
        summary["period_s"],
        summary["median_confidence"],
        color="#2f6f8f",
        marker="o",
        markersize=3,
        linewidth=1.2,
        label="median confidence",
    )
    ax2.set_ylabel("Median confidence")
    ax2.set_ylim(0, 1.05)
    ax1.set_title("Bad-frequency screening rule")
    _save(fig, "kp_bad_frequency_mask.png")


def _processing_chain(functions, sites):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="All-NaN slice encountered")
        notch = functions["notch_powerline"](
            sites,
            mains_hz=50.0,
            n_harm=20,
            tol_hz=0.06,
            recursive=False,
        )
        recovered = functions["recover_low_confidence_frequencies"](
            notch,
            method="composite",
            ci_hi=0.9,
            ci_lo=0.5,
            interpolation="linear",
            recursive=False,
        )
        filtered = functions["hampel_filter_freq"](
            recovered,
            win=3,
            nsig=3.0,
            recursive=False,
        )
        smoothed = functions["smooth_rho_phase"](
            filtered,
            components="offdiag",
            degree=3,
            smooth_rho=True,
            smooth_phase=True,
            recursive=False,
        )
    return smoothed


def _plot_filter_before_after(raw_sites, filtered_sites) -> None:
    stations = STATIONS[:3]
    fig, axes = plt.subplots(
        len(stations), 2, figsize=(11.6, 8.8), constrained_layout=True
    )
    for row, station in enumerate(stations):
        raw = _get_site(raw_sites, station)
        flt = _get_site(filtered_sites, station)
        for col, comp in enumerate([(0, 1), (1, 0)]):
            ax = axes[row, col]
            fr0, rho0, ph0 = _rho_phase(raw, comp)
            fr1, rho1, ph1 = _rho_phase(flt, comp)
            per0 = 1.0 / fr0
            per1 = 1.0 / fr1
            label = "xy" if comp == (0, 1) else "yx"
            ax.plot(
                per0,
                rho0,
                color="#9a9a9a",
                marker="o",
                markersize=3,
                linewidth=0.9,
                label="raw",
            )
            ax.plot(
                per1,
                rho1,
                color="#2f6f8f",
                marker="s",
                markersize=3,
                linewidth=1.2,
                label="conditioned",
            )
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.invert_xaxis()
            ax.set_title(f"{station} rho {label}")
            ax.set_xlabel("Period (s)")
            ax.set_ylabel("rho (ohm m)")
            _style_axis(ax)
    axes[0, 0].legend(fontsize=8)
    _save(fig, "kp_filter_before_after.png")


def _static_shift(functions, estimate_sites, apply_sites):
    factors = functions["estimate_ss_ama"](
        estimate_sites,
        sort_by="name",
        half_window=3,
        max_skew=None,
        recursive=False,
        api=True,
    ).to_pandas(copy=True)
    factors["fac_z_reviewed"] = factors["fac_z"].clip(lower=0.35, upper=2.85)
    applied = factors[["station", "fac_z_reviewed"]].rename(
        columns={"fac_z_reviewed": "fac_z"}
    )
    shifted = functions["apply_ss_factors"](
        apply_sites,
        applied,
        key="fac_z",
        inplace=False,
        recursive=False,
    )
    return factors, shifted


def _plot_static_shift(
    factors: pd.DataFrame, before_sites, after_sites
) -> None:
    fig, ax = plt.subplots(figsize=(10.8, 4.2))
    x = np.arange(len(factors))
    ax.plot(
        x,
        factors["fac_z"],
        color="#9a9a9a",
        marker="o",
        linewidth=1.0,
        label="raw AMA fac_z",
    )
    ax.plot(
        x,
        factors["fac_z_reviewed"],
        color="#c85745",
        marker="s",
        linewidth=1.2,
        label="reviewed/clipped fac_z",
    )
    ax.axhline(1.0, color="#27323a", linestyle="--", linewidth=0.9)
    ax.set_yscale("log")
    ax.set_xticks(x[::2])
    ax.set_xticklabels(factors["station"].iloc[::2], rotation=45, ha="right")
    ax.set_ylabel("Impedance scale factor")
    ax.set_title("Static-shift factors must be reviewed before application")
    ax.legend(fontsize=8)
    _style_axis(ax)
    _save(fig, "kp_static_shift_factors.png")

    fig, axes = plt.subplots(
        1, 2, figsize=(11.2, 4.4), constrained_layout=True
    )
    for ax, station in zip(axes, ["kap103", "kap115"]):
        raw = _get_site(before_sites, station)
        corr = _get_site(after_sites, station)
        for comp, color, label in [
            ((0, 1), "#2f6f8f", "xy"),
            ((1, 0), "#c85745", "yx"),
        ]:
            fr0, rho0, _ = _rho_phase(raw, comp)
            fr1, rho1, _ = _rho_phase(corr, comp)
            ax.plot(
                1.0 / fr0,
                rho0,
                color=color,
                alpha=0.35,
                marker="o",
                markersize=3,
                label=f"raw {label}",
            )
            ax.plot(
                1.0 / fr1,
                rho1,
                color=color,
                linestyle="--",
                marker="s",
                markersize=3,
                label=f"shifted {label}",
            )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.invert_xaxis()
        ax.set_title(station)
        ax.set_xlabel("Period (s)")
        ax.set_ylabel("rho (ohm m)")
        _style_axis(ax)
    axes[0].legend(fontsize=8, ncol=2)
    _save(fig, "kp_static_shift_before_after.png")


def _dominant_strike(functions, sites) -> tuple[float, pd.DataFrame]:
    detail = functions["analyze_anisotropy"](sites, recursive=False)
    strikes = detail["strike_deg"].to_numpy(dtype=float)
    strikes = strikes[np.isfinite(strikes)]
    doubled = np.deg2rad(2.0 * strikes)
    mean = 0.5 * np.rad2deg(
        np.arctan2(np.sin(doubled).mean(), np.cos(doubled).mean())
    )
    return float(mean), detail


def _plot_strike_rose(detail: pd.DataFrame, strike_deg: float) -> None:
    strikes = detail["strike_deg"].to_numpy(dtype=float)
    strikes = strikes[np.isfinite(strikes)]
    theta = np.deg2rad((strikes + 180.0) % 180.0)
    bins = np.linspace(0, np.pi, 19)
    hist, edges = np.histogram(theta, bins=bins)
    centers = 0.5 * (edges[:-1] + edges[1:])
    width = np.diff(edges)
    fig, ax = plt.subplots(figsize=(6.6, 6.0), subplot_kw={"polar": True})
    ax.bar(
        centers,
        hist,
        width=width,
        bottom=0.0,
        color="#2f6f8f",
        edgecolor="#27323a",
        alpha=0.78,
    )
    ax.plot(
        [np.deg2rad((strike_deg + 180) % 180)] * 2,
        [0, max(hist) * 1.1],
        color="#c85745",
        linewidth=2.0,
    )
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_title(f"Swift strike rose, dominant {strike_deg:.1f} deg")
    _save(fig, "kp_strike_rose.png")


def _plot_phase_tensor_grid(functions, sites) -> pd.DataFrame:
    pt = functions["build_phase_tensor_table"](sites, recursive=False)
    stations = list(dict.fromkeys(pt["station"]))
    periods = np.array(sorted(pt["period"].unique()))
    selected = periods[
        np.linspace(0, len(periods) - 1, min(8, len(periods))).astype(int)
    ]
    fig, ax = plt.subplots(figsize=(12.0, 5.8))
    max_s1 = np.nanpercentile(pt["s1"], 90)
    for ix, station in enumerate(stations):
        sdf = pt[pt["station"] == station]
        for period in selected:
            row = sdf.iloc[(sdf["period"] - period).abs().argsort()[:1]]
            if row.empty:
                continue
            r = row.iloc[0]
            y = np.where(selected == period)[0][0]
            width = 0.55 * float(r["s1"]) / max(max_s1, 1e-9)
            height = max(
                0.08,
                width * max(float(r["s2"]) / max(float(r["s1"]), 1e-9), 0.08),
            )
            ell = Ellipse(
                (ix, y),
                width=width,
                height=height,
                angle=float(r["theta"]),
                facecolor=plt.cm.RdBu_r(
                    np.clip((float(r["beta"]) + 20) / 40, 0, 1)
                ),
                edgecolor="#27323a",
                linewidth=0.4,
                alpha=0.92,
            )
            ax.add_patch(ell)
    ax.set_xlim(-0.8, len(stations) - 0.2)
    ax.set_ylim(-0.7, len(selected) - 0.3)
    ax.set_xticks(np.arange(0, len(stations), 2))
    ax.set_xticklabels(stations[::2], rotation=45, ha="right")
    ax.set_yticks(np.arange(len(selected)))
    ax.set_yticklabels([f"{p:.3g}" for p in selected])
    ax.set_ylabel("Period (s)")
    ax.set_title("Phase tensor ellipse grid, color by beta")
    _style_axis(ax)
    _save(fig, "kp_phase_tensor_grid.png")
    return pt


def _rotate_sites(sites, angle_deg: float):
    rotated = copy.deepcopy(sites)
    for site in rotated:
        ed = _edi(site)
        ed.Z.rotate(angle_deg)
        tip = getattr(ed, "Tip", None)
        if tip is not None:
            tip.rotate(angle_deg)
    return rotated


def _plot_rotation(before_sites, rotated_sites, angle_deg: float) -> None:
    fig, axes = plt.subplots(
        1, 2, figsize=(11.0, 4.4), constrained_layout=True
    )
    for ax, station in zip(axes, ["kap112", "kap136"]):
        before = _get_site(before_sites, station)
        after = _get_site(rotated_sites, station)
        for comp, color, label in [
            ((0, 0), "#7c4d79", "xx"),
            ((1, 1), "#3f8f61", "yy"),
            ((0, 1), "#2f6f8f", "xy"),
            ((1, 0), "#c85745", "yx"),
        ]:
            fr0, rho0, _ = _rho_phase(before, comp)
            fr1, rho1, _ = _rho_phase(after, comp)
            style = "--" if label in {"xx", "yy"} else "-"
            ax.plot(1.0 / fr0, rho0, color=color, alpha=0.28, linewidth=0.9)
            ax.plot(
                1.0 / fr1,
                rho1,
                color=color,
                linestyle=style,
                linewidth=1.3,
                label=f"rot {label}",
            )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.invert_xaxis()
        ax.set_title(f"{station}, rotated {angle_deg:.1f} deg")
        ax.set_xlabel("Period (s)")
        ax.set_ylabel("rho (ohm m)")
        _style_axis(ax)
    axes[0].legend(fontsize=8, ncol=2)
    _save(fig, "kp_rotation_before_after.png")


def _decision_table(
    qc: pd.DataFrame, factors: pd.DataFrame, strike_deg: float
) -> None:
    rows = [
        ["Tipper present?", "yes, all stations have tipper rows"],
        [
            "Suppress weak frequencies?",
            "recover/mask rows below confidence 0.5",
        ],
        [
            "Power-line step?",
            "run 50 Hz diagnostic notch; skip if no affected row",
        ],
        ["Static shift?", "apply reviewed factors, not raw extremes"],
        ["Rotation angle?", f"use dominant strike {strike_deg:.1f} deg"],
        ["Before inversion?", "export one conditioned, rotated EDI folder"],
    ]
    fig, ax = plt.subplots(figsize=(10.8, 3.4))
    ax.axis("off")
    table = ax.table(
        cellText=rows,
        colLabels=["processing question", "decision"],
        cellLoc="left",
        colLoc="left",
        loc="center",
        colWidths=[0.42, 0.58],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9.2)
    table.scale(1.0, 1.5)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#39434d")
        cell.set_linewidth(0.5)
        if row == 0:
            cell.set_facecolor("#2f6f8f")
            cell.get_text().set_color("white")
            cell.get_text().set_weight("bold")
        elif col == 1:
            cell.set_facecolor("#eef4ed")
        else:
            cell.set_facecolor("#fbfbf7")
    ax.set_title("KP MT conditioning decisions", pad=12)
    _save(fig, "kp_processing_decision_table.png")


def main() -> int:
    functions = _import_pycsamt()
    survey = functions["read_edis"](
        DATA_DIR,
        recursive=False,
        strict=False,
        on_dup="replace",
        progress=False,
    )
    sites = survey.collection
    inventory = survey.summary().to_pandas(copy=True)
    qc = functions["build_qc_table"](
        sites,
        include_skew=True,
        recursive=False,
        api=True,
    ).to_pandas(copy=True)
    station_ci = functions["station_confidence_table"](
        sites,
        method="composite",
        recursive=False,
        api=True,
    ).to_pandas(copy=True)

    _plot_raw_tensor(sites)
    _plot_tipper(sites)
    freq_ci = _plot_qc_frequency(functions, sites)
    _plot_bad_frequency_mask(freq_ci)
    conditioned = _processing_chain(functions, sites)
    _plot_filter_before_after(sites, conditioned)
    factors, shifted = _static_shift(functions, sites, sites)
    _plot_static_shift(factors, sites, shifted)
    strike_deg, strike_detail = _dominant_strike(functions, shifted)
    _plot_strike_rose(strike_detail, strike_deg)
    pt = _plot_phase_tensor_grid(functions, shifted)
    rotated = _rotate_sites(shifted, strike_deg)
    _plot_rotation(shifted, rotated, strike_deg)
    _decision_table(qc, factors, strike_deg)

    weak_rows = int((freq_ci["confidence"] < 0.5).sum())
    print("survey_summary:")
    print(survey.summary())
    print("inventory:")
    print(
        inventory[["station", "n_freq", "tipper", "spectra"]]
        .head(6)
        .to_string(index=False)
    )
    print("qc_summary:")
    print(
        qc[["station", "n_freq", "n_tip", "frac_ok", "snr_med", "skew_med"]]
        .head(6)
        .to_string(index=False, float_format=lambda value: f"{value:.3f}")
    )
    print("confidence_summary:")
    print(
        station_ci[
            [
                "station",
                "confidence",
                "coverage",
                "offdiag",
                "diagonal",
                "spatial",
            ]
        ]
        .head(6)
        .to_string(index=False, float_format=lambda value: f"{value:.3f}")
    )
    print("frequency_screen:")
    print(
        f"rows={len(freq_ci)} weak_rows={weak_rows} weak_fraction={weak_rows / len(freq_ci):.3f}"
    )
    print("static_shift_factors:")
    print(
        factors[["station", "fac_z", "fac_z_reviewed", "n_used"]]
        .head(8)
        .to_string(index=False, float_format=lambda value: f"{value:.3g}")
    )
    print("strike_rotation:")
    print(f"dominant_strike_deg={strike_deg:.2f}")
    print(f"phase_tensor_rows={len(pt)}")
    print("decision:")
    print(
        "order=load -> raw curves/tipper -> QC -> recover/mask -> notch/filter -> static shift -> strike/PT -> rotate"
    )
    print(f"images: {IMAGE_DIR.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
