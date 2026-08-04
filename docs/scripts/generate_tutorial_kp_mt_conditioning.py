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
TOPO_CSV = ROOT / "data" / "MT" / "kap03_topography_open_meteo.csv"
RESULT_DIR = ROOT / "results" / "kap03_mt_tutorial"
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
    from pycsamt.api import read_edis
    from pycsamt.emtools import (
        apply_ss_factors,
        build_qc_table,
        frequency_confidence_table,
        drop_low_confidence_frequencies,
        drop_freqs_manual,
        hampel_filter_freq,
        station_confidence_table,
        estimate_ss_ama,
        ensure_sites,
        estimate_strike_consensus,
        plot_frequency_confidence_psection,
        plot_frequency_edit_decisions,
        plot_dimensionality_psection,
        plot_ellipticity_psection,
        plot_induction_map,
        plot_induction_multiperiod_map,
        plot_induction_section,
        plot_phase_tensor_psection,
        plot_phase_tensor_strip_grid,
        plot_skew_ellipt_density,
        plot_response_tipper,
        plot_tipper_hodograms,
        plot_strike_analysis,
        plot_strike_rose,
    )
    from pycsamt.emtools.anisotropy import analyze_anisotropy
    from pycsamt.emtools.tensor import build_phase_tensor_table

    return locals()


def _inject_coordinates_topography(sites) -> pd.DataFrame:
    """Promote EDI DEFINEMEAS coordinates and cached DEM elevations."""
    from pycsamt.site.utils import set_coords

    topo = pd.read_csv(TOPO_CSV)
    lookup = topo.set_index("station")
    rows = []
    for site in sites:
        ed = _edi(site)
        dm = ed.get_section("definemeas")
        row = lookup.loc[ed.station]
        lat = float(dm.reflat)
        lon = float(dm.reflong)
        if abs(lat - float(row.latitude)) > 1e-5 or abs(lon - float(row.longitude)) > 1e-5:
            raise ValueError(f"{ed.station}: cached coordinates disagree with DEFINEMEAS")
        set_coords(ed, lat=lat, lon=lon, elev=float(row.elevation_m), inplace=True)
        rows.append((ed.station, lat, lon, float(row.elevation_m)))
    return pd.DataFrame(rows, columns=["station", "latitude", "longitude", "elevation_m"])


def _plot_topography(coords: pd.DataFrame) -> None:
    lat = np.deg2rad(coords.latitude.to_numpy())
    lon = np.deg2rad(coords.longitude.to_numpy())
    dlat, dlon = np.diff(lat), np.diff(lon)
    a = np.sin(dlat / 2) ** 2 + np.cos(lat[:-1]) * np.cos(lat[1:]) * np.sin(dlon / 2) ** 2
    segment_km = 2 * 6371.0088 * np.arcsin(np.sqrt(a))
    distance = np.r_[0.0, np.cumsum(segment_km)]
    coords["chainage_km"] = distance

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8), constrained_layout=True)
    sc = axes[0].scatter(coords.longitude, coords.latitude, c=coords.elevation_m,
                         cmap="terrain", edgecolor="k", s=38)
    axes[0].plot(coords.longitude, coords.latitude, "k-", lw=0.7, alpha=0.5)
    for row in coords.iloc[::3].itertuples():
        axes[0].annotate(row.station, (row.longitude, row.latitude), fontsize=7,
                         xytext=(3, 3), textcoords="offset points")
    axes[0].set(xlabel="Longitude (degrees east)", ylabel="Latitude (degrees north)",
                title="KAP03 station geometry and DEM elevation")
    fig.colorbar(sc, ax=axes[0], label="Open-Meteo elevation (m)")
    axes[1].plot(distance, coords.elevation_m, "o-", ms=4, color="#2f6f8f")
    for row in coords.iloc[::3].itertuples():
        axes[1].annotate(row.station, (row.chainage_km, row.elevation_m),
                         fontsize=7, rotation=55, xytext=(0, 6), textcoords="offset points")
    axes[1].set(xlabel="Cumulative geodesic chainage (km)", ylabel="Elevation (m)",
                title="Terrain used for inversion meshing")
    for ax in axes:
        _style_axis(ax)
    _save(fig, "kp_coordinates_topography.png")


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
    """Use the public response API for raw rho, phase, and tipper."""
    from pycsamt.emtools import plot_response_tipper

    fig = plot_response_tipper(
        sites, stations=STATIONS, components=("xy", "yx"), raw=True,
        ncols_groups=2, show_error_bars=True, show_tipper_error_bars=True,
        recursive=False,
    )
    _save(fig, "kp_raw_response_tipper.png")


def _plot_tipper(sites) -> None:
    """Use the public tipper-component API on three representative sites."""
    from pycsamt.emtools import plot_tipper_components

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.2), constrained_layout=True)
    for ax, station in zip(axes, ["kap103", "kap136", "kap169"]):
        plot_tipper_components(
            [_get_site(sites, station)], kind=("real", "imag"),
            axis="period", recursive=False, ax=ax,
        )
        ax.set_title(station)
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
    ax = functions["plot_frequency_confidence_psection"](
        sites, method="composite", ci_hi=0.9, ci_lo=0.5,
        station_label_step=2, recursive=False,
    )
    _save(ax.figure, "kp_qc_frequency_confidence.png")
    return freq_ci


def _plot_bad_frequency_mask(functions, before_sites, after_sites) -> None:
    """Use the public edit-decision plot to expose every dropped row."""
    ax = functions["plot_frequency_edit_decisions"](
        before_sites, after_sites, method="composite", ci_hi=0.9,
        ci_lo=0.5, station_label_step=2,
    )
    fig = ax.figure if hasattr(ax, "figure") else ax
    _save(fig, "kp_bad_frequency_mask.png")


def _processing_chain(functions, sites):
    """Drop demonstrably weak rows, then repair only isolated spikes."""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="All-NaN slice encountered")
        dropped = functions["drop_low_confidence_frequencies"](
            sites,
            method="composite",
            threshold=0.5,
            also="both",
            recursive=False,
        )
        supported = functions["drop_freqs_manual"](
            dropped, drop_freqs=(999.0,), tol_rel=0.005,
            inplace=False, recursive=False,
        )
        filtered = functions["hampel_filter_freq"](
            supported,
            win=2,
            nsig=3.0,
            on="both",
            domain="magphase",
            recursive=False,
        )
    return filtered


def _plot_dimensionality_and_induction(functions, sites) -> None:
    """Use public dimensionality and induction-vector diagnostics."""
    ax = functions["plot_dimensionality_psection"](
        sites, skew_th=3.0, ellipt_th=0.2, recursive=False
    )
    _save(ax.figure, "kp_dimensionality_psection.png")

    ax = functions["plot_ellipticity_psection"](
        sites, agg="median", recursive=False
    )
    _save(ax.figure, "kp_ellipticity_psection.png")

    ax = functions["plot_induction_section"](
        sites, component="abs", n_periods=18,
        title="KAP03 induction-vector amplitude", recursive=False,
    )
    _save(ax.figure, "kp_induction_section.png")

    ax = functions["plot_induction_map"](
        sites, period=653.2, convention="park", show_real=True,
        show_imag=True, station_labels=True,
        title="KAP03 induction vectors near 653 s", recursive=False,
    )
    _save(ax.figure, "kp_induction_map_653s.png")

    topo = pd.read_csv(TOPO_CSV)
    gx = np.linspace(topo.longitude.min(), topo.longitude.max(), 180)
    gy = np.linspace(topo.latitude.min(), topo.latitude.max(), 140)
    xx, yy = np.meshgrid(gx, gy)
    dx = xx[..., None] - topo.longitude.to_numpy()[None, None, :]
    dy = yy[..., None] - topo.latitude.to_numpy()[None, None, :]
    weights = 1.0 / np.maximum(dx * dx + dy * dy, 1e-10)
    dem = np.sum(weights * topo.elevation_m.to_numpy(), axis=2) / np.sum(weights, axis=2)
    extent = (gx.min(), gx.max(), gy.min(), gy.max())
    fig, _ = functions["plot_induction_multiperiod_map"](
        sites, periods=(30.0, 200.0, 1000.0, 8000.0),
        convention="park", background=dem, background_extent=extent,
        background_cmap="terrain", background_clim=(topo.elevation_m.min(), topo.elevation_m.max()),
        station_labels=False, title="KAP03 real induction vectors across period",
        xlabel="Longitude (degrees east)", ylabel="Latitude (degrees north)",
        recursive=False,
    )
    _save(fig, "kp_induction_multiperiod_map.png")

    fig = functions["plot_tipper_hodograms"](
        sites, station="kap136", n_bands=4, normalize=True,
        ms=2.5, lw=1.15, unit_circle=True, figsize=(8.2, 4.2), recursive=False,
    )
    fig.suptitle("KAP136 complex tipper hodograms by period band", fontsize=11)
    _save(fig, "kp_tipper_hodograms_kap136.png")


def _plot_three_station_raw_corrected(raw_sites, corrected_sites) -> None:
    """Use the public response API for matched raw and corrected panels."""
    from pycsamt.emtools import plot_response_tipper

    stations = ["kap103", "kap136", "kap169"]
    for data, raw, filename in (
        (raw_sites, True, "kp_three_station_raw.png"),
        (corrected_sites, False, "kp_three_station_corrected.png"),
    ):
        fig = plot_response_tipper(
            data, stations=stations, components=("xy", "yx"), raw=raw,
            ncols_groups=3, show_error_bars=True,
            show_tipper_error_bars=True, recursive=False,
        )
        _save(fig, filename)


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
    detail = functions["estimate_strike_consensus"](sites, recursive=False)
    angles = detail["ang"].to_numpy(dtype=float)
    weights = 1.0 / np.maximum(detail["iqr"].to_numpy(dtype=float), 1e-6)
    doubled = np.deg2rad(2.0 * angles)
    axial = 0.5 * np.rad2deg(
        np.arctan2(np.sum(weights * np.sin(doubled)),
                   np.sum(weights * np.cos(doubled)))
    ) % 180.0
    signed = ((axial + 90.0) % 180.0) - 90.0
    return float(signed), detail


def _plot_strike_rose(functions, sites) -> None:
    """Use public rose and Z/PT/tipper comparison APIs."""
    fig = functions["plot_strike_rose"](
        sites, method="consensus", bins=36,
        suptitle="KAP03 consensus geoelectric strike", recursive=False,
    )
    _save(fig, "kp_strike_rose.png")
    fig = functions["plot_strike_analysis"](
        sites, method="consensus", bins=36,
        suptitle="KAP03 impedance strike, phase-tensor azimuth, and tipper",
        recursive=False,
    )
    _save(fig, "kp_strike_analysis.png")


def _plot_phase_tensor_grid(functions, sites) -> pd.DataFrame:
    pt = functions["build_phase_tensor_table"](sites, recursive=False)
    fig = plt.figure(figsize=(15.0, 6.2))
    grid = fig.add_gridspec(1, 3, width_ratios=(1.15, 1.15, 0.8), wspace=0.38)
    ax = fig.add_subplot(grid[0, :2])
    functions["plot_phase_tensor_psection"](
        sites, axis_y="logperiod", period_up=False, c_by="beta",
        normalise_by="shape", min_aspect=0.12, clim=(-3.0, 3.0),
        color_mode="segmented",
        ellipse_kws={"edgecolor": "#202020", "linewidth": 0.55},
        cb_kws={"size": "3.2%", "pad": 0.08, "ticksize": 8},
        title="KAP03 phase-tensor pseudosection", recursive=False, ax=ax,
    )
    density_ax = fig.add_subplot(grid[0, 2])
    functions["plot_skew_ellipt_density"](
        sites, gridsize=24, recursive=False, ax=density_ax,
    )
    density_ax.axvline(3.0, color="#b2182b", ls="--", lw=1.0)
    density_ax.axhline(0.2, color="#2166ac", ls=":", lw=1.0)
    density_ax.set_title("Skew–ellipticity distribution", fontsize=10)
    fig.subplots_adjust(top=0.90, bottom=0.12, left=0.06, right=0.97)
    _save(fig, "kp_phase_tensor_psection_api.png")

    strip_fig = functions["plot_phase_tensor_strip_grid"](
        sites,
        profiles={"KAP03 selected stations": ["kap103", "kap121", "kap148", "kap175"]},
        c_by="beta", clim=(-3.0, 3.0), normalise_by="cell",
        min_aspect=0.12, edgecolor="#202020", linewidth=0.55,
        suptitle="KAP03 phase-tensor strips at representative stations",
        panel_size=(8.0, 1.25), recursive=False,
    )
    _save(strip_fig, "kp_phase_tensor_strip_grid.png")
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
    ax.invert_yaxis()
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


def _export_and_validate(functions, sites):
    """Write conditioned EDIs and prove coordinates, tipper, and inventory survive."""
    from pycsamt.site.export import write_sites

    edi_out = RESULT_DIR / "edi_conditioned_rotated"
    written = write_sites(
        sites, edi_out, exist_ok=True, manifest_csv=edi_out / "manifest.csv"
    )
    reloaded = functions["ensure_sites"](
        edi_out, recursive=False, strict=True
    ).ordered()
    if len(written) != len(sites) or len(reloaded) != len(sites):
        raise RuntimeError("conditioned EDI export/reload inventory mismatch")
    if not all(getattr(_edi(site), "Tip", None) is not None for site in reloaded):
        raise RuntimeError("tipper was lost during EDI round trip")
    elevations = np.asarray(
        [float(_edi(site).get_section("head").elev) for site in reloaded]
    )
    return written, reloaded, elevations


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
            "drop rows below confidence 0.5 from Z and tipper",
        ],
        [
            "Power-line step?",
            "omit: measured maximum is only 0.04 Hz",
        ],
        ["Static shift?", "reject automatic trial; preserve amplitudes"],
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
    coords = _inject_coordinates_topography(sites)
    _plot_topography(coords)
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
    conditioned = _processing_chain(functions, sites)
    _plot_bad_frequency_mask(functions, sites, conditioned)
    _plot_dimensionality_and_induction(functions, conditioned)
    _plot_filter_before_after(sites, conditioned)
    factors, shifted_trial = _static_shift(functions, sites, sites)
    _plot_static_shift(factors, sites, shifted_trial)
    strike_deg, strike_detail = _dominant_strike(functions, conditioned)
    _plot_strike_rose(functions, conditioned)
    pt = _plot_phase_tensor_grid(functions, conditioned)
    rotated = _rotate_sites(conditioned, strike_deg)
    _plot_rotation(conditioned, rotated, strike_deg)
    _plot_three_station_raw_corrected(sites, rotated)
    _decision_table(qc, factors, strike_deg)

    written, reloaded, elevations = _export_and_validate(functions, rotated)

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
    print("coordinates_topography:")
    print(
        f"stations={len(coords)} elevation_m={coords.elevation_m.min():.0f}-"
        f"{coords.elevation_m.max():.0f} chainage_km={coords.chainage_km.iloc[-1]:.1f}"
    )
    print(f"conditioned_edi_written_reloaded={len(written)}/{len(reloaded)}")
    print("decision:")
    print(
        "order=load/coordinates/topography -> raw curves/tipper -> QC -> drop weak rows -> Hampel -> reject static-shift trial -> strike/PT -> rotate Z+tipper -> export/reload"
    )
    print(f"images: {IMAGE_DIR.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
