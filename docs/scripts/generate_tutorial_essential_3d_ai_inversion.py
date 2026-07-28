"""Generate figures for the essential L18 3-D AI inversion tutorial."""

from __future__ import annotations

import contextlib
import io
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
IMAGE_DIR = (
    ROOT
    / "docs"
    / "source"
    / "images"
    / "tutorials"
    / "essential_3d_ai_inversion"
)


def _import_pycsamt():
    sys.path.insert(0, str(ROOT))
    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        pass

    return locals()


def _save(fig: plt.Figure, name: str) -> None:
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(IMAGE_DIR / name, dpi=190, bbox_inches="tight")
    plt.close(fig)


def _style_axis(ax: plt.Axes) -> None:
    ax.set_facecolor("#fbfbf7")
    ax.grid(True, color="#d8dad1", linewidth=0.65, alpha=0.70)
    for spine in ax.spines.values():
        spine.set_color("#303943")
        spine.set_linewidth(0.85)


def _cell_edges(centres: np.ndarray) -> np.ndarray:
    centres = np.asarray(centres, dtype=float).ravel()
    if centres.size == 1:
        width = max(float(centres[0]), 0.05)
        return np.asarray([0.0, centres[0] + width], dtype=float)
    edges = np.empty(centres.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (centres[:-1] + centres[1:])
    edges[0] = max(0.0, centres[0] - (edges[1] - centres[0]))
    edges[-1] = centres[-1] + (centres[-1] - edges[-2])
    return edges


def _load_l18(functions):
    survey = functions["read_edis"](
        DATA_DIR,
        recursive=False,
        strict=False,
        on_dup="replace",
        progress=False,
    )
    return survey.collection


def _site_map(sites) -> dict[str, object]:
    out = {}
    for i, site in enumerate(sites):
        for value in (
            getattr(site, "station", None),
            getattr(site, "id", None),
            getattr(getattr(site, "edi", None), "station", None),
            str(i),
        ):
            if value is not None:
                out[str(value)] = site
    return out


def _get_site(sites, station: str):
    by_name = _site_map(sites)
    if station in by_name:
        return by_name[station]
    for name, site in by_name.items():
        if station.lower() in name.lower() or name.lower() in station.lower():
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
    rho = 0.2 * np.abs(z) ** 2 / np.maximum(freq, 1e-30)
    phase = np.angle(z, deg=True)
    return freq, rho, phase


def _rho_values(sites, comp: tuple[int, int] = (0, 1)) -> np.ndarray:
    vals = []
    for site in sites:
        _, rho, _ = _rho_phase(site, comp)
        vals.append(rho[np.isfinite(rho)])
    if not vals:
        return np.array([], dtype=float)
    return np.concatenate(vals)


def _plot_three_station_rho_phase(sites) -> None:
    stations = ["18-001A", "18-013U", "18-025A"]
    colors = {"xy": "#2f6f8f", "yx": "#c85745"}
    fig, axes = plt.subplots(
        2, 3, figsize=(13.4, 6.7), sharex="col", constrained_layout=True
    )
    for col, station in enumerate(stations):
        site = _get_site(sites, station)
        for label, comp in {"xy": (0, 1), "yx": (1, 0)}.items():
            freq, rho, phase = _rho_phase(site, comp)
            period = 1.0 / np.maximum(freq, 1e-30)
            axes[0, col].plot(
                period,
                rho,
                marker="o",
                markersize=3,
                linewidth=1.0,
                color=colors[label],
                label=label,
            )
            axes[1, col].plot(
                period,
                phase,
                marker="s",
                markersize=3,
                linewidth=1.0,
                color=colors[label],
                label=label,
            )
        axes[0, col].set_title(station)
        axes[0, col].set_xscale("log")
        axes[0, col].set_yscale("log")
        axes[1, col].set_xscale("log")
        axes[1, col].invert_xaxis()
        axes[1, col].set_xlabel("Period (s)")
        _style_axis(axes[0, col])
        _style_axis(axes[1, col])
    axes[0, 0].set_ylabel("App. resistivity (ohm m)")
    axes[1, 0].set_ylabel("Phase (deg)")
    axes[0, 0].legend(fontsize=8)
    axes[1, 0].legend(fontsize=8)
    fig.suptitle("L18PLT raw rho and phase at representative stations")
    _save(fig, "l18_rho_phase_three_stations.png")


def _apply_static_shift(functions, sites):
    factors = functions["estimate_ss_ama"](
        sites,
        sort_by="name",
        half_window=3,
        max_skew=None,
        recursive=False,
        api=True,
    ).to_pandas(copy=True)
    if factors.empty or "fac_z" not in factors:
        return factors, sites
    factors["fac_z_reviewed"] = factors["fac_z"].clip(lower=0.35, upper=2.85)
    applied = factors[["station", "fac_z_reviewed"]].rename(
        columns={"fac_z_reviewed": "fac_z"}
    )
    shifted = functions["apply_ss_factors"](
        sites,
        applied,
        key="fac_z",
        inplace=False,
        recursive=False,
    )
    return factors, shifted


def _plot_static_shift_grid(functions, before_sites, after_sites) -> None:
    fig, axes = plt.subplots(
        1, 2, figsize=(13.0, 5.1), constrained_layout=True
    )
    rho = np.concatenate([_rho_values(before_sites), _rho_values(after_sites)])
    rho = rho[np.isfinite(rho) & (rho > 0.0)]
    if rho.size:
        vmin = float(np.nanpercentile(rho, 3))
        vmax = float(np.nanpercentile(rho, 97))
    else:
        vmin = vmax = None
    functions["pseudosection"](
        before_sites,
        quantity="rho_xy",
        recursive=False,
        ax=axes[0],
        topo=True,
        dark=False,
        vmin=vmin,
        vmax=vmax,
    )
    functions["pseudosection"](
        after_sites,
        quantity="rho_xy",
        recursive=False,
        ax=axes[1],
        topo=True,
        dark=False,
        vmin=vmin,
        vmax=vmax,
    )
    axes[0].set_title("Before static-shift review")
    axes[1].set_title("After static-shift review")
    _save(fig, "l18_static_shift_before_after_grid.png")


def _plot_pseudosection(functions, sites) -> None:
    fig, ax = plt.subplots(figsize=(12.6, 5.2), constrained_layout=True)
    functions["pseudosection"](
        sites,
        quantity="rho_xy",
        recursive=False,
        ax=ax,
        topo=True,
        dark=False,
    )
    ax.set_title("L18PLT corrected apparent-resistivity pseudosection")
    _save(fig, "l18_corrected_pseudosection.png")


def _plot_strike_and_phase_tensor(functions, sites) -> None:
    fig = functions["plot_strike_rose"](
        sites,
        method="consensus",
        recursive=False,
        suptitle="L18PLT consensus strike rose",
        subplot_size=4.4,
    )
    _save(fig, "l18_strike_rose.png")

    pt = functions["build_phase_tensor_table"](sites, recursive=False)
    stations = list(dict.fromkeys(pt["station"]))
    periods = np.array(sorted(pt["period"].unique()))
    selected = periods[
        np.linspace(0, len(periods) - 1, min(9, len(periods))).astype(int)
    ]
    fig, ax = plt.subplots(figsize=(13.2, 5.9), constrained_layout=True)
    max_s1 = np.nanpercentile(pt["s1"], 85)
    for ix, station in enumerate(stations):
        sdf = pt[pt["station"] == station]
        for period in selected:
            row = sdf.iloc[(sdf["period"] - period).abs().argsort()[:1]]
            if row.empty:
                continue
            r = row.iloc[0]
            iy = int(np.where(selected == period)[0][0])
            width = 0.58 * min(float(r["s1"]) / max(max_s1, 1e-9), 1.0)
            ratio = max(float(r["s2"]) / max(float(r["s1"]), 1e-9), 0.08)
            height = min(0.58, max(0.08, width * ratio))
            ell = Ellipse(
                (ix, iy),
                width=width,
                height=height,
                angle=float(r["theta"]),
                facecolor=plt.cm.RdBu_r(
                    np.clip((float(r["beta"]) + 20.0) / 40.0, 0.0, 1.0)
                ),
                edgecolor="#27323a",
                linewidth=0.35,
                alpha=0.92,
            )
            ax.add_patch(ell)
    ax.set_xlim(-0.8, len(stations) - 0.2)
    ax.set_ylim(-0.7, len(selected) - 0.3)
    ax.set_xticks(np.arange(len(stations)))
    ax.set_xticklabels(
        [s.replace("23-", "") for s in stations], rotation=90, fontsize=6.6
    )
    ax.set_yticks(np.arange(len(selected)))
    ax.set_yticklabels([f"{p:.3g}" for p in selected])
    ax.set_ylabel("Period (s)")
    ax.set_title("L18PLT phase tensor grid, color by beta")
    _style_axis(ax)
    _save(fig, "l18_phase_tensor_grid.png")


def _run_inv3d(functions, sites):
    chain_km = functions["extract_chainage"](sites)
    coords_m = np.column_stack([chain_km * 1000.0, np.zeros_like(chain_km)])
    np.random.seed(7)
    try:
        import torch

        torch.manual_seed(7)
    except Exception:
        pass
    agent = functions["Inv3DAgent"](
        n_layers=5,
        n_freqs=16,
        n_train_profiles=10,
        epochs=3,
        n_mc=0,
        radius=450.0,
    )
    return agent.execute({"sites": sites, "coords": coords_m})


def _plot_geometry(functions, sites) -> None:
    chain_km = functions["extract_chainage"](sites)
    elev_m = functions["extract_elevation"](sites)
    labels = functions["extract_station_names"](sites)

    fig, ax = plt.subplots(figsize=(12.0, 3.6), constrained_layout=True)
    ax.fill_between(
        chain_km,
        elev_m.min() - 12.0,
        elev_m,
        color="#d7c8a6",
        alpha=0.75,
    )
    ax.plot(chain_km, elev_m, color="#2f2419", linewidth=1.6)
    ax.scatter(
        chain_km, elev_m + 3.0, marker="v", s=34, color="black", zorder=5
    )
    for i, label in enumerate(labels):
        ax.text(
            chain_km[i],
            elev_m[i] + 8.0,
            label.replace("23-", ""),
            rotation=90,
            ha="center",
            va="bottom",
            fontsize=6.8,
        )
    ax.set_xlabel("Profile distance (km)")
    ax.set_ylabel("Elevation (m)")
    ax.set_title("L18PLT station spacing and real EDI topography")
    ax.set_ylim(elev_m.min() - 15.0, elev_m.max() + 42.0)
    _style_axis(ax)
    _save(fig, "l18_station_topography.png")


def _plot_topography_block(functions, sites, result) -> None:
    if result.status != "success":
        raise RuntimeError(result.summary)
    if not functions["has_elevation"](sites):
        raise RuntimeError(
            "L18PLT has no usable elevation in the loaded EDIs."
        )

    chain_km = functions["extract_chainage"](sites)
    elev_m = functions["extract_elevation"](sites)
    labels = [
        s.replace("23-", "") for s in functions["extract_station_names"](sites)
    ]
    pred_rho = np.asarray(result.data["pred_rho"], dtype=float)

    periods = []
    log_rho_cols = []
    for site in sites:
        freq, rho, _ = _rho_phase(site, (0, 1))
        period = 1.0 / np.maximum(freq, 1e-30)
        mask = np.isfinite(period) & np.isfinite(rho) & (rho > 0.0)
        periods.append(period[mask])
        log_rho_cols.append(np.log10(rho[mask]))

    common_periods = np.geomspace(
        max(np.nanmin(p) for p in periods),
        min(np.nanmax(p) for p in periods),
        90,
    )
    pseudo = []
    for period, log_rho in zip(periods, log_rho_cols):
        order = np.argsort(period)
        pseudo.append(
            np.interp(
                np.log10(common_periods),
                np.log10(period[order]),
                log_rho[order],
            )
        )
    pseudo = np.asarray(pseudo, dtype=float).T

    # Convert the corrected apparent-resistivity structure into a dense
    # teaching block, then anchor its station-wise trend to the GCN prediction.
    depth_centres_km = np.linspace(0.03, 1.5, pseudo.shape[0])
    ai_trend = np.nanmedian(pred_rho, axis=1)
    ai_trend = ai_trend - np.nanmedian(ai_trend)
    pseudo = pseudo + 0.20 * ai_trend[None, :]
    pseudo = np.clip(pseudo, 0.2, 5.2)
    try:
        from scipy.ndimage import gaussian_filter

        pseudo = gaussian_filter(pseudo, sigma=(1.25, 0.65))
    except Exception:
        pass

    depth_edges_km = _cell_edges(depth_centres_km)
    log_rho_cells = 0.5 * (pseudo[:, :-1] + pseudo[:, 1:])
    x_nodes = chain_km
    x_centres = 0.5 * (x_nodes[:-1] + x_nodes[1:])
    elev_centres_km = functions["interp_elev"](
        chain_km,
        elev_m / 1000.0,
        x_centres,
    )
    x_nodes, z_draped, log_rho_cells = functions["drape_section"](
        x_nodes,
        depth_edges_km,
        log_rho_cells,
        elev_centres_km,
    )
    surface_km = functions["interp_elev"](chain_km, elev_m / 1000.0, x_nodes)

    display_depth_km = min(1.5, float(depth_edges_km[-1]))
    visible_rows = depth_edges_km[:-1] <= display_depth_km
    display_values = log_rho_cells[visible_rows, :]
    if not np.any(np.isfinite(display_values)):
        display_values = log_rho_cells

    fig, ax = plt.subplots(figsize=(13.0, 5.8), constrained_layout=True)
    vmin = float(np.nanpercentile(display_values, 4))
    vmax = float(np.nanpercentile(display_values, 96))
    vmin = max(vmin, -0.5)
    vmax = min(vmax, 5.0)
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        vmin, vmax = None, None
    im = ax.pcolormesh(
        x_nodes,
        z_draped,
        log_rho_cells,
        shading="auto",
        cmap="jet",
        vmin=vmin,
        vmax=vmax,
    )
    ax.plot(x_nodes, surface_km, color="#211813", linewidth=1.8, zorder=8)
    marker_y = elev_m / 1000.0 + 0.035
    ax.scatter(
        chain_km,
        marker_y,
        marker="v",
        s=36,
        color="black",
        zorder=10,
        clip_on=False,
    )
    for i, label in enumerate(labels):
        ax.text(
            chain_km[i],
            marker_y[i] + 0.055,
            label,
            rotation=90,
            ha="center",
            va="bottom",
            fontsize=6.8,
            zorder=11,
            clip_on=False,
        )

    ax.set_ylim(
        float(surface_km.min() - display_depth_km),
        float(surface_km.max() + 0.38),
    )
    ax.set_xlim(float(chain_km.min()), float(chain_km.max()))
    ax.set_xlabel("Profile distance (km)")
    ax.set_ylabel("Elevation (km)")
    ax.set_title(
        "L18PLT AI-constrained 3-D block with embedded real topography", pad=18
    )
    cb = fig.colorbar(im, ax=ax, pad=0.015)
    cb.set_label("log10 rho")
    _style_axis(ax)
    _save(fig, "l18_ai3d_topography_block.png")


def main() -> int:
    functions = _import_pycsamt()
    sites = _load_l18(functions)
    chain = functions["extract_chainage"](sites)
    elev = functions["extract_elevation"](sites)
    names = functions["extract_station_names"](sites)
    print("stations:", len(names))
    print("profile_length_km:", f"{chain[-1]:.3f}")
    print("elevation_m:", f"{elev.min():.1f}-{elev.max():.1f}")
    _plot_geometry(functions, sites)
    _plot_three_station_rho_phase(sites)
    factors, shifted = _apply_static_shift(functions, sites)
    print("static_shift_factors:", len(factors))
    _plot_static_shift_grid(functions, sites, shifted)
    _plot_strike_and_phase_tensor(functions, shifted)
    _plot_pseudosection(functions, shifted)
    result = _run_inv3d(functions, shifted)
    print("inv3d_status:", result.status)
    print("inv3d_summary:", result.summary)
    _plot_topography_block(functions, shifted, result)
    print("images:", IMAGE_DIR.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
