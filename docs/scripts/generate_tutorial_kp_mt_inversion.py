"""Prepare KAP03 ModEM inputs and optional MT3D-AI inversion artifacts."""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
EDI_DIR = ROOT / "results" / "kap03_mt_tutorial" / "edi_conditioned_rotated"
RUN_DIR = ROOT / "results" / "kap03_mt_tutorial" / "inversion"
MODEM_DIR = RUN_DIR / "modem3d"
AI_DIR = RUN_DIR / "ai_mt3d"
AI2D_DIR = RUN_DIR / "ai_mt2d_tri"
IMAGE_DIR = ROOT / "docs" / "source" / "images" / "tutorials" / "condition_mt_line_with_tipper_and_rotation"
TOPO_CSV = ROOT / "data" / "MT" / "kap03_topography_open_meteo.csv"


def load_corrected_sites():
    from pycsamt.api import read_edis

    survey = read_edis(EDI_DIR, recursive=False, strict=True, progress=False)
    return survey.collection


def build_modem_inputs(sites):
    """Write and reload the complete classical ModEM input set."""
    from pycsamt.models.modem import InputBuilder, ModEmConfig, ModEmRunner
    from pycsamt.models.modem.model3d import ModEmModel3D

    cfg = ModEmConfig(
        mode="3d",
        component_type="Full_Impedance",
        error_floor_z=0.08,
        freq_min=5.0e-5,
        freq_max=0.04,
        cell_size_h=25_000.0,
        n_padding_xy=5,
        nz=28,
        n_airlayers=6,
        cell_size_v_top=250.0,
        depth_scale=1.25,
        initial_rho=100.0,
        smooth_x=0.3,
        smooth_y=0.3,
        smooth_z=0.2,
        n_smooth_iter=2,
        max_iterations=100,
        target_rms=1.05,
        use_mpi=True,
        n_procs=8,
    )
    builder = InputBuilder(config=cfg)
    files = builder.build(
        sites,
        workdir=MODEM_DIR,
        data_filename="KAP03_ModEM.dat",
        model_filename="KAP03_m0.ws",
        cov_filename="KAP03.cov",
        ctrl_filename="KAP03.inv",
    )
    model = ModEmModel3D.read(files["model"])
    command = ModEmRunner(MODEM_DIR, config=cfg).command(
        files["model"].name,
        files["data"].name,
        files["control"].name,
        covariance=files["covariance"].name,
    )
    manifest = {
        "stations": len(list(sites)),
        "model_shape_nz_ny_nx": list(model.shape),
        "earth_depth_m": float(np.sum(model.z_widths[model.n_air :])),
        "air_layers": int(model.n_air),
        "files": {key: value.name for key, value in files.items()},
        "dry_run_command": command,
        "tipper_note": "Current ModEM builder writes full impedance; retain EDI tipper for independent QC.",
    }
    MODEM_DIR.mkdir(parents=True, exist_ok=True)
    (MODEM_DIR / "build_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))
    return builder, manifest


def plot_modem_mesh(builder):
    """Audit horizontal coverage, padding, air, and earth discretization."""
    model, data = builder.model, builder.data
    x = model.x_nodes / 1000.0
    y = model.y_nodes / 1000.0
    z = model.z_nodes / 1000.0
    coords = np.asarray(list(data.site_coords.values()), dtype=float)
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.2), constrained_layout=True)
    for value in x:
        axes[0].axvline(value, color="0.75", lw=0.35)
    for value in y:
        axes[0].axhline(value, color="0.75", lw=0.35)
    left_x = np.sum(model.x_widths[: builder.config.n_padding_xy]) / 1000.0
    left_y = np.sum(model.y_widths[: builder.config.n_padding_xy]) / 1000.0
    station_x = left_x + (coords[:, 0] - coords[:, 0].min()) / 1000.0
    station_y = left_y + (coords[:, 1] - coords[:, 1].min()) / 1000.0
    axes[0].scatter(station_x, station_y,
                    marker="v", s=28, color="#b2182b", edgecolor="k", zorder=3)
    axes[0].set(xlabel="Local northing mesh coordinate (km)",
                ylabel="Local easting mesh coordinate (km)",
                title=f"ModEM horizontal mesh: {model.nx} × {model.ny} cells")
    for value in z:
        axes[1].axhline(value, color="0.25", lw=0.45)
    axes[1].axhspan(0, z[model.n_air], color="#d9edf7", alpha=0.8, label="air layers")
    axes[1].invert_yaxis()
    axes[1].set(xlim=(0, 1), xticks=[], ylabel="Cumulative mesh coordinate (km)",
                title=f"Vertical mesh: {model.n_air} air + {model.nz-model.n_air} earth layers")
    axes[1].legend(fontsize=8)
    fig.savefig(IMAGE_DIR / "kp_modem3d_mesh_audit.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def build_profile_triangle_mesh():
    """Build an optional topography-following 2-D triangular realization."""
    from pycsamt.api import draw_tri_mesh
    from pycsamt.forward.maxwell.tri_mesh_gen import build_graded_tri_mesh

    topo = pd.read_csv(TOPO_CSV)
    lat = np.deg2rad(topo.latitude.to_numpy())
    lon = np.deg2rad(topo.longitude.to_numpy())
    seg = 2 * 6_371_008.8 * np.arcsin(np.sqrt(
        np.sin(np.diff(lat) / 2) ** 2
        + np.cos(lat[:-1]) * np.cos(lat[1:]) * np.sin(np.diff(lon) / 2) ** 2
    ))
    chain = np.r_[0.0, np.cumsum(seg)]
    surface_depth = topo.elevation_m.max() - topo.elevation_m.to_numpy()
    pad = 50_000.0
    mesh = build_graded_tri_mesh(
        (-pad, float(chain[-1] + pad)),
        (float(surface_depth.min()), 300_000.0),
        chain,
        surface_cell_m=15_000.0,
        growth_rate=1.35,
        max_cell_m=80_000.0,
        min_angle=28.0,
        topo_x_m=chain,
        topo_z_m=surface_depth,
    )
    from matplotlib.ticker import FuncFormatter

    fig, axes = plt.subplots(
        2, 1, figsize=(13.0, 7.2), constrained_layout=True,
        gridspec_kw={"height_ratios": (3.0, 1.15)},
    )
    station_z = np.interp(chain, chain, surface_depth)
    km = FuncFormatter(lambda value, _pos: f"{value / 1000:g}")
    for ax in axes:
        draw_tri_mesh(ax, mesh, preset="diagram")
        ax.plot(chain, surface_depth, color="#2f6f8f", lw=1.4, zorder=4)
        ax.scatter(chain, station_z, marker="v", color="#b2182b", s=26, zorder=5)
        ax.xaxis.set_major_formatter(km)
        ax.yaxis.set_major_formatter(km)
        ax.set_xlabel("Profile distance (km)")
        ax.set_ylabel("Depth below elevation datum (km)")
        ax.invert_yaxis()
    axes[0].set_title(f"Optional 2-D triangular mesh: {mesh.n_triangles:,} elements")
    axes[1].set_ylim(5_000.0, -250.0)
    axes[1].set_title("Shallow 5 km zoom: measured elevation controls the surface boundary")
    fig.savefig(IMAGE_DIR / "kp_optional_triangular_profile_mesh.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"triangle_nodes={mesh.n_nodes} triangle_elements={mesh.n_triangles}")
    return mesh


def profile_geometry():
    """Return KAP03 chainage and depth below the maximum-elevation datum."""
    topo = pd.read_csv(TOPO_CSV)
    lat = np.deg2rad(topo.latitude.to_numpy())
    lon = np.deg2rad(topo.longitude.to_numpy())
    seg = 2 * 6_371_008.8 * np.arcsin(np.sqrt(
        np.sin(np.diff(lat) / 2) ** 2
        + np.cos(lat[:-1]) * np.cos(lat[1:]) * np.sin(np.diff(lon) / 2) ** 2
    ))
    chain = np.r_[0.0, np.cumsum(seg)]
    surface_depth = topo.elevation_m.max() - topo.elevation_m.to_numpy()
    return chain, surface_depth


def run_ai_mt2d_tri(sites, *, n_train_profiles: int = 20,
                    epochs: int = 50, patience: int = 8):
    """Run topography-following triangular MT2D training and field inference.

    Shallower than the first pass at this tutorial (100 km instead of
    250 km -- KAP03's own frequency band, 5e-5 to 0.04 Hz, does not
    resolve structure anywhere near 250 km with useful confidence in the
    first place), but the surface/field cell sizes are now *finer* than
    both earlier passes: ``mesh_target_cell_m=15,000`` matches the
    standalone display mesh built by :func:`build_profile_triangle_mesh`
    rather than the coarser 45,000 m used previously, so the triangles
    actually driving GCN training are no bigger than the ones already
    shown as the "optional" mesh figure.

    ``gcn_adjacency_radius_m`` is now set explicitly to 30,000 m (2x
    ``mesh_target_cell_m``). It was never set in any earlier pass of
    this tutorial, silently leaving ``Inv2DAgent``'s 300 m library
    default in effect -- 15-40x smaller than the actual median distance
    between neighbouring triangle centroids on this mesh (~10 km).
    ``build_adjacency`` degrades to the identity matrix below that
    distance, so the GCN never did any spatial message-passing in any
    earlier run of this branch; every triangle nearest a given station
    predicted the same value regardless of its own depth, which is what
    produced the flat vertical "curtain" per station visible in the
    comparison figure before this fix.
    """
    from pycsamt.agents import Inv2DAgent
    from pycsamt.forward.maxwell.tri_fem2d import TriFEM2DAdapter

    chain, surface_depth = profile_geometry()
    n_sites = len(list(sites))
    spacing = float(chain[-1] / max(n_sites - 1, 1))
    np.random.seed(23)
    try:
        import torch

        torch.manual_seed(23)
    except ImportError:
        pass
    agent = Inv2DAgent(
        physics="mt2d_tri",
        epochs=epochs,
        patience=patience,
        n_freqs=8,
        depth_max=100_000.0,
        n_train_profiles=n_train_profiles,
        n_stations_per_profile=n_sites,
        station_spacing_m=spacing,
        mesh_target_cell_m=15_000.0,
        field_grid_cell_m=7_500.0,
        gcn_adjacency_radius_m=30_000.0,
        correlation_length_x_m=(40_000.0, 180_000.0),
        correlation_length_z_m=(5_000.0, 25_000.0),
        gcn_hidden=(128, 64, 32),
        topo_x_m=chain,
        topo_z_m=surface_depth,
        mare2dem_adapter=TriFEM2DAdapter(),
    )
    AI2D_DIR.mkdir(parents=True, exist_ok=True)
    result = agent.execute({
        "sites": sites,
        "freqs": np.geomspace(5.0e-5, 0.04, 8),
        "output_dir": str(AI2D_DIR),
    })
    if result.status != "success":
        raise RuntimeError(result.summary)
    recovery = result.data.get("mt2d_tri_recovery") or {}
    history = result.data.get("training_history") or {}
    pred = result.data["pred_triangles"]
    station_x = np.arange(n_sites, dtype=float) * spacing
    station_names = np.array([s.station for s in sites])
    np.savez_compressed(
        AI2D_DIR / "mesh_prediction.npz",
        nodes_m=np.asarray(pred["mesh"].nodes_m),
        triangles=np.asarray(pred["mesh"].triangles),
        log10_resistivity=np.asarray(pred["log10_resistivity"]),
        station_x_m=station_x,
        station_z_m=np.interp(station_x, chain, surface_depth),
        station_names=station_names,
    )
    manifest = {
        "status": result.status,
        "epochs_requested": epochs,
        "early_stopping_patience": patience,
        "training_realizations": n_train_profiles,
        "epochs_completed": int(result.data.get("epochs_completed", 0)),
        "best_validation_loss": float(result.data.get("best_validation_loss", np.nan)),
        "held_out_recovery": recovery,
        "training_history": history,
    }
    (AI2D_DIR / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, default=float), encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, default=float))
    return result


def run_ai_mt3d(sites, *, n_train_profiles: int = 20, epochs: int = 50, patience: int = 8):
    """Run the structured-mesh MT3D training branch on request.

    Shallower than the first pass, for the same reason as
    :func:`run_ai_mt2d_tri`: 100 km depth instead of 250 km. The lever
    that actually thins the *plotted* MT3D section is ``n_layers`` (10
    depth rows instead of 6; the 26 station columns are already the
    real station count and cannot be finer) -- ``geology_grid_nx_ny``,
    ``geology_grid_nz``, and ``max_mesh_cells`` instead control the
    realism/cost of the *training-data* forward solves and were left
    close to their original values on purpose. Two smoke tests found a
    real 3-D Maxwell solve's cost on this solver scales far worse than
    linearly with cell count: quadrupling ``max_mesh_cells`` to 20,000
    multiplied per-realization time by ~17x (not 4x), and even a modest
    5,000-to-8,000 bump (with ``geology_grid_nx_ny=4``,
    ``geology_grid_nz=5``) multiplied it by ~5.8x -- a 200-realization
    run at that setting would take about 5.4 hours, not the ~2-2.5
    hours a naive N^2 estimate suggested. ``max_mesh_cells=6,000``
    (a 20% bump) with the original ``geology_grid_nx_ny=3``,
    ``geology_grid_nz=4`` keeps training-data cost close to the
    original run's while still giving every realization a slightly
    finer solver core.
    """
    from pycsamt.agents import Inv3DAgent

    np.random.seed(17)
    try:
        import torch

        torch.manual_seed(17)
    except ImportError:
        pass
    agent = Inv3DAgent(
        physics="mt3d",
        n_layers=10,
        freqs=np.geomspace(5.0e-5, 0.04, 8),
        depth_max=100_000.0,
        n_train_profiles=n_train_profiles,
        epochs=epochs,
        patience=patience,
        radius=120_000.0,
        hidden=(128, 64, 32),
        dropout=0.1,
        n_mc=20,
        correlation_length_x_m=(40_000.0, 180_000.0),
        correlation_length_y_m=(40_000.0, 150_000.0),
        correlation_length_z_m=(5_000.0, 25_000.0),
        geology_grid_nx_ny=3,
        geology_grid_nz=4,
        mesh_safety_factor=8.0,
        max_mesh_cells=6_000,
    )
    result = agent.execute({
        "sites": sites,
        "topography": True,
        "output_dir": str(AI_DIR),
    })
    if result.status != "success":
        raise RuntimeError(result.summary)
    recovery = result.data.get("mt3d_recovery") or {}
    station_names = np.array([s.station for s in sites])
    manifest = {
        "status": result.status,
        "epochs_requested": epochs,
        "early_stopping_patience": patience,
        "training_realizations": n_train_profiles,
        "rms_global": float(result.data["rms_global"]),
        "held_out_recovery": recovery,
        "epochs_completed": int(result.data.get("epochs_completed", 0)),
        "best_validation_loss": float(result.data.get("best_validation_loss", np.nan)),
        "training_history": result.data.get("training_history") or {},
    }
    AI_DIR.mkdir(parents=True, exist_ok=True)
    (AI_DIR / "run_manifest.json").write_text(json.dumps(manifest, indent=2, default=float), encoding="utf-8")
    np.savez_compressed(
        AI_DIR / "mesh_prediction.npz",
        log10_resistivity=np.asarray(result.data["pred_rho"]),
        depths_km=np.asarray(result.data["depths_km"]),
        chainage_km=np.asarray(result.data.get("station_chainage_km")),
        elevation_m=np.asarray(result.data.get("station_elevation_m")),
        station_names=station_names,
    )
    print(json.dumps(manifest, indent=2, default=float))
    return result


def _add_station_labels(ax, x_km, z_km, names):
    """Thin station-name labels above markers (km-scale axes).

    Same :class:`~pycsamt.api.station.StationAxisStyle` thinning logic
    :func:`pycsamt.agents.inv2d_agent.Inv2DAgent`'s own ``mt2d_tri`` figure
    uses, reserving real data-space headroom above the plotted depth range
    rather than relying on ``clip_on=False`` plus external margins -- the
    latter is what silently failed for that figure's first fix attempt.
    """
    from pycsamt.api.station import StationAxisStyle

    y0, y1 = ax.get_ylim()
    span = abs(y1 - y0)
    gap = 0.02 * span
    room = 0.24 * span
    top_val = min(y0, y1) - room
    bottom_val = max(y0, y1)
    if y0 > y1:
        ax.set_ylim(bottom_val, top_val)
    else:
        ax.set_ylim(top_val, bottom_val)
    visible = StationAxisStyle().label_indices(
        list(names), figwidth_in=ax.figure.get_figwidth()
    )
    for i in visible:
        ax.annotate(
            str(names[i]),
            (x_km[i], z_km[i] - gap),
            ha="center",
            va="bottom",
            fontsize=7,
            rotation=90,
            zorder=6,
        )


def plot_ai_comparison(mt2d_result=None, mt3d_result=None):
    """Make the requested 2 x 3 inversion/validation comparison."""
    rows = [
        ("MT2D triangular AI", AI2D_DIR / "mesh_prediction.npz",
         AI2D_DIR / "run_manifest.json", mt2d_result),
        ("MT3D structured AI", AI_DIR / "mesh_prediction.npz",
         AI_DIR / "run_manifest.json", mt3d_result),
    ]
    if not all(mesh_path.exists() and manifest_path.exists()
               for _, mesh_path, manifest_path, _ in rows):
        return
    fig = plt.figure(figsize=(16.5, 9.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, width_ratios=(1.0, 1.0, 0.82))
    model_axes = []
    for row, (title, mesh_path, manifest_path, result) in enumerate(rows):
        ax_model = fig.add_subplot(gs[row, :2])
        model_axes.append(ax_model)
        mesh_data = np.load(mesh_path)
        if row == 0:
            import matplotlib.tri as mtri

            nodes = mesh_data["nodes_m"] / 1000.0
            triangulation = mtri.Triangulation(
                nodes[:, 0], nodes[:, 1], mesh_data["triangles"]
            )
            # Clip the *display* range only -- a handful of triangles
            # nearest a few far-profile stations (real, uncorrected field
            # apparent resistivity/phase well outside anything the
            # synthetic training distribution ever produced) extrapolate
            # to physically nonsensical log10(rho) values in the tens; a
            # few such outliers would otherwise wash out the legitimate
            # depth-varying structure across the rest of the section.
            # ``log10_resistivity`` in the saved .npz is untouched.
            clip_lo, clip_hi = -1.0, 5.0
            log_rho = mesh_data["log10_resistivity"]
            n_clipped = int(np.sum((log_rho < clip_lo) | (log_rho > clip_hi)))
            artist = ax_model.tripcolor(
                triangulation, facecolors=log_rho,
                cmap="turbo_r", shading="flat", edgecolors="0.18",
                linewidth=0.28, vmin=clip_lo, vmax=clip_hi,
            )
            ax_model.plot(mesh_data["station_x_m"] / 1000.0,
                          mesh_data["station_z_m"] / 1000.0,
                          "kv", ms=4, mfc="white", label="MT stations")
            ax_model.set(xlabel="Profile distance (km)",
                         ylabel="Depth below elevation datum (km)")
            ax_model.invert_yaxis()
            ax_model.legend(loc="lower right", fontsize=8)
            if n_clipped:
                ax_model.text(
                    0.99, 0.02,
                    f"{n_clipped}/{log_rho.size} triangles off-scale"
                    f" (clipped to [{clip_lo:g}, {clip_hi:g}])",
                    transform=ax_model.transAxes, ha="right", va="bottom",
                    fontsize=7.5, color="0.25",
                    bbox={"facecolor": "white", "edgecolor": "none",
                          "alpha": 0.75, "pad": 1.5},
                )
            if "station_names" in mesh_data:
                _add_station_labels(
                    ax_model,
                    mesh_data["station_x_m"] / 1000.0,
                    mesh_data["station_z_m"] / 1000.0,
                    mesh_data["station_names"],
                )
        else:
            values = mesh_data["log10_resistivity"].T
            x = mesh_data["chainage_km"]
            z = mesh_data["depths_km"]
            x_edges = _centres_to_edges(x)
            z_edges = _centres_to_edges(z, lower_bound=0.0)
            elevation = mesh_data["elevation_m"]
            surface_depth = (np.nanmax(elevation) - elevation) / 1000.0
            surface_edges = _centres_to_edges(surface_depth)
            x_grid = np.broadcast_to(x_edges[None, :],
                                     (z_edges.size, x_edges.size))
            z_grid = z_edges[:, None] + surface_edges[None, :]
            artist = ax_model.pcolormesh(
                x_grid, z_grid, values, cmap="turbo_r", shading="flat",
                edgecolors=(0.08, 0.08, 0.08, 0.55), linewidth=0.32,
            )
            ax_model.plot(x, surface_depth, "k-", lw=1.0, zorder=4)
            ax_model.plot(x, surface_depth, "kv", ms=4, mfc="white",
                          label="MT stations")
            ax_model.set(xlabel="Profile distance (km)",
                         ylabel="Depth below elevation datum (km)")
            ax_model.invert_yaxis()
            ax_model.legend(loc="lower right", fontsize=8)
            if "station_names" in mesh_data:
                _add_station_labels(ax_model, x, surface_depth, mesh_data["station_names"])
        fig.colorbar(
            artist, ax=ax_model, pad=0.015,
            extend="both" if (row == 0 and n_clipped) else "neither",
            label=r"$\log_{10}\rho$ ($\Omega\,m$)",
        )
        ax_model.text(
            0.01, 0.985, title, transform=ax_model.transAxes,
            ha="left", va="top", fontsize=12, fontweight="bold",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78,
                  "pad": 2.0}, zorder=8,
        )
        audit = json.loads(manifest_path.read_text(encoding="utf-8"))
        history = ((result.data.get("training_history") or {}) if result
                   else (audit.get("training_history") or {}))
        ax_loss = fig.add_subplot(gs[row, 2])
        train = np.asarray(history.get("train_loss", []), dtype=float)
        valid = np.asarray(history.get("val_loss", []), dtype=float)
        if train.size:
            ax_loss.plot(np.arange(1, train.size + 1), train, color="#1f77b4",
                         lw=1.8, label="training")
        if valid.size:
            ax_loss.plot(np.arange(1, valid.size + 1), valid, color="#c43c39",
                         lw=1.8, label="validation")
            best = int(np.nanargmin(valid)) + 1
            ax_loss.axvline(best, color="0.25", ls="--", lw=1.0,
                            label=f"best epoch {best}")
        completed = int(audit.get("epochs_completed", max(train.size, valid.size)))
        requested = int(audit.get("epochs_requested", 50))
        ax_loss.set(xlabel="Epoch", ylabel="Loss",
                    title=f"Validation audit: {completed}/{requested} epochs")
        ax_loss.grid(alpha=0.22)
        if train.size or valid.size:
            ax_loss.legend(fontsize=8)

    # Both model panels compress ~100 km of depth into roughly the same
    # on-screen width as ~1,500 km of profile -- real, close-to-isotropic
    # triangles (independently verified: median width/height ratio 0.99 on
    # the MT2D mesh) therefore render visibly "tall." Rather than let that
    # read as a mesh defect, compute the actual exaggeration from the
    # rendered axes geometry (not an assumed constant) and label it
    # explicitly, the same disclosure convention geological cross-sections
    # use for any vertically exaggerated section.
    fig.canvas.draw()
    for ax_model in model_axes:
        bbox_in = ax_model.get_window_extent().transformed(
            fig.dpi_scale_trans.inverted()
        )
        xlim = ax_model.get_xlim()
        ylim = ax_model.get_ylim()
        data_w = abs(xlim[1] - xlim[0])
        data_h = abs(ylim[1] - ylim[0])
        km_per_in_x = data_w / bbox_in.width
        km_per_in_y = data_h / bbox_in.height
        vert_exag = km_per_in_x / km_per_in_y
        ax_model.text(
            0.99, 0.985, f"vertical exaggeration ≈ {vert_exag:.1f}x",
            transform=ax_model.transAxes, ha="right", va="top", fontsize=8,
            color="0.2",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78,
                  "pad": 1.8}, zorder=8,
        )
    fig.savefig(IMAGE_DIR / "kp_ai_mt2d_mt3d_comparison.png", dpi=190)
    plt.close(fig)


def _centres_to_edges(values, *, lower_bound=None):
    """Convert monotonic cell centres to pcolormesh edges."""
    values = np.asarray(values, dtype=float)
    if values.size == 1:
        half = max(abs(values[0]) * 0.5, 0.5)
        edges = np.array([values[0] - half, values[0] + half])
    else:
        middle = 0.5 * (values[:-1] + values[1:])
        edges = np.r_[values[0] - (middle[0] - values[0]), middle,
                      values[-1] + (values[-1] - middle[-1])]
    if lower_bound is not None:
        edges[0] = lower_bound
    return edges


def capture_ai_outputs():
    """Copy executed AI figures and make a compact quantitative audit."""
    mapping = {
        "inv3d_uncertainty.png": "kp_ai_mt3d_uncertainty.png",
        "inv3d_topography_section.png": "kp_ai_mt3d_topography_section.png",
    }
    for source, target in mapping.items():
        path = AI_DIR / source
        if path.exists():
            shutil.copy2(path, IMAGE_DIR / target)
    manifest_path = AI_DIR / "run_manifest.json"
    if not manifest_path.exists():
        return
    audit = json.loads(manifest_path.read_text(encoding="utf-8"))
    recovery = audit.get("held_out_recovery") or {}
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0), constrained_layout=True)
    axes[0].bar(["Observed\nresponse"], [audit.get("rms_global", np.nan)], color="#2f6f8f")
    axes[0].axhline(1.0, color="0.25", ls="--", lw=1.0)
    axes[0].set(title="Field-response diagnostic", ylabel="Global RMS")
    axes[1].bar(
        ["RMSE", "MAE"],
        [recovery.get("rmse", np.nan), recovery.get("mae", np.nan)],
        color=["#c85745", "#e0a458"],
    )
    axes[1].set(
        title=f"Held-out recovery (n={recovery.get('n_samples', 0)}, R²={recovery.get('r2', np.nan):.2f})",
        ylabel=r"Error in $\log_{10}\rho$",
    )
    fig.savefig(IMAGE_DIR / "kp_ai_mt3d_validation_audit.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-ai", action="store_true", help="Run the costly 50-epoch MT3D training branch.")
    parser.add_argument("--run-ai2d", action="store_true", help="Run triangular-mesh MT2D AI training.")
    parser.add_argument("--n-train-profiles", type=int, default=100)
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--patience", type=int, default=15)
    args = parser.parse_args(argv)
    sites = load_corrected_sites()
    builder, _ = build_modem_inputs(sites)
    plot_modem_mesh(builder)
    build_profile_triangle_mesh()
    result2d = None
    result3d = None
    if args.run_ai2d:
        result2d = run_ai_mt2d_tri(
            sites, n_train_profiles=args.n_train_profiles,
            epochs=args.epochs, patience=args.patience,
        )
    if args.run_ai:
        result3d = run_ai_mt3d(
            sites, n_train_profiles=args.n_train_profiles,
            epochs=args.epochs, patience=args.patience,
        )
    capture_ai_outputs()
    plot_ai_comparison(result2d, result3d)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
