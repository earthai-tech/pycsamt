"""Run and plot the two-line Maxwell 3-D AI workflow used by the porphyry tutorial."""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.agents import Inv3DAgent  # noqa: E402
from pycsamt.ai.geology import (  # noqa: E402
    ElectricalLayer,
    EllipsoidalLens,
    GaussianCorrelation,
    GeologyGrid,
    generate_layered_geology,
    insert_lenses,
)
from pycsamt.api.mesh import PYCSAMT_MESH, cell_edges_from_centres, draw_mesh  # noqa: E402
from pycsamt.emtools import ensure_sites  # noqa: E402
from pycsamt.forward.maxwell import MeshDesign, build_solver_mesh  # noqa: E402
from pycsamt.topo import (  # noqa: E402
    drape_section,
    extract_chainage,
    extract_elevation,
    extract_station_names,
    interp_elev,
)

IMAGE_DIR = ROOT / "docs/source/images/tutorials/map_porphyry_mineralization_from_noisy_amt"
RUN_DIR = ROOT / "runs"

FAST_FREQUENCIES_HZ = np.geomspace(1.0, 1000.0, 8)
DEPTH_MAX_M = 1800.0
LINE_OFFSET_M = 500.0  # stated assumption: no real cross-line coordinate exists


def save(fig: plt.Figure, name: str) -> None:
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(IMAGE_DIR / name, dpi=190, bbox_inches="tight")
    plt.close(fig)


def load_lines():
    lines = {}
    for name in ("L26PLT", "L30PLT"):
        corrected = RUN_DIR / f"{name}_corrected"
        if not corrected.exists():
            raise FileNotFoundError(
                f"{corrected} is missing; run the tutorial's correction/export steps first."
            )
        lines[name] = ensure_sites(corrected, recursive=False, verbose=0).ordered()
    return lines


def build_combined_problem(lines):
    sites26, sites30 = lines["L26PLT"], lines["L30PLT"]
    chain26_m = extract_chainage(sites26) * 1000.0
    chain30_m = extract_chainage(sites30) * 1000.0
    coords = np.column_stack(
        [
            np.concatenate([chain26_m, chain30_m]),
            np.concatenate(
                [np.zeros_like(chain26_m), np.full_like(chain30_m, LINE_OFFSET_M)]
            ),
        ]
    )
    combined_sites = sites26.to_edis() + sites30.to_edis()

    x_span = float(coords[:, 0].max() - coords[:, 0].min())
    grid = GeologyGrid.regular_3d(
        nx=24,
        ny=6,
        nz=18,
        dx_m=x_span / 24.0,
        dy_m=(LINE_OFFSET_M + 300.0) / 6.0,
        dz_m=DEPTH_MAX_M / 18.0,
        x_origin_m=-150.0,
        y_origin_m=-150.0,
    )
    # A seeded display prior: this step is a mesh/geometry illustration only,
    # not the training distribution -- Inv3DAgent(physics="mt3d") builds its
    # own internal GeologyGrid from correlated fields; no agent constructor
    # accepts an externally-built grid.
    correlation = GaussianCorrelation(700.0, 160.0, length_y_m=500.0)
    layers = [
        ElectricalLayer("weathered cover", 45.0, 0.10, correlation),
        ElectricalLayer("altered intrusive", 350.0, 0.16, correlation),
        ElectricalLayer("fresh intrusive", 1600.0, 0.10, correlation),
    ]
    base = generate_layered_geology(
        grid, layers, [320.0, 1050.0], seed=2601,
        interface_relief_std_m=[45.0, 100.0],
        interface_correlation=correlation,
        minimum_thickness_m=120.0,
    )
    target = EllipsoidalLens(
        "combined-line conductive target",
        center_x_m=1500.0, center_y_m=LINE_OFFSET_M / 2.0, center_z_m=760.0,
        radius_x_m=480.0, radius_y_m=350.0, radius_z_m=210.0,
        resistivity_ohm_m=12.0, dip_deg=15.0, transition_fraction=0.20,
    )
    geology = insert_lenses(base, [target])
    solver_model = build_solver_mesh(
        grid,
        resistivity_ohm_m=geology.resistivity_ohm_m,
        frequencies_hz=FAST_FREQUENCIES_HZ,
        design=MeshDesign(
            horizontal_padding_cells=3,
            bottom_padding_cells=4,
            air_layers=3,
            padding_expansion=1.3,
        ),
    )
    return {
        "lines": lines,
        "sites26": sites26,
        "sites30": sites30,
        "chain26_m": chain26_m,
        "chain30_m": chain30_m,
        "coords": coords,
        "combined_sites": combined_sites,
        "grid": grid,
        "solver_model": solver_model,
    }


def plot_geology_and_mesh(problem) -> None:
    mesh_model = problem["solver_model"]
    centres = mesh_model.mesh.cell_centres_m
    x_km = centres["x"] / 1000.0
    y_km = centres["y"] / 1000.0
    z_km = centres["z"] / 1000.0
    x_edges = cell_edges_from_centres(x_km)
    y_edges = cell_edges_from_centres(y_km)
    z_edges = cell_edges_from_centres(z_km)
    log_rho = np.log10(1.0 / mesh_model.conductivity_s_m)
    mesh_style = PYCSAMT_MESH.style_for("review")
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.0), constrained_layout=True)

    xz, _ = draw_mesh(
        axes[0], x_edges, z_edges, log_rho[:, log_rho.shape[1] // 2, :],
        style=mesh_style, cmap="turbo",
    )
    axes[0].invert_yaxis()
    axes[0].set(
        title="Along-line slice (padded mesh)",
        xlabel="x, profile distance (km)", ylabel="depth (km)",
    )

    draw_mesh(
        axes[1], x_edges, y_edges, log_rho[log_rho.shape[0] // 3, :, :],
        style=mesh_style, cmap="turbo",
    )
    axes[1].set(
        title="Horizontal slice, shallow",
        xlabel="x, profile distance (km)", ylabel="y, cross-line (km)",
    )

    draw_mesh(
        axes[2], y_edges, z_edges, log_rho[:, :, log_rho.shape[2] // 2],
        style=mesh_style, cmap="turbo",
    )
    axes[2].invert_yaxis()
    axes[2].set(
        title="Cross-line slice",
        xlabel="y, cross-line (km)", ylabel="depth (km)",
    )
    fig.colorbar(xz, ax=axes, shrink=0.85, label=r"$\log_{10}\rho$ (ohm m)")
    save(fig, "willy_ai3d_geology_maxwell_both_lines.png")


def run_fast_inversion(problem):
    np.random.seed(17)
    try:
        import torch

        torch.manual_seed(17)
    except ImportError:
        pass

    output_dir = RUN_DIR / "willy_ai3d_maxwell"
    output_dir.mkdir(parents=True, exist_ok=True)
    started = time.time()
    agent = Inv3DAgent(
        physics="mt3d",
        n_layers=6,
        freqs=FAST_FREQUENCIES_HZ,
        depth_max=DEPTH_MAX_M,
        n_train_profiles=10,
        epochs=20,
        radius=450.0,
        hidden=(64, 32),
        dropout=0.1,
        n_mc=0,
        geology_grid_nx_ny=4,
        geology_grid_nz=4,
        max_mesh_cells=60_000,
    )
    result = agent.execute(
        {
            "sites": problem["combined_sites"],
            "coords": problem["coords"],
            "topography": True,
            "output_dir": str(output_dir),
        }
    )
    if result.status != "success":
        raise RuntimeError(f"3-D fast run failed: {result.summary}")
    checkpoint = output_dir / "gcn_inverter_3d_maxwell.npz"
    result.data["inverter"].save(checkpoint)
    recovery = result.data.get("mt3d_recovery") or {}
    manifest = {
        "status": result.status,
        "summary": result.summary,
        "seed": 17,
        "n_train_profiles": 10,
        "epochs_requested": 20,
        "n_layers": 6,
        "depth_max_m": DEPTH_MAX_M,
        "line_offset_m": LINE_OFFSET_M,
        "frequencies_hz": FAST_FREQUENCIES_HZ.tolist(),
        "rms_global": float(result.data["rms_global"]),
        "mt3d_recovery": recovery,
        "elapsed_seconds": time.time() - started,
        "checkpoint": checkpoint.name,
    }
    (output_dir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, default=float), encoding="utf-8"
    )
    print(json.dumps(manifest, default=float), flush=True)
    return result


def _cell_edges_km(centres_m: np.ndarray) -> np.ndarray:
    return cell_edges_from_centres(centres_m) / 1000.0


def plot_inversion_with_mesh(problem, result) -> None:
    pred_rho = np.asarray(result.data["pred_rho"], dtype=float)
    depths_km = np.asarray(result.data["depths_km"], dtype=float)
    n26 = len(problem["sites26"])
    pred26, pred30 = pred_rho[:n26], pred_rho[n26:]
    chain26_km = problem["chain26_m"] / 1000.0
    chain30_km = problem["chain30_m"] / 1000.0
    elev26 = extract_elevation(problem["sites26"])
    elev30 = extract_elevation(problem["sites30"])
    labels26 = [n.split("-")[-1] for n in extract_station_names(problem["sites26"])]
    labels30 = [n.split("-")[-1] for n in extract_station_names(problem["sites30"])]
    z_edges_km = cell_edges_from_centres(depths_km)

    fig, axes = plt.subplots(1, 2, figsize=(15.0, 5.6), constrained_layout=True)
    mesh_style = PYCSAMT_MESH.style_for("review")
    fills = []
    for ax, chain_km, elev_m, labels, pred, name in zip(
        axes,
        (chain26_km, chain30_km),
        (elev26, elev30),
        (labels26, labels30),
        (pred26, pred30),
        ("L26PLT", "L30PLT"),
    ):
        x_edges_km = cell_edges_from_centres(chain_km)
        x_centres = 0.5 * (x_edges_km[:-1] + x_edges_km[1:])
        elev_centres_km = interp_elev(chain_km, elev_m / 1000.0, x_centres)
        data = pred.T  # (n_layers, n_sta)
        x_nodes, z_draped, data_draped = drape_section(
            x_edges_km, z_edges_km, data, elev_centres_km
        )
        x_2d = np.broadcast_to(x_nodes[None, :], z_draped.shape)
        fill, _edges = draw_mesh(
            ax, x_2d, z_draped, data_draped, style=mesh_style,
            cmap="turbo", vmin=0.7, vmax=3.5,
        )
        fills.append(fill)
        surface_km = interp_elev(chain_km, elev_m / 1000.0, x_nodes)
        ax.plot(x_nodes, surface_km, color="#211813", linewidth=1.6, zorder=8)
        marker_y = elev_m / 1000.0 + 0.03
        ax.scatter(chain_km, marker_y, marker="v", s=32, color="black", zorder=10)
        for xi, yi, lab in zip(chain_km, marker_y + 0.10, labels):
            ax.text(
                xi, yi, lab, rotation=90, ha="center", va="bottom",
                fontsize=6.6, zorder=11,
            )
        ax.set_ylim(float(surface_km.min() - depths_km[-1]), float(surface_km.max() + 0.55))
        ax.set_xlim(float(chain_km.min()), float(chain_km.max()))
        ax.set_xlabel("Profile distance (km)")
        ax.set_title(f"{name}: Maxwell-trained 3-D AI prediction, mesh + topography")
    axes[0].set_ylabel("Elevation (km)")
    fig.colorbar(fills[0], ax=axes, label=r"$\log_{10}\rho$ (ohm m)", shrink=0.85)
    save(fig, "willy_ai3d_maxwell_predictions_both_lines.png")


def plot_validation(result) -> None:
    recovery = result.data.get("mt3d_recovery") or {}
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.0), constrained_layout=True)
    axes[0].bar(["AI-3D combined"], [float(result.data["rms_global"])], color="#2f6f8f")
    axes[0].axhline(1, color="0.25", ls="--", lw=1)
    axes[0].set(title="Observed-response diagnostic", ylabel="global RMS")
    axes[1].bar(
        ["RMSE", "MAE"],
        [float(recovery.get("rmse", np.nan)), float(recovery.get("mae", np.nan))],
        color=["#c85745", "#e0a458"],
    )
    axes[1].set(
        title=f"Held-out synthetic recovery (n={recovery.get('n_samples', 0)})",
        ylabel=r"error in $\log_{10}\rho$",
    )
    save(fig, "willy_ai3d_maxwell_validation_both_lines.png")


def production_config_reference():
    """Configuration-only reference for a real, many-hour 3-D run.

    Not executed here. Copy this function's body (or run this script with
    ``--production``) in your own environment; a realistic budget for this
    survey is likely several hours to tens of hours on a single CPU core,
    depending on ``n_train_profiles``, ``geology_grid_nx_ny/nz``, and
    ``max_mesh_cells``.
    """
    return Inv3DAgent(
        physics="mt3d",
        n_layers=10,
        freqs=np.geomspace(1.0, 10_000.0, 24),
        depth_max=2200.0,
        n_train_profiles=200,
        epochs=60,
        radius=450.0,
        hidden=(128, 64, 32),
        dropout=0.1,
        n_mc=0,
        correlation_length_x_m=(300.0, 1200.0),
        correlation_length_y_m=(300.0, 900.0),
        correlation_length_z_m=(80.0, 400.0),
        log_resistivity_mean=2.1,
        log_resistivity_std=0.75,
        geology_grid_nx_ny=8,
        geology_grid_nz=8,
        mesh_safety_factor=8.0,
        max_mesh_cells=300_000,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--production",
        action="store_true",
        help="Build the production agent (configuration only; does NOT execute).",
    )
    args = parser.parse_args(argv)

    lines = load_lines()
    problem = build_combined_problem(lines)
    plot_geology_and_mesh(problem)

    if args.production:
        agent = production_config_reference()
        print(
            "production config only, not executed:",
            agent.physics, agent.n_train_profiles, agent.epochs,
        )
        return 0

    result = run_fast_inversion(problem)
    plot_inversion_with_mesh(problem, result)
    plot_validation(result)

    quality = problem["solver_model"].quality
    recovery = result.data.get("mt3d_recovery") or {}
    print(
        "combined",
        "stations", len(problem["combined_sites"]),
        "geology", problem["grid"].shape,
        "mesh", problem["solver_model"].mesh.shape,
        "cells", quality.cell_count,
        "RMS", f"{result.data['rms_global']:.3f}",
        "recovery_RMSE", f"{recovery.get('rmse', float('nan')):.3f}",
        "recovery_n", recovery.get("n_samples", 0),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
