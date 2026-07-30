"""Run and plot the two-line Maxwell-AI workflow used by the porphyry tutorial."""

from __future__ import annotations

import argparse
import json
import random
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.agents import Inv2DAgent  # noqa: E402
from pycsamt.ai.geology import (  # noqa: E402
    ElectricalLayer,
    EllipsoidalLens,
    GaussianCorrelation,
    GeologyGrid,
    generate_layered_geology,
    insert_lenses,
    topography_from_sites,
)
from pycsamt.emtools import ensure_sites  # noqa: E402
from pycsamt.forward.maxwell import MeshDesign, build_solver_mesh  # noqa: E402
from pycsamt.topo import extract_chainage, extract_elevation  # noqa: E402

IMAGE_DIR = ROOT / "docs/source/images/tutorials/map_porphyry_mineralization_from_noisy_amt"
RUN_DIR = ROOT / "runs"
FREQUENCIES_HZ = np.geomspace(1.0, 1000.0, 8)
PRODUCTION_FREQUENCIES_HZ = np.geomspace(1.0, 10_000.0, 24)


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


def build_line_problem(name, sites, *, seed, lens_x_m):
    chain_m = extract_chainage(sites) * 1000.0
    dx = 100.0
    x_origin = -250.0
    nx = int(np.ceil((chain_m[-1] + 500.0) / dx))
    grid = GeologyGrid.regular_2d(
        nx=nx, nz=24, dx_m=dx, dz_m=75.0, x_origin_m=x_origin
    )
    correlation = GaussianCorrelation(700.0, 160.0)
    layers = [
        ElectricalLayer("weathered cover", 45.0, 0.10, correlation),
        ElectricalLayer("altered intrusive", 350.0, 0.16, correlation),
        ElectricalLayer("fresh intrusive", 1600.0, 0.10, correlation),
    ]
    base = generate_layered_geology(
        grid,
        layers,
        [320.0, 1050.0],
        seed=seed,
        interface_relief_std_m=[45.0, 100.0],
        interface_correlation=correlation,
        minimum_thickness_m=120.0,
    )
    target = EllipsoidalLens(
        f"{name} conductive target",
        center_x_m=lens_x_m,
        center_z_m=760.0,
        radius_x_m=480.0,
        radius_z_m=210.0,
        resistivity_ohm_m=12.0,
        dip_deg=18.0 if name == "L26PLT" else -12.0,
        transition_fraction=0.20,
    )
    geology = insert_lenses(base, [target])
    topography = topography_from_sites(
        sites, grid, profile_origin_m=x_origin, interpolation_method="linear"
    )
    solver_model = build_solver_mesh(
        grid,
        resistivity_ohm_m=geology.resistivity_ohm_m,
        frequencies_hz=FREQUENCIES_HZ,
        topography=topography,
        design=MeshDesign(
            horizontal_padding_cells=4,
            bottom_padding_cells=5,
            air_layers=4,
            padding_expansion=1.3,
        ),
    )
    return {
        "name": name,
        "sites": sites,
        "chain_m": chain_m,
        "grid": grid,
        "geology": geology,
        "topography": topography,
        "solver_model": solver_model,
    }


def plot_priors_and_meshes(problems) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13.6, 8.2), constrained_layout=True)
    for row, problem in enumerate(problems):
        grid = problem["grid"]
        topo = problem["topography"]
        geology = problem["geology"]
        mesh_model = problem["solver_model"]
        im = axes[row, 0].pcolormesh(
            grid.x_m / 1000.0,
            grid.z_m / 1000.0,
            np.log10(geology.resistivity_ohm_m),
            cmap="turbo", vmin=0.7, vmax=3.5, shading="auto",
        )
        axes[row, 0].plot(
            grid.x_m / 1000.0, topo.surface_depth_m / 1000.0, "k", lw=1.4
        )
        axes[row, 0].invert_yaxis()
        axes[row, 0].set(
            title=f"{problem['name']}: seeded geology prior",
            xlabel="profile distance (km)", ylabel="depth below datum (km)",
        )
        x = mesh_model.mesh.cell_centres_m["x"] / 1000.0
        z = mesh_model.mesh.cell_centres_m["z"] / 1000.0
        axes[row, 1].pcolormesh(
            x, z, np.log10(1.0 / mesh_model.conductivity_s_m),
            cmap="turbo", vmin=0.7, vmax=8.0, shading="auto",
        )
        axes[row, 1].invert_yaxis()
        axes[row, 1].set(
            title=f"{problem['name']}: padded Maxwell mesh",
            xlabel="profile distance (km)", ylabel="depth below datum (km)",
        )
    fig.colorbar(im, ax=axes[:, 0], shrink=.82, label=r"$\log_{10}\rho$ (ohm m)")
    save(fig, "willy_ai2d_geology_maxwell_both_lines.png")


def run_inversions(
    problems,
    *,
    production: bool = False,
    n_train_profiles: int | None = None,
    epochs: int | None = None,
):
    np.random.seed(17)
    random.seed(17)
    try:
        import torch

        torch.manual_seed(17)
    except ImportError:
        pass
    results = []
    for problem in problems:
        chain = problem["chain_m"]
        spacing = float(np.median(np.diff(chain)))
        frequencies = (
            PRODUCTION_FREQUENCIES_HZ if production else FREQUENCIES_HZ
        )
        profiles = int(
            n_train_profiles
            if n_train_profiles is not None
            else (256 if production else 10)
        )
        epoch_budget = int(
            epochs if epochs is not None else (80 if production else 30)
        )
        line_output = RUN_DIR / (
            f"{problem['name']}_ai2d_maxwell_production"
            if production
            else f"{problem['name']}_ai2d_maxwell"
        )
        line_output.mkdir(parents=True, exist_ok=True)
        started = time.time()
        agent = Inv2DAgent(
            physics="mt2d",
            n_depth=32 if production else 24,
            n_freqs=len(frequencies),
            freqs=frequencies,
            depth_max=2200.0 if production else 1800.0,
            n_train_profiles=profiles,
            n_stations_per_profile=len(problem["sites"]),
            station_spacing_m=spacing,
            epochs=epoch_budget,
            correlation_length_x_m=(250.0, 1600.0) if production else (350.0, 1000.0),
            correlation_length_z_m=(60.0, 450.0) if production else (90.0, 300.0),
            log_resistivity_mean=2.1,
            log_resistivity_std=0.75 if production else 0.5,
            lambda_x=0.01 if production else 0.02,
            lambda_z=0.005 if production else 0.01,
            lambda_tv=0.002 if production else 0.005,
            mesh_safety_factor=8.0 if production else 6.0,
            max_mesh_cells=300_000 if production else 200_000,
        )
        result = agent.execute(
            {
                "sites": problem["sites"],
                "output_dir": str(line_output),
                "topography": True,
            }
        )
        if result.status != "success":
            raise RuntimeError(f"{problem['name']}: {result.summary}")
        checkpoint = line_output / "em_inverter_2d_maxwell.npz"
        result.data["inverter"].save(checkpoint)
        np.savez_compressed(
            line_output / "field_prediction.npz",
            pred_section=np.asarray(result.data["pred_section"]),
            depths_km=np.asarray(result.data["depths_km"]),
            frequency_grid_hz=np.asarray(result.data["frequency_grid_hz"]),
            station_names=np.asarray(result.data["station_names"]),
            chainage_km=np.asarray(problem["chain_m"]) / 1000.0,
            elevation_m=extract_elevation(problem["sites"]),
        )
        recovery = result.data.get("mt2d_recovery") or {}
        manifest = {
            "line": problem["name"],
            "status": result.status,
            "summary": result.summary,
            "production": production,
            "seed": 17,
            "n_train_profiles": profiles,
            "epochs_requested": epoch_budget,
            "n_depth": 32 if production else 24,
            "depth_max_m": 2200.0 if production else 1800.0,
            "frequencies_hz": frequencies.tolist(),
            "rms_global": float(result.data["rms_global"]),
            "mt2d_recovery": recovery,
            "elapsed_seconds": time.time() - started,
            "checkpoint": checkpoint.name,
        }
        (line_output / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2, default=float), encoding="utf-8"
        )
        print(json.dumps(manifest, default=float), flush=True)
        results.append(result)
    return results


def plot_inversions(problems, results) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(13.8, 8.4), constrained_layout=True)
    images = []
    for ax, problem, result in zip(axes, problems, results):
        section = np.asarray(result.data["pred_section"], dtype=float)
        depth_km = np.asarray(result.data["depths_km"], dtype=float)
        chain_km = problem["chain_m"] / 1000.0
        elevation_km = extract_elevation(problem["sites"]) / 1000.0
        reference = float(np.max(elevation_km))
        z = reference - depth_km[:, None] + (elevation_km - reference)[None, :]
        x = np.broadcast_to(chain_km[None, :], section.shape)
        image = ax.contourf(
            x, z, section, levels=np.linspace(.7, 3.5, 29),
            cmap="turbo", vmin=.7, vmax=3.5, extend="both",
        )
        images.append(image)
        ax.plot(chain_km, elevation_km, "k", lw=1.5)
        ax.scatter(chain_km, elevation_km + .025, marker="v", c="k", s=24)
        ax.set(
            title=f"{problem['name']}: direct Maxwell-trained AI prediction",
            xlabel="profile distance (km)", ylabel="elevation (km)",
        )
    fig.colorbar(images[0], ax=axes, shrink=.86, label=r"predicted $\log_{10}\rho$ (ohm m)")
    save(fig, "willy_ai2d_maxwell_predictions_both_lines.png")


def plot_validation(problems, results) -> None:
    names = [p["name"] for p in problems]
    rms = [float(r.data["rms_global"]) for r in results]
    recovery = [r.data.get("mt2d_recovery") or {} for r in results]
    rmse = [float(r.get("rmse", np.nan)) for r in recovery]
    mae = [float(r.get("mae", np.nan)) for r in recovery]
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.2), constrained_layout=True)
    axes[0].bar(names, rms, color=["#2f6f8f", "#c85745"])
    axes[0].axhline(1, color="0.25", ls="--", lw=1)
    axes[0].set(title="Observed-response diagnostic", ylabel="global RMS")
    x = np.arange(2)
    axes[1].bar(x - .18, rmse, .36, label="RMSE")
    axes[1].bar(x + .18, mae, .36, label="MAE")
    axes[1].set_xticks(x, names)
    axes[1].set(title="Held-out synthetic recovery", ylabel=r"error in $\log_{10}\rho$")
    axes[1].legend()
    save(fig, "willy_ai2d_maxwell_validation_both_lines.png")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--production", action="store_true")
    parser.add_argument("--n-train-profiles", type=int)
    parser.add_argument("--epochs", type=int)
    args = parser.parse_args(argv)
    lines = load_lines()
    problems = [
        build_line_problem("L26PLT", lines["L26PLT"], seed=2601, lens_x_m=1550.0),
        build_line_problem("L30PLT", lines["L30PLT"], seed=3001, lens_x_m=1450.0),
    ]
    plot_priors_and_meshes(problems)
    results = run_inversions(
        problems,
        production=args.production,
        n_train_profiles=args.n_train_profiles,
        epochs=args.epochs,
    )
    plot_inversions(problems, results)
    plot_validation(problems, results)
    for problem, result in zip(problems, results):
        quality = problem["solver_model"].quality
        recovery = result.data.get("mt2d_recovery") or {}
        print(
            problem["name"],
            "stations", len(problem["sites"]),
            "geology", problem["grid"].shape,
            "mesh", problem["solver_model"].mesh.shape,
            "cells", quality.cell_count,
            "RMS", f"{result.data['rms_global']:.3f}",
            "recovery_RMSE", f"{recovery.get('rmse', float('nan')):.3f}",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
