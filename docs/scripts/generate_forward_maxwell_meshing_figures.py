"""Generate executed figures for user_guide/forward/maxwell_meshing.rst."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.ai.geology import GeologyGrid, TopographicSurface
from pycsamt.forward.maxwell import MeshDesign, ReceiverSet, build_solver_mesh
from pycsamt.forward.maxwell.mt2d import MT2DAdapter
from pycsamt.forward.maxwell.benchmarks import layered_earth_impedance
from pycsamt.forward.maxwell.tri_mesh_gen import build_graded_tri_mesh

IMAGE_DIR = ROOT / "docs/source/images/user_guide/forward"


def _save(fig: plt.Figure, name: str) -> None:
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(IMAGE_DIR / name, dpi=190, bbox_inches="tight")
    plt.close(fig)


def make_mesh_anatomy() -> tuple[tuple[int, int], bool]:
    grid = GeologyGrid.regular_2d(nx=30, nz=30, dx_m=150, dz_m=40)
    rho = np.full(grid.shape, 400.0)
    rho[:6] = 50.0
    rho[10:18, 10:20] = 15.0

    elevation = 350 + 40 * np.sin(2 * np.pi * grid.x_m / np.ptp(grid.x_m))
    surface = TopographicSurface(grid, elevation, float(elevation.max()), source="synthetic")

    design = MeshDesign(
        horizontal_padding_cells=6, bottom_padding_cells=6, air_layers=6,
        padding_expansion=1.3, air_expansion=1.25,
    )
    model = build_solver_mesh(
        grid, resistivity_ohm_m=rho, frequencies_hz=[200.0, 20.0, 2.0],
        topography=surface, design=design,
    )

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), constrained_layout=True)
    x_edges, z_edges = model.mesh.x_edges_m / 1000, model.mesh.z_edges_m / 1000
    display = np.log10(1.0 / model.conductivity_s_m)
    im = axes[0].pcolormesh(x_edges, z_edges, display, cmap="turbo", shading="flat", vmin=1.0, vmax=4.5)
    cz, cx = model.core_slices
    axes[0].axvline(model.mesh.x_edges_m[cx.start] / 1000, color="white", ls="--", lw=1.2)
    axes[0].axvline(model.mesh.x_edges_m[cx.stop] / 1000, color="white", ls="--", lw=1.2)
    axes[0].axhline(model.mesh.z_edges_m[cz.start] / 1000, color="white", ls="--", lw=1.2)
    axes[0].axhline(model.mesh.z_edges_m[cz.stop] / 1000, color="white", ls="--", lw=1.2)
    axes[0].invert_yaxis()
    axes[0].set(title="Padded mesh: core, air, and padding", xlabel="x (km)", ylabel="z (km)")
    fig.colorbar(im, ax=axes[0], label=r"$\log_{10}\rho$ ($\Omega$ m)")

    x_widths = model.mesh.cell_widths_m["x"]
    axes[1].plot(np.arange(len(x_widths)), x_widths, marker="o", ms=3)
    axes[1].axvspan(cx.start, cx.stop, color="#2a9d8f", alpha=0.15, label="core")
    axes[1].set(title="x cell width by index (padding growth)", xlabel="cell index", ylabel="width (m)")
    axes[1].legend(fontsize=8)
    _save(fig, "maxwell_meshing_anatomy.png")
    return model.mesh.shape, model.quality.acceptable


def make_resolution_convergence() -> tuple[list[float], list[float]]:
    adapter = MT2DAdapter(verbose=False)
    frequencies = [10.0, 1.0]
    receivers = ReceiverSet([[4_000.0, 0.0]], ["S00"])
    analytic = layered_earth_impedance([100.0, 400.0], [500.0], frequencies)

    dz_values = [100.0, 50.0, 25.0, 12.5, 6.25]
    errors: list[float] = []
    cells_per_skin: list[float] = []
    for dz in dz_values:
        nz = int(round(4_000 / dz))
        grid = GeologyGrid.regular_2d(nx=40, nz=nz, dx_m=100.0, dz_m=dz)
        rho = np.full(grid.shape, 400.0)
        rho[: int(round(500 / dz))] = 100.0
        design = MeshDesign(
            horizontal_padding_cells=8, bottom_padding_cells=8, air_layers=6,
            padding_expansion=1.3, air_expansion=1.25,
        )
        model = build_solver_mesh(grid, resistivity_ohm_m=rho, frequencies_hz=frequencies, design=design)
        problem = model.to_problem(frequencies, receivers, mark_air_inactive=False)
        result = adapter.solve(problem)
        z = result.impedance_v_a[0, :, 0]
        errors.append(float(np.linalg.norm(z - analytic) / np.linalg.norm(analytic)))
        cells_per_skin.append(model.quality.cells_per_minimum_skin_depth)

    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    ax.loglog(dz_values, errors, marker="o", label="relative error vs. analytic")
    reference = errors[0] * (np.array(dz_values) / dz_values[0])
    ax.loglog(dz_values, reference, ls="--", color="black", label=r"$O(\Delta z)$ reference")
    ax.invert_xaxis()
    ax.set(title="Near-surface resolution convergence (lateral resolution fixed)",
           xlabel=r"core $\Delta z$ (m)", ylabel="relative impedance error")
    ax.legend(fontsize=9)
    ax.grid(alpha=.25, which="both")
    _save(fig, "maxwell_meshing_resolution_convergence.png")
    return errors, cells_per_skin


def make_3d_mesh_slices() -> tuple[tuple[int, int, int], bool]:
    grid = GeologyGrid.regular_3d(nx=22, ny=18, nz=16, dx_m=150.0, dy_m=150.0, dz_m=80.0)
    rho = np.full(grid.shape, 300.0)
    rho[:3] = 60.0
    rho[6:11, 6:12, 8:16] = 15.0
    design = MeshDesign(
        horizontal_padding_cells=5, bottom_padding_cells=5, air_layers=5,
        padding_expansion=1.3, air_expansion=1.25,
    )
    model = build_solver_mesh(grid, resistivity_ohm_m=rho, frequencies_hz=[100.0, 10.0, 1.0], design=design)
    resistivity = np.log10(1.0 / model.conductivity_s_m)
    cz, cy, cx = model.core_slices
    z_index = cz.start + 8
    y_index = cy.start + 8

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    x_edges, y_edges, z_edges = (model.mesh.x_edges_m / 1000, model.mesh.y_edges_m / 1000,
                                  model.mesh.z_edges_m / 1000)
    im0 = axes[0].pcolormesh(x_edges, y_edges, resistivity[z_index], cmap="turbo", shading="flat", vmin=1, vmax=3.2)
    axes[0].axvline(model.mesh.x_edges_m[cx.start] / 1000, color="white", ls="--", lw=1.1)
    axes[0].axvline(model.mesh.x_edges_m[cx.stop] / 1000, color="white", ls="--", lw=1.1)
    axes[0].axhline(model.mesh.y_edges_m[cy.start] / 1000, color="white", ls="--", lw=1.1)
    axes[0].axhline(model.mesh.y_edges_m[cy.stop] / 1000, color="white", ls="--", lw=1.1)
    axes[0].set(title=f"Horizontal slice at z={model.mesh.cell_centres_m['z'][z_index]:.0f} m",
                xlabel="x (km)", ylabel="y (km)")
    fig.colorbar(im0, ax=axes[0], label=r"$\log_{10}\rho$ ($\Omega$ m)")

    im1 = axes[1].pcolormesh(x_edges, z_edges, resistivity[:, y_index, :], cmap="turbo", shading="flat", vmin=1, vmax=3.2)
    axes[1].axvline(model.mesh.x_edges_m[cx.start] / 1000, color="white", ls="--", lw=1.1)
    axes[1].axvline(model.mesh.x_edges_m[cx.stop] / 1000, color="white", ls="--", lw=1.1)
    axes[1].axhline(model.mesh.z_edges_m[cz.start] / 1000, color="white", ls="--", lw=1.1)
    axes[1].axhline(model.mesh.z_edges_m[cz.stop] / 1000, color="white", ls="--", lw=1.1)
    axes[1].invert_yaxis()
    axes[1].set(title=f"Vertical slice at y={model.mesh.cell_centres_m['y'][y_index]:.0f} m",
                xlabel="x (km)", ylabel="z (km)")
    fig.colorbar(im1, ax=axes[1], label=r"$\log_{10}\rho$ ($\Omega$ m)")
    _save(fig, "maxwell_meshing_3d_slices.png")
    return model.mesh.shape, model.quality.acceptable


def make_tri_grading_law() -> tuple[int, float]:
    station_x = np.linspace(-500, 500, 9)
    surface_cell_m, growth_rate, max_cell_m = 20.0, 1.25, 250.0
    mesh = build_graded_tri_mesh(
        (-1_200, 1_200), (0, 900), station_x,
        surface_cell_m=surface_cell_m, growth_rate=growth_rate, max_cell_m=max_cell_m,
    )
    centroids = mesh.triangle_centroids_m
    station_points = np.column_stack([station_x, np.zeros_like(station_x)])
    distance = np.min(
        np.linalg.norm(centroids[:, None, :] - station_points[None, :, :], axis=2), axis=1,
    )
    edge_equivalent = np.sqrt(mesh.triangle_areas_m2 * 4.0 / np.sqrt(3.0))

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), constrained_layout=True)
    triangulation = plt.matplotlib.tri.Triangulation(mesh.nodes_m[:, 0], mesh.nodes_m[:, 1], mesh.triangles)
    axes[0].triplot(triangulation, lw=0.3, color="#457b9d")
    axes[0].scatter(station_x, np.zeros_like(station_x), marker="v", color="#e76f51", zorder=5, label="stations")
    axes[0].set_ylim(500, -100)
    axes[0].set(title=f"Graded mesh, {mesh.n_triangles} triangles", xlabel="x (m)", ylabel="z (m)")
    axes[0].legend(fontsize=8)

    order = np.argsort(distance)
    axes[1].scatter(distance, edge_equivalent, s=8, alpha=.5, label="triangle equivalent edge length")
    curve_distance = np.linspace(0, distance.max(), 200)
    theoretical = np.minimum(surface_cell_m * growth_rate ** (curve_distance / surface_cell_m), max_cell_m)
    axes[1].plot(curve_distance, theoretical, color="black", ls="--", label="size function")
    axes[1].set(title="Grading law: cell size vs. distance to nearest station",
                xlabel="distance to nearest station (m)", ylabel="edge length (m)")
    axes[1].legend(fontsize=8)
    _save(fig, "maxwell_meshing_tri_grading.png")
    correlation = float(np.corrcoef(theoretical_at(distance, surface_cell_m, growth_rate, max_cell_m), edge_equivalent)[0, 1])
    return mesh.n_triangles, correlation


def theoretical_at(distance, surface_cell_m, growth_rate, max_cell_m):
    return np.minimum(surface_cell_m * growth_rate ** (distance / surface_cell_m), max_cell_m)


def main() -> int:
    anatomy_shape, anatomy_acceptable = make_mesh_anatomy()
    errors, cells_per_skin = make_resolution_convergence()
    slices_shape, slices_acceptable = make_3d_mesh_slices()
    n_triangles, correlation = make_tri_grading_law()
    print("anatomy mesh shape:", anatomy_shape, "acceptable:", anatomy_acceptable)
    print("convergence errors:", [round(e, 4) for e in errors])
    print("convergence cells_per_skin (constant):", round(cells_per_skin[0], 2))
    print("3d mesh shape:", slices_shape, "acceptable:", slices_acceptable)
    print("tri mesh triangles:", n_triangles, "size-function correlation:", round(correlation, 3))
    print("figures generated: 4")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
