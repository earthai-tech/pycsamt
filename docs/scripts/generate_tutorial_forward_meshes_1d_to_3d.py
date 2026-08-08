"""Generate executed figures for tutorials/forward_model_1d_to_3d.rst."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.api import read_edis
from pycsamt.site import Sites
from pycsamt.topo.extract import extract_elevation
from pycsamt.emtools.source_effects import correct_near_field
from pycsamt.site.edit import select_freq_all
from pycsamt.forward import LayeredModel, MT1DForward
from pycsamt.forward.plot import plot_model_1d
from pycsamt.ai.geology import GeologyGrid
from pycsamt.forward.maxwell import MeshDesign, ReceiverSet, build_solver_mesh, skin_depth_m
from pycsamt.forward.maxwell.contracts_tri import TriProblem
from pycsamt.forward.maxwell.tri_fem2d import TriFEM2DAdapter
from pycsamt.forward.maxwell.tri_mesh_gen import build_graded_tri_mesh
from pycsamt.forward.maxwell.modem3d import ModEm3DAdapter

IMAGE_DIR = ROOT / "docs/source/images/tutorials/forward_model_1d_to_3d"
MU0 = 4.0e-7 * np.pi

RHO1, RHO2, Z1 = 400.0, 3000.0, 300.0
STATION_X = [50.0 * i for i in range(10)]
FULL_FREQS = np.array([8196.722, 4098.361, 2049.18, 1023.541, 512.8206, 255.7545])
SOLVER_FREQS = np.array([512.8206, 255.7545])


def _save(fig: plt.Figure, name: str) -> None:
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(IMAGE_DIR / name, dpi=190, bbox_inches="tight")
    plt.close(fig)


def _load_line():
    survey = read_edis(str(ROOT / "data/CSAMT"), recursive=False, strict=False, progress=False)
    sites = Sites(survey.collection).ordered("station")
    corrected = correct_near_field(sites, source_offset=1000.0, inplace=False)
    usable = list(select_freq_all(corrected, fmin=255.0))
    elevation = extract_elevation(sites)
    return sites, usable, elevation


def make_1d_sounding() -> dict:
    _, usable, _ = _load_line()
    s0 = usable[0]
    real_rho = s0.rho[:, 0, 1]

    model = LayeredModel(resistivity=[RHO1, RHO2], thickness=[Z1])
    response = MT1DForward(FULL_FREQS).run(model)

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6), constrained_layout=True)
    plot_model_1d(model, depth_max=1000.0, ax=axes[0], title="Earth model")

    axes[1].loglog(FULL_FREQS, real_rho, "o-", color="#e76f51", label="real (corrected), csa000")
    axes[1].loglog(FULL_FREQS, response.rho_a, "s--", color="#457b9d", label="1-D model, 400/3000 ohm m")
    axes[1].set(title="Apparent resistivity", xlabel="frequency (Hz)", ylabel=r"$\rho_a$ ($\Omega$ m)")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=.25, which="both")

    axes[2].semilogx(FULL_FREQS, np.full_like(FULL_FREQS, np.nan), alpha=0)  # keep x-axis aligned
    axes[2].semilogx(FULL_FREQS, response.phase, "s--", color="#457b9d", label="1-D model")
    axes[2].set(title="Phase", xlabel="frequency (Hz)", ylabel="phase (degree)")
    axes[2].grid(alpha=.25)
    axes[2].legend(fontsize=8)
    _save(fig, "sounding_1d.png")
    return {"real_rho": real_rho.tolist(), "model_rho": response.rho_a.tolist()}


def make_2d_triangular() -> dict:
    sites, usable, elevation = _load_line()
    topo_z = elevation.max() - elevation

    surface_cell_m = 0.04 * float(skin_depth_m(RHO1, SOLVER_FREQS.max()))
    half_width = 8.0 * float(skin_depth_m(RHO2, SOLVER_FREQS.min()))
    mesh = build_graded_tri_mesh(
        (-half_width, 450.0 + half_width), (0.0, half_width), STATION_X,
        surface_cell_m=surface_cell_m, growth_rate=1.3, max_cell_m=half_width / 6.0,
        topo_x_m=STATION_X, topo_z_m=topo_z,
    )
    centroids = mesh.triangle_centroids_m
    local_terrain = np.interp(centroids[:, 0], STATION_X, topo_z)
    depth_below_terrain = centroids[:, 1] - local_terrain
    resistivity = np.where(depth_below_terrain < Z1, RHO1, RHO2)
    conductivity = 1.0 / resistivity

    station_names = [f"csa{int(x):03d}" for x in STATION_X]
    receivers = ReceiverSet(np.column_stack([STATION_X, topo_z]), station_names)
    problem = TriProblem(mesh, conductivity, list(SOLVER_FREQS), receivers)
    result = TriFEM2DAdapter().solve(problem)
    zxy = list(problem.components).index("zxy")
    rhoa_2d = np.abs(result.impedance_v_a[:, :, zxy]) ** 2 / (MU0 * 2 * np.pi * SOLVER_FREQS[None, :])

    real_rho = np.array([s.rho[:, 0, 1][[4, 5]] for s in usable])  # 512.8, 255.8 Hz columns

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), constrained_layout=True)
    triangulation = plt.matplotlib.tri.Triangulation(mesh.nodes_m[:, 0], mesh.nodes_m[:, 1], mesh.triangles)
    im = axes[0].tripcolor(triangulation, facecolors=resistivity, cmap="turbo_r",
                            vmin=RHO1, vmax=RHO2, edgecolors="white", linewidth=0.1)
    axes[0].plot(STATION_X, topo_z, color="black", lw=1.5)
    axes[0].scatter(STATION_X, topo_z, marker="v", color="black", s=34, zorder=5)
    axes[0].set_xlim(-80, 530)
    axes[0].set_ylim(700, -80)
    axes[0].set(title="Real topography-draped triangular mesh (zoomed to the line)",
                xlabel="x (m)", ylabel="z (m)")
    fig.colorbar(im, ax=axes[0], label=r"resistivity ($\Omega$ m)")

    axes[1].plot(STATION_X, real_rho[:, 0], "o-", color="#e76f51", label="real, 512.8 Hz")
    axes[1].plot(STATION_X, real_rho[:, 1], "o--", color="#e76f51", alpha=.6, label="real, 255.8 Hz")
    axes[1].plot(STATION_X, rhoa_2d[:, 0], "s-", color="#457b9d", label="2-D model, 512.8 Hz")
    axes[1].plot(STATION_X, rhoa_2d[:, 1], "s--", color="#457b9d", alpha=.6, label="2-D model, 255.8 Hz")
    axes[1].set(title="Apparent resistivity at each real station", xlabel="station x (m)",
                ylabel=r"$\rho_a$ ($\Omega$ m)")
    axes[1].legend(fontsize=7)
    axes[1].grid(alpha=.25)
    _save(fig, "triangular_2d.png")
    return {
        "n_triangles": mesh.n_triangles,
        "rhoa_2d": rhoa_2d.tolist(),
        "real_rho": real_rho.tolist(),
        "station_names": station_names,
    }


def _build_3d_model():
    grid = GeologyGrid.regular_3d(nx=10, ny=3, nz=14, dx_m=50.0, dy_m=200.0, dz_m=40.0)
    rho_3d = np.where(np.arange(14)[:, None, None] * 40.0 < Z1, RHO1, RHO2)
    rho_3d = np.broadcast_to(rho_3d, grid.shape).copy()
    design = MeshDesign(horizontal_padding_cells=6, bottom_padding_cells=6, air_layers=0, padding_expansion=1.3)
    model_3d = build_solver_mesh(grid, resistivity_ohm_m=rho_3d, frequencies_hz=list(SOLVER_FREQS), design=design)
    receivers_3d = ReceiverSet([[x, 0.0, 0.0] for x in STATION_X], [f"csa{int(x):03d}" for x in STATION_X])
    problem_3d = model_3d.to_problem(list(SOLVER_FREQS), receivers_3d, mark_air_inactive=False)
    result_3d = ModEm3DAdapter().solve(problem_3d)
    zxy_3d = list(problem_3d.components).index("zxy")
    rhoa_3d = np.abs(result_3d.impedance_v_a[:, :, zxy_3d]) ** 2 / (MU0 * 2 * np.pi * SOLVER_FREQS[None, :])
    return model_3d, problem_3d, rhoa_3d


def make_3d_mesh() -> dict:
    model_3d, _, _ = _build_3d_model()
    resistivity = 1.0 / model_3d.conductivity_s_m
    cz, cy, cx = model_3d.core_slices
    y_index = (cy.start + cy.stop) // 2

    x_edges = model_3d.mesh.x_edges_m
    y_edges = model_3d.mesh.y_edges_m
    z_edges = model_3d.mesh.z_edges_m
    z_line = model_3d.mesh.cell_centres_m["z"][cz]

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), constrained_layout=True)

    im0 = axes[0].pcolormesh(x_edges, z_edges, resistivity[:, y_index, :], cmap="turbo_r",
                              vmin=RHO1, vmax=RHO2, shading="flat")
    axes[0].axvline(x_edges[cx.start], color="white", ls="--", lw=1.1)
    axes[0].axvline(x_edges[cx.stop], color="white", ls="--", lw=1.1)
    axes[0].axhline(z_edges[cz.start], color="white", ls="--", lw=1.1)
    axes[0].axhline(z_edges[cz.stop], color="white", ls="--", lw=1.1)
    axes[0].scatter(STATION_X, [0.0] * 10, marker="v", color="black", s=32, zorder=5)
    axes[0].set_xlim(x_edges[cx.start] - 400, x_edges[cx.stop] + 400)
    axes[0].set_ylim(z_edges[cz.stop] + 400, -100)
    axes[0].set(title="Vertical slice through the real profile (y=0)", xlabel="x (m)", ylabel="z (m)")
    fig.colorbar(im0, ax=axes[0], label=r"resistivity ($\Omega$ m)")

    z_index = cz.start + 2
    im1 = axes[1].pcolormesh(x_edges, y_edges, resistivity[z_index], cmap="turbo_r",
                              vmin=RHO1, vmax=RHO2, shading="flat")
    axes[1].axvline(x_edges[cx.start], color="white", ls="--", lw=1.1)
    axes[1].axvline(x_edges[cx.stop], color="white", ls="--", lw=1.1)
    axes[1].axhline(y_edges[cy.start], color="white", ls="--", lw=1.1)
    axes[1].axhline(y_edges[cy.stop], color="white", ls="--", lw=1.1)
    axes[1].scatter(STATION_X, [0.0] * 10, marker="v", color="black", s=32, zorder=5)
    axes[1].set(title=f"Horizontal slice at z={float(z_line[2]):.0f} m", xlabel="x (m)", ylabel="y (m)")
    fig.colorbar(im1, ax=axes[1], label=r"resistivity ($\Omega$ m)")

    _save(fig, "mesh_3d.png")
    return {
        "mesh_shape": model_3d.mesh.shape,
        "quality_acceptable": model_3d.quality.acceptable,
        "y_extent_m": (float(y_edges[0]), float(y_edges[-1])),
    }


def make_dimension_comparison() -> dict:
    _, usable, elevation = _load_line()
    real_rho = np.array([s.rho[:, 0, 1][[4, 5]] for s in usable])

    model = LayeredModel(resistivity=[RHO1, RHO2], thickness=[Z1])
    response_1d = MT1DForward(SOLVER_FREQS).run(model)

    model_3d, _, rhoa_3d = _build_3d_model()

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6), constrained_layout=True)
    for fi, (freq, ax) in enumerate(zip(SOLVER_FREQS, axes)):
        ax.plot(STATION_X, real_rho[:, fi], "o-", color="#e76f51", label="real field data")
        ax.axhline(response_1d.rho_a[fi], color="#2a9d8f", ls=":", lw=2, label="1-D analytic (no lateral structure)")
        ax.axhline(rhoa_3d[0, fi], color="#6a4c93", ls="-.", lw=2, label="3-D volume (y-invariant)")
        ax.plot(STATION_X, [np.nan] * 10, color="#457b9d")  # legend placeholder
        ax.set(title=f"{freq:.1f} Hz", xlabel="station x (m)", ylabel=r"$\rho_a$ ($\Omega$ m)")
        ax.grid(alpha=.25)
    axes[0].legend(fontsize=7)
    _save(fig, "dimension_comparison.png")
    return {
        "rhoa_1d": response_1d.rho_a.tolist(),
        "rhoa_3d_station0": rhoa_3d[0].tolist(),
        "rhoa_3d_station5": rhoa_3d[5].tolist(),
        "mesh_3d_shape": model_3d.mesh.shape,
    }


def main() -> int:
    one = make_1d_sounding()
    two = make_2d_triangular()
    mesh_3d = make_3d_mesh()
    three = make_dimension_comparison()
    print("1-D real csa000 rho:", [round(v, 1) for v in one["real_rho"]])
    print("1-D model rho:", [round(v, 1) for v in one["model_rho"]])
    print("2-D n_triangles:", two["n_triangles"])
    print("3-D mesh shape:", mesh_3d["mesh_shape"])
    print("3-D mesh quality acceptable:", mesh_3d["quality_acceptable"])
    print("3-D mesh y extent (m):", tuple(round(v, 1) for v in mesh_3d["y_extent_m"]))
    print("3-D rhoa (station0):", [round(v, 1) for v in three["rhoa_3d_station0"]])
    print("3-D rhoa (station5):", [round(v, 1) for v in three["rhoa_3d_station5"]])
    print("figures generated: 4")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
