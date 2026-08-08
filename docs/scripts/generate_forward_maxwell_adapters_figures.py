"""Generate executed figures for user_guide/forward/maxwell_adapters.rst."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.forward.maxwell import MaxwellMesh, MaxwellProblem, ReceiverSet
from pycsamt.forward.maxwell.contracts_tri import TriProblem
from pycsamt.forward.maxwell.mt2d import MT2DAdapter
from pycsamt.forward.maxwell.tri_fem2d import TriFEM2DAdapter
from pycsamt.forward.maxwell.tri_mesh_gen import build_graded_tri_mesh

IMAGE_DIR = ROOT / "docs/source/images/user_guide/forward"
_MU0 = 4.0e-7 * np.pi


def _save(fig: plt.Figure, name: str) -> None:
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(IMAGE_DIR / name, dpi=190, bbox_inches="tight")
    plt.close(fig)


def _rhoa_phase(impedance: np.ndarray, frequencies_hz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    rhoa = np.abs(impedance) ** 2 / (_MU0 * 2 * np.pi * frequencies_hz)
    phase = np.rad2deg(np.angle(impedance))
    return rhoa, phase


def make_mt2d_response() -> tuple[tuple[int, int, int], bool]:
    mesh = MaxwellMesh(np.linspace(0, 10_000, 41), np.linspace(0, 5_000, 31))
    resistivity = np.full(mesh.shape, 100.0)
    resistivity[10:16, 15:26] = 5.0
    frequencies = np.geomspace(1000, 1, 12)
    station_x = np.linspace(2000, 8000, 7)
    receivers = ReceiverSet(
        [[float(x), 0.0] for x in station_x],
        [f"S{i:02d}" for i in range(len(station_x))],
    )
    problem = MaxwellProblem(
        mesh, 1.0 / resistivity, frequencies, receivers, ("zxy", "zyx")
    )
    result = MT2DAdapter(verbose=False).solve(problem)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3), constrained_layout=True)
    for si, x in enumerate(station_x):
        rhoa, phase = _rhoa_phase(result.impedance_v_a[si, :, 0], result.frequencies_hz)
        colour = plt.cm.viridis(si / (len(station_x) - 1))
        axes[0].loglog(result.frequencies_hz, rhoa, color=colour, label=f"x={x:.0f} m")
        axes[1].semilogx(result.frequencies_hz, phase, color=colour)
    axes[0].axhline(100.0, color="black", ls="--", lw=1, label="background 100 ohm m")
    axes[0].set(title=r"$Z_{xy}$ apparent resistivity", xlabel="frequency (Hz)",
                ylabel=r"apparent resistivity ($\Omega$ m)")
    axes[1].set(title=r"$Z_{xy}$ phase", xlabel="frequency (Hz)", ylabel="phase (degree)")
    axes[0].legend(fontsize=7, ncol=2)
    for ax in axes:
        ax.grid(alpha=.25, which="both")
    _save(fig, "maxwell_adapters_mt2d_response.png")
    return result.shape, bool(np.all(result.diagnostics.converged))


def make_trifem2d_mesh_response() -> tuple[int, tuple[int, int, int], bool]:
    station_x = np.linspace(-800, 800, 9)
    topo_x = np.linspace(-1200, 1200, 9)
    topo_z = -60 * np.exp(-((topo_x / 500) ** 2))
    mesh = build_graded_tri_mesh(
        (-1200, 1200), (0, 2500), station_x,
        surface_cell_m=30, growth_rate=1.35,
        topo_x_m=topo_x, topo_z_m=topo_z,
    )
    centroid_z = mesh.nodes_m[mesh.triangles, 1].mean(axis=1)
    resistivity_ohm_m = np.where(centroid_z < 150.0, 20.0, 300.0)
    receiver_z = np.interp(station_x, topo_x, topo_z)
    receivers = ReceiverSet(
        np.column_stack([station_x, receiver_z]),
        [f"S{i:02d}" for i in range(len(station_x))],
    )
    problem = TriProblem(mesh, 1.0 / resistivity_ohm_m, [200.0, 20.0], receivers)
    result = TriFEM2DAdapter().solve(problem)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
    triangulation = plt.matplotlib.tri.Triangulation(
        mesh.nodes_m[:, 0], mesh.nodes_m[:, 1], mesh.triangles
    )
    im = axes[0].tripcolor(triangulation, facecolors=resistivity_ohm_m,
                            cmap="turbo_r", vmin=20, vmax=300, edgecolors="white",
                            linewidth=0.15)
    axes[0].plot(topo_x, topo_z, color="black", lw=1.5)
    axes[0].scatter(station_x, receiver_z, marker="v", color="black", s=36, zorder=5)
    axes[0].set_ylim(600, -100)
    axes[0].set(title="Graded triangular mesh on real topography",
                xlabel="x (m)", ylabel="z (m)")
    fig.colorbar(im, ax=axes[0], label=r"resistivity ($\Omega$ m)")

    for ci, label in enumerate(("zxy", "zyx")):
        rhoa, _ = _rhoa_phase(result.impedance_v_a[:, 0, ci], result.frequencies_hz[0])
        axes[1].plot(station_x, rhoa, marker="o", label=f"{label}, 200 Hz")
    for ci, label in enumerate(("zxy", "zyx")):
        rhoa, _ = _rhoa_phase(result.impedance_v_a[:, 1, ci], result.frequencies_hz[1])
        axes[1].plot(station_x, rhoa, marker="s", ls="--", label=f"{label}, 20 Hz")
    axes[1].set(title="Apparent resistivity at each receiver",
                xlabel="station x (m)", ylabel=r"apparent resistivity ($\Omega$ m)")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=.25)
    _save(fig, "maxwell_adapters_trifem2d_mesh_response.png")
    return mesh.n_triangles, result.shape, bool(np.all(result.diagnostics.converged))


def main() -> int:
    mt2d_shape, mt2d_converged = make_mt2d_response()
    n_triangles, tri_shape, tri_converged = make_trifem2d_mesh_response()
    print("mt2d result shape:", mt2d_shape, "converged:", mt2d_converged)
    print("trifem2d triangles:", n_triangles, "result shape:", tri_shape,
          "converged:", tri_converged)
    print("figures generated: 2")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
