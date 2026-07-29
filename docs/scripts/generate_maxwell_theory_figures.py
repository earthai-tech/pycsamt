"""Generate executed figures for theory/maxwell_forward.rst."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.ai.geology import GeologyGrid, TopographicSurface  # noqa: E402
from pycsamt.forward.maxwell import (  # noqa: E402
    MeshDesign,
    build_solver_mesh,
    half_space_impedance,
    layered_earth_impedance,
    skin_depth_m,
)

IMAGE_DIR = ROOT / "docs/source/images/theory"


def _save(fig: plt.Figure, name: str) -> None:
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(IMAGE_DIR / name, dpi=190, bbox_inches="tight")
    plt.close(fig)


def make_mesh_anatomy() -> tuple[int, float, int]:
    grid = GeologyGrid.regular_2d(nx=28, nz=18, dx_m=125, dz_m=80)
    z, x = np.mgrid[:18, :28]
    rho = np.full(grid.shape, 600.0)
    rho[z < 4] = 55.0
    rho[((x - 17) / 5) ** 2 + ((z - 10) / 3) ** 2 <= 1] = 9.0
    elevation = 440 + 42 * np.sin(2 * np.pi * grid.x_m / np.ptp(grid.x_m))
    topo = TopographicSurface(grid, elevation, float(elevation.max()), source="synthetic")
    model = build_solver_mesh(
        grid,
        resistivity_ohm_m=rho,
        frequencies_hz=np.geomspace(.5, 1000, 14),
        topography=topo,
        design=MeshDesign(horizontal_padding_cells=5, bottom_padding_cells=5,
                          air_layers=5, padding_expansion=1.3, air_expansion=1.2),
    )
    display = np.log10(1 / model.conductivity_s_m)
    xe, ze = model.mesh.x_edges_m / 1000, model.mesh.z_edges_m / 1000
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.7), constrained_layout=True)
    im = axes[0].pcolormesh(xe, ze, display, cmap="turbo", shading="flat", vmin=.8, vmax=8)
    cs = model.core_slices
    axes[0].axvline(model.mesh.x_edges_m[cs[1].start] / 1000, color="white", ls="--")
    axes[0].axvline(model.mesh.x_edges_m[cs[1].stop] / 1000, color="white", ls="--")
    axes[0].axhline(model.mesh.z_edges_m[cs[0].start] / 1000, color="white", ls="--")
    axes[0].axhline(model.mesh.z_edges_m[cs[0].stop] / 1000, color="white", ls="--")
    axes[0].invert_yaxis(); axes[0].set(title="Core, air, and geometric padding", xlabel="x (km)", ylabel="z (km)")
    fig.colorbar(im, ax=axes[0], label=r"$\log_{10}\rho$ ($\Omega$ m)")
    frequencies = np.geomspace(.1, 10000, 200)
    for value in (10, 100, 1000):
        axes[1].loglog(frequencies, skin_depth_m(value, frequencies) / 1000,
                       label=fr"$\rho={value}\ \Omega$ m")
    axes[1].axhline(.08 * 4, color="black", ls="--", label="4 core cells")
    axes[1].set(title="Skin depth is a mesh scale, not a guarantee", xlabel="frequency (Hz)", ylabel="skin depth (km)")
    axes[1].legend(fontsize=8)
    _save(fig, "maxwell_mesh_anatomy.png")
    return model.quality.cell_count, model.quality.cells_per_minimum_skin_depth, len(model.quality.warnings)


def make_analytic_responses() -> tuple[float, float]:
    f = np.geomspace(.01, 10000, 160)
    half = half_space_impedance(100, f)
    layered = layered_earth_impedance([30, 600, 10], [250, 700], f)
    mu0 = 4e-7 * np.pi
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.2), constrained_layout=True)
    for values, label in ((half, "100 ohm m half-space"), (layered, "30/600/10 ohm m layers")):
        rhoa = np.abs(values) ** 2 / (mu0 * 2 * np.pi * f)
        phase = np.rad2deg(np.angle(values))
        axes[0].loglog(f, rhoa, label=label)
        axes[1].semilogx(f, phase, label=label)
    axes[0].set(title="Amplitude response", xlabel="frequency (Hz)", ylabel=r"apparent resistivity ($\Omega$ m)")
    axes[1].set(title="Phase response", xlabel="frequency (Hz)", ylabel="phase (degree)")
    for ax in axes: ax.legend(fontsize=8); ax.grid(alpha=.25)
    _save(fig, "maxwell_analytic_benchmarks.png")
    return float(np.min(np.abs(layered))), float(np.max(np.abs(layered)))


def make_adapter_lifecycle() -> None:
    labels = ["geology", "mesh model", "problem", "capability\npreflight",
              "backend solve", "postflight", "forward result"]
    colors = ["#2a9d8f", "#2a9d8f", "#457b9d", "#e9c46a", "#e76f51", "#e9c46a", "#457b9d"]
    fig, ax = plt.subplots(figsize=(14, 3.3), constrained_layout=True)
    ax.set_xlim(-.6, 6.6); ax.set_ylim(-1.0, 1.1); ax.axis("off")
    for i, (label, color) in enumerate(zip(labels, colors)):
        ax.text(i, 0, label, ha="center", va="center", fontsize=10,
                bbox=dict(boxstyle="round,pad=.55", fc=color, ec="white", lw=1.5), color="white" if color != "#e9c46a" else "black")
        if i < len(labels) - 1:
            ax.annotate("", (i + .68, 0), (i + .32, 0), arrowprops=dict(arrowstyle="->", lw=1.8))
    ax.annotate("", (3, -.62), (3, -.24),
                arrowprops=dict(arrowstyle="->", color="#b22222"))
    ax.text(3, -.8, "reject incompatible", ha="center", color="#b22222")
    ax.annotate("", (5, -.62), (5, -.24),
                arrowprops=dict(arrowstyle="->", color="#b22222"))
    ax.text(5, -.8, "reject invalid/unconverged", ha="center", color="#b22222")
    ax.set_title("The adapter is a validation boundary, not merely a function wrapper", fontsize=14)
    _save(fig, "maxwell_adapter_lifecycle.png")


def main() -> int:
    cells, cps, warnings = make_mesh_anatomy()
    zmin, zmax = make_analytic_responses()
    make_adapter_lifecycle()
    print("solver mesh cells:", cells)
    print("cells per minimum skin depth:", f"{cps:.2f}")
    print("mesh quality warnings:", warnings)
    print("layered impedance magnitude range:", f"{zmin:.6g}", f"{zmax:.6g}")
    print("figures generated: 3")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
