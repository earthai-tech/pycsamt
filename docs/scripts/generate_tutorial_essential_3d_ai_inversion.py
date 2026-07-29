"""Generate reproducible geology and Maxwell-mesh figures for the 3-D tutorial."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.ai.geology import (  # noqa: E402
    ElectricalLayer,
    EllipsoidalLens,
    GaussianCorrelation,
    GeologyGrid,
    TopographicSurface,
    generate_layered_geology,
    insert_lenses,
)
from pycsamt.forward.maxwell import MeshDesign, MT3DAdapter, build_solver_mesh  # noqa: E402

IMAGE_DIR = ROOT / "docs/source/images/tutorials/essential_3d_ai_inversion"


def save(fig: plt.Figure, name: str) -> None:
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(IMAGE_DIR / name, dpi=190, bbox_inches="tight")
    plt.close(fig)


def build_model():
    grid = GeologyGrid.regular_3d(
        nx=36, ny=24, nz=24, dx_m=200, dy_m=200, dz_m=100,
        x_origin_m=-3600, y_origin_m=-2400,
    )
    corr = GaussianCorrelation(1200, 180, length_y_m=800, azimuth_deg=25)
    layers = [
        ElectricalLayer("weathered cover", 35, log10_std=0.10, heterogeneity=corr),
        ElectricalLayer("resistive host", 900, log10_std=0.12, heterogeneity=corr),
        ElectricalLayer("deep basement", 2200, log10_std=0.08, heterogeneity=corr),
    ]
    base = generate_layered_geology(
        grid, layers, [450, 1350], seed=41,
        interface_relief_std_m=[70, 130], interface_correlation=corr,
        minimum_thickness_m=150,
    )
    lens = EllipsoidalLens(
        "dipping conductor", center_x_m=400, center_y_m=-100,
        center_z_m=900, radius_x_m=1100, radius_y_m=650,
        radius_z_m=300, resistivity_ohm_m=8, azimuth_deg=30,
        dip_deg=18, transition_fraction=0.20,
    )
    geology = insert_lenses(base, [lens])
    xx, yy = np.meshgrid(grid.x_m, grid.y_m)
    elevation = (
        620 + 85 * np.exp(-((xx + 900) / 1500) ** 2 - ((yy - 300) / 1100) ** 2)
        - 45 * np.exp(-((xx - 1500) / 900) ** 2 - ((yy + 500) / 700) ** 2)
        + 18 * np.sin(xx / 1200)
    )
    topo = TopographicSurface(
        grid, elevation, float(np.max(elevation)), source="deterministic tutorial surface"
    )
    return grid, geology, topo


def plot_geology(grid, geology, topo) -> None:
    rho = np.log10(geology.resistivity_ohm_m)
    iz = int(np.argmin(abs(grid.z_m - 900)))
    iy = len(grid.y_m) // 2
    ix = len(grid.x_m) // 2
    fig, axes = plt.subplots(1, 3, figsize=(14.4, 4.5), constrained_layout=True)
    kw = dict(cmap="turbo", vmin=0.7, vmax=3.5, shading="auto")
    im = axes[0].pcolormesh(grid.x_m / 1000, grid.y_m / 1000, rho[iz], **kw)
    axes[0].set(title=f"Horizontal slice: z = {grid.z_m[iz]:.0f} m", xlabel="x (km)", ylabel="y (km)")
    axes[1].pcolormesh(grid.x_m / 1000, grid.z_m / 1000, rho[:, iy, :], **kw)
    axes[1].plot(grid.x_m / 1000, topo.surface_depth_m[iy] / 1000, color="black", lw=1.5)
    axes[1].invert_yaxis(); axes[1].set(title="Along-line section: y = 0", xlabel="x (km)", ylabel="depth below datum (km)")
    axes[2].pcolormesh(grid.y_m / 1000, grid.z_m / 1000, rho[:, :, ix], **kw)
    axes[2].plot(grid.y_m / 1000, topo.surface_depth_m[:, ix] / 1000, color="black", lw=1.5)
    axes[2].invert_yaxis(); axes[2].set(title="Cross-line section: x = 0", xlabel="y (km)", ylabel="depth below datum (km)")
    fig.colorbar(im, ax=axes, shrink=.82, pad=.02, label=r"$\log_{10}\rho$ (ohm m)")
    save(fig, "essential3d_geology_volume_slices.png")


def plot_mesh(grid, geology, topo):
    model = build_solver_mesh(
        grid, resistivity_ohm_m=geology.resistivity_ohm_m,
        frequencies_hz=[100, 10, 1, .1], topography=topo,
        design=MeshDesign(horizontal_padding_cells=4, bottom_padding_cells=5, air_layers=4),
    )
    sigma = model.conductivity_s_m
    iy = sigma.shape[1] // 2
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.8), constrained_layout=True)
    x = model.mesh.cell_centres_m["x"] / 1000
    z = model.mesh.cell_centres_m["z"] / 1000
    im = axes[0].pcolormesh(x, z, np.log10(1 / sigma[:, iy, :]), cmap="turbo", shading="auto", vmin=.7, vmax=8)
    axes[0].invert_yaxis(); axes[0].set(title="Padded Maxwell mesh: central y slice", xlabel="x (km)", ylabel="depth below datum (km)")
    axes[1].imshow(model.earth_mask[:, iy, :], aspect="auto", origin="upper", cmap="cividis")
    axes[1].set(title="Active earth (yellow) and air (blue)", xlabel="x-cell index", ylabel="z-cell index")
    fig.colorbar(im, ax=axes[0], label=r"$\log_{10}\rho$ (ohm m)")
    save(fig, "essential3d_maxwell_mesh.png")
    return model


def plot_gate(model):
    cap = MT3DAdapter().capabilities
    checks = [
        ("3-D geometry", 3 in cap.dimensions),
        ("cell ceiling", model.quality.cell_count <= cap.maximum_cells),
        ("nonuniform padding", cap.supports_nonuniform_mesh),
        ("inactive/air cells", cap.supports_inactive_cells),
        ("topography", cap.supports_topography),
        ("production validation", False),
    ]
    fig, ax = plt.subplots(figsize=(9.4, 3.8), constrained_layout=True)
    y = np.arange(len(checks)); ok = np.array([v for _, v in checks])
    ax.barh(y, np.ones_like(y), color=np.where(ok, "#2f8f72", "#c6534f"))
    ax.set_yticks(y, [name for name, _ in checks]); ax.set_xlim(0, 1); ax.set_xticks([]); ax.invert_yaxis()
    for yi, value in zip(y, ok): ax.text(.5, yi, "PASS" if value else "STOP", ha="center", va="center", color="white", weight="bold")
    ax.set_title("Bundled mt3d backend gate for this topographic mesh")
    save(fig, "essential3d_backend_gate.png")


def main() -> int:
    grid, geology, topo = build_model()
    plot_geology(grid, geology, topo)
    model = plot_mesh(grid, geology, topo)
    plot_gate(model)
    print("geology shape:", grid.shape)
    print("resistivity range (ohm m):", f"{geology.resistivity_ohm_m.min():.1f}", f"{geology.resistivity_ohm_m.max():.1f}")
    print("topographic relief (m):", f"{topo.relief_m:.1f}")
    print("Maxwell mesh shape:", model.mesh.shape)
    print("Maxwell cells:", model.quality.cell_count)
    print("mesh warnings:", len(model.quality.warnings))
    print("mt3d maximum cells:", MT3DAdapter().capabilities.maximum_cells)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
