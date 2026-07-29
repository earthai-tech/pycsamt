"""Generate figures for the AI inversion tutorial.

Every figure here comes from actually running the real pipeline
(:func:`~pycsamt.ai.domain_gap.audit_survey`,
:class:`~pycsamt.agents.Inv2DAgent`, :class:`~pycsamt.agents.Inv3DAgent`,
:func:`~pycsamt.ai.validation.flag_out_of_distribution`) against the K2
CSAMT line. K2 is a local-only survey, not bundled with the repository --
this script is a no-op wherever ``K2_DIR`` does not exist on disk.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".tmp" / "matplotlib"))
os.environ.setdefault("PYCSAMT_DATA", str(ROOT / ".tmp" / "docs-data"))

K2_DIR = ROOT / "k2_corrected"
IMAGE_DIR = (
    ROOT
    / "docs"
    / "source"
    / "images"
    / "tutorials"
    / "ai_inversion_from_corrected_edis"
)

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

_K2_FULL_FREQS_HZ = np.array(
    [15.8, 20.0, 25.1, 31.6, 39.8, 50.1, 63.1, 79.4, 100.0, 126.0,
     158.0, 200.0, 251.0, 316.0, 398.0, 501.0, 631.0, 794.0, 1000.0,
     1260.0, 1580.0, 2000.0, 2510.0, 3160.0, 3980.0, 5010.0, 6310.0,
     7940.0, 10000.0]
)


def _save(fig: plt.Figure, name: str) -> None:
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(IMAGE_DIR / name, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _k2_clean_sites():
    """Load the K2 survey, ordered, with exact-duplicate stations
    dropped (four station pairs share identical EDI coordinates).
    """
    from pycsamt.emtools._core import ensure_sites

    sites = ensure_sites(K2_DIR, recursive=False, verbose=0).ordered()
    lat = np.array([s.coords[0] for s in sites])
    lon = np.array([s.coords[1] for s in sites])
    names = [s.name for s in sites]
    keep = [names[0]]
    last = (lat[0], lon[0])
    for i in range(1, len(names)):
        if (lat[i], lon[i]) == last:
            continue
        keep.append(names[i])
        last = (lat[i], lon[i])
    return sites.select(names=keep).ordered()


def make_geology_topography_mesh() -> None:
    """Build the tutorial prior, terrain surface, and padded solver mesh."""
    from pycsamt.ai.geology import (
        ElectricalLayer,
        GaussianCorrelation,
        GeologyGrid,
        generate_layered_geology,
        interpolate_topography,
    )
    from pycsamt.forward.maxwell.mesh import MeshDesign, build_solver_mesh

    grid = GeologyGrid.regular_2d(nx=21, nz=20, dx_m=80.8, dz_m=100.0)
    layers = (
        ElectricalLayer("conductive cover", 35.0, log10_std=0.12,
                        heterogeneity=GaussianCorrelation(350, 100)),
        ElectricalLayer("weathered host", 140.0, log10_std=0.16,
                        heterogeneity=GaussianCorrelation(500, 160)),
        ElectricalLayer("resistive basement", 1200.0, log10_std=0.10,
                        heterogeneity=GaussianCorrelation(700, 220)),
    )
    geology = generate_layered_geology(
        grid, layers, [380.0, 1050.0], seed=17,
        interface_relief_std_m=[45.0, 90.0],
        interface_correlation=GaussianCorrelation(550, 150),
        minimum_thickness_m=120.0,
    )
    control_x = np.linspace(grid.x_m[0], grid.x_m[-1], 9)
    elevation = 285 + 42 * np.sin(control_x / 310) + 18 * np.cos(control_x / 145)
    topography = interpolate_topography(
        grid, control_x, elevation, source="corrected EDI elevations",
        interpolation_method="cubic",
    )
    solver = build_solver_mesh(
        grid, resistivity_ohm_m=geology.resistivity_ohm_m,
        frequencies_hz=[15.8, 10_000.0], topography=topography,
        design=MeshDesign(horizontal_padding_cells=5, bottom_padding_cells=5,
                          air_layers=5),
    )

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.7))
    extent = (0, grid.x_m[-1] / 1000, grid.z_m[-1] / 1000, 0)
    image = axes[0].imshow(np.log10(geology.resistivity_ohm_m), extent=extent,
                           aspect="auto", cmap="viridis_r", vmin=1, vmax=3.3)
    for interface in geology.interface_depth_m:
        axes[0].plot(grid.x_m / 1000, interface / 1000, "w-", lw=1)
    axes[0].set(title="Geological target before meshing", xlabel="Distance (km)",
                ylabel="Depth below reference (km)")
    fig.colorbar(image, ax=axes[0], label=r"$\log_{10}\rho$ [$\Omega\cdot$m]")

    topo_image = axes[1].imshow(topography.earth_mask(), extent=extent,
                                aspect="auto", cmap="Greys", vmin=0, vmax=1)
    axes[1].plot(grid.x_m / 1000, topography.surface_depth_m / 1000,
                 color="#f15a29", lw=2, label="interpolated terrain")
    axes[1].scatter(control_x / 1000,
                    topography.surface_depth_m[np.searchsorted(grid.x_m, control_x)] / 1000,
                    marker="v", color="#2563eb", s=24, label="EDI elevations")
    axes[1].set(title=f"Terrain mask: relief={topography.relief_m:.0f} m",
                xlabel="Distance (km)", ylabel="Depth below reference (km)")
    axes[1].legend(facecolor="white", framealpha=0.9, fontsize=8,
                   loc="lower left")
    fig.colorbar(topo_image, ax=axes[1], ticks=[0, 1], label="air (0) / earth (1)")

    x_edges = solver.mesh.x_edges_m / 1000
    z_edges = solver.mesh.z_edges_m / 1000
    mesh_values = np.log10(1.0 / solver.conductivity_s_m)
    mesh_image = axes[2].pcolormesh(x_edges, z_edges, mesh_values,
                                    shading="flat", cmap="viridis_r", vmin=1, vmax=8)
    axes[2].set(title=f"Padded Maxwell mesh: {solver.quality.cell_count:,} cells",
                xlabel="Padded distance (km)", ylabel="Depth (km)")
    axes[2].invert_yaxis()
    axes[2].axvline(grid.x_m[0] / 1000, color="white", ls="--", lw=0.9)
    axes[2].axvline(grid.x_m[-1] / 1000, color="white", ls="--", lw=0.9,
                    label="geological core")
    axes[2].legend(frameon=False, fontsize=8)
    fig.colorbar(mesh_image, ax=axes[2], label=r"$\log_{10}\rho$ [$\Omega\cdot$m]")
    fig.suptitle("Prior geometry becomes a physics-facing air/earth mesh", fontsize=13)
    fig.tight_layout()
    _save(fig, "geology_topography_mesh.png")


def make_maxwell_training_pair() -> None:
    """Generate one solved 2-D training pair and expose model/response structure."""
    from pycsamt.ai.geology import GeologyGrid
    from pycsamt.ai.training.dataset2d import Maxwell2DDatasetConfig, generate_2d_maxwell_dataset

    grid = GeologyGrid.regular_2d(nx=21, nz=20, dx_m=80.8, dz_m=100.0)
    frequencies = np.array([15.8, 31.6, 63.1, 126, 251, 501, 1000, 2000, 3980, 7940])
    config = Maxwell2DDatasetConfig(
        dataset_id="corrected-edi-tutorial-figure", grid=grid,
        correlation_length_x_m=(200.0, 600.0),
        correlation_length_z_m=(50.0, 200.0), frequencies_hz=frequencies,
        station_x_m=grid.x_m, n_realizations=1, seed=23,
        log_resistivity_mean=2.0, log_resistivity_std=0.5,
        components=("zxy",), validation_fraction=0, test_fraction=0,
    )
    sample = generate_2d_maxwell_dataset(config).samples[0]
    impedance = sample.survey.impedance[:, :, 0]
    mu0 = 4e-7 * np.pi
    apparent = np.abs(impedance) ** 2 / (2 * np.pi * frequencies[None, :] * mu0)
    phase = np.angle(impedance, deg=True)

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.5))
    image = axes[0].imshow(np.log10(sample.resistivity_ohm_m), aspect="auto",
                           extent=(0, 1.616, 2.0, 0), cmap="viridis_r")
    axes[0].set(title="Target earth model", xlabel="Distance (km)", ylabel="Depth (km)")
    fig.colorbar(image, ax=axes[0], label=r"$\log_{10}\rho$ [$\Omega\cdot$m]")
    for ax, values, title, label, cmap in (
        (axes[1], np.log10(apparent).T, "Maxwell input: apparent resistivity",
         r"$\log_{10}\rho_a$ [$\Omega\cdot$m]", "viridis_r"),
        (axes[2], phase.T, "Maxwell input: impedance phase", "Phase (degrees)", "twilight"),
    ):
        response_image = ax.imshow(values, origin="lower", aspect="auto",
                                   extent=(0, 1.616, np.log10(frequencies[0]),
                                           np.log10(frequencies[-1])), cmap=cmap)
        ax.set(xlabel="Station distance (km)", ylabel=r"$\log_{10}$ frequency (Hz)",
               title=title)
        fig.colorbar(response_image, ax=ax, label=label)
    fig.suptitle(f"One supervised pair: {sample.mesh_cells:,} mesh cells, "
                 f"solver residual={sample.relative_residual:.1e}", fontsize=12.5)
    fig.tight_layout()
    _save(fig, "maxwell_training_pair.png")


def make_training_convergence_smoke() -> None:
    """Run a small bundled-data fit so tutorial readers can audit learning dynamics."""
    from pycsamt.agents import Inv2DAgent
    from pycsamt.emtools._core import ensure_sites

    sites = ensure_sites(ROOT / "data" / "AMT" / "WILLY_data" / "L18PLT",
                         recursive=True, verbose=0).ordered()
    names = [site.name for site in sites][:8]
    subset = sites.select(names=names).ordered()
    frequencies = np.logspace(0, 3, 10)
    result = Inv2DAgent(
        physics="mt1d", n_depth=12, n_stations_per_profile=8,
        n_train_profiles=24, epochs=8, depth_max=1200,
        api_key=None,
    ).execute({"sites": subset, "freqs": frequencies})
    if result.status != "success":
        raise RuntimeError(result.error)
    history = result.data["inverter"]._history
    train = np.asarray(history["train_loss"])
    valid = np.asarray(history["val_loss"])
    epochs = np.arange(1, len(train) + 1)
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.2))
    axes[0].plot(epochs, train, "o-", label="training", color="#2563eb")
    axes[0].plot(epochs, valid, "s-", label="validation", color="#f15a29")
    axes[0].set(xlabel="Epoch", ylabel="Loss", title="Learning curves are evidence, not decoration")
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False)
    gap = valid - train
    axes[1].bar(epochs, gap, color=np.where(gap > 0, "#dc2626", "#16a34a"))
    axes[1].axhline(0, color="#111827", lw=1)
    axes[1].set(xlabel="Epoch", ylabel="Validation - training loss",
                title=f"Generalization gap; field RMS={result.data['rms_global']:.2f}")
    axes[1].grid(alpha=0.25, axis="y")
    fig.suptitle("Bundled WILLY smoke run (24 profiles, 8 epochs): dynamics, not acceptance", fontsize=12)
    fig.tight_layout()
    _save(fig, "training_convergence_smoke.png")


def make_inv2d_topography() -> None:
    """Real physics="mt2d" Inv2DAgent run on all 82 unique K2 stations."""
    from pycsamt.agents import Inv2DAgent

    clean = _k2_clean_sites()
    freqs = _K2_FULL_FREQS_HZ[::3]
    agent = Inv2DAgent(
        physics="mt2d",
        n_depth=20,
        n_stations_per_profile=len(clean),
        n_train_profiles=24,
        epochs=30,
        depth_max=2000.0,
        station_spacing_m=20.2,
        correlation_length_x_m=(200.0, 600.0),
        correlation_length_z_m=(50.0, 200.0),
        log_resistivity_mean=2.0,
        log_resistivity_std=0.5,
        lambda_x=0.01,
        lambda_z=0.005,
        lambda_tv=0.002,
        mesh_safety_factor=4.0,
    )
    result = agent.execute({"sites": clean, "freqs": freqs, "topography": True})
    if result.status == "failed":
        raise RuntimeError(result.error)
    fig = result["figures"]["topography_section"]
    _save(fig, "inv2d_topography_section.png")
    print(
        "K2 2-D:",
        {"rms_global": result.data["rms_global"],
         "recovery": result.data.get("mt2d_recovery"),
         "log10_rho_range": (
             float(np.nanmin(result.data["pred_section"])),
             float(np.nanmax(result.data["pred_section"])),
         )},
    )
    history = result.data["inverter"]._history
    train = np.asarray(history["train_loss"])
    valid = np.asarray(history["val_loss"])
    epoch = np.arange(1, len(train) + 1)
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.2))
    axes[0].plot(epoch, train, "o-", color="#2563eb", label="training")
    axes[0].plot(epoch, valid, "s-", color="#f15a29", label="validation")
    axes[0].set(xlabel="Epoch", ylabel="Loss", title="K2 Maxwell-trained U-Net history")
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False)
    recovery = result.data.get("mt2d_recovery") or {}
    metric_names = ["rmse", "mae", "r2"]
    metric_values = [float(recovery.get(name, np.nan)) for name in metric_names]
    axes[1].bar(metric_names, metric_values,
                color=["#2563eb", "#60a5fa", "#dc2626"])
    axes[1].axhline(0, color="#111827", lw=1)
    axes[1].set(ylabel="Metric value", title="Held-out synthetic recovery")
    axes[1].grid(alpha=0.25, axis="y")
    fig.suptitle("Training fit and recovery skill must be read together", fontsize=12.5)
    fig.tight_layout()
    _save(fig, "inv2d_training_and_recovery.png")


def make_survey_audit() -> None:
    """Visualize the K2 geometry, frequency coverage, and dimensionality gate."""
    from pycsamt.ai.domain_gap import audit_survey

    sites = _k2_clean_sites()
    report = audit_survey(sites, verbose=0)
    names = [site.name for site in sites]
    elevations = np.array([site.coords[2] for site in sites], dtype=float)
    counts = np.array([len(site.freq) for site in sites], dtype=int)
    dim = report.dimensionality
    fractions = np.array([dim.frac_1d, dim.frac_2d, dim.frac_3d])
    lat = np.array([site.coords[0] for site in sites], dtype=float)
    lon = np.array([site.coords[1] for site in sites], dtype=float)
    lat0 = np.radians(np.nanmean(lat))
    dx = np.diff(lon) * 111_320.0 * np.cos(lat0)
    dy = np.diff(lat) * 110_574.0
    spacing = np.sqrt(dx**2 + dy**2)

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.4))
    index = np.arange(len(sites))
    axes[0, 0].plot(index, elevations, color="#2563eb", lw=1.8)
    axes[0, 0].fill_between(index, elevations.min(), elevations,
                            color="#bfdbfe", alpha=0.7)
    axes[0, 0].set(xlabel="Ordered station index", ylabel="Elevation (m)",
                   title="Terrain carried by corrected EDI headers")
    axes[0, 0].grid(alpha=0.2)
    axes[0, 1].bar(index, counts, color=np.where(counts == counts.max(), "#16a34a", "#f59e0b"))
    axes[0, 1].set(xlabel="Ordered station index", ylabel="Available frequencies",
                   title="Frequency support is not uniform")
    axes[0, 1].grid(alpha=0.2, axis="y")
    axes[1, 0].bar(["1-D", "2-D", "3-D"], fractions,
                   color=["#16a34a", "#2563eb", "#dc2626"])
    axes[1, 0].set(ylim=(0, 1), ylabel="Fraction of station-period samples",
                   title="Tensor dimensionality controls model choice")
    for i, value in enumerate(fractions):
        axes[1, 0].text(i, value + 0.025, f"{value:.1%}", ha="center")
    axes[1, 0].grid(alpha=0.2, axis="y")
    axes[1, 1].hist(spacing[np.isfinite(spacing)], bins=18, color="#64748b",
                    edgecolor="white")
    axes[1, 1].axvline(report.station_spacing_m["median"], color="#f15a29",
                       ls="--", label=f"median={report.station_spacing_m['median']:.1f} m")
    axes[1, 1].set(xlabel="Adjacent spacing (m)", ylabel="Count",
                   title="Geometry audit precedes meshing")
    axes[1, 1].legend(frameon=False)
    fig.suptitle("Corrected does not mean inversion-ready: audit K2 before training", fontsize=13)
    fig.tight_layout()
    _save(fig, "survey_audit.png")


def make_inv3d_graph_diagnostic() -> None:
    """Real Inv3DAgent (GCN) run on the full 82-station clean K2 line."""
    from pycsamt.agents import Inv3DAgent

    clean = _k2_clean_sites()
    lat = np.array([s.coords[0] for s in clean])
    lon = np.array([s.coords[1] for s in clean])
    elev = np.array([s.coords[2] for s in clean])
    lat0 = np.radians(lat[0])
    x_m = (lon - lon[0]) * 111_320.0 * np.cos(lat0)
    y_m = (lat - lat[0]) * 110_574.0
    chain_km = np.concatenate(
        [[0.0], np.cumsum(np.sqrt(np.diff(x_m) ** 2 + np.diff(y_m) ** 2))]
    ) / 1000.0
    agent = Inv3DAgent(
        n_layers=6, epochs=30, n_train_profiles=150, n_mc=20,
        radius=300.0, depth_max=2000.0,
    )
    result = agent.execute({
        "sites": clean,
        "freqs": _K2_FULL_FREQS_HZ,
        "topography": {"elevation_m": elev, "chainage_km": chain_km},
    })
    if result.status == "failed":
        raise RuntimeError(result.error)
    raw_log_rho = np.asarray(result.data["pred_rho"])
    physical = np.isfinite(raw_log_rho) & (raw_log_rho >= 0.0) & (raw_log_rho <= 5.0)
    bounded_log_rho = np.where(physical, raw_log_rho, np.nan)
    from pycsamt.topo import plot_topo_section

    rejected_model = {
        "pred_rho": bounded_log_rho,
        "depths_km": result.data["depths_km"],
        "station_names": result.data["station_names"],
        "rms_global": result.data["rms_global"],
    }
    ax = plot_topo_section(
        rejected_model,
        elevation=elev,
        chainage=chain_km,
        station_names=result.data["station_names"],
        station_x=chain_km,
        topo_source="array",
        model_unit="km",
        depth_max=2.0,
        title="Rejected GCN candidate: only physically bounded cells are shown",
    )
    ax.text(
        0.5, 0.96,
        f"REJECTED — {100 * (1 - physical.mean()):.1f}% of cells outside "
        r"$1$–$10^5\ \Omega\,\mathrm{m}$ display bounds",
        transform=ax.transAxes, ha="center", va="top", color="white",
        fontsize=10, fontweight="bold",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "#b91c1c",
              "edgecolor": "white", "alpha": 0.94},
    )
    # Deliberately do not publish a subsurface section from a rejected
    # candidate.  The raw values are retained below for gate diagnostics.
    plt.close(ax.get_figure())
    print(
        "K2 graph candidate:",
        {"rms_global": result.data["rms_global"],
         "log10_rho_range": (float(np.nanmin(raw_log_rho)),
                              float(np.nanmax(raw_log_rho))),
         "uncertainty_range": (
             float(np.nanmin(result.data["pred_uncertainty"])),
             float(np.nanmax(result.data["pred_uncertainty"])),
         ),
         "fraction_outside_display_bounds": float(1 - physical.mean())},
    )
    adjacency = np.asarray(result.data["adjacency"])
    uncertainty = np.asarray(result.data["pred_uncertainty"])
    neighbours = (adjacency > 0).sum(axis=1) - 1
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.3))
    axes[0].imshow(adjacency, origin="lower", aspect="auto", cmap="Blues")
    axes[0].set(xlabel="Station node", ylabel="Station node",
                title="Radius graph adjacency")
    axes[1].plot(np.arange(len(neighbours)), neighbours, color="#2563eb")
    axes[1].axhline(neighbours.mean(), color="#f15a29", ls="--",
                    label=f"mean={neighbours.mean():.1f}")
    axes[1].set(xlabel="Station node", ylabel="Connected neighbours",
                title="Spatial context varies at line ends")
    axes[1].legend(frameon=False)
    axes[1].grid(alpha=0.2)
    image = axes[2].imshow(uncertainty.T, origin="upper", aspect="auto",
                           extent=(0, uncertainty.shape[0] - 1, 2.0, 0), cmap="magma")
    axes[2].set(xlabel="Station node", ylabel="Depth (km)",
                title="MC-dropout spread")
    fig.colorbar(image, ax=axes[2], label=r"$\sigma(\log_{10}\rho)$")
    fig.suptitle("The 3-D-labelled workflow is a graph model, not a 3-D Maxwell solve",
                 fontsize=12.5)
    fig.tight_layout()
    _save(fig, "inv3d_graph_uncertainty.png")


def make_ood_scores() -> None:
    """Mahalanobis OOD score of every clean K2 station against the
    per-station synthetic 1-D feature distribution :class:`Inv3DAgent`
    trains on -- the finding that explains the previous figure.
    """
    from pycsamt.agents.ai_inversion import _z_to_features
    from pycsamt.ai.validation import flag_out_of_distribution
    from pycsamt.emtools._core import _get_z_block
    from pycsamt.forward.batch import generate_dataset

    clean = _k2_clean_sites()
    freqs = _K2_FULL_FREQS_HZ
    n = freqs.size

    def summarize(X: np.ndarray) -> np.ndarray:
        rho, pha = X[:, :n], X[:, n : 2 * n]
        return np.stack(
            [rho.mean(1), rho.std(1), pha.mean(1), pha.std(1)], axis=1
        )

    x_obs = []
    for s in clean:
        z_obj, z, fr = _get_z_block(s)
        feat = _z_to_features(z_obj, z, fr, freqs)
        x_obs.append(feat[: 2 * n])
    feat_obs = summarize(np.array(x_obs))

    ds = generate_dataset(
        solver="mt1d", n_samples=300, freqs=freqs, n_layers=6,
        noise_level=0.03, seed=1, n_jobs=1, verbose=False,
    )
    feat_train = summarize(ds.X[:, : 2 * n])

    report = flag_out_of_distribution(
        feat_obs, feat_train, method="mahalanobis", quantile=0.95
    )
    order = np.argsort(report.scores)
    scores = report.scores[order]
    flagged = report.flagged[order]

    fig, ax = plt.subplots(figsize=(11.0, 4.6), constrained_layout=True)
    colors = np.where(flagged, "#dc2626", "#16a34a")
    ax.bar(np.arange(scores.size), scores, color=colors, width=0.85)
    ax.axhline(
        report.threshold, color="#111827", ls="--", lw=1.3,
        label=f"threshold (95th pct. of training self-scores) = "
        f"{report.threshold:.2f}",
    )
    ax.set_yscale("log")
    ax.set_xlabel("K2 station, sorted by OOD score")
    ax.set_ylabel("Mahalanobis distance (log scale)")
    ax.set_xlim(-0.6, scores.size - 0.4)
    ax.legend(loc="upper left", fontsize=9)
    n_flagged = int(flagged.sum())
    ax.set_title(
        f"{n_flagged}/{scores.size} K2 stations flagged out-of-distribution "
        "against the GCN's synthetic 1-D training features",
        fontsize=11,
    )
    _save(fig, "ood_scores.png")


def main() -> int:
    make_geology_topography_mesh()
    make_maxwell_training_pair()
    make_training_convergence_smoke()
    if not K2_DIR.is_dir():
        print(f"skip K2-only figures: {K2_DIR} not found (local-only dataset)")
        print(f"images: {IMAGE_DIR.relative_to(ROOT)}")
        return 0
    make_survey_audit()
    make_inv2d_topography()
    make_inv3d_graph_diagnostic()
    make_ood_scores()
    print(f"images: {IMAGE_DIR.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
