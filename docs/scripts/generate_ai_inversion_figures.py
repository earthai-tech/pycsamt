#!/usr/bin/env python
"""Generate AI inversion figures for the user guide.

The figures are deterministic documentation examples.  They avoid long model
training while still showing the visual checks users should perform around an
AI inversion workflow.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
os.environ.setdefault(
    "MPLCONFIGDIR", str(PROJECT_ROOT / ".tmp" / "matplotlib")
)
os.environ.setdefault("PYCSAMT_DOCS_BUILD", "1")
os.environ.setdefault("PYCSAMT_DATA", str(PROJECT_ROOT / ".tmp" / "docs-data"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from pycsamt.agents import (
    AIInversionAgent,
    EnsembleAgent,
    Inv2DAgent,
    Inv3DAgent,
    PINNInversionAgent,
)
from pycsamt.ai.data.normalization import ComplexZScore
from pycsamt.ai.domain_gap import (
    CorruptionConfig,
    apply_corruption_suite,
    compare_survey_distributions,
    survey_data_from_sites,
)
from pycsamt.ai.data.contracts import SurveyData
from pycsamt.ai.inversion import (
    EMInverter1D,
    PINNInverter2D,
    sites_to_coords_3d,
    sites_to_features_1d,
    sites_to_obs_1d,
    sites_to_obs_2d,
)
from pycsamt.ai.nets import build_adjacency
from pycsamt.ai.experiments import (
    AcceptanceCriterion,
    DatasetReference,
    ExperimentConfig,
    SeedPlan,
)
from pycsamt.ai.training.dataset2d import (
    Maxwell2DDatasetConfig,
    generate_2d_maxwell_dataset,
)
from pycsamt.ai.training.dataset3d import (
    Maxwell3DDatasetConfig,
    generate_3d_maxwell_dataset,
)
from pycsamt.ai.validation import (
    flag_out_of_distribution,
    recovery_report,
    reliability_curve,
    response_residual_report,
)
from pycsamt.forward.maxwell import MaxwellMesh, MaxwellProblem, ReceiverSet
from pycsamt.forward.maxwell.benchmarks import half_space_impedance
from pycsamt.forward.maxwell.mt2d import MT2DAdapter
from pycsamt.ai.plot.convergence import plot_convergence
from pycsamt.ai.plot.diagnostics import plot_uncertainty_bands
from pycsamt.ai.plot.inversion import plot_inversion_result_2d
from pycsamt.emtools import frequency_for_depth
from pycsamt.emtools._core import ensure_sites
from pycsamt.emtools.dimensionality import classify_dimensionality
from pycsamt.forward.batch import generate_dataset

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "source" / "images" / "user_guide" / "ai_inversion"


def _save(fig: plt.Figure, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / name, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _layer_model(
    depths: np.ndarray, rho_layers: list[float], breaks: list[float]
) -> np.ndarray:
    log_rho = np.empty_like(depths, dtype=float)
    edges = [0.0, *breaks, float(depths[-1]) + 1.0]
    for lo, hi, rho in zip(edges[:-1], edges[1:], rho_layers):
        log_rho[(depths >= lo) & (depths < hi)] = np.log10(rho)
    return log_rho


def make_workflow() -> None:
    fig, ax = plt.subplots(figsize=(11.0, 3.3))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    labels = [
        ("Survey QC", "clean EDI, static shift,\nfrequency grid"),
        ("Synthetic models", "sample resistivity,\nthickness, noise"),
        ("Forward responses", "MT/CSAMT/TEM\ntraining pairs"),
        ("Train inverter", "ResNet, CNN1D,\nU-Net, ensemble"),
        ("Predict + validate", "sections, residuals,\nuncertainty"),
    ]
    xs = np.linspace(0.08, 0.92, len(labels))
    colors = ["#eff6ff", "#f0fdf4", "#fff7ed", "#f5f3ff", "#fef2f2"]
    borders = ["#2563eb", "#16a34a", "#ea580c", "#7c3aed", "#dc2626"]

    for i, ((title, body), x, face, edge) in enumerate(
        zip(labels, xs, colors, borders)
    ):
        box = FancyBboxPatch(
            (x - 0.078, 0.35),
            0.156,
            0.34,
            boxstyle="round,pad=0.018,rounding_size=0.025",
            linewidth=1.4,
            edgecolor=edge,
            facecolor=face,
        )
        ax.add_patch(box)
        ax.text(
            x,
            0.61,
            title,
            ha="center",
            va="center",
            fontsize=10.5,
            weight="bold",
            color="#111827",
        )
        ax.text(
            x,
            0.47,
            body,
            ha="center",
            va="center",
            fontsize=8.4,
            color="#374151",
            linespacing=1.25,
        )
        if i < len(labels) - 1:
            ax.add_patch(
                FancyArrowPatch(
                    (x + 0.087, 0.52),
                    (xs[i + 1] - 0.087, 0.52),
                    arrowstyle="-|>",
                    mutation_scale=14,
                    linewidth=1.4,
                    color="#64748b",
                )
            )

    ax.text(
        0.5,
        0.18,
        "The trained network is a fast surrogate: reliability depends on training coverage and post-prediction QC.",
        ha="center",
        va="center",
        fontsize=9.5,
        color="#475569",
    )
    _save(fig, "workflow.png")


def make_training_distribution() -> None:
    rng = np.random.default_rng(2026)
    depths = np.linspace(0, 1500, 180)
    periods = np.logspace(-2, 2, 80)

    fig, axes = plt.subplots(
        1, 2, figsize=(11.0, 4.1), gridspec_kw={"width_ratios": [0.95, 1.25]}
    )
    ax_m, ax_r = axes

    response_stack = []
    for _ in range(90):
        breaks = sorted(rng.uniform([120, 350, 650], [260, 720, 1250]))
        rho_layers = 10 ** rng.uniform(0.4, 3.7, size=4)
        model = _layer_model(depths, list(rho_layers), breaks)
        if _ < 16:
            ax_m.plot(model, depths, lw=0.9, alpha=0.45, color="#2563eb")
        conductance = np.interp(
            np.log10(periods),
            [-2, 2],
            [model[:30].mean(), model[-40:].mean()],
        )
        bowl = 0.22 * np.sin(
            np.log(periods) * rng.uniform(0.7, 1.2) + rng.uniform(-1.2, 1.2)
        )
        response_stack.append(
            conductance + bowl + rng.normal(0, 0.035, size=periods.size)
        )

    response_stack = np.vstack(response_stack)
    p10, p50, p90 = np.percentile(response_stack, [10, 50, 90], axis=0)
    observed = p50 + 0.10 * np.exp(-((np.log10(periods) - 0.35) ** 2) / 0.55)

    ax_m.invert_yaxis()
    ax_m.set_xlabel(r"$\log_{10}(\rho)$  ($\Omega$ m)")
    ax_m.set_ylabel("Depth (m)")
    ax_m.set_title("Sampled earth models", fontsize=10)
    ax_m.grid(alpha=0.25)

    ax_r.fill_between(
        periods,
        p10,
        p90,
        color="#93c5fd",
        alpha=0.35,
        label="training envelope",
    )
    ax_r.plot(periods, p50, color="#2563eb", lw=1.8, label="synthetic median")
    ax_r.plot(
        periods,
        observed,
        color="#dc2626",
        lw=1.6,
        ls="--",
        label="observed survey",
    )
    ax_r.set_xscale("log")
    ax_r.set_xlabel("Period (s)")
    ax_r.set_ylabel(r"Response feature, $\log_{10}(\rho_a)$")
    ax_r.set_title(
        "Observed data should sit inside the synthetic envelope", fontsize=10
    )
    ax_r.legend(frameon=False, fontsize=8)
    ax_r.grid(alpha=0.25, which="both")

    fig.suptitle("Training distribution coverage", fontsize=13, y=1.02)
    fig.tight_layout()
    _save(fig, "training_distribution.png")


def make_convergence() -> None:
    epochs = np.arange(1, 71)
    rng = np.random.default_rng(7)
    histories = []
    for offset in [0.0, 0.08, -0.05]:
        train = (
            0.42 * np.exp(-epochs / 19.5)
            + 0.018
            + rng.normal(0, 0.003, epochs.size)
            + offset * 0.002
        )
        val = (
            0.48 * np.exp(-epochs / 16.0)
            + 0.028
            + 0.000025 * (epochs - 48).clip(min=0) ** 2
        )
        val += rng.normal(0, 0.004, epochs.size) + offset * 0.003
        lr = np.where(epochs < 45, 1e-3, 3e-4)
        histories.append(
            {
                "train_loss": train.clip(1e-4),
                "val_loss": val.clip(1e-4),
                "lr": lr,
            }
        )

    fig = plot_convergence(
        histories,
        smoothing=0.25,
        best_epoch=45,
        title="AI inverter training convergence",
    )
    _save(fig, "training_convergence.png")


def make_predicted_section() -> None:
    depths = np.linspace(0, 1400, 70)
    stations = np.linspace(0, 18, 34)
    xx, zz = np.meshgrid(stations, depths)

    background = 2.35 + 0.45 * np.tanh((zz - 500) / 350)
    conductor = -0.95 * np.exp(
        -((xx - 7.2) ** 2 / 10.0 + (zz - 420) ** 2 / 90000)
    )
    basement = 0.48 * np.exp(
        -((xx - 14.5) ** 2 / 18.0 + (zz - 1000) ** 2 / 180000)
    )
    log_true = background + conductor + basement
    log_pred = (
        log_true
        + 0.08 * np.sin(xx / 2.4) * np.exp(-zz / 1200)
        - 0.05 * np.exp(-((xx - 4) ** 2) / 4)
    )

    train_loss = 0.38 * np.exp(-np.arange(1, 55) / 14) + 0.018
    val_loss = 0.42 * np.exp(-np.arange(1, 55) / 12) + 0.024

    fig = plot_inversion_result_2d(
        log_pred,
        log_true=log_true,
        depths=depths,
        stations=stations,
        train_loss=train_loss,
        val_loss=val_loss,
        fault_positions=[10.2],
        annotations=[
            {
                "text": "conductive zone",
                "xy": (7.1, 420),
                "xytext": (2.5, 250),
                "arrowprops": {
                    "arrowstyle": "->",
                    "color": "#111827",
                    "lw": 0.8,
                },
                "bbox": {
                    "boxstyle": "round,pad=0.25",
                    "fc": "white",
                    "ec": "#cbd5e1",
                    "alpha": 0.85,
                },
            }
        ],
        suptitle="AI inversion result - synthetic profile L22PLT",
        figsize=(12.2, 6.8),
        tick_every=4,
        tick_label_rotation=0,
    )
    _save(fig, "predicted_section.png")


def make_uncertainty() -> None:
    depths = np.linspace(0, 1500, 120)
    mean = (
        2.05
        + 0.35 * np.tanh((depths - 600) / 360)
        - 0.55 * np.exp(-((depths - 360) ** 2) / 70000)
    )
    spread = (
        0.07
        + 0.00012 * depths
        + 0.08 * np.exp(-((depths - 850) ** 2) / 150000)
    )
    truth = (
        mean
        - 0.10 * np.exp(-((depths - 430) ** 2) / 55000)
        + 0.05 * np.sin(depths / 260)
    )

    fig = plot_uncertainty_bands(
        depths,
        mean,
        mean + spread,
        mean - spread,
        y_true=truth,
        title="Station 18-023A - prediction interval",
    )
    ax = fig.axes[0]
    ax.text(
        0.03,
        0.06,
        "Spread increases where training coverage is weak",
        transform=ax.transAxes,
        fontsize=8.5,
        color="#475569",
        bbox={
            "boxstyle": "round,pad=0.3",
            "fc": "white",
            "ec": "#cbd5e1",
            "alpha": 0.9,
        },
    )
    _save(fig, "uncertainty_profile.png")


def make_willy_depth_support() -> None:
    """Plot an observed-data attenuation diagnostic, not an inversion depth."""
    line = PROJECT_ROOT / "data" / "AMT" / "WILLY_data" / "L18PLT"
    observations = sites_to_obs_1d(
        ensure_sites(line, recursive=True, verbose=0), comp="xy"
    )
    frequency = np.asarray(observations[0].freq, dtype=float)
    rho = np.vstack([item.rho_obs for item in observations]).astype(float)
    order = np.argsort(frequency)
    frequency = frequency[order]
    rho = rho[:, order]
    depth_km = 503.0 * np.sqrt(rho / frequency[None, :]) / 1000.0

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4))
    ax_rho, ax_depth = axes
    for curve in rho:
        ax_rho.loglog(frequency, curve, color="#93c5fd", lw=0.7, alpha=0.38)
    ax_rho.loglog(
        frequency,
        np.nanmedian(rho, axis=0),
        color="#1d4ed8",
        lw=2.1,
        label="station median",
    )
    ax_rho.set(
        xlabel="Frequency (Hz)",
        ylabel=r"XY apparent resistivity ($\Omega$ m)",
    )
    ax_rho.set_title("Observed response coverage")
    ax_rho.grid(alpha=0.22, which="both")
    ax_rho.legend(frameon=False)

    p10, p50, p90 = np.nanpercentile(depth_km, [10, 50, 90], axis=0)
    ax_depth.fill_between(
        frequency,
        p10,
        p90,
        color="#fdba74",
        alpha=0.45,
        label="station P10--P90",
    )
    ax_depth.loglog(
        frequency, p50, color="#c2410c", lw=2.1, label="station median"
    )
    ax_depth.axhline(
        2.0,
        color="#111827",
        ls="--",
        lw=1.2,
        label="2 km modelling target",
    )
    ax_depth.set(xlabel="Frequency (Hz)", ylabel="Skin-depth scale (km)")
    ax_depth.set_title("Attenuation scale, not recoverable depth")
    ax_depth.grid(alpha=0.22, which="both")
    ax_depth.legend(frameon=False, fontsize=8)
    fig.suptitle(
        "WILLY L18: observed bandwidth and depth support", fontsize=13
    )
    fig.tight_layout()
    _save(fig, "willy_l18_depth_support.png")


def _willy_features() -> tuple[np.ndarray, np.ndarray, list[str]]:
    line = PROJECT_ROOT / "data" / "AMT" / "WILLY_data" / "L18PLT"
    sites = ensure_sites(line, recursive=True, verbose=0)
    return sites_to_features_1d(
        sites, comp="xy", n_freqs=32, freq_min=1.01, freq_max=10_000.0
    )


def make_data_contracts_normalization() -> None:
    line = PROJECT_ROOT / "data" / "AMT" / "WILLY_data" / "L18PLT"
    sites = ensure_sites(line, recursive=True, verbose=0)
    survey = survey_data_from_sites(sites, recursive=False, verbose=0)
    state = ComplexZScore.fit(survey)
    normalized = state.transform_survey(survey)

    idx = survey.component_index("xy")
    z, _, valid = survey.component_data("xy")
    raw_real = z.real[valid]
    raw_imag = z.imag[valid]

    values = normalized.values[:, :, idx, :]
    mask = normalized.valid[:, :, idx, :]
    norm_real = values[..., 0][mask[..., 0]]
    norm_imag = values[..., 1][mask[..., 1]]

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.2))
    axes[0].hist(
        raw_real, bins=60, color="#1d4ed8", alpha=0.6, label="Re(Z)"
    )
    axes[0].hist(
        raw_imag, bins=60, color="#dc2626", alpha=0.6, label="Im(Z)"
    )
    axes[0].set_title("Raw impedance, all 53 frequencies pooled")
    axes[0].set_xlabel("Z_xy (V/A)")
    axes[0].set_ylabel("Count")
    axes[0].legend(frameon=False, fontsize=9)

    axes[1].hist(
        norm_real, bins=60, color="#1d4ed8", alpha=0.6, label="Re(Z)"
    )
    axes[1].hist(
        norm_imag, bins=60, color="#dc2626", alpha=0.6, label="Im(Z)"
    )
    axes[1].axvline(0.0, color="#374151", lw=1.0, ls="--")
    axes[1].set_title("ComplexZScore.transform output")
    axes[1].set_xlabel("Standardized channel value")
    axes[1].legend(frameon=False, fontsize=9)

    fig.suptitle(
        "WILLY L18, Z_xy: per-frequency z-scoring removes the "
        "cross-frequency scale spread", fontsize=12
    )
    fig.tight_layout()
    _save(fig, "data_contracts_normalization.png")


def make_forward_physics_halfspace_benchmark() -> None:
    dz0, growth, n = 25.0, 1.22, 34
    dz = dz0 * growth ** np.arange(n)
    z_edges = np.concatenate([[0.0], np.cumsum(dz)])
    x_edges = np.linspace(0, 240_000, 25)
    mesh = MaxwellMesh(x_edges, z_edges)
    receivers = ReceiverSet([[120_000.0, 0.0]], ["S00"])
    freqs = np.logspace(-1, 2, 12)

    problem = MaxwellProblem(
        mesh, np.full(mesh.shape, 0.01), freqs, receivers, ("zxy",)
    )
    result = MT2DAdapter().solve(problem)
    zxy = result.impedance_v_a[0, :, 0]
    analytic = half_space_impedance(100.0, freqs)

    mu0 = 4.0e-7 * np.pi
    rho_pred = np.abs(zxy) ** 2 / (2 * np.pi * freqs * mu0)
    rho_true = np.abs(analytic) ** 2 / (2 * np.pi * freqs * mu0)
    phase_pred = np.degrees(np.angle(zxy))
    phase_true = np.degrees(np.angle(analytic))

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.2))
    axes[0].plot(
        freqs, rho_true, "o-", color="#374151", label="analytic", lw=1.6
    )
    axes[0].plot(
        freqs, rho_pred, "x--", color="#dc2626", label="MT2DAdapter", lw=1.4
    )
    axes[0].set_xscale("log")
    axes[0].set_xlabel("Frequency (Hz)")
    axes[0].set_ylabel(r"Apparent resistivity ($\Omega\cdot$m)")
    axes[0].set_ylim(95, 105)
    axes[0].legend(frameon=False, fontsize=9)
    axes[0].grid(alpha=0.25, which="both")

    axes[1].plot(freqs, phase_true, "o-", color="#374151", lw=1.6)
    axes[1].plot(freqs, phase_pred, "x--", color="#dc2626", lw=1.4)
    axes[1].set_xscale("log")
    axes[1].set_xlabel("Frequency (Hz)")
    axes[1].set_ylabel("Phase (degrees)")
    axes[1].set_ylim(40, 50)
    axes[1].grid(alpha=0.25, which="both")

    fig.suptitle(
        "MT2DAdapter vs. the analytic 100 Ohm.m half-space limit "
        "(zxy, single receiver)",
        fontsize=12,
    )
    fig.tight_layout()
    _save(fig, "forward_physics_halfspace_benchmark.png")


def make_forward_physics_mesh_sensitivity() -> None:
    """Measure half-space error as lateral domain width changes."""
    dz = 25.0 * 1.22 ** np.arange(34)
    z_edges = np.concatenate([[0.0], np.cumsum(dz)])
    frequencies = np.logspace(-1, 2, 10)
    analytic = half_space_impedance(100.0, frequencies)
    widths_km = (20, 40, 80, 240)
    errors = []
    rho_curves = []
    mu0 = 4.0e-7 * np.pi

    for width_km in widths_km:
        width_m = width_km * 1000.0
        mesh = MaxwellMesh(np.linspace(0.0, width_m, 25), z_edges)
        receiver = ReceiverSet([[width_m / 2.0, 0.0]], ["S00"])
        problem = MaxwellProblem(
            mesh,
            np.full(mesh.shape, 0.01),
            frequencies,
            receiver,
            ("zxy",),
        )
        predicted = MT2DAdapter(verbose=False).solve(problem).impedance_v_a[0, :, 0]
        errors.append(np.abs(np.abs(predicted) - np.abs(analytic)) / np.abs(analytic))
        rho_curves.append(
            np.abs(predicted) ** 2 / (2.0 * np.pi * frequencies * mu0)
        )

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4))
    colors = plt.cm.viridis(np.linspace(0.08, 0.92, len(widths_km)))
    for width, error, rho, color in zip(widths_km, errors, rho_curves, colors):
        axes[0].plot(frequencies, 100.0 * error, "o-", color=color,
                     label=f"{width} km")
        axes[1].plot(frequencies, rho, "o-", color=color, label=f"{width} km")
    axes[0].axhline(5.0, color="#dc2626", ls="--", lw=1.2,
                    label="5% benchmark limit")
    axes[0].set(
        xscale="log", yscale="log", xlabel="Frequency (Hz)",
        ylabel="Impedance-amplitude error (%)",
        title="Boundary sensitivity is frequency dependent",
    )
    axes[1].axhline(100.0, color="#111827", ls="--", lw=1.2,
                    label="analytic half-space")
    axes[1].set(
        xscale="log", xlabel="Frequency (Hz)",
        ylabel=r"Apparent resistivity ($\Omega\cdot$m)",
        title="The same error biases the observable response",
    )
    for ax in axes:
        ax.grid(alpha=0.25, which="both")
        ax.legend(fontsize=8, ncol=2)
    fig.suptitle(
        "MT2D half-space verification across lateral mesh widths "
        "(fixed vertical mesh)", fontsize=12,
    )
    fig.tight_layout()
    _save(fig, "forward_physics_mesh_sensitivity.png")


def make_geology_prior_diagnostic() -> None:
    """Visualise correlation recovery and compositional geological priors."""
    from pycsamt.ai.geology import (
        ElectricalLayer,
        EllipsoidalLens,
        GaussianCorrelation,
        GeologyGrid,
        directional_variogram,
        generate_gaussian_field,
        generate_layered_geology,
        insert_lenses,
        interpolate_topography,
    )

    grid = GeologyGrid.regular_2d(nx=64, nz=32, dx_m=100, dz_m=50)
    correlation = GaussianCorrelation(length_x_m=700, length_z_m=150)
    example = generate_gaussian_field(grid, correlation, seed=12)

    variograms_x = []
    variograms_z = []
    for seed in range(32):
        realization = generate_gaussian_field(grid, correlation, seed=seed)
        variograms_x.append(
            directional_variogram(realization, "x", max_lag_cells=21).semivariance
        )
        variograms_z.append(
            directional_variogram(realization, "z", max_lag_cells=15).semivariance
        )
    variograms_x = np.asarray(variograms_x)
    variograms_z = np.asarray(variograms_z)
    lag_x = np.arange(1, 22) * 100.0
    lag_z = np.arange(1, 16) * 50.0

    layers = (
        ElectricalLayer(
            "weathered cover", 30.0, log10_std=0.10,
            heterogeneity=GaussianCorrelation(500, 120),
            resistivity_bounds_ohm_m=(8.0, 120.0),
        ),
        ElectricalLayer("sedimentary unit", 180.0),
        ElectricalLayer(
            "basement", 1200.0, log10_std=0.12,
            heterogeneity=GaussianCorrelation(900, 180),
            resistivity_bounds_ohm_m=(300.0, 4000.0),
        ),
    )
    layered = generate_layered_geology(
        grid,
        layers,
        [350.0, 950.0],
        seed=24,
        interface_relief_std_m=[60.0, 110.0],
        interface_correlation=GaussianCorrelation(850, 150),
        minimum_thickness_m=100.0,
    )
    lenses = (
        EllipsoidalLens(
            "conductor", 2100.0, 650.0, 650.0, 170.0, 6.0,
            dip_deg=12.0, transition_fraction=0.25,
        ),
        EllipsoidalLens(
            "resistor", 4700.0, 1150.0, 520.0, 210.0, 3500.0,
            dip_deg=-18.0, transition_fraction=0.20,
        ),
    )
    composed = insert_lenses(layered, lenses, conflict_policy="error")
    topo_x = np.linspace(0.0, 6400.0, 9)
    topo_elevation = 420.0 + 70.0 * np.sin(topo_x / 900.0) + 25.0 * np.cos(
        topo_x / 430.0
    )
    surface = interpolate_topography(
        grid, topo_x, topo_elevation, source="synthetic survey",
        interpolation_method="cubic",
    )

    fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.4))
    extent = [0.0, 6.4, 1.6, 0.0]
    image = axes[0, 0].imshow(
        example.values, extent=extent, aspect="auto", cmap="RdBu_r",
        vmin=-2.5, vmax=2.5,
    )
    axes[0, 0].set(
        title="One standardized anisotropic Gaussian field",
        xlabel="Profile distance (km)", ylabel="Depth (km)",
    )
    fig.colorbar(image, ax=axes[0, 0], label="Standard deviations")

    ax = axes[0, 1]
    for lag, values, length, label, color in (
        (lag_x, variograms_x, 700.0, "x direction", "#2563eb"),
        (lag_z, variograms_z, 150.0, "z direction", "#dc2626"),
    ):
        normalized_lag = lag / length
        median = np.median(values, axis=0)
        lower, upper = np.percentile(values, [10, 90], axis=0)
        ax.plot(normalized_lag, median, "o-", color=color, label=f"{label}: median")
        ax.fill_between(normalized_lag, lower, upper, color=color, alpha=0.16)
    theory_lag = np.linspace(0.0, 3.1, 180)
    ax.plot(
        theory_lag, 1.0 - np.exp(-0.5 * theory_lag**2),
        "k--", lw=1.5, label="requested Gaussian model",
    )
    ax.set(
        title="32-realization directional variogram audit",
        xlabel="Lag / requested correlation length", ylabel="Semivariance",
        xlim=(0, 3.1), ylim=(0, 1.65),
    )
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)

    for ax, values, title in (
        (axes[1, 0], layered.resistivity_ohm_m, "Correlated interfaces + within-unit variability"),
        (axes[1, 1], composed.resistivity_ohm_m, "Lenses + interpolated topographic mask"),
    ):
        display = np.log10(values)
        if ax is axes[1, 1]:
            display = np.where(surface.earth_mask(), display, np.nan)
        image = ax.imshow(
            display, extent=extent, aspect="auto", cmap="viridis_r",
            vmin=0.6, vmax=3.7,
        )
        for interface in layered.interface_depth_m:
            ax.plot(grid.x_m / 1000.0, interface / 1000.0,
                    color="white", lw=0.9, alpha=0.75)
        if ax is axes[1, 1]:
            ax.plot(
                grid.x_m / 1000.0, surface.surface_depth_m / 1000.0,
                color="#111827", lw=1.4, label="terrain",
            )
            ax.legend(fontsize=8, loc="lower left")
        ax.set(title=title, xlabel="Profile distance (km)", ylabel="Depth (km)")
        fig.colorbar(image, ax=ax, label=r"$\log_{10}\rho$ [$\Omega\cdot$m]")

    fig.suptitle(
        "A geological prior is audited statistically, then composed geometrically",
        fontsize=13,
    )
    fig.tight_layout()
    _save(fig, "geology_prior_diagnostic.png")


def make_geology_topographic_surface() -> None:
    """Exercise terrain interpolation, masks, local depth, and slope."""
    from pycsamt.ai.geology import GeologyGrid, interpolate_topography

    grid = GeologyGrid.regular_2d(nx=80, nz=36, dx_m=100, dz_m=35)
    sample_x = np.array([0, 700, 1450, 2300, 3200, 4100, 5050, 6100, 7000, 8000])
    sample_z = np.array([418, 455, 438, 510, 487, 552, 515, 575, 548, 590], dtype=float)
    surfaces = {
        method: interpolate_topography(
            grid,
            sample_x,
            sample_z,
            interpolation_method=method,
            source="surveyed benchmarks",
            station_names=tuple(f"T{i:02d}" for i in range(sample_x.size)),
        )
        for method in ("nearest", "linear", "cubic")
    }
    surface = surfaces["cubic"]

    fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.0))
    ax = axes[0, 0]
    for method, item in surfaces.items():
        ax.plot(grid.x_m / 1000, item.elevation_m, lw=1.8, label=method)
    ax.scatter(sample_x / 1000, sample_z, c="black", marker="v", zorder=4,
               label="input samples")
    ax.set(title="Interpolation is an explicit modelling choice",
           xlabel="Profile distance (km)", ylabel="Elevation (m)")
    ax.grid(alpha=0.25)
    ax.legend(ncol=2, fontsize=8)

    local_depth = surface.local_depth_m()
    image = axes[0, 1].imshow(
        local_depth,
        extent=(0, 8.0, 1.26, 0),
        aspect="auto",
        cmap="coolwarm",
        vmin=-160,
        vmax=1100,
    )
    axes[0, 1].contour(
        grid.x_m / 1000,
        grid.z_m / 1000,
        local_depth,
        levels=[0],
        colors="black",
        linewidths=1.4,
    )
    axes[0, 1].set(title="Signed depth below local terrain",
                   xlabel="Profile distance (km)", ylabel="Grid depth (km)")
    fig.colorbar(image, ax=axes[0, 1], label="local depth (m)")

    mask_image = axes[1, 0].imshow(
        surface.earth_mask(),
        extent=(0, 8.0, 1.26, 0),
        aspect="auto",
        cmap="Greys",
        vmin=0,
        vmax=1,
    )
    axes[1, 0].set(title="Physics-facing mask: white air, black earth",
                   xlabel="Profile distance (km)", ylabel="Grid depth (km)")
    fig.colorbar(mask_image, ax=axes[1, 0], ticks=[0, 1], label="earth mask")

    slope = surface.slope_degrees()
    axes[1, 1].plot(grid.x_m / 1000, slope, color="#dc2626", lw=1.8)
    axes[1, 1].fill_between(grid.x_m / 1000, 0, slope, color="#fecaca")
    axes[1, 1].set(title="Slope is derived on the raster grid",
                   xlabel="Profile distance (km)", ylabel="Slope (degrees)")
    axes[1, 1].grid(alpha=0.25)
    fig.suptitle(
        f"TopographicSurface audit: relief={surface.relief_m:.1f} m, "
        f"air fraction={surface.summary()['air_cell_fraction']:.3f}",
        fontsize=13,
    )
    fig.tight_layout()
    _save(fig, "geology_topographic_surface.png")


def make_geology_3d_composition() -> None:
    """Build and inspect a correlated layered 3-D volume with a body."""
    from pycsamt.ai.geology import (
        ElectricalLayer,
        EllipsoidalLens,
        GaussianCorrelation,
        GeologyGrid,
        generate_layered_geology,
        insert_lenses,
        interpolate_topography,
    )

    grid = GeologyGrid.regular_3d(
        nx=42, ny=30, nz=24, dx_m=150, dy_m=150, dz_m=60
    )
    correlation = GaussianCorrelation(900, 240, length_y_m=600, azimuth_deg=30)
    layers = (
        ElectricalLayer("cover", 35, log10_std=0.08, heterogeneity=correlation),
        ElectricalLayer("host", 420),
        ElectricalLayer("basement", 1800, log10_std=0.10, heterogeneity=correlation),
    )
    layered = generate_layered_geology(
        grid,
        layers,
        [360, 900],
        seed=31,
        interface_relief_std_m=[70, 120],
        interface_correlation=correlation,
        minimum_thickness_m=120,
    )
    body = EllipsoidalLens(
        "dipping conductor",
        center_x_m=3300,
        center_y_m=2250,
        center_z_m=690,
        radius_x_m=1150,
        radius_y_m=650,
        radius_z_m=260,
        resistivity_ohm_m=8,
        azimuth_deg=35,
        dip_deg=18,
        transition_fraction=0.22,
    )
    volume = insert_lenses(layered, [body], conflict_policy="error")
    samples_xy = np.array([
        [0, 0], [6300, 0], [0, 4500], [6300, 4500], [3150, 2250],
        [1575, 1125], [4725, 3375],
    ])
    samples_elevation = 510 + 0.018 * samples_xy[:, 0] - 0.012 * samples_xy[:, 1]
    topography = interpolate_topography(
        grid, samples_xy, samples_elevation, source="projected control points"
    )

    iy = grid.shape[1] // 2
    ix = grid.shape[2] // 2
    iz = int(np.argmin(np.abs(grid.z_m - body.center_z_m)))
    log_rho = np.log10(volume.resistivity_ohm_m)
    earth = topography.earth_mask()
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 4.8))
    views = (
        (np.where(earth[:, iy, :], log_rho[:, iy, :], np.nan),
         [0, 6.3, 1.44, 0], "Vertical x-z slice", "x (km)", "depth (km)"),
        (np.where(earth[:, :, ix], log_rho[:, :, ix], np.nan),
         [0, 4.5, 1.44, 0], "Vertical y-z slice", "y (km)", "depth (km)"),
        (log_rho[iz], [0, 6.3, 4.5, 0],
         f"Horizontal slice at {grid.z_m[iz]:.0f} m", "x (km)", "y (km)"),
    )
    for ax, (values, extent, title, xlabel, ylabel) in zip(axes, views):
        image = ax.imshow(values, extent=extent, aspect="auto", cmap="viridis_r",
                          vmin=0.7, vmax=3.4)
        ax.set(title=title, xlabel=xlabel, ylabel=ylabel)
    colorbar_axis = fig.add_axes([0.925, 0.19, 0.015, 0.60])
    fig.colorbar(
        image,
        cax=colorbar_axis,
        label=r"$\log_{10}\rho$ [$\Omega\cdot$m]",
    )
    fig.suptitle(
        "One immutable 3-D prior: correlated stratigraphy, oriented lens, and terrain mask",
        fontsize=12.5,
    )
    fig.subplots_adjust(left=0.06, right=0.90, bottom=0.13, top=0.82, wspace=0.32)
    _save(fig, "geology_3d_composition.png")


def make_inference_calibration_diagnostic() -> None:
    """Compare calibrated and overconfident predictive distributions."""
    from scipy.stats import norm

    from pycsamt.ai.validation import reliability_curve

    rng = np.random.default_rng(2026)
    n_station, n_layer = 120, 6
    depth = np.arange(1, n_layer + 1)
    mean = 1.55 + 0.27 * depth[None, :] + rng.normal(
        0.0, 0.08, size=(n_station, 1)
    )
    true_sigma = 0.07 + 0.035 * depth[None, :]
    true_sigma = np.broadcast_to(true_sigma, mean.shape)
    truth = mean + rng.normal(size=mean.shape) * true_sigma
    calibrated_std = true_sigma.copy()
    overconfident_std = 0.45 * true_sigma
    levels = np.linspace(0.10, 0.99, 18)
    calibrated = reliability_curve(
        truth, mean, calibrated_std, levels=levels, kind="l1"
    )
    overconfident = reliability_curve(
        truth, mean, overconfident_std, levels=levels, kind="l1"
    )

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.3))
    axes[0].plot(levels, levels, "k--", label="ideal")
    axes[0].plot(
        calibrated.levels, calibrated.coverage, "o-", label="calibrated scale"
    )
    axes[0].plot(
        overconfident.levels,
        overconfident.coverage,
        "s-",
        label="45% uncertainty scale",
    )
    axes[0].set(
        xlabel="Nominal coverage",
        ylabel="Empirical coverage",
        title="Reliability is measured against truth",
        xlim=(0, 1),
        ylim=(0, 1),
    )
    axes[0].grid(alpha=0.25)
    axes[0].legend(fontsize=8)

    station = 7
    z90 = norm.ppf(0.95)
    axes[1].plot(mean[station], depth, "ko-", label="predictive mean")
    axes[1].fill_betweenx(
        depth,
        mean[station] - z90 * calibrated_std[station],
        mean[station] + z90 * calibrated_std[station],
        alpha=0.25,
        label="calibrated 90% interval",
    )
    axes[1].scatter(truth[station], depth, c="#dc2626", zorder=3, label="truth")
    axes[1].invert_yaxis()
    axes[1].set(
        xlabel=r"$\log_{10}\rho$ [$\Omega\cdot$m]",
        ylabel="Layer index",
        title="One synthetic held-out station",
    )
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)

    labels = ["calibrated\nscale", "overconfident\nscale"]
    calibration_error = [
        calibrated.calibration.value,
        overconfident.calibration.value,
    ]
    sharpness = [calibrated.sharpness, overconfident.sharpness]
    x = np.arange(2)
    axes[2].bar(x - 0.18, calibration_error, 0.36, label="mean |coverage error|")
    axes[2].bar(x + 0.18, sharpness, 0.36, label="sharpness (mean std)")
    axes[2].set_xticks(x, labels)
    axes[2].set(
        ylabel="Diagnostic value",
        title="Sharper is not better when miscalibrated",
    )
    axes[2].grid(alpha=0.25, axis="y")
    axes[2].legend(fontsize=8)
    fig.suptitle(
        "Inference uncertainty audit on 720 held-out layer parameters",
        fontsize=13,
    )
    fig.tight_layout()
    _save(fig, "inference_calibration_diagnostic.png")


def make_losses_penalty_anatomy() -> None:
    """Visualize robust, probabilistic, and response-space loss behavior."""
    from pycsamt.ai.losses import (
        gaussian_nll_loss,
        model_huber_loss,
        model_l1_loss,
        model_l2_loss,
        response_residual_loss,
    )

    residual = np.linspace(-4.0, 4.0, 401)
    zeros = np.zeros_like(residual)
    l1 = np.array([model_l1_loss([r], [0]).value for r in residual])
    l2 = np.array([model_l2_loss([r], [0]).value for r in residual])
    huber = np.array([
        model_huber_loss([r], [0], delta=1.0).value for r in residual
    ])

    sigma = np.logspace(-2.2, 0.7, 260)
    nll = {}
    for error in (0.1, 0.5, 1.0):
        nll[error] = np.array([
            gaussian_nll_loss([error], [0.0], [np.log(s**2)]).value
            for s in sigma
        ])

    frequency = np.logspace(-1, 3, 18)
    observed = (1.0 + 0.5j) * np.sqrt(frequency / frequency.max())
    predicted = observed + (0.035 + 0.025j)
    errors = np.linspace(0.015, 0.12, frequency.size)
    raw_cell = np.abs(predicted - observed) ** 2
    normalized_cell = np.abs((predicted - observed) / errors) ** 2
    raw = response_residual_loss(predicted, observed, kind="l2")
    normalized = response_residual_loss(
        predicted, observed, errors=errors, kind="l2"
    )

    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.2))
    axes[0, 0].plot(residual, l1, label="L1")
    axes[0, 0].plot(residual, l2, label="L2")
    axes[0, 0].plot(residual, huber, label=r"Huber $\delta=1$")
    axes[0, 0].set(
        xlabel="Residual", ylabel="Per-cell penalty",
        title="Large residuals dominate L2",
    )
    axes[0, 0].set_ylim(0, 8)
    axes[0, 0].grid(alpha=0.25)
    axes[0, 0].legend()

    axes[0, 1].plot(residual, np.sign(residual), label="L1 influence")
    axes[0, 1].plot(residual, 2 * residual, label="L2 influence")
    axes[0, 1].plot(
        residual, np.clip(residual, -1, 1), label="Huber influence"
    )
    axes[0, 1].set(
        xlabel="Residual", ylabel="Derivative with respect to residual",
        title="Influence controls optimization pressure",
    )
    axes[0, 1].set_ylim(-4.5, 4.5)
    axes[0, 1].grid(alpha=0.25)
    axes[0, 1].legend()

    for error, values in nll.items():
        axes[1, 0].plot(sigma, values, label=f"|residual|={error}")
        axes[1, 0].axvline(error, color="grey", lw=0.7, alpha=0.4)
    axes[1, 0].set_xscale("log")
    axes[1, 0].set(
        xlabel="Predicted standard deviation",
        ylabel="Gaussian NLL",
        title="NLL penalizes both over- and under-confidence",
    )
    axes[1, 0].set_ylim(-2.0, 15.0)
    axes[1, 0].grid(alpha=0.25, which="both")
    axes[1, 0].legend(fontsize=8)

    axes[1, 1].plot(frequency, raw_cell, "o-", label="raw contribution")
    axes[1, 1].plot(
        frequency, normalized_cell, "s-", label="error-normalized contribution"
    )
    axes[1, 1].set_xscale("log")
    axes[1, 1].set_yscale("log")
    axes[1, 1].set(
        xlabel="Frequency (Hz)", ylabel="Per-cell L2 contribution",
        title=f"Same residual: mean raw={raw.value:.3g}, normalized={normalized.value:.2f}",
    )
    axes[1, 1].grid(alpha=0.25, which="both")
    axes[1, 1].legend(fontsize=8)
    fig.suptitle("Loss choice changes which errors drive an inversion", fontsize=13)
    fig.tight_layout()
    _save(fig, "losses_penalty_anatomy.png")


def make_losses_regularization_tradeoff() -> None:
    """Demonstrate data-fit versus spatial regularization on a known model."""
    from scipy.ndimage import gaussian_filter

    from pycsamt.ai.losses import (
        boundary_condition_loss,
        model_l2_loss,
        total_variation_loss,
    )

    rng = np.random.default_rng(44)
    nz, nx = 44, 76
    truth = np.full((nz, nx), 2.65)
    truth[:9] = 1.45
    truth[9:25] = 2.15
    z, x = np.mgrid[:nz, :nx]
    conductor = ((x - 38) / 15) ** 2 + ((z - 24) / 7) ** 2 <= 1
    truth[conductor] = 0.75
    noisy = truth + rng.normal(0, 0.30, truth.shape)
    air = z < (3 + 1.8 * np.sin(x / 8.0))
    noisy[air] = rng.normal(0.2, 0.35, air.sum())
    valid_earth = ~air

    scales = np.linspace(0.0, 4.0, 25)
    candidates = [noisy if scale == 0 else gaussian_filter(noisy, scale) for scale in scales]
    model_loss = np.array([
        model_l2_loss(item, truth, valid=valid_earth).value for item in candidates
    ])
    tv_loss = np.array([
        total_variation_loss(item, valid=valid_earth).value for item in candidates
    ])
    boundary_loss = np.array([
        boundary_condition_loss(item, boundary_mask=air, target=0.0).value
        for item in candidates
    ])
    lambda_tv, lambda_boundary = 0.35, 0.20
    objective = model_loss + lambda_tv * tv_loss + lambda_boundary * boundary_loss
    selected = int(np.argmin(objective))

    fig, axes = plt.subplots(2, 3, figsize=(13.2, 8.0))
    extent = [0, 7.6, 2.2, 0]
    for ax, values, title in (
        (axes[0, 0], truth, "Known blocky truth"),
        (axes[0, 1], noisy, "Noisy prediction"),
        (axes[0, 2], candidates[selected],
         f"Selected smoothing scale = {scales[selected]:.2f}"),
    ):
        image = ax.imshow(values, extent=extent, aspect="auto", cmap="turbo",
                          vmin=0.3, vmax=3.1)
        ax.contour(
            x / 10.0,
            z * (2.2 / nz),
            air,
            levels=[0.5],
            colors="white",
            linewidths=0.8,
        )
        ax.set(xlabel="Profile distance (km)", ylabel="Depth (km)", title=title)
    colorbar_axis = fig.add_axes([0.925, 0.57, 0.014, 0.29])
    fig.colorbar(image, cax=colorbar_axis,
                 label=r"$\log_{10}\rho$ [$\Omega\cdot$m]")

    axes[1, 0].plot(scales, model_loss, "o-", label=r"$L_{model}$")
    axes[1, 0].plot(scales, tv_loss, "s-", label=r"$L_{TV}$")
    axes[1, 0].plot(scales, boundary_loss, "^-", label=r"$L_{boundary}$")
    axes[1, 0].set(
        xlabel="Gaussian smoothing scale (cells)", ylabel="Mean loss",
        title="Terms prefer different models",
    )
    axes[1, 0].grid(alpha=0.25)
    axes[1, 0].legend(fontsize=8)

    axes[1, 1].plot(scales, objective, "o-", color="#7c3aed")
    axes[1, 1].axvline(scales[selected], color="#111827", ls="--")
    axes[1, 1].set(
        xlabel="Gaussian smoothing scale (cells)", ylabel="Combined objective",
        title=rf"$L_m+{lambda_tv}L_{{TV}}+{lambda_boundary}L_b$",
    )
    axes[1, 1].grid(alpha=0.25)

    scatter = axes[1, 2].scatter(
        tv_loss, model_loss, c=scales, cmap="viridis", s=55
    )
    axes[1, 2].scatter(tv_loss[selected], model_loss[selected], marker="*",
                       s=180, c="#dc2626", label="selected")
    axes[1, 2].set(
        xlabel=r"$L_{TV}$", ylabel=r"$L_{model}$",
        title="Regularity has a model-fit cost",
    )
    axes[1, 2].grid(alpha=0.25)
    axes[1, 2].legend(fontsize=8)
    fig.colorbar(scatter, ax=axes[1, 2], label="smoothing scale")
    fig.suptitle("A regularization weight selects a trade-off, not a true model", fontsize=13)
    fig.subplots_adjust(left=0.06, right=0.90, bottom=0.08, top=0.90,
                        hspace=0.34, wspace=0.31)
    _save(fig, "losses_regularization_tradeoff.png")


def make_hybrid_physics_refinement_audit() -> None:
    """Run the differentiable pseudo-2-D refinement on a controlled model."""
    from pycsamt.ai.inversion._pinn_ops import fit_2d_joint
    from pycsamt.forward.em1d import MT1DForward
    from pycsamt.forward.synthetic import LayeredModel

    frequencies = np.logspace(-1, 2, 18)
    n_station, n_layer = 7, 5
    station = np.linspace(-1.0, 1.0, n_station)
    depth_pattern = np.array([1.35, 2.05, 2.65, 3.05, 3.25])
    true_log_rho = np.tile(depth_pattern, (n_station, 1))
    true_log_rho[:, 1] -= 0.55 * np.exp(-(station / 0.48) ** 2)
    true_log_rho[:, 2] += 0.35 * station
    true_log_thick = np.tile(np.log10([120.0, 210.0, 360.0, 600.0]),
                             (n_station, 1))

    observed_log_rho = np.empty((n_station, frequencies.size))
    observed_phase = np.empty_like(observed_log_rho)
    for index in range(n_station):
        model = LayeredModel(
            resistivity=10.0 ** true_log_rho[index],
            thickness=10.0 ** true_log_thick[index],
        )
        response = MT1DForward(frequencies).run(model)
        observed_log_rho[index] = np.log10(response.rho_a)
        observed_phase[index] = response.phase

    initial_log_rho = true_log_rho + 0.30
    initial_log_rho[:, 1] += 0.28
    initial_log_rho[:, 3] -= 0.18 * station
    result = fit_2d_joint(
        observed_log_rho,
        observed_phase,
        frequencies,
        n_layers=n_layer,
        depth_max=1800.0,
        lam_z=0.002,
        lam_x=0.003,
        lr=0.025,
        epochs=140,
        device="cpu",
        log_every=0,
        init_log_rho=initial_log_rho,
        init_log_thick=true_log_thick,
        backend="torch",
    )
    refined = result["log_rho"]
    initial_rmse = np.sqrt(np.mean((initial_log_rho - true_log_rho) ** 2))
    refined_rmse = np.sqrt(np.mean((refined - true_log_rho) ** 2))

    def response_rms(log_rho, log_thick):
        residuals = []
        for index in range(n_station):
            model = LayeredModel(
                resistivity=10.0 ** log_rho[index],
                thickness=10.0 ** log_thick[index],
            )
            response = MT1DForward(frequencies).run(model)
            residuals.extend(
                (np.log10(response.rho_a) - observed_log_rho[index]) ** 2
                + ((response.phase - observed_phase[index]) / 90.0) ** 2
            )
        return float(np.sqrt(np.mean(residuals)))

    initial_response_rms = response_rms(initial_log_rho, true_log_thick)
    refined_response_rms = response_rms(refined, result["log_thick"])

    fig, axes = plt.subplots(1, 4, figsize=(13.4, 5.4))
    extent = [1, n_station, n_layer, 0]
    for ax, values, title in (
        (axes[0], true_log_rho.T, "Known layered truth"),
        (axes[1], initial_log_rho.T,
         f"Warm start\nmodel RMSE = {initial_rmse:.3f}\nresponse RMS = {initial_response_rms:.3f}"),
        (axes[2], refined.T,
         f"Physics-refined\nmodel RMSE = {refined_rmse:.3f}\nresponse RMS = {refined_response_rms:.3f}"),
    ):
        image = ax.imshow(values, extent=extent, aspect="auto", cmap="viridis",
                          vmin=0.8, vmax=3.6)
        ax.set(xlabel="Station", ylabel="Layer index", title=title)
        ax.set_xticks(range(1, n_station + 1))
    fig.colorbar(
        image,
        ax=axes[:3],
        orientation="horizontal",
        fraction=0.055,
        pad=0.22,
        label=r"$\log_{10}\rho$ [$\Omega\cdot$m]",
    )
    axes[3].plot(np.arange(1, 141), result["history"], color="#7c3aed", lw=1.8)
    axes[3].set_yscale("log")
    axes[3].set(
        xlabel="Adam iteration", ylabel="Total objective", title="Stage-2 convergence"
    )
    axes[3].grid(alpha=0.25, which="both")
    fig.suptitle(
        "Executed hybrid Stage 2: differentiable 1-D responses with lateral coupling",
        fontsize=12.5,
    )
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.27, top=0.72, wspace=0.38)
    _save(fig, "hybrid_physics_refinement_audit.png")


def make_dataset2d_realization_gallery() -> None:
    from pycsamt.ai.geology import GeologyGrid

    grid = GeologyGrid.regular_2d(nx=30, nz=18, dx_m=200, dz_m=120)
    station_x = list(np.linspace(500.0, 5500.0, 9))
    config = Maxwell2DDatasetConfig(
        dataset_id="gallery-demo",
        grid=grid,
        correlation_length_x_m=(600.0, 1500.0),
        correlation_length_z_m=(200.0, 500.0),
        frequencies_hz=[10.0, 3.0, 1.0],
        station_x_m=station_x,
        n_realizations=3,
        seed=4,
        log_resistivity_mean=2.0,
        log_resistivity_std=0.4,
        validation_fraction=0.0,
        test_fraction=0.0,
    )
    dataset = generate_2d_maxwell_dataset(config)

    nz, nx = dataset.samples[0].resistivity_ohm_m.shape
    dz, dx = grid.spacing_m
    x_edges = np.linspace(0, grid.x_m[-1] + dx / 2, nx + 1)
    z_edges = np.linspace(0, grid.z_m[-1] + dz / 2, nz + 1)

    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.0), sharey=True)
    vmin = min(np.log10(s.resistivity_ohm_m).min() for s in dataset.samples)
    vmax = max(np.log10(s.resistivity_ohm_m).max() for s in dataset.samples)
    mesh = None
    for ax, sample in zip(axes, dataset.samples):
        mesh = ax.pcolormesh(
            x_edges,
            z_edges,
            np.log10(sample.resistivity_ohm_m),
            cmap="viridis_r",
            vmin=vmin,
            vmax=vmax,
            shading="flat",
        )
        ax.scatter(
            station_x,
            np.zeros_like(station_x),
            marker="v",
            color="white",
            edgecolor="black",
            s=35,
            zorder=3,
        )
        ax.invert_yaxis()
        ax.set_title(sample.realization_id.split("-")[-1])
        ax.set_xlabel("x (m)")
    axes[0].set_ylabel("depth (m)")
    fig.colorbar(
        mesh, ax=axes, shrink=0.85, pad=0.02, label=r"$\log_{10}(\rho)$ [Ohm.m]"
    )
    fig.suptitle(
        "Three independent realizations from one Maxwell2DDatasetConfig "
        "(triangles: receiver x-positions)",
        fontsize=12,
    )
    _save(fig, "dataset2d_realization_gallery.png")


def make_dataset2d_tri_topo_gallery() -> None:
    """Triangular-mesh section: topography drape, filled cells, electrodes.

    Demonstrates M15's ``draw_tri_mesh`` on a real, topography-following
    triangular mesh and a real correlated resistivity field sampled onto
    its triangle centroids (``pycsamt.ai.training.dataset2d_tri``'s own
    nearest-mapping technique) -- not a forward-solved inversion result
    (``Mare2DEMAdapter`` needs a real compiled MARE2DEM binary, unverified
    in this environment; see its module docstring). The mesh is a real,
    graded Shewchuk-Triangle quality mesh -- fine near the electrodes,
    coarsening with depth -- built with
    :func:`pycsamt.forward.maxwell.tri_mesh_gen.build_graded_tri_mesh`'s
    own ``topo_x_m``/``topo_z_m`` support. An earlier version of this
    figure inlined its own duplicate copy of that same two-pass Triangle
    technique directly in this script (before topography support existed
    in the shared function); it now just calls the shared builder.
    """
    from pycsamt.ai.geology import (
        GaussianCorrelation,
        GeologyGrid,
        generate_gaussian_field,
    )
    from pycsamt.api.mesh import draw_tri_mesh
    from pycsamt.api.station import StationMarkerStyle
    from pycsamt.forward.maxwell.tri_mesh_gen import build_graded_tri_mesh

    x_min, x_max = 0.0, 2000.0
    station_x = np.linspace(200.0, 1800.0, 9)

    def topo(x: np.ndarray) -> np.ndarray:
        # z is positive-down, so a negative value is terrain *above* the
        # flat reference datum -- a ridge, not a depression.
        return -80.0 * np.exp(-(((x - 1100.0) / 500.0) ** 2)) - 20.0 * np.sin(
            x / 600.0
        )

    x_profile = np.linspace(x_min, x_max, 60)
    z_top = topo(x_profile)
    bottom = 700.0

    mesh = build_graded_tri_mesh(
        (x_min, x_max),
        (0.0, bottom),
        station_x,
        surface_cell_m=40.0,
        growth_rate=1.3,
        topo_x_m=x_profile,
        topo_z_m=z_top,
    )
    stations = np.column_stack(
        [station_x, np.interp(station_x, x_profile, z_top)]
    )

    z0 = float(z_top.min())
    nz = 20
    grid = GeologyGrid.regular_2d(
        nx=40,
        nz=nz,
        dx_m=(x_max - x_min) / 40,
        dz_m=(bottom - z0) / nz,
        x_origin_m=x_min,
        z_origin_m=z0,
    )
    field = generate_gaussian_field(
        grid, GaussianCorrelation(length_x_m=350, length_z_m=150), seed=3
    )
    centroids = mesh.triangle_centroids_m
    dx = grid.x_m[1] - grid.x_m[0]
    dz = grid.z_m[1] - grid.z_m[0]
    ix = np.clip(
        np.round((centroids[:, 0] - grid.x_m[0]) / dx).astype(int),
        0,
        len(grid.x_m) - 1,
    )
    iz = np.clip(
        np.round((centroids[:, 1] - grid.z_m[0]) / dz).astype(int),
        0,
        len(grid.z_m) - 1,
    )
    log_rho = 2.0 + 0.5 * field.values[iz, ix]

    fig, ax = plt.subplots(figsize=(9, 5))
    fill, _edges = draw_tri_mesh(ax, mesh, log_rho, preset="review", cmap="viridis_r")
    ax.scatter(
        stations[:, 0],
        stations[:, 1],
        label="Electrodes",
        **StationMarkerStyle().kwargs(),
    )
    ax.plot(x_profile, z_top, color="black", lw=1.2)
    ax.invert_yaxis()
    ax.set_xlabel("x (m)")
    ax.set_ylabel("Depth (m)")
    ax.legend(loc="lower right", frameon=True)
    fig.colorbar(fill, ax=ax, label=r"$\log_{10}(\rho)$ [Ohm.m]")
    ax.set_title(
        f"Triangular mesh, topography-draped ({mesh.n_triangles} cells) "
        "-- true resistivity model, not a solved inversion"
    )
    _save(fig, "dataset2d_tri_topo_gallery.png")


def make_dataset2d_tri_topo_ai_inversion() -> None:
    """Real, trained triangular AI inversion, draped over real topography.

    Companion to :func:`make_dataset2d_tri_topo_gallery`, which only ever
    paints the *true* geological model -- this one is a genuine
    ``Inv2DAgent(physics="mt2d_tri", topo_x_m=..., topo_z_m=...)`` run
    (the same synthetic ridge topography, for visual continuity),
    rendering the actual **predicted** per-triangle resistivity, not the
    target field. Building this required a real fix first: an earlier
    ``tri_fem2d.py`` hardcoded "surface" as literal ``z=0``, which
    rejected any receiver away from that datum -- exactly what a real
    topography-following station line needs. See that module's own
    docstring for the fix and its datum-shift-invariance benchmark.

    Real WILLY ``L18PLT`` stations supply the observed impedance
    features (``Inv2DAgent`` always reads its observed panel from a
    real ``Sites`` collection, even on this fully-synthetic-training
    path); the topography and training geology are synthetic, same
    teaching-scale honesty discipline as the Tongkeng CSAMT tutorial
    (real recovery RMSE/R2 reported plainly, no convergence claimed).
    """
    from pycsamt.agents import Inv2DAgent
    from pycsamt.api.mesh import draw_tri_mesh
    from pycsamt.api.station import StationMarkerStyle
    from pycsamt.forward.maxwell.tri_fem2d import TriFEM2DAdapter

    np.random.seed(7)
    try:
        import torch

        torch.manual_seed(7)
    except ImportError:
        pass

    sites = ensure_sites(
        PROJECT_ROOT / "data" / "AMT" / "WILLY_data" / "L18PLT",
        recursive=True,
        verbose=0,
    )
    n_sta = len(list(sites))
    x_min, x_max = 0.0, 2000.0
    station_spacing_m = (x_max - x_min) / (n_sta - 1)

    def topo(x: np.ndarray) -> np.ndarray:
        return -80.0 * np.exp(-(((x - 1100.0) / 500.0) ** 2)) - 20.0 * np.sin(
            x / 600.0
        )

    x_profile = np.linspace(x_min, x_max, 60)
    z_top = topo(x_profile)

    agent = Inv2DAgent(
        physics="mt2d_tri",
        epochs=15,
        n_freqs=4,
        depth_max=700.0,
        n_train_profiles=40,
        n_stations_per_profile=n_sta,
        station_spacing_m=station_spacing_m,
        mesh_target_cell_m=40.0,
        field_grid_cell_m=20.0,
        correlation_length_x_m=(300.0, 600.0),
        correlation_length_z_m=(80.0, 200.0),
        topo_x_m=x_profile,
        topo_z_m=z_top,
        mare2dem_adapter=TriFEM2DAdapter(),
    )
    result = agent.execute(
        {"sites": sites, "freqs": [100.0, 30.0, 10.0, 3.0]}
    )
    if result.status != "success":
        raise RuntimeError(f"make_dataset2d_tri_topo_ai_inversion: {result.summary}")

    pred = result.data["pred_triangles"]
    mesh = pred["mesh"]
    log_rho = pred["log10_resistivity"]
    recovery = result.data.get("mt2d_tri_recovery") or {}
    station_x = np.arange(n_sta, dtype=float) * station_spacing_m
    station_z = np.interp(station_x, x_profile, z_top)

    fig, ax = plt.subplots(figsize=(9, 5))
    fill, _edges = draw_tri_mesh(ax, mesh, log_rho, preset="review", cmap="viridis_r")
    ax.scatter(
        station_x, station_z, label="Electrodes", **StationMarkerStyle().kwargs()
    )
    ax.plot(x_profile, z_top, color="black", lw=1.2)
    ax.invert_yaxis()
    ax.set_xlabel("x (m)")
    ax.set_ylabel("Depth (m)")
    ax.legend(loc="lower right", frameon=True)
    fig.colorbar(fill, ax=ax, label=r"$\log_{10}(\rho)$ [Ohm.m]")
    rmse = recovery.get("rmse", float("nan"))
    r2 = recovery.get("r2", float("nan"))
    n_recovery = recovery.get("n_samples", 0)
    ax.set_title(
        f"Triangular AI inversion, topography-draped ({mesh.n_triangles} cells)\n"
        f"predicted log10(resistivity) -- recovery RMSE={rmse:.2f}, "
        f"R2={r2:.2f} (n={n_recovery}, teaching-scale, not converged)"
    )
    _save(fig, "dataset2d_tri_topo_ai_inversion.png")


def make_dataset2d_response_anatomy() -> None:
    """Plot one generated model beside its canonical TE/TM responses."""
    from pycsamt.ai.geology import GeologyGrid

    grid = GeologyGrid.regular_2d(nx=12, nz=8, dx_m=300, dz_m=150)
    station_x = np.linspace(450.0, 3150.0, 7)
    frequencies = np.array([30.0, 10.0, 3.0, 1.0])
    config = Maxwell2DDatasetConfig(
        dataset_id="response-anatomy",
        grid=grid,
        correlation_length_x_m=(600.0, 1200.0),
        correlation_length_z_m=(180.0, 450.0),
        frequencies_hz=frequencies,
        station_x_m=station_x,
        n_realizations=1,
        seed=11,
        log_resistivity_mean=2.0,
        log_resistivity_std=0.35,
        mesh_safety_factor=2.0,
        validation_fraction=0.0,
        test_fraction=0.0,
    )
    sample = generate_2d_maxwell_dataset(config).samples[0]
    impedance = sample.survey.impedance
    mu0 = 4.0e-7 * np.pi
    omega = 2.0 * np.pi * frequencies
    rho_a = np.abs(impedance) ** 2 / (omega[None, :, None] * mu0)
    phase = np.degrees(np.angle(impedance))

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.2))
    ax = axes[0, 0]
    dx = grid.spacing_m[1]
    dz = grid.spacing_m[0]
    image = ax.imshow(
        np.log10(sample.resistivity_ohm_m),
        extent=(0, grid.shape[1] * dx, grid.shape[0] * dz, 0),
        aspect="auto",
        cmap="viridis_r",
    )
    ax.scatter(station_x, np.zeros_like(station_x), marker="v", c="white",
               edgecolor="black", s=42, clip_on=False)
    ax.set(title="Known target model", xlabel="x (m)", ylabel="depth (m)")
    fig.colorbar(image, ax=ax, label=r"$\log_{10}(\rho)$ [$\Omega\cdot$m]")

    colors = plt.cm.plasma(np.linspace(0.08, 0.92, len(station_x)))
    for station, color, x_value in zip(range(len(station_x)), colors, station_x):
        axes[0, 1].plot(
            frequencies, rho_a[station, :, 0], "o-", color=color,
            label=f"x={x_value:.0f} m",
        )
        axes[1, 0].plot(
            frequencies, rho_a[station, :, 1], "o-", color=color
        )
    for ax, title in zip(
        (axes[0, 1], axes[1, 0]),
        (r"TE $Z_{xy}$ apparent resistivity", r"TM $Z_{yx}$ apparent resistivity"),
    ):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set(title=title, xlabel="frequency (Hz)", ylabel=r"$\rho_a$ [$\Omega\cdot$m]")
        ax.grid(alpha=0.25, which="both")
    axes[0, 1].legend(frameon=False, fontsize=7, ncol=2)

    axes[1, 1].plot(
        frequencies, np.median(phase[:, :, 0], axis=0), "o-",
        label=r"median TE $Z_{xy}$",
    )
    tm_phase_adjusted = np.degrees(np.angle(-impedance[:, :, 1]))
    axes[1, 1].plot(
        frequencies, np.median(tm_phase_adjusted, axis=0), "s--",
        label=r"median adjusted TM $-Z_{yx}$",
    )
    axes[1, 1].set_xscale("log")
    axes[1, 1].set(
        title="Median phase response",
        xlabel="frequency (Hz)",
        ylabel="phase magnitude (degrees)",
    )
    axes[1, 1].grid(alpha=0.25, which="both")
    axes[1, 1].legend(frameon=False)
    fig.suptitle(
        "One accepted 2-D realization: target and canonical response channels",
        fontsize=13,
    )
    fig.tight_layout()
    _save(fig, "dataset2d_response_anatomy.png")


def make_dataset3d_realization_gallery() -> None:
    """Three independent 3-D realizations, each shown as two orthogonal
    slices, to demonstrate genuine (not extruded-2-D) 3-D structure.
    """
    from pycsamt.ai.geology import GeologyGrid

    grid = GeologyGrid.regular_3d(nx=10, ny=10, nz=8, dx_m=200, dy_m=200, dz_m=120)
    config = Maxwell3DDatasetConfig(
        dataset_id="gallery-demo",
        grid=grid,
        correlation_length_x_m=(400.0, 800.0),
        correlation_length_y_m=(400.0, 800.0),
        correlation_length_z_m=(120.0, 280.0),
        frequencies_hz=[50.0, 20.0],
        station_xy_m=[[700.0, 700.0], [700.0, 1100.0], [1100.0, 700.0], [1100.0, 1100.0]],
        n_realizations=3,
        seed=7,
        log_resistivity_mean=2.0,
        log_resistivity_std=0.4,
        mesh_safety_factor=2.0,
        max_mesh_cells=10_000,
        validation_fraction=0.0,
        test_fraction=0.0,
    )
    dataset = generate_3d_maxwell_dataset(config)
    station_xy = config.station_xy_m

    nz, ny, nx = dataset.samples[0].resistivity_ohm_m.shape
    dz, dy, dx = grid.spacing_m
    x_edges = np.linspace(0, grid.x_m[-1] + dx / 2, nx + 1)
    y_edges = np.linspace(0, grid.y_m[-1] + dy / 2, ny + 1)
    z_edges = np.linspace(0, grid.z_m[-1] + dz / 2, nz + 1)
    z_shallow_idx = 1
    y_mid_idx = ny // 2

    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.6))
    vmin = min(np.log10(s.resistivity_ohm_m).min() for s in dataset.samples)
    vmax = max(np.log10(s.resistivity_ohm_m).max() for s in dataset.samples)
    mesh = None
    for col, sample in enumerate(dataset.samples):
        log_rho = np.log10(sample.resistivity_ohm_m)

        ax_top = axes[0, col]
        mesh = ax_top.pcolormesh(
            x_edges, y_edges, log_rho[z_shallow_idx],
            cmap="viridis_r", vmin=vmin, vmax=vmax, shading="flat",
        )
        ax_top.scatter(
            station_xy[:, 0], station_xy[:, 1], marker="v", color="white",
            edgecolor="black", s=45, zorder=3,
        )
        ax_top.set_title(
            f"{sample.realization_id.split('-')[-1]}  "
            f"(z={grid.z_m[z_shallow_idx]:.0f} m)"
        )
        ax_top.set_xlabel("x (m)")
        ax_top.set_aspect("equal")

        ax_bottom = axes[1, col]
        mesh = ax_bottom.pcolormesh(
            x_edges, z_edges, log_rho[:, y_mid_idx, :],
            cmap="viridis_r", vmin=vmin, vmax=vmax, shading="flat",
        )
        ax_bottom.scatter(
            station_xy[:, 0], np.zeros(station_xy.shape[0]), marker="v",
            color="white", edgecolor="black", s=45, zorder=3,
        )
        ax_bottom.invert_yaxis()
        ax_bottom.set_xlabel("x (m)")
    axes[0, 0].set_ylabel("y (m)")
    axes[1, 0].set_ylabel("depth (m)")
    fig.colorbar(
        mesh, ax=axes, shrink=0.85, pad=0.02, label=r"$\log_{10}(\rho)$ [Ohm.m]"
    )
    fig.suptitle(
        "Three independent realizations from one Maxwell3DDatasetConfig "
        f"(top: horizontal slice at z={grid.z_m[z_shallow_idx]:.0f} m; "
        f"bottom: vertical slice at y={grid.y_m[y_mid_idx]:.0f} m; "
        "triangles: receivers)",
        fontsize=11,
    )
    _save(fig, "dataset3d_realization_gallery.png")


def make_dataset3d_response_anatomy() -> None:
    """Plot one generated 3-D model beside its full-tensor response."""
    from pycsamt.ai.geology import GeologyGrid

    grid = GeologyGrid.regular_3d(nx=6, ny=6, nz=6, dx_m=250, dy_m=250, dz_m=150)
    station_xy = [
        [500.0, 500.0], [500.0, 1000.0], [1000.0, 500.0],
        [1000.0, 1000.0], [750.0, 750.0],
    ]
    frequencies = np.array([2.0, 1.0, 0.5, 0.2])
    config = Maxwell3DDatasetConfig(
        dataset_id="response-anatomy-3d",
        grid=grid,
        correlation_length_x_m=(500.0, 900.0),
        correlation_length_y_m=(500.0, 900.0),
        correlation_length_z_m=(150.0, 350.0),
        frequencies_hz=frequencies,
        station_xy_m=station_xy,
        n_realizations=1,
        seed=13,
        log_resistivity_mean=2.0,
        log_resistivity_std=0.35,
        components=("zxx", "zxy", "zyx", "zyy"),
        mesh_safety_factor=2.0,
        max_mesh_cells=10_000,
        validation_fraction=0.0,
        test_fraction=0.0,
    )
    sample = generate_3d_maxwell_dataset(config).samples[0]
    impedance = sample.survey.impedance
    mu0 = 4.0e-7 * np.pi
    omega = 2.0 * np.pi * frequencies
    rho_a = np.abs(impedance) ** 2 / (omega[None, :, None] * mu0)

    fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.6))
    ax = axes[0, 0]
    dz_, dy_, dx_ = grid.spacing_m
    z_shallow_idx = 1
    x_edges = np.linspace(0, grid.x_m[-1] + dx_ / 2, grid.shape[2] + 1)
    y_edges = np.linspace(0, grid.y_m[-1] + dy_ / 2, grid.shape[1] + 1)
    image = ax.pcolormesh(
        x_edges, y_edges, np.log10(sample.resistivity_ohm_m[z_shallow_idx]),
        cmap="viridis_r", shading="flat",
    )
    station_arr = np.asarray(station_xy)
    ax.scatter(
        station_arr[:, 0], station_arr[:, 1], marker="v", c="white",
        edgecolor="black", s=55, clip_on=False, zorder=3,
    )
    ax.set(
        title=f"Target model, z={grid.z_m[z_shallow_idx]:.0f} m",
        xlabel="x (m)", ylabel="y (m)",
    )
    ax.set_aspect("equal")
    fig.colorbar(image, ax=ax, label=r"$\log_{10}(\rho)$ [$\Omega\cdot$m]")

    colors = plt.cm.plasma(np.linspace(0.08, 0.92, station_arr.shape[0]))
    for s_idx, color in enumerate(colors):
        axes[0, 1].plot(
            frequencies, rho_a[s_idx, :, 1], "o-", color=color,
            label=f"S{s_idx}",
        )
        axes[1, 0].plot(frequencies, rho_a[s_idx, :, 2], "o-", color=color)
    for ax_, title in zip(
        (axes[0, 1], axes[1, 0]),
        (r"$Z_{xy}$ apparent resistivity", r"$Z_{yx}$ apparent resistivity"),
    ):
        ax_.set_xscale("log")
        ax_.set_yscale("log")
        ax_.set(title=title, xlabel="frequency (Hz)", ylabel=r"$\rho_a$ [$\Omega\cdot$m]")
        ax_.grid(alpha=0.25, which="both")
    axes[0, 1].legend(frameon=False, fontsize=7, ncol=2)

    median_zxx = np.median(np.abs(impedance[:, :, 0]), axis=0)
    median_zxy = np.median(np.abs(impedance[:, :, 1]), axis=0)
    median_zyx = np.median(np.abs(impedance[:, :, 2]), axis=0)
    median_zyy = np.median(np.abs(impedance[:, :, 3]), axis=0)
    axes[1, 1].plot(frequencies, median_zxy, "o-", label=r"median $|Z_{xy}|$")
    axes[1, 1].plot(frequencies, median_zyx, "s-", label=r"median $|Z_{yx}|$")
    axes[1, 1].plot(frequencies, median_zxx, "^--", label=r"median $|Z_{xx}|$")
    axes[1, 1].plot(frequencies, median_zyy, "d--", label=r"median $|Z_{yy}|$")
    axes[1, 1].set_xscale("log")
    axes[1, 1].set_yscale("log")
    axes[1, 1].set(
        title="Full tensor magnitude", xlabel="frequency (Hz)",
        ylabel=r"$|Z|$ [V/A]",
    )
    axes[1, 1].grid(alpha=0.25, which="both")
    axes[1, 1].legend(frameon=False, fontsize=8)
    fig.suptitle(
        "One accepted 3-D realization: target volume and full-tensor "
        "response channels",
        fontsize=13,
    )
    fig.tight_layout()
    _save(fig, "dataset3d_response_anatomy.png")


def make_data_preparation_contract() -> None:
    features, frequency, station_names = _willy_features()
    n_freq = frequency.size
    fig, axes = plt.subplots(1, 3, figsize=(12.8, 5.0), sharey=True)
    blocks = [
        (features[:, :n_freq], r"$\log_{10}(\rho_a)$", "viridis"),
        (features[:, n_freq:], "Phase (degrees)", "magma"),
        (
            (
                np.isfinite(features[:, :n_freq])
                & np.isfinite(features[:, n_freq:])
            ).astype(float),
            "Both blocks finite",
            "Greys",
        ),
    ]
    x = np.arange(len(station_names) + 1)
    frequency_edges = np.log10(
        np.r_[frequency, frequency[-1] * frequency[-1] / frequency[-2]]
    )
    for ax, (values, title, cmap) in zip(axes, blocks):
        mesh = ax.pcolormesh(
            x, frequency_edges, values.T, cmap=cmap, shading="auto"
        )
        ax.set_title(title)
        ax.set_xlabel("WILLY L18 station index")
        ax.xaxis.set_label_position("top")
        ax.xaxis.tick_top()
        ax.set_xticks(
            np.arange(0.5, len(station_names), 5),
            labels=np.arange(0, len(station_names), 5),
        )
        fig.colorbar(mesh, ax=ax, shrink=0.78, pad=0.02)
    axes[0].set_ylabel("Frequency (Hz)")
    axes[0].set_yticks(np.arange(0, 5), [f"$10^{tick}$" for tick in range(5)])
    fig.suptitle(
        "Field feature contract: 28 stations x 64 features", fontsize=13
    )
    fig.tight_layout()
    _save(fig, "data_preparation_contract.png")


def make_data_preparation_coverage() -> None:
    field, frequency, _ = _willy_features()
    synthetic = generate_dataset(
        solver="mt1d",
        n_samples=240,
        freqs=frequency,
        n_layers=5,
        rho_range=(1.0, 10_000.0),
        depth_max=2000.0,
        noise_level=0.05,
        noise_type="field",
        include_phase=True,
        seed=42,
        n_jobs=1,
        output=None,
        verbose=False,
    ).X
    n_freq = frequency.size
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.3))
    labels = [(r"$\log_{10}(\rho_a)$", 0), ("Phase (degrees)", n_freq)]
    for ax, (ylabel, start) in zip(axes, labels):
        stop = start + n_freq
        p01, p50, p99 = np.nanpercentile(
            synthetic[:, start:stop], [1, 50, 99], axis=0
        )
        field_median = np.nanmedian(field[:, start:stop], axis=0)
        ax.fill_between(
            frequency,
            p01,
            p99,
            color="#93c5fd",
            alpha=0.45,
            label="synthetic P1--P99",
        )
        ax.plot(
            frequency, p50, color="#1d4ed8", lw=1.7, label="synthetic median"
        )
        ax.plot(
            frequency,
            field_median,
            color="#dc2626",
            lw=2.0,
            label="WILLY median",
        )
        ax.set_xscale("log")
        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.22, which="both")
        ax.legend(frameon=False, fontsize=8)
    axes[0].set_title("Apparent-resistivity coverage")
    axes[1].set_title("Phase coverage")
    fig.suptitle("Generic synthetic prior versus WILLY L18", fontsize=13)
    fig.tight_layout()
    _save(fig, "data_preparation_coverage.png")


def make_model_selection_willy_dimension() -> None:
    line = PROJECT_ROOT / "data" / "AMT" / "WILLY_data" / "L18PLT"
    sites = ensure_sites(line, recursive=True, verbose=0)
    table = classify_dimensionality(sites)
    colors = {0: "#16a34a", 1: "#2563eb", 2: "#dc2626"}
    names = {0: "1-D", 1: "2-D", 2: "3-D"}

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.6))
    ax_feature, ax_period = axes
    for label in (2, 1, 0):
        selected = table["dim"].to_numpy() == label
        ax_feature.scatter(
            table.loc[selected, "beta_abs"],
            table.loc[selected, "ellipt_abs"],
            s=10,
            alpha=0.42,
            color=colors[label],
            label=f"{names[label]} (n={int(selected.sum())})",
        )
    ax_feature.axvline(3.0, color="#111827", ls="--", lw=1.1)
    ax_feature.axhline(0.2, color="#111827", ls="--", lw=1.1)
    ax_feature.set(xlabel=r"$|\beta|$ (degrees)", ylabel="Ellipticity")
    ax_feature.set_xlim(0, 90)
    ax_feature.set_title("Classification in phase-tensor space")
    ax_feature.grid(alpha=0.2)
    ax_feature.legend(frameon=False, fontsize=8)

    period = table["period"].to_numpy(dtype=float)
    edges = np.logspace(
        np.log10(np.nanmin(period)), np.log10(np.nanmax(period)), 9
    )
    centres = np.sqrt(edges[:-1] * edges[1:])
    fractions = np.zeros((3, centres.size))
    for index in range(centres.size):
        in_bin = (period >= edges[index]) & (period <= edges[index + 1])
        if np.any(in_bin):
            labels = table.loc[in_bin, "dim"].to_numpy(dtype=int)
            fractions[:, index] = [
                (labels == value).mean() for value in (0, 1, 2)
            ]
    ax_period.stackplot(
        centres,
        fractions,
        colors=[colors[0], colors[1], colors[2]],
        labels=[names[0], names[1], names[2]],
        alpha=0.88,
    )
    ax_period.set_xscale("log")
    ax_period.set_ylim(0, 1)
    ax_period.set(xlabel="Period (s)", ylabel="Fraction of classified samples")
    ax_period.set_title("Dimensionality occupancy across period")
    ax_period.grid(alpha=0.2, which="both")
    ax_period.legend(frameon=False, fontsize=8, loc="lower right")
    fig.suptitle(
        "WILLY L18: dimension must be justified before architecture",
        fontsize=13,
    )
    fig.tight_layout()
    _save(fig, "model_selection_willy_dimension.png")


def make_model_selection_graph_radius() -> None:
    """Show how radius changes graph topology and normalized GCN weights."""
    coords = np.array(
        [
            [0, 0], [310, 90], [690, -35], [1080, 115], [1530, 30],
            [2050, 160], [190, 760], [610, 680], [1040, 850],
            [1510, 710], [2160, 900], [2860, 390],
        ],
        dtype=float,
    )
    radii = (450.0, 900.0, 1300.0)
    fig = plt.figure(figsize=(13.2, 4.3))
    grid = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 0.045], wspace=0.18)
    axes = [fig.add_subplot(grid[0, i]) for i in range(3)]
    cax = fig.add_subplot(grid[0, 3])
    for ax, radius in zip(axes, radii):
        raw = build_adjacency(
            coords, radius, self_loops=False, normalise=False
        )
        normalized = build_adjacency(coords, radius)
        degree = raw.sum(axis=1).astype(int)
        for i in range(len(coords)):
            for j in range(i + 1, len(coords)):
                if raw[i, j] > 0:
                    ax.plot(
                        coords[[i, j], 0], coords[[i, j], 1],
                        color="#94a3b8", lw=0.8 + 2.8 * normalized[i, j],
                        zorder=1,
                    )
        points = ax.scatter(
            coords[:, 0], coords[:, 1], c=degree, cmap="viridis",
            vmin=0, vmax=max(1, int(degree.max())), s=80,
            edgecolor="white", linewidth=0.8, zorder=2,
        )
        isolated = int(np.sum(degree == 0))
        edges = int(raw.sum() // 2)
        ax.set_title(
            f"r = {radius:.0f} m\n{edges} edges; {isolated} isolated"
        )
        ax.set_aspect("equal", adjustable="box")
        ax.grid(alpha=0.18)
        ax.set_xlabel("Easting offset (m)")
    axes[0].set_ylabel("Northing offset (m)")
    colorbar = fig.colorbar(points, cax=cax)
    colorbar.set_label("Node degree (self-loop excluded)")
    fig.suptitle("Graph radius is a spatial prior, not a cosmetic setting")
    fig.subplots_adjust(left=0.06, right=0.94, bottom=0.16, top=0.78)
    _save(fig, "model_selection_graph_radius.png")


def make_model_selection_tradeoff() -> None:
    """Illustrate gates, Pareto trade-offs, and score-weight sensitivity."""
    names = np.array(["FCN", "CNN1D", "ResNet", "U-Net2D", "GCN"])
    metrics = np.array(
        [
            [1.00, 1.00, 1.00, 1.00],
            [0.91, 0.95, 0.93, 1.25],
            [0.76, 0.80, 0.78, 1.70],
            [0.69, 0.73, 1.31, 3.10],
            [0.72, 0.76, 0.89, 2.45],
        ]
    )
    gate_pass = np.array([True, True, True, False, True])
    rng = np.random.default_rng(42)
    seed_scores = np.clip(
        metrics[:, 1, None] * rng.normal(1.0, 0.045, (5, 12)), 0, None
    )

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.4))
    ax_heat, ax_seed, ax_sensitivity = axes
    image = ax_heat.imshow(metrics, cmap="YlGnBu_r", vmin=0.65, vmax=3.1)
    ax_heat.set_xticks(range(4), ["model", "response", "calibration", "runtime"], rotation=28, ha="right")
    ax_heat.set_yticks(range(5), names)
    for i in range(5):
        for j in range(4):
            ax_heat.text(j, i, f"{metrics[i, j]:.2f}", ha="center", va="center", fontsize=8)
    ax_heat.set_title("Baseline-normalized costs")
    fig.colorbar(image, ax=ax_heat, shrink=0.75, label="lower is better")

    parts = ax_seed.violinplot(seed_scores.T, showmeans=True, showextrema=False)
    for body in parts["bodies"]:
        body.set_facecolor("#60a5fa")
        body.set_edgecolor("#1d4ed8")
        body.set_alpha(0.68)
    parts["cmeans"].set_color("#111827")
    ax_seed.set_xticks(range(1, 6), names, rotation=28, ha="right")
    ax_seed.set_ylabel("Response NRMS / FCN mean")
    ax_seed.set_title("Twelve repeated-seed outcomes")
    ax_seed.grid(axis="y", alpha=0.22)

    response_weights = np.linspace(0.10, 0.70, 121)
    fixed_calibration, fixed_runtime = 0.20, 0.10
    for i, name in enumerate(names):
        model_weights = 1.0 - response_weights - fixed_calibration - fixed_runtime
        scores = (
            model_weights * metrics[i, 0]
            + response_weights * metrics[i, 1]
            + fixed_calibration * metrics[i, 2]
            + fixed_runtime * metrics[i, 3]
        )
        style = "--" if not gate_pass[i] else "-"
        ax_sensitivity.plot(response_weights, scores, ls=style, lw=1.8, label=name)
    ax_sensitivity.axvspan(0.10, 0.70, color="#f8fafc", zorder=-2)
    ax_sensitivity.set(xlabel="Response-error weight", ylabel="Composite score")
    ax_sensitivity.set_title("Scores depend on declared weights")
    ax_sensitivity.grid(alpha=0.22)
    ax_sensitivity.legend(frameon=False, fontsize=8, ncol=2)
    fig.suptitle("A controlled selection illustration: evidence before one score")
    fig.tight_layout()
    _save(fig, "model_selection_tradeoff.png")


def make_training_executed_audit() -> None:
    frequency = np.logspace(np.log10(1.01), 4, 24)
    samples = generate_dataset(
        solver="mt1d",
        n_samples=240,
        freqs=frequency,
        n_layers=5,
        rho_range=(1.0, 10_000.0),
        depth_max=2000.0,
        noise_level=0.05,
        noise_type="field",
        include_phase=True,
        seed=137,
        n_jobs=1,
        output=None,
        verbose=False,
    )
    train, validation, test = samples.split(
        val_frac=0.15, test_frac=0.15, seed=137
    )
    pool_x = np.vstack([train.X, validation.X])
    pool_y = np.vstack([train.y, validation.y])
    inverter = EMInverter1D(
        n_features=48,
        n_layers=5,
        arch="fcn",
        solver="mt1d",
        device="cpu",
        log_thickness=False,
        augment_noise=0.01,
    )
    inverter.fit(
        pool_x,
        pool_y,
        epochs=25,
        batch_size=64,
        lr=1e-3,
        patience=7,
        val_frac=0.15,
        grad_clip=1.0,
        seed=137,
        verbose=False,
    )
    prediction = inverter.predict(test.X)
    mae = np.mean(np.abs(prediction - test.y), axis=0)
    history = inverter._history
    epochs = np.arange(1, len(history["train_loss"]) + 1)

    fig, axes = plt.subplots(
        1, 3, figsize=(12.4, 4.2), gridspec_kw={"width_ratios": [1.35, 1, 1]}
    )
    ax_loss, ax_rho, ax_h = axes
    ax_loss.plot(
        epochs, history["train_loss"], color="#2563eb", label="training"
    )
    ax_loss.plot(
        epochs, history["val_loss"], color="#dc2626", label="validation"
    )
    best = int(inverter._meta["best_epoch"])
    ax_loss.axvline(
        best, color="#111827", ls="--", lw=1.1, label=f"restored epoch {best}"
    )
    ax_loss.set(xlabel="Epoch", ylabel="Normalized masked MSE")
    ax_loss.set_title("Optimization and early stopping")
    ax_loss.grid(alpha=0.22)
    ax_loss.legend(frameon=False, fontsize=8)

    ax_rho.bar(np.arange(1, 6), mae[:5], color="#0f766e")
    ax_rho.set(xlabel="Layer", ylabel=r"MAE in $\log_{10}(\rho)$")
    ax_rho.set_title("External-test resistivity error")
    ax_rho.set_xticks(np.arange(1, 6))
    ax_rho.grid(alpha=0.2, axis="y")

    ax_h.bar(np.arange(1, 5), mae[5:], color="#c2410c")
    ax_h.set(xlabel="Interface thickness", ylabel="MAE (m)")
    ax_h.set_title("External-test thickness error")
    ax_h.set_xticks(np.arange(1, 5))
    ax_h.grid(alpha=0.2, axis="y")
    fig.suptitle(
        "Executed 1-D training audit: small FCN smoke test", fontsize=13
    )
    fig.tight_layout()
    _save(fig, "training_executed_audit.png")


def make_training_augmentation_audit() -> None:
    """Execute the public augmenters on response-like feature curves."""
    from pycsamt.ai.training import (
        AugmentFreqDrop,
        AugmentMixup,
        AugmentNoise,
        AugmentStaticShift,
    )

    frequency = np.logspace(-1, 4, 24)
    log_frequency = np.log10(frequency)
    x = np.stack(
        [
            1.8 + 0.32 * np.tanh(log_frequency - centre)
            + 0.08 * np.sin(1.8 * log_frequency + phase)
            for centre, phase in zip(np.linspace(0.4, 2.0, 8), np.linspace(0, 2, 8))
        ]
    ).astype(np.float32)
    y = np.stack([np.linspace(1.2 + 0.08 * i, 3.1 - 0.04 * i, 5)
                  for i in range(len(x))]).astype(np.float32)

    noisy, _ = AugmentNoise(sigma=0.06)(x, y, rng=np.random.default_rng(11))
    shifted, _ = AugmentStaticShift(
        shift_range=(0.5, 2.0), n_amp_features=24
    )(x, y, rng=np.random.default_rng(12))
    dropped, _ = AugmentFreqDrop(
        drop_rate=0.25, contiguous=True, fill_value=np.nan
    )(x, y, rng=np.random.default_rng(13))
    mixed, y_mixed = AugmentMixup(alpha=0.4)(
        x, y, rng=np.random.default_rng(14)
    )

    fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.2), sharex=True)
    panels = [
        (noisy, "Additive feature noise", "local scatter; target unchanged"),
        (shifted, "Static shift", "whole amplitude curve translated"),
        (dropped, "Contiguous frequency drop", "dead band is explicit, not interpolated"),
        (mixed, "Mixup", "response and target move together"),
    ]
    for ax, (changed, title, subtitle) in zip(axes.flat, panels):
        ax.semilogx(frequency, x[0], color="#0f172a", lw=2.1, label="original")
        ax.semilogx(frequency, changed[0], color="#f15a29", lw=1.8,
                    marker="o", ms=3, label="augmented")
        ax.set(title=f"{title}\n{subtitle}", ylabel=r"feature $\log_{10}\rho_a$")
        ax.grid(alpha=0.22)
        ax.legend(frameon=False, fontsize=8)
    for ax in axes[-1]:
        ax.set_xlabel("Frequency (Hz)")
    axes[1, 1].text(
        0.03, 0.06,
        f"target change (sample 1): {np.linalg.norm(y_mixed[0] - y[0]):.3f}",
        transform=axes[1, 1].transAxes, fontsize=8.5,
        bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": "#94a3b8"},
    )
    fig.suptitle("Executed training augmentations: each encodes a different nuisance model")
    fig.tight_layout()
    _save(fig, "training_augmentation_audit.png")


def make_training_trainer_controls() -> None:
    """Run EMTrainer and expose validation control, LR, timing, and recovery."""
    import torch
    from torch import nn
    from torch.utils.data import TensorDataset

    from pycsamt.ai.training import EMTrainer

    torch.manual_seed(29)
    rng = np.random.default_rng(29)
    x = rng.normal(size=(480, 12)).astype(np.float32)
    weights = rng.normal(scale=0.35, size=(12, 4)).astype(np.float32)
    y = (x @ weights + 0.18 * np.sin(x[:, :4])).astype(np.float32)
    y += rng.normal(scale=0.08, size=y.shape).astype(np.float32)
    # A noisier validation acquisition creates an honest irreducible floor,
    # making the scheduler and best-weight restoration visible in a short run.
    y[360:] += rng.normal(scale=0.30, size=y[360:].shape).astype(np.float32)
    train = TensorDataset(torch.from_numpy(x[:360]), torch.from_numpy(y[:360]))
    validation = TensorDataset(torch.from_numpy(x[360:]), torch.from_numpy(y[360:]))
    model = nn.Sequential(nn.Linear(12, 24), nn.ReLU(), nn.Linear(24, 4))
    trainer = EMTrainer(
        model, lr=3e-3, weight_decay=1e-5, patience=30,
        min_delta=1e-5, batch_size=48, device="cpu", grad_clip=1.0,
        verbose=False,
    ).fit(train, validation, epochs=120)
    history = trainer.history
    epoch = np.arange(1, len(history["train_loss"]) + 1)
    with torch.no_grad():
        predicted = trainer.model(torch.from_numpy(x[360:])).numpy()
    per_target_rmse = np.sqrt(np.mean((predicted - y[360:]) ** 2, axis=0))

    fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.0))
    axes[0, 0].plot(epoch, history["train_loss"], color="#2563eb", label="training")
    axes[0, 0].plot(epoch, history["val_loss"], color="#dc2626", label="validation")
    axes[0, 0].axvline(trainer.best_epoch, color="#111827", ls="--",
                       label=f"restored epoch {trainer.best_epoch}")
    axes[0, 0].set(xlabel="Epoch", ylabel="Masked MSE", title="Validation selects weights")
    axes[0, 0].legend(frameon=False, fontsize=8)
    axes[0, 0].grid(alpha=0.22)
    axes[0, 1].step(epoch, history["lr"], where="post", color="#7c3aed")
    axes[0, 1].set(xlabel="Epoch", ylabel="Learning rate", title="Plateau scheduler state",
                   yscale="log")
    axes[0, 1].grid(alpha=0.22)
    axes[1, 0].plot(epoch, 1000 * np.asarray(history["epoch_time"]),
                    color="#0f766e", marker="o", ms=3)
    axes[1, 0].axhline(1000 * np.median(history["epoch_time"]), color="#f15a29",
                       ls="--", label="median")
    axes[1, 0].set(xlabel="Epoch", ylabel="CPU time (ms)", title="Runtime is part of the record")
    axes[1, 0].legend(frameon=False, fontsize=8)
    axes[1, 0].grid(alpha=0.22)
    axes[1, 1].bar(np.arange(1, 5), per_target_rmse, color="#60a5fa",
                   edgecolor="#1d4ed8")
    axes[1, 1].set(xlabel="Target parameter", ylabel="Validation RMSE",
                   title="Aggregate loss can hide target differences")
    axes[1, 1].set_xticks(np.arange(1, 5))
    axes[1, 1].grid(alpha=0.22, axis="y")
    fig.suptitle("Executed EMTrainer control audit")
    fig.tight_layout()
    _save(fig, "training_trainer_controls.png")
    print(
        "training controls:",
        {"epochs": len(epoch), "best_epoch": trainer.best_epoch,
         "best_val_loss": trainer.best_val_loss,
         "final_lr": history["lr"][-1],
         "target_rmse": per_target_rmse.tolist()},
    )


def make_hybrid_paired_diagnostic() -> None:
    """Render the executed Stage-1/Stage-2 audit used by hybrid.rst."""
    stage1_residual = np.array(
        [
            [1.4, -0.8, 0.5, 1.2],
            [0.9, -1.1, 0.7, 1.0],
            [1.6, -0.6, 0.4, 1.3],
        ]
    )
    stage2_residual = np.array(
        [
            [0.8, -0.5, 0.4, 0.7],
            [0.6, -0.7, 0.5, 0.6],
            [0.9, -0.4, 0.3, 0.8],
        ]
    )
    stage1_model = np.array(
        [
            [2.0, 2.5, 3.0, 1.7],
            [2.1, 2.6, 2.9, 1.8],
            [2.2, 2.4, 3.1, 1.7],
        ]
    )
    stage2_model = np.array(
        [
            [2.1, 2.45, 2.95, 1.72],
            [2.15, 2.55, 2.95, 1.82],
            [2.35, 2.35, 3.0, 1.75],
        ]
    )
    station_gain = np.sqrt(np.mean(stage1_residual**2, axis=1)) - np.sqrt(
        np.mean(stage2_residual**2, axis=1)
    )
    model_move = np.sqrt(np.mean((stage2_model - stage1_model) ** 2, axis=1))

    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.1))
    for ax, values, title in zip(
        axes[:2],
        (stage1_residual, stage2_residual),
        ("Stage 1 normalized residual", "Stage 2 normalized residual"),
    ):
        image = ax.imshow(
            values.T, cmap="RdBu_r", vmin=-1.6, vmax=1.6, aspect="auto"
        )
        ax.set(xlabel="Station", ylabel="Frequency-bin index", title=title)
        ax.xaxis.set_label_position("top")
        ax.xaxis.tick_top()
        ax.set_xticks(range(3), labels=["S1", "S2", "S3"])
        ax.set_yticks(range(4), labels=["F1", "F2", "F3", "F4"])
    color_axis = fig.add_axes([0.10, 0.07, 0.53, 0.025])
    fig.colorbar(
        image,
        cax=color_axis,
        orientation="horizontal",
        label="Residual / adopted error",
    )

    axes[2].scatter(model_move, station_gain, s=65, color="#0f766e")
    for index, (move, gain) in enumerate(
        zip(model_move, station_gain), start=1
    ):
        axes[2].annotate(
            f"S{index}",
            (move, gain),
            xytext=(5, 4),
            textcoords="offset points",
        )
    axes[2].axhline(0, color="#111827", lw=1)
    axes[2].set(
        xlabel=r"Model move $D_i$ in log-parameter space",
        ylabel="NRMS gain",
        title="Fit gain must justify model change",
    )
    axes[2].grid(alpha=0.22)
    fig.suptitle("Executed hybrid paired diagnostic", fontsize=13)
    fig.subplots_adjust(
        left=0.06, right=0.97, bottom=0.24, top=0.78, wspace=0.42
    )
    _save(fig, "hybrid_paired_diagnostic.png")


def make_uncertainty_coverage_reliability() -> None:
    """Independent-test split-conformal reliability method check."""
    rng = np.random.default_rng(31415)
    n_calibration, n_test, n_parameters = 1200, 4000, 9
    scales = np.r_[np.linspace(0.08, 0.16, 5), np.linspace(35, 90, 4)]
    calibration_error = (
        rng.standard_t(6, (n_calibration, n_parameters)) * scales
    )
    test_error = rng.standard_t(6, (n_test, n_parameters)) * scales
    groups = {
        "resistivity dims (0–4)": np.arange(5),
        "thickness dims (5–8)": np.arange(5, 9),
        "joint (all 9 dims)": np.arange(9),
    }
    nominal = np.array([0.50, 0.70, 0.80, 0.90, 0.95])
    colors = ["#2563eb", "#b91c1c", "#4b5563"]

    fig, ax = plt.subplots(figsize=(7.0, 5.7))
    ax.plot([0, 1], [0, 1], "--", color="#111827", label="nominal = actual")
    for (label, indices), color in zip(groups.items(), colors):
        calibration_score = np.max(
            np.abs(calibration_error[:, indices]) / scales[indices], axis=1
        )
        test_score = np.max(
            np.abs(test_error[:, indices]) / scales[indices], axis=1
        )
        actual = []
        for target in nominal:
            rank = int(np.ceil((n_calibration + 1) * target))
            rank = min(max(rank, 1), n_calibration)
            threshold = np.partition(calibration_score, rank - 1)[rank - 1]
            actual.append(np.mean(test_score <= threshold))
        ax.plot(nominal, actual, marker="o", lw=2, color=color, label=label)
    ax.set(
        xlabel=r"Nominal simultaneous coverage $(1-\alpha)$",
        ylabel="Empirical coverage on independent test set",
        title="Split-conformal reliability by parameter group",
        xlim=(0.45, 1.0),
        ylim=(0.45, 1.0),
    )
    ax.grid(alpha=0.22)
    ax.legend(frameon=True, fontsize=8, loc="lower right")
    fig.tight_layout()
    _save(fig, "uncertainty_coverage_reliability.png")


def make_uncertainty_calibration_regimes() -> None:
    """Execute reliability_curve for sharp, calibrated, and broad scales."""
    from pycsamt.ai.validation import reliability_curve

    rng = np.random.default_rng(8128)
    truth = rng.normal(size=(5000, 5))
    mean = np.zeros_like(truth)
    levels = np.array([0.50, 0.68, 0.80, 0.90, 0.95, 0.99])
    regimes = {
        "too narrow (sigma=0.55)": np.full_like(truth, 0.55),
        "calibrated (sigma=1.00)": np.full_like(truth, 1.00),
        "too broad (sigma=1.80)": np.full_like(truth, 1.80),
    }
    colors = ["#dc2626", "#16a34a", "#2563eb"]
    reports = {
        name: reliability_curve(truth, mean, std, levels=levels, kind="l1")
        for name, std in regimes.items()
    }

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.3))
    ax_rel, ax_trade, ax_res = axes
    ax_rel.plot([0, 1], [0, 1], "--", color="#111827", label="ideal")
    for (name, report), color in zip(reports.items(), colors):
        ax_rel.plot(report.levels, report.coverage, "o-", color=color, label=name)
    ax_rel.set(xlabel="Nominal coverage", ylabel="Empirical coverage",
               title="Reliability distinguishes scale errors", xlim=(0.45, 1.0),
               ylim=(0.35, 1.02))
    ax_rel.legend(frameon=False, fontsize=8)
    ax_rel.grid(alpha=0.22)

    names = list(reports)
    calibration_error = [reports[name].calibration.value for name in names]
    sharpness = [reports[name].sharpness for name in names]
    xpos = np.arange(3)
    width = 0.36
    ax_trade.bar(xpos - width / 2, calibration_error, width, color="#f15a29",
                 label="mean |coverage - nominal|")
    ax_trade.bar(xpos + width / 2, sharpness, width, color="#60a5fa",
                 label="sharpness (mean sigma)")
    ax_trade.set_xticks(xpos, ["narrow", "calibrated", "broad"], rotation=18)
    ax_trade.set(ylabel="Metric value", title="Sharpness alone rewards overconfidence")
    ax_trade.legend(frameon=False, fontsize=8)
    ax_trade.grid(alpha=0.22, axis="y")

    absolute_error = np.abs(truth - mean).ravel()
    for (name, std), color in zip(regimes.items(), colors):
        z = np.sort(absolute_error / std.ravel())
        probability = np.arange(1, len(z) + 1) / len(z)
        ax_res.plot(z, probability, color=color, label=name)
    ax_res.axvline(1.96, color="#111827", ls="--", label="Gaussian 95% z=1.96")
    ax_res.set(xlabel=r"Absolute standardized residual $|y-\mu|/\sigma$",
               ylabel="Empirical cumulative fraction",
               title="Why the same errors imply different coverage", xlim=(0, 4))
    ax_res.legend(frameon=False, fontsize=7.5)
    ax_res.grid(alpha=0.22)
    fig.suptitle("Executed Gaussian uncertainty audit: calibration and sharpness are paired")
    fig.tight_layout()
    _save(fig, "uncertainty_calibration_regimes.png")
    print(
        "uncertainty regimes:",
        {name: {"mace": round(report.calibration.value, 4),
                "sharpness": round(report.sharpness, 2),
                "coverage_95": round(float(report.coverage[-2]), 4)}
         for name, report in reports.items()},
    )


def make_uncertainty_depth_propagation() -> None:
    """Propagate correlated thickness draws before taking depth quantiles."""
    rng = np.random.default_rng(2718)
    mean_log_h = np.log10([110.0, 180.0, 290.0, 430.0])
    sigma = np.array([0.10, 0.12, 0.14, 0.16])
    correlation = np.array(
        [[1.00, -0.55, 0.20, 0.00],
         [-0.55, 1.00, -0.45, 0.15],
         [0.20, -0.45, 1.00, -0.35],
         [0.00, 0.15, -0.35, 1.00]]
    )
    covariance = correlation * np.outer(sigma, sigma)
    log_h = rng.multivariate_normal(mean_log_h, covariance, size=6000)
    thickness = 10.0 ** log_h
    interface_depth = np.cumsum(thickness, axis=1)
    probability = np.array([0.05, 0.50, 0.95])
    correct = np.quantile(interface_depth, probability, axis=0)
    wrong = np.cumsum(np.quantile(thickness, probability, axis=0), axis=1)

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.4))
    image = axes[0].imshow(correlation, vmin=-1, vmax=1, cmap="coolwarm")
    axes[0].set_xticks(range(4), ["h1", "h2", "h3", "h4"])
    axes[0].set_yticks(range(4), ["h1", "h2", "h3", "h4"])
    axes[0].set_title("Declared log-thickness correlation")
    fig.colorbar(image, ax=axes[0], label="correlation")

    interface = np.arange(1, 5)
    axes[1].fill_between(interface, correct[0], correct[2], color="#93c5fd",
                         alpha=0.65, label="draws -> cumsum -> quantiles")
    axes[1].plot(interface, correct[1], "o-", color="#1d4ed8", label="median")
    axes[1].plot(interface, wrong[0], "--", color="#dc2626")
    axes[1].plot(interface, wrong[2], "--", color="#dc2626",
                 label="quantiles -> cumsum (wrong)")
    axes[1].set(xlabel="Interface", ylabel="Cumulative depth (m)",
                title="Nonlinear propagation changes the band")
    axes[1].set_xticks(interface)
    axes[1].legend(frameon=False, fontsize=7.5)
    axes[1].grid(alpha=0.22)

    axes[2].hist(interface_depth[:, -1], bins=45, density=True,
                 color="#bfdbfe", edgecolor="white")
    for value, label, color in zip(correct[:, -1], ["5%", "median", "95%"],
                                   ["#2563eb", "#111827", "#2563eb"]):
        axes[2].axvline(value, color=color, ls="--" if label != "median" else "-",
                        label=f"{label}: {value:.0f} m")
    axes[2].set(xlabel="Fourth-interface depth (m)", ylabel="Density",
                title="Preserve the propagated draw distribution")
    axes[2].legend(frameon=False, fontsize=8)
    axes[2].grid(alpha=0.18)
    fig.suptitle("Executed propagation of correlated layer-thickness uncertainty")
    fig.tight_layout()
    _save(fig, "uncertainty_depth_propagation.png")
    print(
        "depth propagation:",
        {"correct_final_m": np.round(correct[:, -1], 1).tolist(),
         "wrong_final_m": np.round(wrong[:, -1], 1).tolist()},
    )


def make_agents_execution_contract() -> None:
    """Show the result-state gate around an agent execution."""
    fig, ax = plt.subplots(figsize=(12.0, 4.0))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    nodes = [
        (
            0.08,
            0.56,
            "Reviewed input",
            "QC · order · components",
            "#dbeafe",
            "#2563eb",
        ),
        (
            0.31,
            0.56,
            "Agent execution",
            "load · train/load · predict",
            "#ede9fe",
            "#7c3aed",
        ),
        (0.55, 0.78, "success", "workflow returned", "#dcfce7", "#16a34a"),
        (
            0.55,
            0.50,
            "needs_review",
            "warning or optional failure",
            "#ffedd5",
            "#ea580c",
        ),
        (0.55, 0.22, "failed", "required step failed", "#fee2e2", "#dc2626"),
        (
            0.82,
            0.64,
            "Scientific gate",
            "residuals · uncertainty · baseline",
            "#f1f5f9",
            "#475569",
        ),
        (
            0.82,
            0.28,
            "Stop and repair",
            "inspect error + fix hint",
            "#fef2f2",
            "#991b1b",
        ),
    ]
    for x, y, title, body, face, edge in nodes:
        box = FancyBboxPatch(
            (x - 0.09, y - 0.105),
            0.18,
            0.21,
            boxstyle="round,pad=0.012,rounding_size=0.018",
            facecolor=face,
            edgecolor=edge,
            linewidth=1.4,
        )
        ax.add_patch(box)
        ax.text(
            x,
            y + 0.035,
            title,
            ha="center",
            va="center",
            fontsize=10,
            fontweight="bold",
            color=edge,
        )
        ax.text(x, y - 0.035, body, ha="center", va="center", fontsize=8)
    arrows = [
        ((0.17, 0.56), (0.22, 0.56)),
        ((0.40, 0.58), (0.46, 0.76)),
        ((0.40, 0.56), (0.46, 0.50)),
        ((0.40, 0.53), (0.46, 0.24)),
        ((0.64, 0.76), (0.72, 0.66)),
        ((0.64, 0.50), (0.72, 0.62)),
        ((0.64, 0.22), (0.72, 0.28)),
    ]
    for start, end in arrows:
        ax.add_patch(
            FancyArrowPatch(
                start,
                end,
                arrowstyle="-|>",
                mutation_scale=12,
                lw=1.2,
                color="#64748b",
            )
        )
    ax.text(
        0.82,
        0.92,
        "Agent status is a software outcome, not acceptance",
        ha="center",
        fontsize=11,
        fontweight="bold",
    )
    fig.tight_layout()
    _save(fig, "agents_execution_contract.png")


def make_agents_executed_1d_audit() -> None:
    """Run a small WILLY AIInversionAgent example and audit its outputs."""
    try:
        import torch

        torch.manual_seed(42)
    except ImportError:
        pass
    np.random.seed(42)
    result = AIInversionAgent(
        arch="fcn",
        n_layers=5,
        n_train_samples=240,
        epochs=12,
    ).execute(
        {"path": PROJECT_ROOT / "data" / "AMT" / "WILLY_data" / "L18PLT"}
    )
    if result.status == "failed":
        raise RuntimeError(result.error)
    history = result["train_history"]
    predictions = result["predictions"]
    station_names = list(predictions)
    section = np.column_stack([predictions[name] for name in station_names])
    rms = np.array([result["rms_per_station"][name] for name in station_names])

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(14.0, 4.5),
        gridspec_kw={"width_ratios": [1.0, 1.35, 1.2]},
    )
    epoch = np.arange(1, len(history["train_loss"]) + 1)
    axes[0].plot(
        epoch, history["train_loss"], label="training", color="#2563eb"
    )
    axes[0].plot(
        epoch, history["val_loss"], label="validation", color="#dc2626"
    )
    axes[0].set(
        xlabel="Epoch",
        ylabel="Normalized loss",
        title="Synthetic training history",
    )
    axes[0].grid(alpha=0.22)
    axes[0].legend(frameon=False, fontsize=8)

    image = axes[1].imshow(
        section, aspect="auto", origin="upper", cmap="turbo"
    )
    axes[1].set(
        xlabel="Station index",
        ylabel="Layer index",
        title="Executed five-layer prediction",
    )
    axes[1].xaxis.set_label_position("top")
    axes[1].xaxis.tick_top()
    axes[1].set_xticks(
        np.arange(0, len(station_names), 5),
        labels=np.arange(0, len(station_names), 5),
    )
    fig.colorbar(
        image, ax=axes[1], shrink=0.78, label=r"$\log_{10}\rho$ ($\Omega$ m)"
    )

    order = np.argsort(rms)[::-1]
    axes[2].bar(np.arange(len(rms)), rms[order], color="#0f766e")
    axes[2].axhline(
        result["rms_global"],
        color="#b91c1c",
        ls="--",
        lw=1.4,
        label=f"station mean = {result['rms_global']:.3f}",
    )
    axes[2].set(
        xlabel="Stations ranked worst to best",
        ylabel=r"Log$_{10}$ apparent-resistivity RMS",
        title="Response fit is not spatially uniform",
    )
    axes[2].grid(alpha=0.22, axis="y")
    axes[2].legend(frameon=False, fontsize=8)
    fig.suptitle(
        "Executed WILLY AIInversionAgent smoke audit: success ≠ acceptance",
        fontsize=13,
    )
    fig.tight_layout()
    _save(fig, "agents_executed_1d_audit.png")


def make_agents_real_data_sections(names: set[str] | None = None) -> None:
    """Regenerate legacy-named sections from current WILLY agents."""
    line = PROJECT_ROOT / "data" / "AMT" / "WILLY_data" / "L18PLT"
    sites = ensure_sites(line, recursive=True, verbose=0)
    runs = (
        (
            "agents_ai1d_section.png",
            AIInversionAgent(
                arch="resnet", n_layers=5,
                n_train_samples=10_000, epochs=100,
            ),
            "ai_section",
        ),
        (
            "agents_ensemble_section.png",
            EnsembleAgent(
                n_estimators=5, arch="resnet", n_layers=5,
                n_train_samples=2_000, epochs=30, calibrate=True,
            ),
            "uncertainty_section",
        ),
    )
    for filename, agent, figure_key in runs:
        if names is not None and filename not in names:
            continue
        result = agent.execute({"sites": sites})
        if result.status == "failed":
            raise RuntimeError(f"{filename}: {result.error}")
        figure = result["figures"].get(figure_key)
        if figure is None:
            raise RuntimeError(f"{filename}: missing figure {figure_key!r}")
        _save(figure, filename)
        print(filename, result.summary)


def make_agents_architecture_comparison() -> None:
    """Execute matched CNN1D, ResNet, and FCN agent smoke runs."""
    line = PROJECT_ROOT / "data" / "AMT" / "WILLY_data" / "L18PLT"
    architectures = ("cnn1d", "resnet", "fcn")
    results = []
    for index, architecture in enumerate(architectures):
        seed = 510 + index
        np.random.seed(seed)
        try:
            import torch

            torch.manual_seed(seed)
        except ImportError:
            pass
        result = AIInversionAgent(
            arch=architecture,
            n_layers=5,
            n_train_samples=240,
            epochs=12,
        ).execute({"path": line})
        if result.status == "failed":
            raise RuntimeError(f"{architecture}: {result.error}")
        results.append(result)

    sections = [
        np.column_stack(
            [result["predictions"][name] for name in result["predictions"]]
        )
        for result in results
    ]
    finite = np.concatenate(
        [section[np.isfinite(section)] for section in sections]
    )
    color_min, color_max = np.percentile(finite, [2, 98])
    fig, axes = plt.subplots(
        3,
        3,
        figsize=(14.3, 11.0),
        gridspec_kw={"width_ratios": [1.0, 1.25, 1.1]},
    )
    last_image = None
    for row, (architecture, result, section) in enumerate(
        zip(architectures, results, sections)
    ):
        history = result["train_history"]
        epoch = np.arange(1, len(history["train_loss"]) + 1)
        axes[row, 0].plot(
            epoch, history["train_loss"], color="#2563eb", label="training"
        )
        axes[row, 0].plot(
            epoch, history["val_loss"], color="#dc2626", label="validation"
        )
        axes[row, 0].set(
            xlabel="Epoch",
            ylabel="Normalized loss",
            title=f"{architecture.upper()} optimization",
        )
        axes[row, 0].grid(alpha=0.22)
        axes[row, 0].legend(frameon=False, fontsize=7)

        last_image = axes[row, 1].imshow(
            section,
            aspect="auto",
            origin="upper",
            cmap="turbo",
            vmin=color_min,
            vmax=color_max,
        )
        axes[row, 1].set(
            xlabel="Station index",
            ylabel="Layer index",
            title=f"{architecture.upper()} five-layer model",
        )
        axes[row, 1].xaxis.set_label_position("top")
        axes[row, 1].xaxis.tick_top()
        axes[row, 1].set_xticks(
            np.arange(0, section.shape[1], 5),
            labels=np.arange(0, section.shape[1], 5),
        )

        station_names = list(result["predictions"])
        rms = np.array(
            [result["rms_per_station"][name] for name in station_names]
        )
        axes[row, 2].bar(np.arange(len(rms)), rms, color="#0f766e")
        axes[row, 2].axhline(
            result["rms_global"],
            color="#b91c1c",
            ls="--",
            lw=1.2,
            label=f"mean={result['rms_global']:.3f}",
        )
        axes[row, 2].set(
            xlabel="Station index",
            ylabel=r"Log$_{10}$ $\rho_a$ RMS",
            title=f"{architecture.upper()} response audit",
        )
        axes[row, 2].grid(alpha=0.22, axis="y")
        axes[row, 2].legend(frameon=False, fontsize=7)
    color_axis = fig.add_axes([0.405, 0.035, 0.25, 0.014])
    fig.colorbar(
        last_image,
        cax=color_axis,
        orientation="horizontal",
        label=r"Shared $\log_{10}\rho$ ($\Omega$ m) scale",
    )
    fig.suptitle(
        "Matched WILLY architecture smoke tests: structure and response fit",
        fontsize=14,
        y=0.995,
    )
    fig.subplots_adjust(
        left=0.06, right=0.98, bottom=0.085, top=0.90, hspace=0.52, wspace=0.32
    )
    _save(fig, "agents_architecture_comparison.png")


def make_agents_inv3d_mt1d_baseline() -> None:
    """The historical Inv3DAgent(physics="mt1d") baseline: tiled independent
    1-D columns shared through a GCN over real WILLY station geometry. Kept
    to explain agents_inv3d_section.png/agents_inv3d_uncertainty.png as a
    contrast with the genuine 3-D Maxwell physics="mt3d" route documented
    alongside it -- not the new 3-D Maxwell workflow itself.
    """
    try:
        import torch

        torch.manual_seed(73)
    except ImportError:
        pass
    np.random.seed(73)
    sites = ensure_sites(
        PROJECT_ROOT / "data" / "AMT" / "WILLY_data" / "L18PLT",
        recursive=True,
        verbose=0,
    )
    # Use the WILLY stations' own measured band (~1 Hz-10 kHz) rather than
    # the class default (1e-4-1e3 Hz): that legacy default reaches skin
    # depths of hundreds of kilometres -- far below anything this AMT
    # data actually constrains -- and was the root cause of the
    # mantle-depth section this figure used to show.
    _, measured_frequency, _ = sites_to_features_1d(
        sites, comp="xy", n_freqs=64,
    )
    agent = Inv3DAgent(
        physics="mt1d",
        n_layers=10,
        freqs=np.geomspace(
            float(measured_frequency.min()),
            float(measured_frequency.max()),
            32,
        ),
        n_train_profiles=300,
        epochs=80,
        radius=250.0,
        hidden=(256, 128, 64),
        dropout=0.1,
        n_mc=50,
    )
    result = agent.execute({"sites": sites})
    if result.status == "failed":
        raise RuntimeError(result.error)
    # The resistivity_section itself is intentionally not saved: it is a
    # tiled mt1d baseline whose value is in the printed contract below and
    # in agents.rst's worked source, not in a captured "reference" image.
    print("agents_inv3d_uncertainty.png", ascii(result.summary))
    _save(
        result["figures"]["uncertainty_map"], "agents_inv3d_uncertainty.png"
    )


def make_agents_inv3d_willy_2km() -> None:
    """Execute genuine 3-D Maxwell training on the real WILLY geometry."""
    try:
        import torch

        torch.manual_seed(73)
    except ImportError:
        pass
    np.random.seed(73)
    sites = ensure_sites(
        PROJECT_ROOT / "data" / "AMT" / "WILLY_data" / "L18PLT",
        recursive=True,
        verbose=0,
    )
    dense_features, measured_frequency, _ = sites_to_features_1d(
        sites,
        comp="xy",
        n_freqs=64,
    )
    rho_reference = float(
        10.0 ** np.nanmedian(dense_features[:, : measured_frequency.size])
    )
    depth_frequency = float(frequency_for_depth(2000.0, rho_reference))
    frequency_min = float(measured_frequency.min())
    frequency_max = min(float(measured_frequency.max()), depth_frequency)
    _, frequency, _ = sites_to_features_1d(
        sites,
        comp="xy",
        n_freqs=3,
        freq_min=frequency_min,
        freq_max=frequency_max,
    )
    result = Inv3DAgent(
        n_layers=5,
        freqs=frequency,
        depth_max=2000.0,
        n_train_profiles=20,
        epochs=40,
        radius=250.0,
        hidden=(64, 32),
        dropout=0.1,
        n_mc=0,
        physics="mt3d",
        geology_grid_nx_ny=4,
        geology_grid_nz=4,
        max_mesh_cells=60_000,
    ).execute({"sites": sites, "topography": True})
    if result.status == "failed":
        raise RuntimeError(result.error)
    print(
        "agents_inv3d_willy_2km.png",
        ascii(result.summary),
        "recovery=",
        ascii(result.get("mt3d_recovery")),
    )
    fig = result["figures"]["resistivity_section"]
    ax = fig.axes[0]
    ax.set_title(
        f"WILLY MT3D-trained GCN — {frequency[0]:.2f}–"
        f"{frequency[-1]:.2f} Hz, 2 km model",
        fontsize=10,
        fontweight="bold",
    )
    ax.xaxis.set_label_position("top")
    ax.xaxis.tick_top()
    _save(fig, "agents_inv3d_willy_2km_section.png")
    # topography_section is intentionally not saved here: see the
    # topography discussion in agents.rst, which points readers at the
    # code-dropdown above instead of a fixed smoke-test capture.


def make_agents_inv2d_willy_topography() -> None:
    """Execute a compact WILLY Inv2DAgent terrain-output contract."""
    try:
        import torch

        torch.manual_seed(43)
    except ImportError:
        pass
    np.random.seed(43)
    sites = ensure_sites(
        PROJECT_ROOT / "data" / "AMT" / "WILLY_data" / "L18PLT",
        recursive=True,
        verbose=0,
    )
    dense_features, measured_frequency, _ = sites_to_features_1d(
        sites,
        comp="xy",
        n_freqs=64,
    )
    rho_reference = float(
        10.0 ** np.nanmedian(dense_features[:, : measured_frequency.size])
    )
    frequency_min = max(
        float(measured_frequency.min()),
        float(frequency_for_depth(2000.0, rho_reference)),
    )
    _, frequency, _ = sites_to_features_1d(
        sites,
        comp="xy",
        n_freqs=16,
        freq_min=frequency_min,
        freq_max=float(measured_frequency.max()),
    )
    result = Inv2DAgent(
        n_depth=8,
        freqs=frequency,
        depth_max=2000.0,
        n_train_profiles=8,
        n_stations_per_profile=12,
        epochs=2,
    ).execute({"sites": sites, "topography": True})
    if result.status == "failed":
        raise RuntimeError(result.error)
    fig = result["figures"]["topography_section"]
    fig.set_constrained_layout(False)
    fig.set_size_inches(12.5, 6.2)
    fig.subplots_adjust(left=0.09, right=0.88, bottom=0.14, top=0.70)
    _save(fig, "agents_inv2d_willy_topography.png")


def make_agents_inv2d_synthetic_demo() -> None:
    """EMInverter2D on a fully synthetic target -- no real survey involved.

    Unlike the WILLY examples above, nothing here is constrained by a real
    line's station count, frequency band, or wall-clock budget: the "true"
    model is a hand-built, smoothed compact conductor (matched to the
    training distribution's own correlation length so it is not an
    out-of-distribution shape), forward-modelled with the same
    :class:`~pycsamt.forward.maxwell.mt2d.MT2DAdapter` the training data
    generator uses internally, then inverted by a freshly trained
    :class:`~pycsamt.ai.inversion.inv2d.EMInverter2D`. This is the
    EMInverter2D/Maxwell2DDatasetConfig science-API layer beneath
    ``Inv2DAgent(physics="mt2d")``, not a captured ``Inv2DAgent.execute()``
    result -- ``Inv2DAgent`` always reads its observed panel from a real
    ``Sites`` collection, even on the synthetic-training ``"mt2d"`` path.
    """
    from scipy.ndimage import gaussian_filter

    from pycsamt.ai.geology import GeologyGrid
    from pycsamt.ai.inversion.inv2d import EMInverter2D
    from pycsamt.ai.training.dataset2d import _solver_mesh_and_conductivity
    from pycsamt.api.mesh import PYCSAMT_MESH

    try:
        import torch

        torch.manual_seed(7)
    except ImportError:
        pass
    np.random.seed(7)

    n_sta, n_depth = 24, 24
    dx_m, depth_max_m = 100.0, 1200.0
    grid = GeologyGrid.regular_2d(
        nx=n_sta, nz=n_depth, dx_m=dx_m, dz_m=depth_max_m / n_depth
    )
    freqs = np.geomspace(2.0, 400.0, 12)

    # Hand-built true model: a background with a compact conductor blob,
    # Gaussian-smoothed so its spatial scale matches the correlation
    # lengths below rather than presenting the network with a sharp box
    # it never saw an analogue of during training.
    log_true_rho = np.full((n_depth, n_sta), 2.0)  # log10(100 ohm.m)
    z_core = slice(int(n_depth * 0.35), int(n_depth * 0.65))
    x_core = slice(int(n_sta * 0.35), int(n_sta * 0.65))
    log_true_rho[z_core, x_core] = np.log10(5.0)  # 5 ohm.m core
    log_true_rho = gaussian_filter(log_true_rho, sigma=(1.8, 1.8))
    true_rho = 10.0**log_true_rho

    mesh_true, conductivity_true = _solver_mesh_and_conductivity(
        grid, true_rho, freqs, safety_factor=6.0, max_mesh_cells=20_000
    )
    receivers = ReceiverSet(
        [[float(x), 0.0] for x in grid.x_m],
        [f"S{i:02d}" for i in range(n_sta)],
    )
    true_problem = MaxwellProblem(
        mesh_true, conductivity_true, freqs, receivers, components=("zxy",)
    )
    true_result = MT2DAdapter().solve(true_problem)
    if not true_result.success:
        raise RuntimeError("synthetic true-model forward solve failed.")
    true_survey = SurveyData(
        impedance=true_result.impedance_v_a,
        frequencies_hz=true_result.frequencies_hz,
        station_names=true_result.receiver_names,
        components=true_result.components,
        coordinates_m=[[float(x), 0.0, 0.0] for x in grid.x_m],
        valid=true_result.valid,
    )

    mu0 = 4.0e-7 * np.pi

    def _samples_to_arrays(samples):
        n_freq = len(freqs)
        X = np.empty((len(samples), 2, n_freq, n_sta), dtype=np.float32)
        y = np.empty((len(samples), n_depth, n_sta), dtype=np.float32)
        for i, sample in enumerate(samples):
            zxy = sample.survey.impedance[:, :, 0]
            f = sample.survey.frequencies_hz[None, :]
            rho_a = np.abs(zxy) ** 2 / (2.0 * np.pi * f * mu0)
            phase = np.degrees(np.angle(zxy))
            X[i, 0] = np.log10(np.clip(rho_a, 1e-12, None)).T
            X[i, 1] = phase.T
            y[i] = np.log10(sample.resistivity_ohm_m)
        return X, y

    def _survey_to_x(survey):
        zxy = survey.impedance[:, :, 0]
        f = survey.frequencies_hz[None, :]
        rho_a = np.abs(zxy) ** 2 / (2.0 * np.pi * f * mu0)
        phase = np.degrees(np.angle(zxy))
        x = np.empty((1, 2, len(freqs), n_sta), dtype=np.float32)
        x[0, 0] = np.log10(np.clip(rho_a, 1e-12, None)).T
        x[0, 1] = phase.T
        return x

    config = Maxwell2DDatasetConfig(
        dataset_id="agents-inv2d-synthetic-demo",
        grid=grid,
        correlation_length_x_m=(300.0, 900.0),
        correlation_length_z_m=(100.0, 300.0),
        frequencies_hz=freqs,
        station_x_m=grid.x_m,
        n_realizations=100,
        seed=0,
        log_resistivity_mean=2.0,
        log_resistivity_std=0.4,
        components=("zxy",),
        mesh_safety_factor=4.0,
        max_mesh_cells=20_000,
    )
    dataset = generate_2d_maxwell_dataset(config)
    X_train, y_train = _samples_to_arrays(dataset.select("train"))

    inverter = EMInverter2D(
        n_components=2, n_depth=n_depth, n_stations=n_sta, n_freqs=len(freqs)
    )
    inverter.fit(X_train, y_train, epochs=80, patience=15, verbose=False, seed=7)
    history = inverter._history

    X_true = _survey_to_x(true_survey)
    log_pred = inverter.predict(X_true)[0]
    log_true = np.log10(true_rho)
    rmse = float(np.sqrt(np.nanmean((log_pred - log_true) ** 2)))
    print(
        "agents_inv2d_synthetic_demo.png",
        f"n_train={len(dataset.select('train'))}",
        f"epochs_run={len(history['train_loss'])}",
        f"rmse={rmse:.3f}",
    )

    fig = plot_inversion_result_2d(
        log_pred,
        log_true=log_true,
        depths=None,
        stations=grid.x_m / 1000.0,
        depth_max=depth_max_m,
        show_mesh=True,
        mesh_style=PYCSAMT_MESH.style_for("review"),
        train_loss=np.array(history["train_loss"]),
        val_loss=np.array(history["val_loss"]),
        suptitle="EMInverter2D on a fully synthetic target",
    )
    _save(fig, "agents_inv2d_synthetic_demo.png")


def make_pinn2d_willy_topography_audit() -> None:
    """Execute and display a deliberately short, auditable WILLY PINN run."""
    try:
        import torch

        torch.manual_seed(241)
    except ImportError:
        pass
    np.random.seed(241)
    line = PROJECT_ROOT / "data" / "AMT" / "WILLY_data" / "L18PLT"
    sites = ensure_sites(line, recursive=True, verbose=0)
    inverter = PINNInverter2D(
        line,
        n_layers=8,
        depth_max=2000.0,
        n_freqs=24,
        mode="te",
        smoothness_weight=0.01,
        lateral_weight=0.005,
        epochs=20,
        lr=1e-2,
        device="cpu",
        verbose=0,
    ).fit(verbose=False)
    loss = inverter.convergence_curve()["loss"].to_numpy()
    section = inverter.resistivity_section(as_log10=True)
    thickness = inverter.thickness_section()

    station_list = list(sites)
    latitude = np.array([site.coords[0] for site in station_list], dtype=float)
    longitude = np.array(
        [site.coords[1] for site in station_list], dtype=float
    )
    elevation = np.array(
        [site.coords[2] for site in station_list], dtype=float
    )
    radius = 6_371_000.0
    lat_rad = np.deg2rad(latitude)
    lon_rad = np.deg2rad(longitude)
    dlat = np.diff(lat_rad)
    dlon = np.diff(lon_rad)
    hav = np.sin(dlat / 2) ** 2 + (
        np.cos(lat_rad[:-1]) * np.cos(lat_rad[1:]) * np.sin(dlon / 2) ** 2
    )
    segment = 2 * radius * np.arcsin(np.sqrt(np.clip(hav, 0, 1)))
    chainage = np.r_[0.0, np.cumsum(segment)] / 1000.0

    depth = np.linspace(0.0, 2000.0, 201)
    resampled = np.empty((len(depth), section.shape[1]))
    for station in range(section.shape[1]):
        interfaces = np.cumsum(thickness[:, station])
        layer = np.searchsorted(interfaces, depth, side="right")
        layer = np.clip(layer, 0, section.shape[0] - 1)
        resampled[:, station] = section[layer, station]

    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.4))
    axes[0].plot(np.arange(1, len(loss) + 1), loss, color="#b91c1c", lw=2)
    axes[0].axhline(
        loss[0], color="#111827", ls="--", lw=1, label="initial recorded loss"
    )
    axes[0].set(
        xlabel="Epoch",
        ylabel="Total objective",
        title="Executed convergence audit",
    )
    axes[0].grid(alpha=0.22)
    axes[0].legend(frameon=False, fontsize=8)

    image = axes[1].pcolormesh(
        chainage,
        depth,
        resampled,
        shading="nearest",
        cmap="turbo",
        vmin=0.0,
        vmax=4.0,
    )
    axes[1].invert_yaxis()
    axes[1].set(
        xlabel="Profile chainage (km)",
        ylabel="Depth below station (m)",
        title="Flat-datum section",
    )

    x_grid = np.broadcast_to(chainage, resampled.shape)
    y_grid = elevation[None, :] - depth[:, None]
    axes[2].contourf(
        x_grid,
        y_grid,
        resampled,
        levels=np.linspace(0, 4, 33),
        cmap="turbo",
        vmin=0.0,
        vmax=4.0,
        extend="both",
    )
    axes[2].plot(chainage, elevation, color="#111827", lw=1.4)
    axes[2].scatter(chainage, elevation, marker="v", s=18, color="#111827")
    axes[2].set(
        xlabel="Profile chainage (km)",
        ylabel="Elevation (m)",
        title="Same model draped on EDI elevations",
    )
    color_axis = fig.add_axes([0.955, 0.19, 0.012, 0.60])
    colorbar = fig.colorbar(image, cax=color_axis)
    colorbar.set_label(r"$\log_{10}\rho$ ($\Omega$ m)")
    fig.suptitle(
        "WILLY L18 PINN2D smoke test — rejected because loss increased",
        fontsize=13,
    )
    fig.subplots_adjust(
        left=0.055, right=0.935, bottom=0.16, top=0.82, wspace=0.3
    )
    _save(fig, "pinn2d_willy_topography_audit.png")


def make_pinn2d_input_diagnostic() -> None:
    """Inspect WILLY frequency support and TE/TM disagreement before fitting."""
    line = PROJECT_ROOT / "data" / "AMT" / "WILLY_data" / "L18PLT"
    observations = sites_to_obs_2d(line, comp_te="xy", comp_tm="yx")
    frequency = np.logspace(
        np.log10(min(o.freq.min() for o in observations)),
        np.log10(max(o.freq.max() for o in observations)),
        48,
    )

    def interpolate(obs, values):
        order = np.argsort(obs.freq)
        result = np.interp(
            np.log10(frequency), np.log10(obs.freq[order]), values[order]
        )
        outside = ((frequency < obs.freq.min() * (1.0 - 1e-12)) |
                   (frequency > obs.freq.max() * (1.0 + 1e-12)))
        result[outside] = np.nan
        return result

    te_rho = np.vstack(
        [interpolate(o, np.log10(o.rho_te)) for o in observations]
    )
    tm_rho = np.vstack(
        [interpolate(o, np.log10(o.rho_tm)) for o in observations]
    )
    te_phase = np.vstack([interpolate(o, o.phase_te) for o in observations])
    tm_phase = np.vstack([interpolate(o, o.phase_tm) for o in observations])
    valid = np.isfinite(te_rho) & np.isfinite(te_phase)

    fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.0))
    ax_coverage, ax_rho, ax_phase, ax_count = axes.ravel()
    coverage = ax_coverage.pcolormesh(
        frequency, np.arange(len(observations)), valid,
        cmap="Blues", shading="nearest", vmin=0, vmax=1,
    )
    ax_coverage.set_xscale("log")
    ax_coverage.set(
        xlabel="Frequency (Hz)", ylabel="Station index",
        title="Common-grid support (blue = valid)",
    )

    ax_rho.plot(
        frequency, np.nanmedian(te_rho, axis=0), color="#2563eb",
        lw=2, label="TE / xy",
    )
    ax_rho.plot(
        frequency, np.nanmedian(tm_rho, axis=0), color="#f15a29",
        lw=2, label="TM / |yx|",
    )
    ax_rho.set_xscale("log")
    ax_rho.set(
        xlabel="Frequency (Hz)", ylabel=r"Median $\log_{10}\rho_a$",
        title="Modes are not interchangeable",
    )
    ax_rho.grid(alpha=0.22, which="both")
    ax_rho.legend(frameon=False)

    phase_difference = te_phase - tm_phase
    phase_image = ax_phase.pcolormesh(
        frequency, np.arange(len(observations)),
        phase_difference, cmap="coolwarm", shading="nearest", vmin=-180, vmax=180,
    )
    ax_phase.set_xscale("log")
    ax_phase.set(
        xlabel="Frequency (Hz)", ylabel="Station index",
        title="TE minus TM phase (degrees)",
    )
    fig.colorbar(phase_image, ax=ax_phase, shrink=0.86, label="Phase difference (°)")

    phase_rms = np.sqrt(np.nanmean(phase_difference**2, axis=1))
    rho_rms = np.sqrt(np.nanmean((te_rho - tm_rho) ** 2, axis=1))
    station_index = np.arange(len(observations))
    ax_count.bar(station_index, phase_rms, color="#3e65b0", alpha=0.78)
    ax_count.set(
        xlabel="Station index", ylabel="TE–TM phase RMS (degrees)",
        title="Mode disagreement by station",
    )
    ax_rho_rms = ax_count.twinx()
    ax_rho_rms.plot(station_index, rho_rms, color="#f15a29", marker="o", ms=3)
    ax_rho_rms.set_ylabel(r"TE–TM $\log_{10}\rho_a$ RMS", color="#f15a29")
    ax_count.grid(axis="y", alpha=0.22)
    fig.suptitle("WILLY L18 inputs to PINNInverter2D — inspect before optimization")
    fig.tight_layout()
    _save(fig, "pinn2d_input_diagnostic.png")


def make_pinn2d_regularization_anatomy() -> None:
    """Visualize the separate effects of vertical and lateral roughness."""
    stations, layers = 30, 12
    x = np.linspace(0, 1, stations)
    truth = np.full((stations, layers), 2.5)
    truth[:, :3] = 2.0
    truth[(x > 0.28) & (x < 0.68), 4:8] = 0.8
    truth[x > 0.76, 6:10] = 3.5
    rng = np.random.default_rng(27)
    proxy = truth + rng.normal(0, 0.42, truth.shape)

    dz = np.eye(layers, k=1) - np.eye(layers)
    dz = dz[:-1]
    dx = np.eye(stations, k=1) - np.eye(stations)
    dx = dx[:-1]
    identity = np.eye(stations * layers)
    pz = np.kron(np.eye(stations), dz.T @ dz)
    px = np.kron(dx.T @ dx, np.eye(layers))

    settings = [(0.0, 0.0), (3.0, 0.0), (0.0, 3.0), (1.2, 1.2)]
    estimates = []
    diagnostics = []
    vector = proxy.ravel()
    for lam_z, lam_x in settings:
        estimate = np.linalg.solve(identity + lam_z * pz + lam_x * px, vector)
        estimate = estimate.reshape(stations, layers)
        estimates.append(estimate)
        data = np.mean((estimate - proxy) ** 2)
        vertical = np.mean(np.diff(estimate, axis=1) ** 2)
        lateral = np.mean(np.diff(estimate, axis=0) ** 2)
        diagnostics.append((data, vertical + lateral))

    fig, axes = plt.subplots(2, 3, figsize=(12.8, 7.7))
    panels = [truth, proxy, estimates[1], estimates[2], estimates[3]]
    titles = [
        "Known log-resistivity", "Noisy data-fit proxy",
        r"Vertical only: $\lambda_z=3$", r"Lateral only: $\lambda_x=3$",
        r"Balanced: $\lambda_z=\lambda_x=1.2$",
    ]
    for ax, values, title in zip(axes.ravel()[:5], panels, titles):
        image = ax.imshow(
            values.T, origin="upper", aspect="auto", cmap="turbo",
            vmin=0.5, vmax=3.8,
        )
        ax.set(title=title, xlabel="Station index", ylabel="Layer index")

    ax_trade = axes.ravel()[5]
    for (data, roughness), (lam_z, lam_x) in zip(diagnostics, settings):
        ax_trade.scatter(data, roughness, s=65)
        ax_trade.annotate(
            f"({lam_z:g}, {lam_x:g})", (data, roughness),
            xytext=(5, 5), textcoords="offset points", fontsize=8,
        )
    ax_trade.set(
        xlabel="Data-proxy MSE", ylabel="Vertical + lateral roughness",
        title=r"Trade-off labels: $(\lambda_z,\lambda_x)$",
    )
    ax_trade.grid(alpha=0.22)
    color_axis = fig.add_axes([0.92, 0.18, 0.012, 0.64])
    fig.colorbar(image, cax=color_axis, label=r"$\log_{10}\rho$ ($\Omega$ m)")
    fig.suptitle("Regularization changes what the pseudo-2-D model is allowed to express")
    fig.subplots_adjust(left=0.07, right=0.89, bottom=0.09, top=0.88, hspace=0.36, wspace=0.30)
    _save(fig, "pinn2d_regularization_anatomy.png")


def make_reporting_validation_dashboard() -> None:
    """Build one dashboard entirely from structured validation reports."""
    rng = np.random.default_rng(73)

    depth = np.linspace(0, 1800, 24)
    x = np.linspace(0, 1, 32)
    truth = 2.5 + 0.8 * np.exp(-((depth[:, None] - 850) / 260) ** 2)
    truth = truth - 1.5 * np.exp(-((x[None, :] - 0.55) / 0.16) ** 2) * np.exp(
        -((depth[:, None] - 650) / 330) ** 2
    )
    prediction = truth + rng.normal(0, 0.08 + depth[:, None] / 9000, truth.shape)
    recovery = recovery_report(prediction, truth, depth_axis=0)

    n_station, n_frequency, n_component = 12, 18, 2
    observed = rng.normal(size=(n_station, n_frequency, n_component)) + 1j * rng.normal(
        size=(n_station, n_frequency, n_component)
    )
    error = np.full(observed.shape, 0.12)
    bias = np.linspace(0.02, 0.22, n_station)[:, None, None]
    predicted = observed + bias + 1j * rng.normal(0, 0.06, observed.shape)
    residual = response_residual_report(
        predicted, observed, errors=error, kind="l2",
        station_names=[f"S{i:02d}" for i in range(n_station)],
        frequencies_hz=np.logspace(3, -1, n_frequency),
        components=("zxy", "zyx"),
    )

    reference = rng.normal(0, 1, (160, 2))
    field = np.vstack([rng.normal(0, 0.8, (8, 2)), [[3.4, 3.0], [-3.2, 2.8]]])
    ood = flag_out_of_distribution(
        field, reference, method="mahalanobis", quantile=0.95
    )

    y_true = rng.normal(size=500)
    y_mean = y_true + rng.normal(0, 0.45, y_true.shape)
    y_std = np.full(y_true.shape, 0.38)
    calibration = reliability_curve(y_true, y_mean, y_std)

    fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.0))
    ax_depth, ax_residual, ax_ood, ax_calibration = axes.ravel()
    ax_depth.plot(recovery.depth_rmse, depth, color="#f15a29", lw=2, label="RMSE")
    ax_depth.plot(recovery.depth_mae, depth, color="#3e65b0", lw=2, label="MAE")
    ax_depth.invert_yaxis()
    ax_depth.set(
        xlabel=r"Log-resistivity error ($\log_{10}\Omega$ m)",
        ylabel="Depth (m)", title="Synthetic recovery by depth",
    )
    ax_depth.grid(alpha=0.22)
    ax_depth.legend(frameon=False)

    ax_residual.bar(
        np.arange(n_station), residual.by_station,
        color="#3e65b0", alpha=0.82,
    )
    ax_residual.set(
        xlabel="Station index", ylabel="Mean squared normalized residual",
        title="Response residuals retain axis structure",
    )
    ax_residual.grid(axis="y", alpha=0.22)

    colors = np.where(ood.flagged, "#f15a29", "#3e65b0")
    ax_ood.bar(np.arange(len(ood.scores)), ood.scores, color=colors)
    ax_ood.axhline(
        ood.threshold, color="#111827", ls="--", lw=1.4,
        label=f"95% reference threshold = {ood.threshold:.2f}",
    )
    ax_ood.set(
        xlabel="Field sample", ylabel="Mahalanobis score",
        title="Domain gate preserves rejected rows",
    )
    ax_ood.legend(frameon=False, fontsize=8)
    ax_ood.grid(axis="y", alpha=0.22)

    ax_calibration.plot([0, 1], [0, 1], color="#111827", ls="--", label="ideal")
    ax_calibration.plot(
        calibration.levels, calibration.coverage,
        color="#f15a29", marker="o", lw=2, label="empirical",
    )
    ax_calibration.set(
        xlim=(0.45, 1.0), ylim=(0.45, 1.0),
        xlabel="Nominal coverage", ylabel="Empirical coverage",
        title=f"Calibration and sharpness ({calibration.sharpness:.2f})",
    )
    ax_calibration.grid(alpha=0.22)
    ax_calibration.legend(frameon=False)

    fig.suptitle("Machine-readable evidence behind an AI inversion report")
    fig.tight_layout()
    _save(fig, "reporting_validation_dashboard.png")


def make_scientific_validation_anatomy() -> None:
    """Compare the four independent validation views on one synthetic audit."""
    from scipy.ndimage import gaussian_filter

    rng = np.random.default_rng(314)
    nz, nx = 28, 48
    z = np.linspace(0.0, 1.0, nz)[:, None]
    x = np.linspace(-1.0, 1.0, nx)[None, :]
    truth = 2.25 + 0.55 * z
    truth = truth - 1.15 * np.exp(-((x + 0.20) / 0.24) ** 2
                                  - ((z - 0.58) / 0.16) ** 2)
    truth = truth + 0.45 * np.exp(-((x - 0.55) / 0.20) ** 2
                                  - ((z - 0.30) / 0.12) ** 2)
    prediction = gaussian_filter(truth, sigma=(1.1, 1.6))
    prediction += rng.normal(0.0, 0.035, truth.shape)
    prediction[z[:, 0] > 0.78] += 0.18
    recovery = recovery_report(prediction, truth, ssim_window=7)

    n_station, n_frequency = 14, 20
    frequency = np.logspace(-1, 3, n_frequency)
    observed = np.ones((n_station, n_frequency, 2), dtype=complex) * (70 + 45j)
    predicted = observed.copy()
    predicted[8:12, 6:14, 1] += 9 + 7j
    predicted[:, :3, :] += 3 + 2j
    error = np.full(observed.shape, 5.0)
    residual = response_residual_report(
        predicted, observed, errors=error, kind="l2",
        station_names=[f"S{i + 1:02d}" for i in range(n_station)],
        frequencies_hz=frequency, components=("zxy", "zyx"),
    )
    _finish_scientific_validation_anatomy(
        rng, nz, n_station, frequency, observed, predicted, error,
        recovery, residual,
    )


def make_validation_aggregation_trap() -> None:
    """Show why equal global RMSE does not imply equal earth recovery."""
    from scipy.ndimage import gaussian_filter

    from pycsamt.ai.validation import recovery_report

    rng = np.random.default_rng(90210)
    nz, nx = 42, 72
    z = np.linspace(0, 1, nz)[:, None]
    x = np.linspace(-1, 1, nx)[None, :]
    truth = 2.7 + 0.35 * z + 0.10 * np.sin(2 * np.pi * x)
    conductor = ((x + 0.18) / 0.34) ** 2 + ((z - 0.43) / 0.17) ** 2 <= 1
    truth = np.where(conductor, 1.35, truth)

    shifted_mask = ((x + 0.04) / 0.34) ** 2 + ((z - 0.50) / 0.17) ** 2 <= 1
    localized = np.where(shifted_mask, 1.35, 2.7 + 0.35 * z + 0.10 * np.sin(2 * np.pi * x))
    localized += gaussian_filter(rng.normal(scale=0.015, size=(nz, nx)), 1.2)
    localized_rmse = float(np.sqrt(np.mean((localized - truth) ** 2)))

    diffuse_noise = gaussian_filter(rng.normal(size=(nz, nx)), sigma=1.6)
    diffuse_noise -= diffuse_noise.mean()
    diffuse_noise *= localized_rmse / np.sqrt(np.mean(diffuse_noise**2))
    diffuse = truth + diffuse_noise

    report_local = recovery_report(localized, truth)
    report_diffuse = recovery_report(diffuse, truth)
    depth = np.linspace(0, 1800, nz)

    fig, axes = plt.subplots(2, 3, figsize=(13.0, 7.2))
    extent = (0, 3.5, 1.8, 0)
    panels = [truth, localized, diffuse]
    titles = ["Known truth", "Shifted conductor", "Diffuse correlated error"]
    for ax, panel, title in zip(axes[0], panels, titles):
        image = ax.imshow(panel, aspect="auto", extent=extent, cmap="turbo",
                          vmin=1.2, vmax=3.2)
        ax.set(title=title, xlabel="Distance (km)", ylabel="Depth (km)")
    fig.colorbar(image, ax=axes[0, 2], label=r"$\log_{10}\rho$ ($\Omega\cdot$m)",
                 fraction=0.046, pad=0.04)

    error_local = localized - truth
    error_diffuse = diffuse - truth
    vmax = max(np.max(np.abs(error_local)), np.max(np.abs(error_diffuse)))
    axes[1, 0].imshow(error_local, aspect="auto", extent=extent, cmap="coolwarm",
                      vmin=-vmax, vmax=vmax)
    axes[1, 0].set(title="Localized signed error", xlabel="Distance (km)",
                   ylabel="Depth (km)")
    axes[1, 1].imshow(error_diffuse, aspect="auto", extent=extent, cmap="coolwarm",
                      vmin=-vmax, vmax=vmax)
    axes[1, 1].set(title="Diffuse signed error", xlabel="Distance (km)",
                   ylabel="Depth (km)")
    axes[1, 2].plot(report_local.depth_rmse, depth, color="#dc2626",
                    label=f"shifted, SSIM={report_local.ssim:.3f}")
    axes[1, 2].plot(report_diffuse.depth_rmse, depth, color="#2563eb",
                    label=f"diffuse, SSIM={report_diffuse.ssim:.3f}")
    axes[1, 2].invert_yaxis()
    axes[1, 2].set(xlabel=r"Depth-row RMSE in $\log_{10}\rho$", ylabel="Depth (m)",
                   title="The error location changes the verdict")
    axes[1, 2].legend(frameon=False, fontsize=8)
    axes[1, 2].grid(alpha=0.22)
    fig.suptitle("Equal global RMSE, unequal geological recovery")
    fig.subplots_adjust(left=0.06, right=0.94, bottom=0.08, top=0.88,
                        wspace=0.28, hspace=0.30)
    _save(fig, "validation_aggregation_trap.png")
    print(
        "validation aggregation:",
        {"shifted": {"rmse": report_local.rmse, "mae": report_local.mae,
                     "ssim": report_local.ssim},
         "diffuse": {"rmse": report_diffuse.rmse, "mae": report_diffuse.mae,
                     "ssim": report_diffuse.ssim}},
    )


def make_validation_paired_bootstrap() -> None:
    """Compare row and parent-survey bootstrap uncertainty for a baseline delta."""
    rng = np.random.default_rng(440)
    n_survey, rows_per_survey = 18, 8
    survey_effect = rng.normal(0.0, 0.055, n_survey)
    improvement = 0.035 + survey_effect
    row_delta = np.repeat(improvement, rows_per_survey)
    row_delta += rng.normal(0.0, 0.012, row_delta.size)
    survey_id = np.repeat(np.arange(n_survey), rows_per_survey)
    survey_delta = np.array([row_delta[survey_id == i].mean() for i in range(n_survey)])

    n_boot = 8000
    row_boot = row_delta[rng.integers(0, len(row_delta), (n_boot, len(row_delta)))].mean(1)
    group_boot = survey_delta[
        rng.integers(0, n_survey, (n_boot, n_survey))
    ].mean(1)
    row_ci = np.quantile(row_boot, [0.025, 0.975])
    group_ci = np.quantile(group_boot, [0.025, 0.975])

    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.2))
    axes[0].bar(np.arange(1, n_survey + 1), survey_delta,
                color=np.where(survey_delta >= 0, "#2563eb", "#dc2626"))
    axes[0].axhline(0, color="#111827", lw=1)
    axes[0].set(xlabel="Parent survey", ylabel="Baseline RMSE - AI RMSE",
                title="Improvement is heterogeneous")
    axes[0].grid(alpha=0.22, axis="y")

    axes[1].hist(row_boot, bins=55, density=True, alpha=0.62, color="#f59e0b",
                 label="resample 144 rows")
    axes[1].hist(group_boot, bins=55, density=True, alpha=0.55, color="#2563eb",
                 label="resample 18 surveys")
    axes[1].axvline(0, color="#111827", ls="--")
    axes[1].set(xlabel="Mean paired RMSE improvement", ylabel="Bootstrap density",
                title="Rows understate sampling uncertainty")
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].grid(alpha=0.2)

    means = [row_boot.mean(), group_boot.mean()]
    lows = [row_ci[0], group_ci[0]]
    highs = [row_ci[1], group_ci[1]]
    axes[2].errorbar([0, 1], means,
                     yerr=[np.array(means) - lows, highs - np.array(means)],
                     fmt="o", color="#1d4ed8", capsize=6, lw=2)
    axes[2].axhline(0, color="#dc2626", ls="--", label="no improvement")
    axes[2].set_xticks([0, 1], ["row bootstrap", "survey bootstrap"], rotation=15)
    axes[2].set(ylabel="Mean improvement with 95% interval",
                title="Inference follows the independent unit")
    axes[2].legend(frameon=False, fontsize=8)
    axes[2].grid(alpha=0.22, axis="y")
    fig.suptitle("Paired baseline validation must preserve parent-survey dependence")
    fig.tight_layout()
    _save(fig, "validation_paired_bootstrap.png")
    print(
        "validation bootstrap:",
        {"mean": float(group_boot.mean()),
         "row_ci": np.round(row_ci, 4).tolist(),
         "survey_ci": np.round(group_ci, 4).tolist(),
         "surveys_improved": int(np.sum(survey_delta > 0)),
         "n_surveys": n_survey},
    )


def _finish_scientific_validation_anatomy(
    rng, nz, n_station, frequency, observed, predicted, error,
    recovery, residual,
) -> None:
    residual_cells = np.mean(np.abs((predicted - observed) / error) ** 2, axis=2).T

    true_u = rng.normal(size=1800)
    levels = np.linspace(0.1, 0.99, 16)
    calibrated = reliability_curve(true_u, np.zeros_like(true_u),
                                   np.ones_like(true_u), levels=levels)
    narrow = reliability_curve(true_u, np.zeros_like(true_u),
                               np.full_like(true_u, 0.55), levels=levels)

    reference = rng.multivariate_normal([0.0, 0.0], [[1.0, 0.55], [0.55, 0.8]], 260)
    field = np.vstack([rng.multivariate_normal([0.1, -0.1], [[0.8, 0.4], [0.4, 0.7]], 24),
                       [[3.2, 3.0], [-3.1, 1.8], [2.8, -2.5]]])
    ood = flag_out_of_distribution(field, reference, quantile=0.975)

    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.6))
    ax = axes[0, 0]
    depth = np.linspace(0, 2.1, nz)
    ax.plot(recovery.depth_rmse, depth, "o-", color="#2563eb", label="RMSE")
    ax.plot(recovery.depth_mae, depth, "s-", color="#f15a29", label="MAE")
    ax.invert_yaxis()
    ax.set(xlabel=r"Error in $\log_{10}\rho$", ylabel="Depth (km)",
           title=f"Recovery: RMSE={recovery.rmse:.3f}, SSIM={recovery.ssim:.3f}")
    ax.grid(alpha=0.22)
    ax.legend(frameon=False)

    ax = axes[0, 1]
    image = ax.imshow(residual_cells, origin="lower", aspect="auto",
                      extent=(0.5, n_station + 0.5, np.log10(frequency[0]),
                              np.log10(frequency[-1])), cmap="magma")
    ax.set(xlabel="Station", ylabel=r"$\log_{10}$ frequency (Hz)",
           title=f"Response: mean normalized squared residual={residual.overall.value:.2f}")
    fig.colorbar(image, ax=ax, label="component-mean squared residual")

    ax = axes[1, 0]
    ax.plot(levels, levels, "k--", label="ideal")
    ax.plot(calibrated.levels, calibrated.coverage, "o-", color="#2563eb",
            label=f"calibrated (sharpness={calibrated.sharpness:.2f})")
    ax.plot(narrow.levels, narrow.coverage, "s-", color="#dc2626",
            label=f"overconfident (sharpness={narrow.sharpness:.2f})")
    ax.set(xlim=(0, 1), ylim=(0, 1), xlabel="Nominal coverage",
           ylabel="Empirical coverage", title="Uncertainty: coverage before sharpness")
    ax.grid(alpha=0.22)
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 1]
    ax.scatter(reference[:, 0], reference[:, 1], s=14, alpha=0.25,
               color="#64748b", label="training reference")
    colors = np.where(ood.flagged, "#dc2626", "#16a34a")
    ax.scatter(field[:, 0], field[:, 1], s=45, c=colors, edgecolor="white",
               linewidth=0.5, label="field samples")
    ax.set(xlabel="Standardized feature 1", ylabel="Standardized feature 2",
           title=f"Domain support: {ood.flagged.sum()}/{len(field)} flagged")
    ax.grid(alpha=0.22)
    ax.legend(frameon=False, fontsize=8)

    fig.suptitle("Four diagnostics answer four different scientific questions", fontsize=13)
    fig.tight_layout()
    _save(fig, "scientific_validation_anatomy.png")
    print(
        "validation anatomy:",
        {"recovery_rmse": recovery.rmse, "recovery_ssim": recovery.ssim,
         "response_loss": residual.overall.value,
         "calibration_mace": calibrated.calibration.value,
         "overconfident_mace": narrow.calibration.value,
         "ood_flagged": int(ood.flagged.sum()), "ood_total": len(field)},
    )


def make_validation_stress_envelope() -> None:
    """Execute all validation reports across one declared stress axis."""
    from scipy.ndimage import gaussian_filter

    rng = np.random.default_rng(1618)
    severity = np.linspace(0.0, 1.5, 13)
    nz, nx = 24, 36
    z = np.linspace(0, 1, nz)[:, None]
    x = np.linspace(-1, 1, nx)[None, :]
    truth = 2.4 + 0.45 * z - 1.0 * np.exp(
        -((x + 0.15) / 0.24) ** 2 - ((z - 0.52) / 0.17) ** 2
    )
    observed = np.ones((10, 18, 2), dtype=complex) * (60.0 + 35.0j)
    errors = np.full(observed.shape, 5.0)
    reference = rng.normal(size=(500, 4))

    recovery_rmse, response_loss, coverage, flagged = [], [], [], []
    for value in severity:
        prediction = gaussian_filter(truth, sigma=(0.35 + value, 0.5 + value))
        prediction = prediction + value * 0.13 * z
        recovery_rmse.append(recovery_report(prediction, truth).rmse)

        response_prediction = observed + value * (4.0 + 2.5j)
        report = response_residual_report(
            response_prediction, observed, errors=errors, kind="l2"
        )
        response_loss.append(report.overall.value)

        true_u = rng.normal(size=2400)
        mean_u = true_u + rng.normal(scale=0.25 + 0.24 * value, size=true_u.shape)
        std_u = np.full_like(true_u, 0.32)
        curve = reliability_curve(true_u, mean_u, std_u, levels=[0.90])
        coverage.append(curve.coverage[0])

        field = rng.normal(loc=value * 0.95, size=(120, 4))
        domain = flag_out_of_distribution(
            field, reference, method="mahalanobis", quantile=0.975
        )
        flagged.append(domain.fraction_flagged)

    recovery_rmse = np.asarray(recovery_rmse)
    response_loss = np.asarray(response_loss)
    coverage = np.asarray(coverage)
    flagged = np.asarray(flagged)
    thresholds = {"recovery": 0.10, "response": 1.0,
                  "coverage_low": 0.85, "ood_fraction": 0.10}

    fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.5), sharex=True)
    panels = [
        (recovery_rmse, thresholds["recovery"], "Synthetic recovery RMSE",
         r"RMSE in $\log_{10}\rho$", False),
        (response_loss, thresholds["response"], "Normalized response loss",
         "Mean squared normalized residual", False),
        (coverage, thresholds["coverage_low"], "Nominal 90% interval coverage",
         "Empirical coverage", True),
        (flagged, thresholds["ood_fraction"], "Domain-support failures",
         "Fraction flagged OOD", False),
    ]
    for ax, (values, gate, title, ylabel, pass_above) in zip(axes.flat, panels):
        ax.plot(severity, values, "o-", color="#2563eb", lw=2)
        ax.axhline(gate, color="#dc2626", ls="--", label=f"gate={gate:g}")
        failed = values < gate if pass_above else values > gate
        ax.scatter(severity[failed], values[failed], color="#dc2626", zorder=4,
                   label="failed")
        ax.set(title=title, ylabel=ylabel)
        ax.grid(alpha=0.22)
        ax.legend(frameon=False, fontsize=8)
    for ax in axes[-1]:
        ax.set_xlabel("Declared stress severity")
    fig.suptitle("Executed challenge sweep: the operating envelope ends at the first mandatory gate")
    fig.tight_layout()
    _save(fig, "validation_stress_envelope.png")

    first_failures = {}
    for name, values, gate, pass_above in [
        ("recovery", recovery_rmse, thresholds["recovery"], False),
        ("response", response_loss, thresholds["response"], False),
        ("coverage", coverage, thresholds["coverage_low"], True),
        ("ood", flagged, thresholds["ood_fraction"], False),
    ]:
        failed = values < gate if pass_above else values > gate
        first_failures[name] = (
            None if not np.any(failed) else float(severity[np.flatnonzero(failed)[0]])
        )
    print("validation stress first failures:", first_failures)


def make_validation_gate_dashboard() -> None:
    frequency = np.logspace(np.log10(1.01), 4, 24)
    samples = generate_dataset(
        solver="mt1d",
        n_samples=240,
        freqs=frequency,
        n_layers=5,
        rho_range=(1.0, 10_000.0),
        depth_max=2000.0,
        noise_level=0.05,
        noise_type="field",
        include_phase=True,
        seed=137,
        n_jobs=1,
        output=None,
        verbose=False,
    )
    train, validation, test = samples.split(
        val_frac=0.15, test_frac=0.15, seed=137
    )
    pool_x = np.vstack([train.X, validation.X])
    pool_y = np.vstack([train.y, validation.y])
    inverter = EMInverter1D(
        n_features=48,
        n_layers=5,
        arch="fcn",
        solver="mt1d",
        device="cpu",
        log_thickness=False,
        augment_noise=0.01,
    )
    inverter.fit(
        pool_x,
        pool_y,
        epochs=25,
        batch_size=64,
        lr=1e-3,
        patience=7,
        val_frac=0.15,
        grad_clip=1.0,
        seed=137,
        verbose=False,
    )
    mae = np.mean(np.abs(inverter.predict(test.X) - test.y), axis=0)

    line = PROJECT_ROOT / "data" / "AMT" / "WILLY_data" / "L18PLT"
    sites = ensure_sites(line, recursive=True, verbose=0)
    field, _, station_names = sites_to_features_1d(
        sites, comp="xy", n_freqs=24, freq_min=1.01, freq_max=10_000.0
    )
    low, high = np.percentile(pool_x, [1, 99], axis=0)
    outside = np.mean((field < low) | (field > high), axis=1)
    dimension = classify_dimensionality(sites)["dim"].to_numpy(dtype=int)
    dim_fraction = np.array(
        [(dimension == value).mean() for value in (0, 1, 2)]
    )

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.3))
    ax_error, ax_domain, ax_dim = axes
    normalized_error = np.r_[mae[:5] / 0.3, mae[5:] / 100.0]
    colors = [
        "#16a34a" if value <= 1 else "#dc2626" for value in normalized_error
    ]
    ax_error.bar(np.arange(1, 10), normalized_error, color=colors)
    ax_error.axhline(1.0, color="#111827", ls="--", lw=1.1, label="gate")
    ax_error.set(
        xlabel="R1--R5, then H1--H4",
        ylabel="Error / allowed error",
        title="Synthetic external-test gates",
    )
    ax_error.legend(frameon=False, fontsize=8)
    ax_error.grid(alpha=0.2, axis="y")

    station_index = np.arange(len(station_names))
    ax_domain.bar(station_index, outside, color="#dc2626")
    ax_domain.axhline(
        0.10, color="#111827", ls="--", lw=1.1, label="review threshold"
    )
    ax_domain.set(
        xlabel="WILLY L18 station index",
        ylabel="Fraction outside P1--P99",
        title="Field-domain gate",
    )
    ax_domain.legend(frameon=False, fontsize=8)
    ax_domain.grid(alpha=0.2, axis="y")

    ax_dim.bar(
        ["1-D", "2-D", "3-D"],
        dim_fraction,
        color=["#16a34a", "#2563eb", "#dc2626"],
    )
    ax_dim.set_ylim(0, 1)
    ax_dim.set(ylabel="Fraction of samples", title="WILLY tensor evidence")
    ax_dim.grid(alpha=0.2, axis="y")
    ax_dim.text(
        2,
        dim_fraction[2] + 0.025,
        f"{dim_fraction[2]:.1%}",
        ha="center",
        fontsize=9,
    )
    fig.suptitle(
        "Validation integrates independent gates: decision = rejected",
        fontsize=13,
    )
    fig.tight_layout()
    _save(fig, "validation_gate_dashboard.png")


def make_gcn_3d_context() -> None:
    """The real spatial graph GCNInverter3D/Inv3DAgent build from the real
    WILLY L18 station geometry, at the same adjacency radius used by the
    executed physics="mt3d" example, coloured by each station's own real
    measured apparent resistivity -- not a fabricated synthetic layout.
    """
    sites = ensure_sites(
        PROJECT_ROOT / "data" / "AMT" / "WILLY_data" / "L18PLT",
        recursive=True,
        verbose=0,
    )
    coords = sites_to_coords_3d(sites, recursive=False, verbose=0)
    dense_features, measured_frequency, _ = sites_to_features_1d(
        sites, comp="xy", n_freqs=64,
    )
    # Each station's own median log10 apparent resistivity across the
    # measured band -- real field data, robust to the ragged per-station
    # frequency coverage sites_to_features_1d leaves as NaN at the edges.
    station_log_rho = np.nanmedian(
        dense_features[:, : measured_frequency.size], axis=1
    )

    radius = 250.0  # matches Inv3DAgent(physics="mt3d")'s own WILLY example
    adjacency_raw = build_adjacency(
        coords, radius, self_loops=False, normalise=False
    )
    adjacency_norm = build_adjacency(coords, radius)
    degree = adjacency_raw.sum(axis=1).astype(int)
    x, y = coords[:, 0], coords[:, 1]
    n_stations = coords.shape[0]

    fig = plt.figure(figsize=(10.2, 4.6))
    ax = fig.add_subplot(111, projection="3d")
    for i in range(n_stations):
        for j in range(i + 1, n_stations):
            if adjacency_raw[i, j] > 0:
                ax.plot(
                    [x[i], x[j]],
                    [y[i], y[j]],
                    [station_log_rho[i], station_log_rho[j]],
                    color="#94a3b8",
                    alpha=0.4,
                    lw=0.6 + 2.6 * adjacency_norm[i, j],
                )

    sc = ax.scatter(
        x,
        y,
        station_log_rho,
        c=station_log_rho,
        cmap="viridis_r",
        s=60,
        edgecolor="#111827",
        linewidth=0.35,
        depthshade=False,
    )
    ax.set_xlabel("Easting offset (m)")
    ax.set_ylabel("Northing offset (m)")
    ax.set_zlabel(r"median $\log_{10}\rho_a$ (measured)")
    ax.set_title(
        f"Real GCN spatial graph: WILLY L18, r={radius:.0f} m "
        f"({int(adjacency_raw.sum() // 2)} edges, "
        f"{int(np.sum(degree == 0))} isolated)",
        pad=12,
    )
    ax.view_init(elev=16, azim=-70)
    cb = fig.colorbar(sc, ax=ax, shrink=0.74, pad=0.08)
    cb.set_label(r"median $\log_{10}\rho_a$ [$\Omega\cdot$m]")
    fig.tight_layout()
    _save(fig, "gcn_3d_context.png")


def make_domain_gap_corruption_diagnostic() -> None:
    """Show signal corruption, missingness, and distribution diagnostics."""
    rng = np.random.default_rng(8)
    n_station, n_frequency = 18, 28
    frequencies = np.logspace(3.5, -1.0, n_frequency)
    station_x = np.linspace(0.0, 5100.0, n_station)
    log_f = np.log10(frequencies)[None, :]
    lateral = np.sin(station_x[:, None] / 850.0)
    log_mag = 1.85 - 0.22 * log_f + 0.18 * lateral
    phase = 43.0 + 8.0 * np.tanh(log_f - 1.0) + 5.0 * lateral
    impedance = 10.0**log_mag * np.exp(1j * np.deg2rad(phase))
    clean = SurveyData(
        impedance[:, :, None],
        frequencies,
        [f"S{i:02d}" for i in range(n_station)],
        ["xy"],
        np.column_stack([station_x, np.zeros(n_station)]),
    )
    config = CorruptionConfig(
        noise_level_range=(0.04, 0.12),
        error_floor_fraction=0.03,
        static_shift_log10_sigma=0.12,
        distortion_gain_log10_sigma=0.06,
        distortion_twist_deg_sigma=4.0,
        station_dropout_rate=0.06,
        frequency_dropout_rate=0.04,
        random_dropout_rate=0.06,
        outlier_rate=0.025,
        coordinate_sigma_m=8.0,
    )
    corrupted, _ = apply_corruption_suite(clean, config=config, seed=19)

    field_like_z = np.array(corrupted.impedance)
    valid = np.isfinite(field_like_z)
    field_like_z[valid] *= np.exp(
        1j * np.deg2rad(rng.normal(6.0, 4.0, valid.sum()))
    )
    field_like = SurveyData(
        field_like_z,
        frequencies,
        clean.station_names,
        clean.components,
        clean.coordinates_m,
        impedance_error=corrupted.impedance_error,
    )
    report = compare_survey_distributions(corrupted, field_like)

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.2))
    extent = [station_x[0] / 1000, station_x[-1] / 1000,
              np.log10(frequencies[-1]), np.log10(frequencies[0])]
    panels = [
        (np.log10(np.abs(clean.impedance[:, :, 0])).T, "Clean log$_{10}|Z|$"),
        (np.log10(np.abs(corrupted.impedance[:, :, 0])).T,
         "Corrupted log$_{10}|Z|$ (white = dropout)"),
    ]
    for ax, (values, title) in zip(axes[0], panels):
        image = ax.imshow(
            values,
            origin="lower",
            aspect="auto",
            extent=extent,
            cmap="viridis",
            vmin=1.0,
            vmax=2.8,
        )
        ax.set_title(title)
        ax.set_xlabel("Profile distance (km)")
        ax.set_ylabel("log$_{10}$ frequency (Hz)")
        fig.colorbar(image, ax=ax, label="log$_{10}|Z|$")

    feature_specs = [
        ("log_impedance_magnitude", r"$\log_{10}|Z|$"),
        ("phase_deg", "Phase (degrees)"),
    ]
    for ax, (feature, label) in zip(axes[1], feature_specs):
        if feature == "log_impedance_magnitude":
            sim = np.log10(np.abs(corrupted.impedance[corrupted.valid]))
            fld = np.log10(np.abs(field_like.impedance[field_like.valid]))
        else:
            sim = np.angle(corrupted.impedance[corrupted.valid], deg=True)
            fld = np.angle(field_like.impedance[field_like.valid], deg=True)
        for values, name, color in (
            (sim, "corrupted synthetic", "#2563eb"),
            (fld, "field-like target", "#dc2626"),
        ):
            values = np.sort(values)
            ax.step(values, np.arange(1, values.size + 1) / values.size,
                    where="post", label=name, color=color, lw=1.8)
        comparison = report.comparisons[feature]
        ax.set_title(
            f"Empirical CDF: KS = {comparison.ks_statistic:.3f}, "
            f"mean difference = {comparison.mean_difference:.2f}"
        )
        ax.set_xlabel(label)
        ax.set_ylabel("Cumulative probability")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)

    fig.suptitle("A domain-gap audit separates corrupted appearance from statistical match")
    fig.tight_layout()
    _save(fig, "domain_gap_corruption_diagnostic.png")


def make_experiment_gate_diagnostic() -> None:
    """Visualise pinned provenance and multi-seed gate outcomes."""
    reference = DatasetReference("demo-v1", "a" * 64, "b" * 64, "c" * 64)
    criteria = (
        AcceptanceCriterion("test.nrms", "<=", 2.0),
        AcceptanceCriterion("test.rmse", "<=", 0.50),
        AcceptanceCriterion("test.coverage", ">=", 0.90),
    )
    config = ExperimentConfig(
        "multi-seed-demo",
        "learning_2d",
        reference,
        SeedPlan(42, "docs/experiments"),
        {"architecture": "unet"},
        {"epochs": 80},
        {"solver": "mt2d"},
        criteria,
    )
    run_ids = np.arange(1, 13)
    nrms = np.array([1.52, 1.78, 1.91, 2.08, 1.67, 1.86,
                     1.73, 2.16, 1.95, 1.81, 1.69, 1.88])
    rmse = np.array([0.39, 0.44, 0.47, 0.46, 0.51, 0.43,
                     0.41, 0.48, 0.45, 0.42, 0.38, 0.46])
    coverage = np.array([0.93, 0.91, 0.90, 0.92, 0.94, 0.89,
                         0.91, 0.93, 0.88, 0.92, 0.95, 0.91])
    evaluations = [
        config.evaluate_gate(
            {"test.nrms": a, "test.rmse": b, "test.coverage": c}
        )
        for a, b, c in zip(nrms, rmse, coverage)
    ]
    passed = np.array([item.passed for item in evaluations])

    fig = plt.figure(figsize=(11.2, 7.1))
    grid = fig.add_gridspec(2, 1, height_ratios=[0.75, 1.8], hspace=0.32)
    ax_flow = fig.add_subplot(grid[0])
    ax_flow.set_axis_off()
    labels = [
        ("Dataset", "manifest + split +\nnormalizer hashes"),
        ("Configuration", "model + training +\nphysics + seed plan"),
        ("Execution", "derived run seed +\nobserved metrics"),
        ("Decision", "all predeclared\ncriteria must pass"),
    ]
    xs = np.linspace(0.1, 0.9, len(labels))
    for index, (x, (title, body)) in enumerate(zip(xs, labels)):
        box = FancyBboxPatch(
            (x - 0.095, 0.28), 0.19, 0.47,
            boxstyle="round,pad=0.015", facecolor="#eff6ff",
            edgecolor="#2563eb", linewidth=1.3,
        )
        ax_flow.add_patch(box)
        ax_flow.text(x, 0.62, title, ha="center", weight="bold")
        ax_flow.text(x, 0.43, body, ha="center", va="center", fontsize=8.5)
        if index < len(labels) - 1:
            ax_flow.add_patch(FancyArrowPatch(
                (x + 0.1, 0.51), (xs[index + 1] - 0.1, 0.51),
                arrowstyle="-|>", mutation_scale=13, color="#475569",
            ))
    ax_flow.set_xlim(0, 1)
    ax_flow.set_ylim(0, 1)

    subgrid = grid[1].subgridspec(3, 1, hspace=0.12)
    series = [
        (nrms, 2.0, "NRMS", "lower"),
        (rmse, 0.50, "Recovery RMSE", "lower"),
        (coverage, 0.90, "Coverage", "higher"),
    ]
    for row, (values, threshold, label, direction) in enumerate(series):
        ax = fig.add_subplot(subgrid[row])
        colors = np.where(passed, "#16a34a", "#dc2626")
        ax.scatter(run_ids, values, c=colors, s=58, zorder=3,
                   edgecolor="white", linewidth=0.6)
        ax.axhline(threshold, color="#111827", ls="--", lw=1.2,
                   label=f"threshold = {threshold:.2f}")
        if direction == "lower":
            ax.axhspan(threshold, max(values.max(), threshold) * 1.04,
                       color="#fee2e2", alpha=0.45)
        else:
            ax.axhspan(min(values.min(), threshold) * 0.99, threshold,
                       color="#fee2e2", alpha=0.45)
        ax.set_ylabel(label)
        ax.set_xlim(0.4, 12.6)
        ax.grid(axis="x", alpha=0.2)
        ax.legend(loc="upper right", fontsize=8)
        if row < 2:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel("Independent network-seed run")
            ax.set_xticks(run_ids)
    fig.suptitle(
        f"One frozen experiment, twelve seeds: {passed.sum()} complete gate passes",
        fontsize=13,
    )
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.08, top=0.91)
    _save(fig, "experiment_gate_diagnostic.png")


def make_roadmap_capability_matrix() -> None:
    """Show implementation maturity without conflating it with release readiness."""
    milestones = [f"M{i}" for i in range(11)]
    levels = [1, 1, 1, 1, 2, 2, 1, 0, -1, 2, -1]
    labels = ["Implemented", "Verified", "Wired", "Gated", "Released"]
    colors = ["#dbeafe", "#bfdbfe", "#93c5fd", "#60a5fa", "#2563eb"]

    fig, (ax, ax_key) = plt.subplots(
        2, 1, figsize=(11.2, 5.8), gridspec_kw={"height_ratios": [4.4, 1.0]}
    )
    for row, label in enumerate(labels):
        ax.axhspan(row - 0.5, row + 0.5, color=colors[row], alpha=0.30)
        ax.axhline(row + 0.5, color="white", lw=1.5)
    ax.bar(
        milestones, np.asarray(levels) + 1, width=0.68,
        color="#2563eb", edgecolor="#1e3a8a", linewidth=0.8,
    )
    for index, level in enumerate(levels):
        status = "Planned" if level < 0 else labels[level]
        ax.text(index, max(level + 0.72, 0.18), status, ha="center", va="bottom",
                fontsize=8, rotation=35)
    ax.set_ylim(0, 5.55)
    ax.set_yticks(np.arange(5) + 0.5, labels)
    ax.set_ylabel("Highest demonstrated readiness level")
    ax.set_title("AI inversion roadmap: evidence reached by each milestone")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", alpha=0.12)

    ax_key.set_axis_off()
    ax_key.text(0.01, 0.72, "A milestone advances only when evidence accumulates:",
                weight="bold", color="#0f172a")
    stages = ["implementation", "tests/benchmarks", "workflow wiring",
              "automatic gate", "blind release"]
    xs = np.linspace(0.08, 0.92, len(stages))
    for i, (x, stage) in enumerate(zip(xs, stages)):
        ax_key.text(x, 0.28, stage, ha="center", va="center", fontsize=8.5,
                    bbox={"boxstyle": "round,pad=0.32", "fc": colors[i],
                          "ec": "#2563eb", "lw": 0.8})
        if i < len(stages) - 1:
            ax_key.add_patch(FancyArrowPatch(
                (x + 0.075, 0.28), (xs[i + 1] - 0.075, 0.28),
                arrowstyle="-|>", mutation_scale=11, color="#475569",
                transform=ax_key.transAxes,
            ))
    fig.subplots_adjust(left=0.13, right=0.98, bottom=0.06, top=0.90, hspace=0.22)
    _save(fig, "roadmap_capability_matrix.png")


def main() -> None:
    make_workflow()
    make_training_distribution()
    make_convergence()
    make_predicted_section()
    make_uncertainty()
    make_willy_depth_support()
    make_data_contracts_normalization()
    make_forward_physics_halfspace_benchmark()
    make_forward_physics_mesh_sensitivity()
    make_geology_prior_diagnostic()
    make_geology_topographic_surface()
    make_geology_3d_composition()
    make_inference_calibration_diagnostic()
    make_losses_penalty_anatomy()
    make_losses_regularization_tradeoff()
    make_hybrid_physics_refinement_audit()
    make_dataset2d_realization_gallery()
    make_dataset2d_tri_topo_gallery()
    make_dataset2d_tri_topo_ai_inversion()
    make_dataset2d_response_anatomy()
    make_dataset3d_realization_gallery()
    make_dataset3d_response_anatomy()
    make_data_preparation_contract()
    make_data_preparation_coverage()
    make_model_selection_willy_dimension()
    make_model_selection_graph_radius()
    make_model_selection_tradeoff()
    make_training_executed_audit()
    make_training_augmentation_audit()
    make_training_trainer_controls()
    make_hybrid_paired_diagnostic()
    make_uncertainty_coverage_reliability()
    make_uncertainty_calibration_regimes()
    make_uncertainty_depth_propagation()
    make_agents_execution_contract()
    make_agents_real_data_sections()
    make_agents_executed_1d_audit()
    make_agents_architecture_comparison()
    # make_agents_inv2d_willy_topography() is intentionally not called here:
    # its only output, agents_inv2d_willy_topography.png, was dropped from
    # agents.rst in favour of a code-dropdown pointing at this function, so
    # running it during every doc build would just burn compute on an image
    # nothing references. The function itself stays so the dropdown has a
    # real, importable source to show and run.
    make_agents_inv3d_mt1d_baseline()
    make_agents_inv3d_willy_2km()
    make_agents_inv2d_synthetic_demo()
    make_pinn2d_willy_topography_audit()
    make_pinn2d_input_diagnostic()
    make_pinn2d_regularization_anatomy()
    make_reporting_validation_dashboard()
    make_scientific_validation_anatomy()
    make_validation_aggregation_trap()
    make_validation_paired_bootstrap()
    make_validation_stress_envelope()
    make_validation_gate_dashboard()
    make_gcn_3d_context()
    make_domain_gap_corruption_diagnostic()
    make_experiment_gate_diagnostic()
    make_roadmap_capability_matrix()
    print(f"Wrote AI inversion documentation figures to {OUT}")


if __name__ == "__main__":
    main()
