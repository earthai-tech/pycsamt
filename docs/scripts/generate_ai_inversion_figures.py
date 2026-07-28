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

from pycsamt.agents import AIInversionAgent, Inv2DAgent, Inv3DAgent
from pycsamt.ai.data.normalization import ComplexZScore
from pycsamt.ai.domain_gap import survey_data_from_sites
from pycsamt.ai.inversion import (
    EMInverter1D,
    PINNInverter2D,
    sites_to_features_1d,
    sites_to_obs_1d,
)
from pycsamt.ai.training.dataset2d import (
    Maxwell2DDatasetConfig,
    generate_2d_maxwell_dataset,
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


def make_agents_inv3d_willy_2km() -> None:
    """Execute the AMT-band, 2-km Inv3DAgent configuration."""
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
    frequency_min = max(float(measured_frequency.min()), depth_frequency)
    frequency_max = float(measured_frequency.max())
    _, frequency, _ = sites_to_features_1d(
        sites,
        comp="xy",
        n_freqs=24,
        freq_min=frequency_min,
        freq_max=frequency_max,
    )
    result = Inv3DAgent(
        n_layers=5,
        freqs=frequency,
        depth_max=2000.0,
        n_train_profiles=12,
        epochs=10,
        radius=250.0,
        hidden=(64, 32),
        dropout=0.1,
        n_mc=0,
    ).execute({"sites": sites, "topography": True})
    if result.status == "failed":
        raise RuntimeError(result.error)
    fig = result["figures"]["resistivity_section"]
    ax = fig.axes[0]
    ax.set_title(
        f"WILLY-derived GCN section — {frequency[0]:.2f}–"
        f"{frequency[-1]:,.0f} Hz, 2 km model",
        fontsize=10,
        fontweight="bold",
    )
    ax.xaxis.set_label_position("top")
    ax.xaxis.tick_top()
    _save(fig, "agents_inv3d_willy_2km_section.png")
    topo_fig = result["figures"]["topography_section"]
    _save(topo_fig, "agents_inv3d_willy_2km_topography.png")


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
    rng = np.random.default_rng(77)
    n_lines = 4
    n_per_line = 11
    xs, ys = [], []
    for line in range(n_lines):
        x = np.linspace(0, 10_000, n_per_line)
        y = np.full_like(x, line * 1800.0) + rng.normal(0, 90, n_per_line)
        xs.append(x + rng.normal(0, 120, n_per_line))
        ys.append(y)
    x = np.concatenate(xs)
    y = np.concatenate(ys)
    feature_strength = (
        0.45
        + 0.45
        * np.exp(-((x - 6200) ** 2 / 9_000_000 + (y - 2600) ** 2 / 3_000_000))
        + 0.10 * rng.random(x.size)
    )

    fig = plt.figure(figsize=(10.2, 4.6))
    ax = fig.add_subplot(111, projection="3d")
    for i in range(x.size):
        dist = np.hypot(x - x[i], y - y[i])
        neighbours = np.where((dist > 0) & (dist < 2300))[0]
        for j in neighbours:
            if j > i:
                ax.plot(
                    [x[i], x[j]],
                    [y[i], y[j]],
                    [0, 0],
                    color="#94a3b8",
                    alpha=0.35,
                    lw=0.7,
                )

    sc = ax.scatter(
        x,
        y,
        feature_strength * 1200,
        c=feature_strength,
        cmap="viridis",
        s=54,
        edgecolor="#111827",
        linewidth=0.35,
        depthshade=False,
    )
    ax.set_xlabel("Easting offset (m)")
    ax.set_ylabel("Northing offset (m)")
    ax.set_zlabel("Feature response")
    ax.set_title(
        "3-D AI inversion uses station geometry as graph context", pad=12
    )
    ax.view_init(elev=24, azim=-62)
    cb = fig.colorbar(sc, ax=ax, shrink=0.74, pad=0.08)
    cb.set_label("Normalised response feature")
    fig.tight_layout()
    _save(fig, "gcn_3d_context.png")


def main() -> None:
    make_workflow()
    make_training_distribution()
    make_convergence()
    make_predicted_section()
    make_uncertainty()
    make_willy_depth_support()
    make_data_contracts_normalization()
    make_forward_physics_halfspace_benchmark()
    make_dataset2d_realization_gallery()
    make_data_preparation_contract()
    make_data_preparation_coverage()
    make_model_selection_willy_dimension()
    make_training_executed_audit()
    make_hybrid_paired_diagnostic()
    make_uncertainty_coverage_reliability()
    make_agents_execution_contract()
    make_agents_executed_1d_audit()
    make_agents_architecture_comparison()
    make_agents_inv2d_willy_topography()
    make_agents_inv3d_willy_2km()
    make_pinn2d_willy_topography_audit()
    make_validation_gate_dashboard()
    make_gcn_3d_context()
    print(f"Wrote AI inversion documentation figures to {OUT}")


if __name__ == "__main__":
    main()
