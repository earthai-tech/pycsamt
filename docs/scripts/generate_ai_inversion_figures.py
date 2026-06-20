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
os.environ.setdefault("MPLCONFIGDIR", "/tmp/pycsamt-matplotlib")
os.environ.setdefault("PYCSAMT_DOCS_BUILD", "1")
os.environ.setdefault("PYCSAMT_DATA", "/tmp/pycsamt-docs-data")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from pycsamt.ai.plot.convergence import plot_convergence
from pycsamt.ai.plot.diagnostics import plot_uncertainty_bands
from pycsamt.ai.plot.inversion import plot_inversion_result_2d


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "source" / "images" / "user_guide" / "ai_inversion"


def _save(fig: plt.Figure, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / name, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _layer_model(depths: np.ndarray, rho_layers: list[float], breaks: list[float]) -> np.ndarray:
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

    for i, ((title, body), x, face, edge) in enumerate(zip(labels, xs, colors, borders)):
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
        ax.text(x, 0.61, title, ha="center", va="center", fontsize=10.5, weight="bold", color="#111827")
        ax.text(x, 0.47, body, ha="center", va="center", fontsize=8.4, color="#374151", linespacing=1.25)
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

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.1), gridspec_kw={"width_ratios": [0.95, 1.25]})
    ax_m, ax_r = axes

    response_stack = []
    for _ in range(90):
        breaks = sorted(rng.uniform([120, 350, 650], [260, 720, 1250]))
        rho_layers = 10 ** rng.uniform(0.4, 3.7, size=4)
        model = _layer_model(depths, list(rho_layers), breaks)
        if _ < 16:
            ax_m.plot(model, depths, lw=0.9, alpha=0.45, color="#2563eb")
        conductance = np.interp(np.log10(periods), [-2, 2], [model[:30].mean(), model[-40:].mean()])
        bowl = 0.22 * np.sin(np.log(periods) * rng.uniform(0.7, 1.2) + rng.uniform(-1.2, 1.2))
        response_stack.append(conductance + bowl + rng.normal(0, 0.035, size=periods.size))

    response_stack = np.vstack(response_stack)
    p10, p50, p90 = np.percentile(response_stack, [10, 50, 90], axis=0)
    observed = p50 + 0.10 * np.exp(-((np.log10(periods) - 0.35) ** 2) / 0.55)

    ax_m.invert_yaxis()
    ax_m.set_xlabel(r"$\log_{10}(\rho)$  ($\Omega$ m)")
    ax_m.set_ylabel("Depth (m)")
    ax_m.set_title("Sampled earth models", fontsize=10)
    ax_m.grid(alpha=0.25)

    ax_r.fill_between(periods, p10, p90, color="#93c5fd", alpha=0.35, label="training envelope")
    ax_r.plot(periods, p50, color="#2563eb", lw=1.8, label="synthetic median")
    ax_r.plot(periods, observed, color="#dc2626", lw=1.6, ls="--", label="observed survey")
    ax_r.set_xscale("log")
    ax_r.set_xlabel("Period (s)")
    ax_r.set_ylabel(r"Response feature, $\log_{10}(\rho_a)$")
    ax_r.set_title("Observed data should sit inside the synthetic envelope", fontsize=10)
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
        train = 0.42 * np.exp(-epochs / 19.5) + 0.018 + rng.normal(0, 0.003, epochs.size) + offset * 0.002
        val = 0.48 * np.exp(-epochs / 16.0) + 0.028 + 0.000025 * (epochs - 48).clip(min=0) ** 2
        val += rng.normal(0, 0.004, epochs.size) + offset * 0.003
        lr = np.where(epochs < 45, 1e-3, 3e-4)
        histories.append({"train_loss": train.clip(1e-4), "val_loss": val.clip(1e-4), "lr": lr})

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
    conductor = -0.95 * np.exp(-((xx - 7.2) ** 2 / 10.0 + (zz - 420) ** 2 / 90000))
    basement = 0.48 * np.exp(-((xx - 14.5) ** 2 / 18.0 + (zz - 1000) ** 2 / 180000))
    log_true = background + conductor + basement
    log_pred = log_true + 0.08 * np.sin(xx / 2.4) * np.exp(-zz / 1200) - 0.05 * np.exp(-((xx - 4) ** 2) / 4)

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
                "arrowprops": {"arrowstyle": "->", "color": "#111827", "lw": 0.8},
                "bbox": {"boxstyle": "round,pad=0.25", "fc": "white", "ec": "#cbd5e1", "alpha": 0.85},
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
    mean = 2.05 + 0.35 * np.tanh((depths - 600) / 360) - 0.55 * np.exp(-((depths - 360) ** 2) / 70000)
    spread = 0.07 + 0.00012 * depths + 0.08 * np.exp(-((depths - 850) ** 2) / 150000)
    truth = mean - 0.10 * np.exp(-((depths - 430) ** 2) / 55000) + 0.05 * np.sin(depths / 260)

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
        bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": "#cbd5e1", "alpha": 0.9},
    )
    _save(fig, "uncertainty_profile.png")


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
        + 0.45 * np.exp(-((x - 6200) ** 2 / 9_000_000 + (y - 2600) ** 2 / 3_000_000))
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
    ax.set_title("3-D AI inversion uses station geometry as graph context", pad=12)
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
    make_gcn_3d_context()
    print(f"Wrote AI inversion documentation figures to {OUT}")


if __name__ == "__main__":
    main()
