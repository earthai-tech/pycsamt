"""Generate the executed figures for theory/ai_inversion.rst."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.ai.validation import recovery_report  # noqa: E402
from pycsamt.ai.geology import (  # noqa: E402
    ElectricalLayer,
    EllipsoidalLens,
    GaussianCorrelation,
    GeologyGrid,
    TopographicSurface,
    generate_layered_geology,
    insert_lenses,
)
from pycsamt.ai.losses import (  # noqa: E402
    model_huber_loss,
    model_l1_loss,
    model_l2_loss,
    total_variation_loss,
)
from pycsamt.forward.maxwell import skin_depth_m  # noqa: E402

IMAGE_DIR = ROOT / "docs/source/images/theory"


def _save(fig: plt.Figure, name: str) -> None:
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(IMAGE_DIR / name, dpi=190, bbox_inches="tight")
    plt.close(fig)


def make_ai_nonuniqueness() -> tuple[float, float]:
    """Show why two different depth models can have similar responses."""
    depth_m = np.linspace(25.0, 1975.0, 40)
    frequency_hz = np.geomspace(1.0, 1000.0, 14)
    reference_rho = 100.0
    delta = skin_depth_m(reference_rho, frequency_hz)
    kernel = np.exp(-depth_m[None, :] / delta[:, None])
    kernel /= kernel.sum(axis=1, keepdims=True)

    truth = np.full(depth_m.size, 2.3)
    truth[(depth_m >= 550) & (depth_m <= 850)] = 1.0
    alternative = np.full(depth_m.size, 2.3)
    alternative[(depth_m >= 430) & (depth_m <= 1050)] = 1.55
    response_true = kernel @ truth
    response_alt = kernel @ alternative
    response_rmse = float(np.sqrt(np.mean((response_true - response_alt) ** 2)))
    model_rmse = float(np.sqrt(np.mean((truth - alternative) ** 2)))

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4), constrained_layout=True)
    axes[0].plot(truth, depth_m / 1000, lw=2, label="compact conductor")
    axes[0].plot(alternative, depth_m / 1000, lw=2, ls="--", label="broad conductor")
    axes[0].invert_yaxis(); axes[0].set(xlabel=r"$\log_{10}\rho$", ylabel="depth (km)", title="Different earth models")
    axes[0].legend(fontsize=8)
    im = axes[1].imshow(kernel, aspect="auto", cmap="magma", origin="upper",
                        extent=[depth_m[0]/1000, depth_m[-1]/1000,
                                np.log10(frequency_hz[-1]), np.log10(frequency_hz[0])])
    axes[1].set(xlabel="depth (km)", ylabel=r"$\log_{10} f$ (Hz)", title="Depth-sensitivity kernel")
    fig.colorbar(im, ax=axes[1], label="normalized sensitivity")
    axes[2].plot(response_true, frequency_hz, "o-", label="compact")
    axes[2].plot(response_alt, frequency_hz, "s--", label="broad")
    axes[2].set_yscale("log"); axes[2].set(xlabel="response proxy", ylabel="frequency (Hz)", title="Similar observable responses")
    axes[2].legend(fontsize=8)
    _save(fig, "ai_inversion_nonuniqueness.png")
    return model_rmse, response_rmse


def make_ai_domain_gap() -> tuple[float, float]:
    """Contrast in-distribution and out-of-distribution response features."""
    rng = np.random.default_rng(27)
    training = rng.multivariate_normal([2.0, 45.0], [[0.12, 0.8], [0.8, 18]], 500)
    validation = rng.multivariate_normal([2.03, 44.7], [[0.13, 0.7], [0.7, 19]], 160)
    field = rng.multivariate_normal([2.55, 58.0], [[0.20, 1.5], [1.5, 28]], 130)
    centre = training.mean(axis=0)
    scale = training.std(axis=0, ddof=1)
    val_score = np.linalg.norm((validation - centre) / scale, axis=1)
    field_score = np.linalg.norm((field - centre) / scale, axis=1)

    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.3), constrained_layout=True)
    axes[0].scatter(training[:, 0], training[:, 1], s=10, alpha=.25, label="synthetic train")
    axes[0].scatter(validation[:, 0], validation[:, 1], s=15, alpha=.45, label="synthetic test")
    axes[0].scatter(field[:, 0], field[:, 1], s=17, alpha=.55, label="field survey")
    axes[0].set(xlabel=r"median $\log_{10}\rho_a$", ylabel="median phase (degree)", title="Feature support")
    axes[0].legend(fontsize=8)
    bins = np.linspace(0, 8, 35)
    axes[1].hist(val_score, bins=bins, alpha=.7, density=True, label="synthetic test")
    axes[1].hist(field_score, bins=bins, alpha=.7, density=True, label="field")
    axes[1].axvline(3, color="black", ls="--", label="example OOD gate")
    axes[1].set(xlabel="standardized distance from training support", ylabel="density", title="Domain-gap diagnostic")
    axes[1].legend(fontsize=8)
    _save(fig, "ai_inversion_domain_gap.png")
    return float(np.median(val_score)), float(np.median(field_score))


def make_ai_validation() -> tuple[float, float, float]:
    """Compare scalar recovery with depth-resolved error and uncertainty."""
    rng = np.random.default_rng(8)
    nz, nx = 28, 42
    z, x = np.mgrid[:nz, :nx]
    truth = 2.5 - 1.45 * np.exp(-((x - 24) / 7.0) ** 2 - ((z - 14) / 4.0) ** 2)
    bias = 0.35 * (z / (nz - 1))
    prediction = 2.5 - 1.05 * np.exp(-((x - 22) / 9.0) ** 2 - ((z - 14) / 5.5) ** 2) + bias
    prediction += rng.normal(0, 0.035, truth.shape)
    report = recovery_report(prediction, truth, compute_ssim=True)
    uncertainty = 0.04 + 0.20 * (z / (nz - 1)) + 0.10 * np.exp(-((x - 24) / 8) ** 2)
    absolute_error = np.abs(prediction - truth)
    coverage = float(np.mean(absolute_error <= 1.96 * uncertainty))

    fig, axes = plt.subplots(1, 4, figsize=(15.4, 3.9), constrained_layout=True)
    kw = dict(cmap="turbo", vmin=.8, vmax=2.9, origin="upper", aspect="auto")
    im = axes[0].imshow(truth, **kw); axes[0].set_title("Known truth")
    axes[1].imshow(prediction, **kw); axes[1].set_title("AI recovery")
    er = axes[2].imshow(absolute_error, cmap="magma", vmin=0, vmax=.7, origin="upper", aspect="auto")
    axes[2].set_title("Absolute error")
    depth = np.arange(nz)
    axes[3].plot(report.depth_rmse, depth, label="depth RMSE")
    axes[3].plot(np.mean(1.96 * uncertainty, axis=1), depth, label="95% half-width")
    axes[3].invert_yaxis(); axes[3].set(xlabel=r"$\log_{10}\rho$", ylabel="depth-cell index", title="Error versus uncertainty")
    axes[3].legend(fontsize=8)
    fig.colorbar(im, ax=axes[:2], shrink=.8, label=r"$\log_{10}\rho$")
    fig.colorbar(er, ax=axes[2], shrink=.8, label="absolute error")
    _save(fig, "ai_inversion_validation.png")
    return report.rmse, report.r2, coverage


def make_ai_geology_prior() -> tuple[float, int]:
    """Compose layers, lenses, and an explicit terrain mask."""
    grid = GeologyGrid.regular_2d(nx=64, nz=34, dx_m=75, dz_m=50)
    base = generate_layered_geology(
        grid,
        [ElectricalLayer("cover", 45), ElectricalLayer("host", 900)],
        [550],
        seed=41,
        interface_relief_std_m=65,
        interface_correlation=GaussianCorrelation(650, 180),
    )
    bodies = [
        EllipsoidalLens("ore", 2900, 900, 700, 180, 8,
                        dip_deg=18, transition_fraction=.22),
        EllipsoidalLens("resistor", 1350, 1200, 480, 240, 2800,
                        dip_deg=-12, transition_fraction=.18),
    ]
    model = insert_lenses(base, bodies, conflict_policy="last")
    x = grid.x_m
    elevation = 405 + 48 * np.sin(2 * np.pi * x / np.ptp(x)) + 16 * np.cos(5 * np.pi * x / np.ptp(x))
    surface = TopographicSurface(grid, elevation, float(elevation.max()), source="synthetic profile")
    earth = surface.earth_mask()
    display = np.where(earth, np.log10(model.resistivity_ohm_m), np.nan)

    fig, axes = plt.subplots(1, 3, figsize=(14.3, 4.2), constrained_layout=True)
    extent = [grid.x_m[0] / 1000, grid.x_m[-1] / 1000,
              grid.z_m[-1] / 1000, grid.z_m[0] / 1000]
    im0 = axes[0].imshow(np.log10(base.resistivity_ohm_m), cmap="turbo", aspect="auto", extent=extent)
    axes[0].set(title="Correlated layered prior", xlabel="distance (km)", ylabel="depth (km)")
    axes[1].imshow(np.log10(model.resistivity_ohm_m), cmap="turbo", aspect="auto", extent=extent)
    axes[1].contour(grid.x_m / 1000, grid.z_m / 1000, model.lens_index >= 0,
                    levels=[.5], colors="white", linewidths=1)
    axes[1].set(title="Declared lens geometries", xlabel="distance (km)")
    axes[2].imshow(display, cmap="turbo", aspect="auto", extent=extent)
    axes[2].plot(grid.x_m / 1000, surface.surface_depth_m / 1000, color="black", lw=2)
    axes[2].set(title="Topography-aware active earth", xlabel="distance (km)")
    fig.colorbar(im0, ax=axes, shrink=.82, label=r"$\log_{10}\rho$ ($\Omega$ m)")
    _save(fig, "ai_inversion_geology_prior.png")
    return surface.relief_m, int(np.count_nonzero(~earth))


def make_ai_loss_behaviour() -> tuple[float, float, float, float]:
    """Compare robust pointwise losses and spatial penalties."""
    residual = np.linspace(-3, 3, 401)
    l1 = np.abs(residual)
    l2 = residual ** 2
    huber = np.where(np.abs(residual) <= .5, .5 * residual ** 2,
                     .5 * (np.abs(residual) - .25))
    truth = np.r_[np.full(18, 2.5), np.full(12, 1.0), np.full(18, 2.5)]
    smooth = np.convolve(np.pad(truth, 3, mode="edge"), np.ones(7) / 7, mode="valid")
    blocky = truth.copy(); blocky[18:30] = 1.2
    noisy = truth + np.random.default_rng(15).normal(0, .13, truth.size)
    values = [total_variation_loss(v).value for v in (smooth, blocky, noisy)]

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.1), constrained_layout=True)
    axes[0].plot(residual, l2, label="L2")
    axes[0].plot(residual, l1, label="L1")
    axes[0].plot(residual, huber, label=r"Huber, $\delta=0.5$")
    axes[0].set(xlabel="cell residual", ylabel="penalty", title="Pointwise influence of an error")
    axes[0].legend()
    cell = np.arange(truth.size)
    axes[1].plot(cell, truth, color="black", lw=2, label="truth")
    axes[1].plot(cell, smooth, label=f"smooth, TV={values[0]:.3f}")
    axes[1].plot(cell, blocky, label=f"blocky, TV={values[1]:.3f}")
    axes[1].plot(cell, noisy, alpha=.75, label=f"noisy, TV={values[2]:.3f}")
    axes[1].set(xlabel="horizontal cell", ylabel=r"$\log_{10}\rho$", title="Spatial penalty is a structural choice")
    axes[1].legend(fontsize=8)
    _save(fig, "ai_inversion_loss_behaviour.png")
    outlier = np.array([0., .1, -.1, 3.])
    zero = np.zeros_like(outlier)
    return (model_l1_loss(outlier, zero).value,
            model_l2_loss(outlier, zero).value,
            model_huber_loss(outlier, zero, delta=.5).value,
            values[2])


def main() -> int:
    model_rmse, response_rmse = make_ai_nonuniqueness()
    val_ood, field_ood = make_ai_domain_gap()
    recovery_rmse, recovery_r2, coverage = make_ai_validation()
    relief, air_cells = make_ai_geology_prior()
    l1, l2, huber, noisy_tv = make_ai_loss_behaviour()
    print("model RMSE:", f"{model_rmse:.3f}")
    print("response-proxy RMSE:", f"{response_rmse:.3f}")
    print("median OOD score (synthetic test, field):", f"{val_ood:.2f}", f"{field_ood:.2f}")
    print("recovery RMSE / R2:", f"{recovery_rmse:.3f}", f"{recovery_r2:.3f}")
    print("nominal 95% coverage:", f"{coverage:.3f}")
    print("topographic relief / air cells:", f"{relief:.1f}", air_cells)
    print("outlier example L1 / L2 / Huber:", f"{l1:.3f}", f"{l2:.3f}", f"{huber:.3f}")
    print("noisy-model TV:", f"{noisy_tv:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
