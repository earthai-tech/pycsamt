"""Visual validation of WILLY station ordering and L22 AI inversion.

Outputs are written to ``results/willy_site_ordering_validation`` by default.
The script intentionally loads every line through ``ensure_sites`` so all
plots exercise the package's validated spatial ordering policy.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pycsamt.agents.inv3d_agent import Inv3DAgent  # noqa: E402
from pycsamt.emtools._core import ensure_sites  # noqa: E402
from pycsamt.topo import (  # noqa: E402
    extract_chainage,
    extract_station_names,
    plot_topo_section,
)

matplotlib.use("Agg", force=True)


DEFAULT_DATA = ROOT / "data" / "AMT" / "WILLY_DATA"
DEFAULT_OUTPUT = ROOT / "results" / "willy_site_ordering_validation"


def _save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=190, bbox_inches="tight")
    plt.close(fig)
    print(f"saved: {path}")


def _line_directories(data_root: Path) -> list[Path]:
    return sorted(
        path
        for path in data_root.iterdir()
        if path.is_dir() and any(path.glob("*.edi"))
    )


def _write_station_table(lines: dict[str, object], output: Path) -> None:
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "line",
                "order",
                "station",
                "latitude",
                "longitude",
                "elevation_m",
                "chainage_km",
            ]
        )
        for line_name, sites in lines.items():
            chain = extract_chainage(sites)
            for index, (site, chain_km) in enumerate(zip(sites, chain), 1):
                lat, lon, elev = site.coords
                writer.writerow(
                    [
                        line_name,
                        index,
                        site.name,
                        lat,
                        lon,
                        elev,
                        float(chain_km),
                    ]
                )


def plot_all_lines(lines: dict[str, object], output: Path) -> None:
    """Plot every WILLY profile in geographic position and site order."""

    fig, ax = plt.subplots(figsize=(12.8, 8.0), constrained_layout=True)
    colors = plt.cm.tab10(np.linspace(0.0, 1.0, max(len(lines), 1)))
    all_lat = np.asarray(
        [site.coords[0] for sites in lines.values() for site in sites],
        dtype=float,
    )
    all_lon = np.asarray(
        [site.coords[1] for sites in lines.values() for site in sites],
        dtype=float,
    )
    lat0, lon0 = float(np.mean(all_lat)), float(np.mean(all_lon))
    for color, (line_name, sites) in zip(colors, lines.items()):
        lat = np.asarray([site.coords[0] for site in sites], dtype=float)
        lon = np.asarray([site.coords[1] for site in sites], dtype=float)
        east = (lon - lon0) * 111.0 * math.cos(math.radians(lat0))
        north = (lat - lat0) * 111.0
        names = [site.name for site in sites]
        ax.plot(
            east,
            north,
            "-o",
            ms=4.2,
            lw=1.5,
            color=color,
            label=f"{line_name} ({len(sites)})",
        )
        # Mark direction and label endpoints without overcrowding the map.
        ax.annotate(
            "",
            xy=(east[-1], north[-1]),
            xytext=(east[-2], north[-2]),
            arrowprops=dict(arrowstyle="->", color=color, lw=1.8),
        )
        ax.annotate(
            names[0],
            (east[0], north[0]),
            xytext=(5, -12),
            textcoords="offset points",
            fontsize=7,
        )
        ax.annotate(
            names[-1],
            (east[-1], north[-1]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=7,
        )

    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel(f"East offset from {lon0:.5f}° (km)")
    ax.set_ylabel(f"North offset from {lat0:.5f}° (km)")
    ax.set_title("WILLY AMT survey lines — coordinate-derived Sites ordering")
    ax.grid(alpha=0.25)
    ax.legend(loc="best", fontsize=8)
    _save(fig, output / "01_willy_all_station_lines.png")


def _rho_phase(
    site, component: tuple[int, int] = (0, 1)
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    z_obj = site.edi.Z
    freq = np.asarray(z_obj.freq, dtype=float)
    z = np.asarray(z_obj.z, dtype=complex)[:, component[0], component[1]]
    rho = 0.2 * np.abs(z) ** 2 / np.maximum(freq, 1e-30)
    phase = np.angle(z, deg=True)
    return freq, rho, phase


def _l22_grids(
    sites, n_periods: int = 90
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    columns = []
    for site in sites:
        freq, rho, phase = _rho_phase(site)
        period = 1.0 / np.maximum(freq, 1e-30)
        valid = (
            np.isfinite(period)
            & np.isfinite(rho)
            & np.isfinite(phase)
            & (rho > 0.0)
        )
        columns.append((period[valid], rho[valid], phase[valid]))
    period_min = max(float(np.min(col[0])) for col in columns)
    period_max = min(float(np.max(col[0])) for col in columns)
    common = np.geomspace(period_min, period_max, n_periods)
    log_rho_grid = np.empty((n_periods, len(columns)))
    phase_grid = np.empty_like(log_rho_grid)
    for index, (period, rho, phase) in enumerate(columns):
        order = np.argsort(period)
        xp = np.log10(period[order])
        x = np.log10(common)
        log_rho_grid[:, index] = np.interp(x, xp, np.log10(rho[order]))
        phase_grid[:, index] = np.interp(x, xp, phase[order])
    return common, log_rho_grid, phase_grid


def plot_l22_pseudosections(sites, output: Path) -> None:
    """Plot ordered L22 Zxy apparent resistivity and phase sections."""

    periods, log_rho, phase = _l22_grids(sites)
    chain = extract_chainage(sites)
    names = extract_station_names(sites)
    fig, axes = plt.subplots(
        2, 1, figsize=(14.0, 9.0), sharex=True, constrained_layout=True
    )
    x, y = np.meshgrid(chain, np.log10(periods))
    im_rho = axes[0].pcolormesh(x, y, log_rho, shading="auto", cmap="turbo")
    im_phi = axes[1].pcolormesh(
        x, y, phase, shading="auto", cmap="twilight_shifted"
    )
    fig.colorbar(
        im_rho, ax=axes[0], label=r"$\log_{10}\rho_{a,xy}$ ($\Omega$ m)"
    )
    fig.colorbar(im_phi, ax=axes[1], label=r"Phase $\phi_{xy}$ (degrees)")
    axes[0].set_title("L22PLT Zxy apparent-resistivity pseudosection")
    axes[1].set_title("L22PLT Zxy phase pseudosection")
    for ax in axes:
        ax.set_ylabel(r"$\log_{10}$ period (s)")
        ax.invert_yaxis()
        ax.grid(alpha=0.15)
    axes[1].set_xlabel("Profile chainage (km)")
    top = axes[0].secondary_xaxis("top")
    top.set_xticks(chain)
    top.set_xticklabels(names, rotation=90, fontsize=6.5)
    top.set_xlabel("Sites order established by ensure_sites(order_by='auto')")
    _save(fig, output / "02_l22_rho_phase_pseudosections.png")


def run_ai_inversion(
    sites,
    output: Path,
    *,
    epochs: int,
    train_profiles: int,
    n_layers: int,
) -> None:
    """Run a compact real GCN AI inversion and plot it with EDI topography."""

    chain = extract_chainage(sites)
    coords_m = np.column_stack((chain * 1000.0, np.zeros_like(chain)))
    np.random.seed(7)
    try:
        import torch

        torch.manual_seed(7)
    except ImportError:
        pass

    agent_output = output / "ai_agent_artifacts"
    agent = Inv3DAgent(
        n_layers=n_layers,
        n_freqs=16,
        n_train_profiles=train_profiles,
        epochs=epochs,
        n_mc=0,
        radius=450.0,
    )
    result = agent.execute(
        {"sites": sites, "coords": coords_m, "output_dir": str(agent_output)}
    )
    if result.status != "success":
        raise RuntimeError(result.summary)

    ax, section = plot_topo_section(
        result,
        sites=sites,
        topo_source="sites",
        kind="imshow",
        depth_max=min(25.0, float(np.max(result.data["depths_km"]))),
        model_unit="km",
        cmap="turbo",
        smooth_sigma=(0.35, 0.55),
        show_stations=True,
        show_station_names=True,
        title=(
            f"L22PLT GCN AI inversion with EDI topography — "
            f"{len(sites)} stations, RMS={result.data['rms_global']:.3g}"
        ),
        figsize=(14.0, 6.8),
        return_data=True,
    )
    _save(ax.figure, output / "03_l22_ai_inversion_with_topography.png")

    summary = {
        "status": result.status,
        "summary": result.summary,
        "elapsed_seconds": result.elapsed_seconds,
        "warnings": result.warnings,
        "n_stations": len(sites),
        "station_names": result.data["station_names"],
        "n_layers": n_layers,
        "epochs": epochs,
        "train_profiles": train_profiles,
        "rms_global": float(result.data["rms_global"]),
        "depths_km": np.asarray(result.data["depths_km"]).tolist(),
        "topography_source": section.topo_source,
        "ordering": sites.ordering,
    }
    (output / "l22_ai_inversion_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--epochs", type=int, default=3)
    parser.add_argument("--train-profiles", type=int, default=10)
    parser.add_argument("--layers", type=int, default=5)
    parser.add_argument("--skip-ai", action="store_true")
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    lines = {
        directory.name: ensure_sites(directory, recursive=False)
        for directory in _line_directories(args.data_root)
    }
    if "L22PLT" not in lines:
        raise FileNotFoundError(f"L22PLT not found below {args.data_root}")
    for name, sites in lines.items():
        print(f"{name}: {len(sites)} sites; {sites.ordering}")

    _write_station_table(lines, args.output / "willy_ordered_stations.csv")
    plot_all_lines(lines, args.output)
    plot_l22_pseudosections(lines["L22PLT"], args.output)
    if not args.skip_ai:
        run_ai_inversion(
            lines["L22PLT"],
            args.output,
            epochs=args.epochs,
            train_profiles=args.train_profiles,
            n_layers=args.layers,
        )
    print(f"validation complete: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
