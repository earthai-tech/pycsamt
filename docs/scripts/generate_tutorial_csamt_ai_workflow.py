"""Run and plot the Tongkeng triangular-mesh Maxwell-AI workflow used by the
CSAMT groundwater-geology tutorial.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.agents import Inv2DAgent  # noqa: E402
from pycsamt.emtools import ensure_sites  # noqa: E402
from pycsamt.forward.maxwell.tri_fem2d import TriFEM2DAdapter  # noqa: E402
from pycsamt.site.edit import select_freq_all  # noqa: E402
from pycsamt.topo.extract import extract_elevation  # noqa: E402

IMAGE_DIR = ROOT / "docs/source/images/tutorials/map_groundwater_geology_from_csamt"
RUN_DIR = ROOT / "runs"
RAW_DIR = ROOT / "data" / "CSAMT"
CORRECTED_DIR = RUN_DIR / "csamt_corrected_edis"
FREQUENCIES_HZ = [8196.722, 4098.361, 2049.18, 1023.541, 512.8206, 255.7545]
STATION_SPACING_M = 50.0


def load_topography():
    """Real per-station elevation, extracted from the *raw* CSAMT delivery.

    Must come from the original EDIs, not the corrected/exported copies
    used by ``load_line()``: ``EDIExportAgent`` rewrites the ``>INFO``
    block without this delivery's non-standard ``LATITUDE=``/
    ``LONGITUDE=``/``ELEVATION=`` keys, so elevation does not survive the
    export round trip -- confirmed empirically (the exported files carry
    an explicit ``ELEV=0.0`` and no elevation-bearing ``>INFO`` fields at
    all). Elevation is a static per-station property untouched by any of
    the frequency-domain corrections in between, so extracting it once
    here, independently of the corrected-EDI pipeline, is safe.
    """
    raw_sites = ensure_sites(RAW_DIR, recursive=False, verbose=0).ordered("station")
    station_x = np.array([STATION_SPACING_M * i for i in range(len(list(raw_sites)))])
    elevation = extract_elevation(raw_sites)
    topo_z = elevation.max() - elevation
    return station_x, topo_z


def save(fig: plt.Figure, name: str) -> None:
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(IMAGE_DIR / name, dpi=190, bbox_inches="tight")
    plt.close(fig)


def load_line():
    if not CORRECTED_DIR.exists():
        raise FileNotFoundError(
            f"{CORRECTED_DIR} is missing; run the tutorial's "
            "'Drop Weak Frequencies And Export Corrected EDIs' step first."
        )
    sites = ensure_sites(CORRECTED_DIR, recursive=False, verbose=0).ordered(
        "station"
    )
    return select_freq_all(sites, fmin=255.0)


def run_inversion(
    sites,
    *,
    production: bool = False,
    n_train_profiles: int | None = None,
    epochs: int | None = None,
):
    np.random.seed(0)
    try:
        import torch

        torch.manual_seed(0)
    except ImportError:
        pass

    profiles = int(
        n_train_profiles if n_train_profiles is not None else (300 if production else 40)
    )
    epoch_budget = int(epochs if epochs is not None else (60 if production else 15))
    output_dir = RUN_DIR / (
        "csamt_ai2d_tri_production" if production else "csamt_ai2d_tri"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    topo_x_m, topo_z_m = load_topography()

    started = time.time()
    agent = Inv2DAgent(
        physics="mt2d_tri",
        epochs=epoch_budget,
        n_freqs=len(FREQUENCIES_HZ),
        depth_max=300.0,
        n_train_profiles=profiles,
        n_stations_per_profile=len(list(sites)),
        station_spacing_m=STATION_SPACING_M,
        mesh_target_cell_m=20.0 if production else 30.0,
        field_grid_cell_m=10.0 if production else 15.0,
        correlation_length_x_m=(60.0, 200.0),
        correlation_length_z_m=(20.0, 80.0),
        gcn_hidden=(128, 64, 32) if production else (64, 32, 16),
        topo_x_m=topo_x_m,
        topo_z_m=topo_z_m,
        mare2dem_adapter=TriFEM2DAdapter(),
    )
    result = agent.execute(
        {
            "sites": sites,
            "freqs": FREQUENCIES_HZ,
            "output_dir": str(output_dir),
        }
    )
    if result.status != "success":
        raise RuntimeError(f"csamt: {result.summary}")

    pred = result.data["pred_triangles"]
    np.savez_compressed(
        output_dir / "field_prediction.npz",
        log10_resistivity=np.asarray(pred["log10_resistivity"]),
        n_triangles=pred["mesh"].n_triangles,
        frequencies_hz=np.asarray(FREQUENCIES_HZ),
    )
    recovery = result.data.get("mt2d_tri_recovery") or {}
    manifest = {
        "status": result.status,
        "summary": result.summary,
        "production": production,
        "seed": 0,
        "n_train_profiles": profiles,
        "epochs_requested": epoch_budget,
        "n_stations": len(list(sites)),
        "n_triangles": pred["mesh"].n_triangles,
        "frequencies_hz": FREQUENCIES_HZ,
        "mt2d_tri_recovery": recovery,
        "elapsed_seconds": time.time() - started,
    }
    (output_dir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, default=float), encoding="utf-8"
    )
    print(json.dumps(manifest, default=float), flush=True)
    return result


def plot_result(result, *, production: bool) -> None:
    fig_path = result.data.get("figure_paths", {}).get("inv2d_tri_section")
    if not fig_path:
        return
    import shutil

    name = (
        "csamt_ai2d_tri_section_production.png"
        if production
        else "csamt_ai2d_tri_section.png"
    )
    shutil.copy(fig_path, IMAGE_DIR / name)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--production", action="store_true")
    parser.add_argument("--n-train-profiles", type=int)
    parser.add_argument("--epochs", type=int)
    args = parser.parse_args(argv)

    sites = load_line()
    result = run_inversion(
        sites,
        production=args.production,
        n_train_profiles=args.n_train_profiles,
        epochs=args.epochs,
    )
    plot_result(result, production=args.production)

    pred = result.data["pred_triangles"]
    recovery = result.data.get("mt2d_tri_recovery") or {}
    print(
        "Tongkeng CSAMT",
        "stations", len(list(sites)),
        "triangles", pred["mesh"].n_triangles,
        "recovery_RMSE", f"{recovery.get('rmse', float('nan')):.3f}",
        "recovery_n", recovery.get("n_samples", 0),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
