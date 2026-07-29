"""Run the experimental graph candidate without publishing it as 3-D inversion.

The command writes a machine-readable gate report into a quarantine directory.
It intentionally never promotes or renders a subsurface section: the current
Inv3DAgent uses graph learning with tiled 1-D MT physics, not a 3-D Maxwell
forward operator.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.agents import Inv3DAgent
from pycsamt.emtools import ensure_sites


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("edi_directory", type=Path)
    parser.add_argument("--output", type=Path, default=Path("runs/inv3d_candidate"))
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--profiles", type=int, default=150)
    parser.add_argument("--radius", type=float, default=300.0)
    parser.add_argument("--depth-max", type=float, default=2000.0)
    parser.add_argument("--n-mc", type=int, default=20)
    parser.add_argument("--rms-max", type=float, default=2.0)
    parser.add_argument("--rho-min", type=float, default=1.0)
    parser.add_argument("--rho-max", type=float, default=1.0e5)
    parser.add_argument("--max-outside-fraction", type=float, default=0.0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    sites = ensure_sites(args.edi_directory, recursive=True, verbose=0).ordered()
    agent = Inv3DAgent(
        n_layers=6,
        epochs=args.epochs,
        n_train_profiles=args.profiles,
        n_mc=args.n_mc,
        radius=args.radius,
        depth_max=args.depth_max,
    )
    result = agent.execute({"sites": sites, "topography": True})
    if result.status == "failed":
        raise RuntimeError(result.error)

    log10_rho = np.asarray(result.data["pred_rho"], dtype=float)
    uncertainty = np.asarray(result.data["pred_uncertainty"], dtype=float)
    lo, hi = np.log10([args.rho_min, args.rho_max])
    in_bounds = np.isfinite(log10_rho) & (log10_rho >= lo) & (log10_rho <= hi)
    outside_fraction = float(1.0 - in_bounds.mean())
    rms = float(result.data["rms_global"])

    numerical_gates = {
        "finite_prediction": bool(np.isfinite(log10_rho).all()),
        "finite_uncertainty": bool(np.isfinite(uncertainty).all()),
        "rms": bool(np.isfinite(rms) and rms <= args.rms_max),
        "physical_bounds": bool(outside_fraction <= args.max_outside_fraction),
    }
    report = {
        "execution_status": result.status,
        "scientific_release": "rejected",
        "reason": (
            "The current candidate uses a station graph and tiled 1-D MT "
            "forward responses; it is not validated 3-D Maxwell inversion."
        ),
        "forward_physics": "tiled_mt1d_graph",
        "n_stations": len(sites),
        "measurements": {
            "rms_global": rms,
            "log10_rho_min": float(np.nanmin(log10_rho)),
            "log10_rho_max": float(np.nanmax(log10_rho)),
            "outside_physical_bounds_fraction": outside_fraction,
            "uncertainty_min": float(np.nanmin(uncertainty)),
            "uncertainty_max": float(np.nanmax(uncertainty)),
        },
        "thresholds": {
            "rms_max": args.rms_max,
            "rho_min_ohm_m": args.rho_min,
            "rho_max_ohm_m": args.rho_max,
            "max_outside_fraction": args.max_outside_fraction,
        },
        "numerical_gates": numerical_gates,
        "all_numerical_gates_pass": all(numerical_gates.values()),
        "physics_gate": False,
    }
    report_path = args.output / "candidate_gate.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"candidate directory: {args.output.resolve()}")
    print(f"gate report: {report_path.resolve()}")
    print("scientific release: rejected")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
