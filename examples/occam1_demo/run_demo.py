"""Build, invert, and plot the three bundled Occam1D EDI soundings.

Run from any directory with::

    python examples/occam1_demo/run_demo.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

from pycsamt.models.occam1d import (
    Occam1DConfig,
    Occam1DInputBuilder,
    Occam1DInversion,
    Occam1DSolverPolicy,
)

matplotlib.use("Agg")


ROOT = Path(__file__).resolve().parents[2]
EDI_DIRECTORY = ROOT / "data" / "3edis"
DEFAULT_OUTPUT = Path(__file__).resolve().parent / "output"


def _arguments():
    parser = argparse.ArgumentParser(
        description="Run the native pyCSAMT Occam1D EDI demonstration."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Demo output root (default: examples/occam1_demo/output).",
    )
    parser.add_argument(
        "--max-iterations",
        type=int,
        default=20,
        help="Absolute nonlinear iteration budget per station.",
    )
    parser.add_argument(
        "--target-rms",
        type=float,
        default=1.0,
        help="Target normalized RMS per station.",
    )
    parser.add_argument(
        "--station",
        action="append",
        default=None,
        help="EDI stem to run; repeat for several. Default: all bundled EDIs.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Disable iterative standard-output progress.",
    )
    return parser.parse_args()


def _sources(selected):
    paths = sorted(EDI_DIRECTORY.glob("*.edi"))
    if selected:
        wanted = {name.casefold() for name in selected}
        paths = [path for path in paths if path.stem.casefold() in wanted]
    if not paths:
        raise FileNotFoundError(
            f"No selected EDI files were found in {EDI_DIRECTORY}."
        )
    return paths


def run_station(source, output, args):
    """Run one EDI through build, inversion, export, and plotting."""
    station_root = output / source.stem
    inversion_directory = station_root / "occam1d-inversion"
    text_directory = station_root / "model-text"
    image_directory = station_root / "model-image"
    config = Occam1DConfig(
        mode="determinant",
        error_floor_rho=0.05,
        error_floor_phase=2.0,
        n_layers=30,
        first_thickness=5.0,
        depth_max=5000.0,
        starting_resistivity=100.0,
        target_misfit=args.target_rms,
        max_iterations=args.max_iterations,
        roughness_type=1,
        stepsize_cut_count=8,
        lagrange_start=5.0,
        lagrange_min=1e-10,
        lagrange_max=1e10,
    )
    builder = Occam1DInputBuilder(
        source,
        inversion_directory,
        config,
        verbose=not args.quiet,
    ).build(overwrite=True)
    policy = Occam1DSolverPolicy(
        condition_limit=1e12,
        ill_condition_action="damp",
        svd_failure_action="fallback",
        max_retries=3,
    )
    inversion = Occam1DInversion(
        builder.data,
        builder.model,
        config=builder.config,
        startup=builder.startup,
        log_resistivity_bounds=(0.0, 6.0),
        lagrange_bounds=(config.lagrange_min, config.lagrange_max),
        solver_policy=policy,
        verbose=0 if args.quiet else 1,
    )
    result = inversion.run()
    restart = inversion.restart(result)
    restart_path = restart.write(inversion_directory / "native-restart.json")
    text_paths = inversion.export_result(text_directory, result)
    image_paths = inversion.save_main_images(
        image_directory,
        result,
        dpi=180,
    )
    print(inversion.result_summary(result))
    return {
        "source": str(source.relative_to(ROOT)),
        "station": builder.data.station,
        "status": result.convergence.value,
        "converged": result.converged,
        "iterations": result.n_iterations,
        "initial_rms": result.initial.rms,
        "final_rms": result.final.rms,
        "target_rms": result.target_rms,
        "restart": str(restart_path.relative_to(output)),
        "native_inputs": {
            key: str(Path(value).relative_to(output))
            for key, value in builder.paths.to_dict().items()
            if key != "workdir"
        },
        "model_text": {
            key: str(path.relative_to(output))
            for key, path in text_paths.items()
        },
        "model_image": {
            key: str(path.relative_to(output))
            for key, path in image_paths.items()
        },
    }


def main():
    """Run selected bundled soundings and write the demo summary."""
    args = _arguments()
    output = args.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    records = [
        run_station(source, output, args)
        for source in _sources(args.station)
    ]
    summary = output / "demo-summary.json"
    summary.write_text(
        json.dumps(
            {
                "schema": "pycsamt.occam1d.demo/v1",
                "source_directory": str(EDI_DIRECTORY.relative_to(ROOT)),
                "stations": records,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf8",
    )
    print(f"Demo completed: {summary}")


if __name__ == "__main__":
    main()
