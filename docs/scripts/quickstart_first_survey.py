"""Run the complete first-survey quickstart and save its fingerprint plot."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA = ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pycsamt.api import read_edis  # noqa: E402
from pycsamt.emtools import ensure_sites, plot_survey_fingerprint  # noqa: E402


def run(data_dir: Path, output_dir: Path) -> Path:
    """Load *data_dir*, print its inventory, and write one diagnostic plot."""
    survey = read_edis(
        data_dir,
        recursive=False,
        strict=True,
        on_dup="replace",
        progress=False,
    )
    if survey.errors():
        raise RuntimeError(f"EDI parser errors: {survey.errors()}")

    inventory = survey.summary().to_pandas(copy=True)
    print(f"Loaded {survey.n_sites} stations")
    print(inventory[["station", "n_freq", "tipper"]].head())

    sites = ensure_sites(
        survey.collection,
        recursive=False,
        order_by="input",
        strict=True,
    )
    fig = plot_survey_fingerprint(
        sites,
        quantities=["skew", "ellipt", "s1"],
        render="imshow",
        plot_kws={"interpolation": "bilinear"},
        station_grid=True,
        period_range=(1e-4, 1.0),
        recursive=False,
        title="L18PLT quick survey fingerprint",
        figsize=(11.2, 7.6),
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    target = output_dir / "survey_fingerprint.png"
    fig.savefig(target, dpi=180, bbox_inches="tight")
    print(f"Saved diagnostic: {target.as_posix()}")
    return target


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results") / "quickstart",
    )
    args = parser.parse_args()
    run(args.data_dir, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
