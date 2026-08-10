"""Generate figures and sample outputs for the read-EDI tutorial."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
ONE_EDI = DATA_DIR / "18-001A.edi"
IMAGE_DIR = (
    ROOT / "docs" / "source" / "images" / "tutorials" / "read_edi_survey"
)


def _import_pycsamt():
    """Import the public APIs used by this standalone generator."""
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))

    from pycsamt.api import read_edi, read_edis
    from pycsamt.emtools import (
        plot_survey_fingerprint,
        plot_survey_inventory_overview,
    )

    return {
        "plot_survey_fingerprint": plot_survey_fingerprint,
        "plot_survey_inventory_overview": plot_survey_inventory_overview,
        "read_edi": read_edi,
        "read_edis": read_edis,
    }


def _style_axis(ax: plt.Axes) -> None:
    ax.set_facecolor("#fbfbf7")
    ax.grid(True, color="#d8dad1", linewidth=0.7, alpha=0.65)
    for spine in ax.spines.values():
        spine.set_color("#39434d")
        spine.set_linewidth(0.8)


def _save(fig: plt.Figure, name: str) -> None:
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(IMAGE_DIR / name, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _inventory_plot(functions, survey) -> None:
    labels = [name.replace("23-", "") for name in survey.stations]
    fig = functions["plot_survey_inventory_overview"](
        survey.collection,
        station_labels=labels,
        recursive=False,
        title="L18PLT acquisition inventory and period coverage",
        figsize=(11.2, 7.2),
    )
    _save(fig, "survey_inventory_overview.png")


def _path_plot(inventory) -> None:
    names = inventory["path"].map(lambda value: Path(value).name)
    lengths = names.str.len().to_numpy(dtype=int)
    fig, ax = plt.subplots(figsize=(9.0, 3.8))
    ax.plot(
        np.arange(1, len(names) + 1),
        lengths,
        color="#7c4d79",
        marker="s",
        markersize=4,
        linewidth=1.2,
    )
    ax.set_xlabel("Loaded file index")
    ax.set_ylabel("Filename length")
    ax.set_title("Source filenames are regular across the line")
    _style_axis(ax)
    _save(fig, "source_filename_check.png")


def _fingerprint(functions, survey) -> None:
    fig = functions["plot_survey_fingerprint"](
        survey.collection,
        quantities=["skew", "ellipt", "s1"],
        render="imshow",
        plot_kws={"interpolation": "bilinear"},
        station_grid=True,
        period_range=(1e-4, 1.0),
        recursive=False,
        title="L18PLT quick survey fingerprint",
        figsize=(11.2, 7.6),
    )
    _save(fig, "survey_fingerprint.png")


def main() -> int:
    functions = _import_pycsamt()
    survey = functions["read_edis"](
        DATA_DIR,
        recursive=False,
        strict=False,
        on_dup="replace",
        progress=False,
    )
    summary = survey.summary()
    inventory = summary.to_pandas(copy=True)
    one = functions["read_edi"](ONE_EDI)
    single = functions["read_edis"](ONE_EDI, recursive=False, progress=False)

    _inventory_plot(functions, survey)
    _path_plot(inventory)
    # Some processing-oriented plotting paths normalize station metadata on
    # their input collection.  Use a separate load so figure generation cannot
    # alter the survey used for captured inventory and selection output.
    plot_survey = functions["read_edis"](
        DATA_DIR,
        recursive=False,
        strict=False,
        on_dup="replace",
        progress=False,
    )
    _fingerprint(functions, plot_survey)

    print("survey_text:")
    print(survey)
    print("survey_properties:")
    print(f"n_sites={survey.n_sites}")
    print(f"stations={survey.stations[:5]}")
    print(f"paths={[Path(p).name for p in survey.paths[:5]]}")
    print("summary:")
    print(summary)
    print("inventory_head:")
    compact = inventory.copy()
    compact["path"] = compact["path"].map(lambda value: Path(value).name)
    print(compact.head(5).to_string(index=False))
    print(f"columns={inventory.columns.tolist()}")
    print("single_file:")
    print(f"read_edi station={one.station}")
    print(f"read_edis n_sites={single.n_sites}")
    print(single.summary())
    print("selection:")
    first = survey[0]
    selected = survey.get_site(survey.stations[0])
    print(f"first={first.station} {Path(first.path).name}")
    print(f"selected={selected.station} {Path(selected.path).name}")
    print(f"missing={survey.get_site('S999')}")
    print("errors:")
    print(f"{len(survey.errors())} read error(s)")
    print("focused_table:")
    focused = survey.summary(
        fields=["station", "n_freq", "tipper", "path"]
    ).to_pandas(copy=True)
    focused["path"] = focused["path"].map(lambda value: Path(value).name)
    print(focused.head(3).to_string(index=False))
    print(f"images: {IMAGE_DIR.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
