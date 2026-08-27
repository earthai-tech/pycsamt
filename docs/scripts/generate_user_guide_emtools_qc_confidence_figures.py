"""Generate confidence figures for the EMTools QC user guide.

Run from the repository root::

    python docs/scripts/generate_user_guide_emtools_qc_confidence_figures.py \
        --save-dir docs/source/images/user_guide/emtools/confidence

Use ``--save-dir DIR`` to save both figures instead of opening them.
The single-line example uses L18PLT. The contour example uses all five
WILLY_DATA profiles, because interpolation from one line is not a valid
two-dimensional confidence map.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from pycsamt.emtools import (
    ensure_sites,
    export_confidence_map,
    plot_confidence_before_after,
    plot_confidence_component_map,
    plot_confidence_coverage_curve,
    plot_confidence_distribution,
    plot_confidence_grid_map,
    plot_confidence_heatmap,
    plot_confidence_map,
    plot_confidence_method_comparison,
    plot_confidence_rank,
    plot_confidence_risk_map,
    station_confidence_table,
)


def _repository_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _line_labels(sites) -> dict[str, str]:
    table = station_confidence_table(sites, method="composite", api=False)
    return {
        str(station): f"L{str(station).split('-', 1)[0]}"
        for station in table["station"]
    }


def make_figures():
    """Return route, filled-contour, and labelled-isoline map axes."""
    data_root = _repository_root() / "data" / "AMT" / "WILLY_DATA"
    line18 = data_root / "L18PLT"

    route_ax = plot_confidence_map(
        line18,
        method="composite",
        mode="route",
        station_labels=True,
        station_label_step=4,
        show_confidence_values=True,
        confidence_value_step=2,
        confidence_value_fmt="{:.2f}",
        confidence_value_fontsize=6.5,
        recursive=False,
    )
    route_ax.set_title("L18 confidence along the surveyed route")

    all_lines = ensure_sites(data_root, recursive=True)
    contour_ax = plot_confidence_map(
        all_lines,
        method="composite",
        mode="contour",
        line_labels=_line_labels(all_lines),
        boundary_levels=[0.50, 0.85, 0.90, 1.00],
        show_contour_lines=False,
        show_threshold_contours=True,
        threshold_linewidth=1.6,
    )
    contour_ax.set_title("WILLY_DATA confidence map — boundary contours only")

    line_ax = plot_confidence_map(
        all_lines,
        method="composite",
        mode="contour",
        line_labels=_line_labels(all_lines),
        show_contour_lines=True,
        contour_line_levels=[0.55, 0.60, 0.65, 0.70, 0.75, 0.80],
        contour_line_colors="#303030",
        contour_linewidths=0.8,
        contour_linestyles="solid",
        contour_labels=True,
        contour_label_fmt="%.2f",
        contour_label_fontsize=6.5,
        threshold_line_color="black",
        threshold_linewidth=1.6,
        threshold_linestyle="--",
    )
    line_ax.set_title("WILLY_DATA confidence contours with labelled isolines")
    return route_ax, contour_ax, line_ax


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--save-dir", type=Path)
    parser.add_argument(
        "--export-dir",
        type=Path,
        help="Write station CSV and Surfer DSAA confidence grid.",
    )
    args = parser.parse_args()
    axes = make_figures()
    data_root = _repository_root() / "data" / "AMT" / "WILLY_DATA"
    all_lines = ensure_sites(data_root, recursive=True)
    component_fig = plot_confidence_component_map(
        all_lines,
        method="composite",
        line_labels=_line_labels(all_lines),
        ncols=4,
        marker_size=32.0,
    )
    comparison_fig = plot_confidence_method_comparison(
        all_lines,
        line_labels=_line_labels(all_lines),
        marker_size=34.0,
        map_aspect="auto",
    )
    risk_ax = plot_confidence_risk_map(
        all_lines,
        method="composite",
        line_labels=_line_labels(all_lines),
        mode="contour",
        risk_levels=[0.0, 0.05, 0.10, 0.15, 0.25, 0.35, 0.50, 0.75, 1.0],
        station_labels=False,
        map_aspect="auto",
    )
    heatmap_ax = plot_confidence_heatmap(
        all_lines,
        method="composite",
        line_labels=_line_labels(all_lines),
        station_order="route",
        annotate=False,
        station_label_step=4,
    )
    distribution_fig = plot_confidence_distribution(
        all_lines,
        method="both",
        line_labels=_line_labels(all_lines),
        bins=18,
    )
    rank_fig = plot_confidence_rank(
        all_lines,
        method="composite",
        line_labels=_line_labels(all_lines),
        order="worst",
        annotate_stations=False,
        show_iqr=True,
    )
    grid_ax = plot_confidence_grid_map(
        all_lines,
        method="composite",
        line_labels=_line_labels(all_lines),
        grid_shape=(220, 180),
        interpolation="linear",
        show_grid_edges=False,
        map_aspect="auto",
    )
    coverage_fig = plot_confidence_coverage_curve(
        all_lines,
        method="composite",
        line_labels=_line_labels(all_lines),
    )
    before_after_fig = plot_confidence_before_after(
        all_lines,
        before_method="presence",
        after_method="composite",
        before_label="Presence",
        after_label="Composite",
        line_labels=_line_labels(all_lines),
        show_station_labels=False,
    )
    if args.save_dir is None:
        plt.show()
        return
    args.save_dir.mkdir(parents=True, exist_ok=True)
    names = (
        "confidence_route",
        "confidence_contour",
        "confidence_contour_lines",
    )
    for name, ax in zip(names, axes):
        ax.figure.savefig(
            args.save_dir / f"{name}.png",
            dpi=180,
            bbox_inches="tight",
        )
    component_fig.savefig(
        args.save_dir / "confidence_component_map.png",
        dpi=220,
        bbox_inches="tight",
    )
    comparison_fig.savefig(
        args.save_dir / "confidence_method_comparison.png",
        dpi=220,
        bbox_inches="tight",
    )
    risk_ax.figure.savefig(
        args.save_dir / "confidence_risk_map.png",
        dpi=220,
        bbox_inches="tight",
    )
    heatmap_ax.figure.savefig(
        args.save_dir / "confidence_heatmap.png",
        dpi=220,
        bbox_inches="tight",
    )
    distribution_fig.savefig(
        args.save_dir / "confidence_distribution.png",
        dpi=220,
        bbox_inches="tight",
    )
    rank_fig.savefig(
        args.save_dir / "confidence_rank.png",
        dpi=220,
        bbox_inches="tight",
    )
    grid_ax.figure.savefig(
        args.save_dir / "confidence_grid_map.png",
        dpi=220,
        bbox_inches="tight",
    )
    coverage_fig.savefig(
        args.save_dir / "confidence_coverage_curve.png",
        dpi=220,
        bbox_inches="tight",
    )
    before_after_fig.savefig(
        args.save_dir / "confidence_before_after.png",
        dpi=220,
        bbox_inches="tight",
    )
    if args.export_dir is not None:
        args.export_dir.mkdir(parents=True, exist_ok=True)
        export_confidence_map(
            all_lines,
            method="composite",
            line_labels=_line_labels(all_lines),
            csv_path=args.export_dir / "confidence_map.csv",
            surfer_path=args.export_dir / "confidence_map.grd",
            grid_shape=(200, 200),
        )


if __name__ == "__main__":
    main()
