"""Generate static figures for the MapView user guide."""

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.map import MapView, StationMapOptions, save_png


IMAGES = ROOT / "docs/source/images/user_guide/map"
EDI = ROOT / "data/AMT/WILLY_DATA"


def make_elevation_contour_map() -> None:
    """Render station elevations as contours while retaining field support."""
    view = MapView.from_folder(EDI, detect="folder", recursive=True)
    figure = view.station(
        options=StationMapOptions(
            backend="matplotlib",
            overlay="elevation",
            elevation_mode="contours",
            contour_mode="filled+lines",
            contour_levels=16,
            contour_opacity=0.82,
            contour_interp="linear",
            contour_smooth=0.5,
            contour_grid_res=220,
            cmap="terrain",
            marker_size=4,
            show_labels=True,
            label_fontsize=5.2,
            label_rotation=28.0,
            show_profiles=True,
            title="WILLY_DATA elevation contours and survey support",
        )
    )
    figure.set_size_inches(12.0, 7.4)
    axis = figure.axes[0]
    axis.grid(color="#94a3b8", linestyle=":", linewidth=0.65, alpha=0.5)
    axis.set_axisbelow(True)
    save_png(
        figure,
        IMAGES / "user-guide-map-mapview-01.png",
        dpi=190,
        bbox_inches="tight",
    )


if __name__ == "__main__":
    IMAGES.mkdir(parents=True, exist_ok=True)
    make_elevation_contour_map()
