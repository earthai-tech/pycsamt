"""Generate the static figures for the Map View boreholes guide.

Uses the real Baohuashan AMT block (5 lines, 128 stations) and the
georeferenced ZK2203 borehole on line L22.

Run from the repository root::

    python docs/scripts/generate_map_boreholes_figures.py
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.app._borehole import (  # noqa: E402
    scene_borehole_traces,
    strip_log_figure,
)
from pycsamt.format.borehole import DisplayRadiusPolicy, read_pcbh  # noqa: E402
from pycsamt.map import MapView, load_lines  # noqa: E402
from pycsamt.map.borehole_align import (  # noqa: E402
    align_boreholes_to_scene,
    surface_from_sections,
)
from pycsamt.map.geometry import survey_frame, survey_uv  # noqa: E402

IMAGES = ROOT / "docs/source/images/user_guide/map"
PCBH = ROOT / "PAGEO-paper/boreholes/baohuashan_l22.pcbh.json"
EDI = ROOT / "PAGEO-paper/data/processed_edi"


def make_strip_log() -> None:
    document = read_pcbh(PCBH)
    figure = strip_log_figure(document, theme="light")
    figure.update_layout(
        width=640,
        height=760,
        title="Baohuashan L22 drill holes — logged lithology",
    )
    figure.write_image(str(IMAGES / "map_boreholes_strip_log.png"), scale=2)


def make_3d_scene() -> None:
    document = read_pcbh(PCBH)
    data = load_lines(EDI, detect="folder", recursive=True)
    view = MapView(data)
    figure = view.map3d(
        mode="fence",
        quantity="resistivity",
        depth_range=(0.0, 1500.0),
        topography=True,
        show_terrain=True,
        opacity=0.35,
    )

    ids, lats, lons, lines, elevs = [], [], [], [], []
    for s in data.stations:
        if s.latitude is None or s.longitude is None:
            continue
        ids.append(str(s.id))
        lats.append(float(s.latitude))
        lons.append(float(s.longitude))
        lines.append(s.line or "line")
        elevs.append(float(s.elevation) if s.elevation is not None else 0.0)
    frame = survey_frame(lats, lons, lines)
    uv = survey_uv(ids, lats, lons, lines)
    surface = surface_from_sections([([uv[i][0] for i in ids], elevs)])
    alignment = align_boreholes_to_scene(
        document,
        frame,
        datum="surface",
        surface=surface,
        radius_policy=DisplayRadiusPolicy(mode="fixed", fixed_radius=28.0),
    )
    for trace in scene_borehole_traces(alignment, as_tubes=True, opacity=1.0):
        figure.add_trace(trace)
    figure.update_layout(
        width=900,
        height=620,
        title="ZK221 / ZK222 / ZK2203 on line L22 of the Baohuashan block",
        scene_camera={"eye": {"x": 1.5, "y": -1.6, "z": 0.9}},
    )
    figure.write_image(str(IMAGES / "map_boreholes_3d_scene.png"), scale=2)


def main() -> None:
    IMAGES.mkdir(parents=True, exist_ok=True)
    make_strip_log()
    print("wrote map_boreholes_strip_log.png")
    make_3d_scene()
    print("wrote map_boreholes_3d_scene.png")


if __name__ == "__main__":
    main()
