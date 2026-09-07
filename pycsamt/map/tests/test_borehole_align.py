"""Tests for placing PCBH boreholes in the 3-D map scene."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.format.borehole import (
    Collar,
    CoordinateReferenceSystem,
    LogInterval,
    PCBHBorehole,
    PCBHDocument,
    SurveyStation,
    Trajectory,
    VocabularyEntry,
)
from pycsamt.map.borehole_align import (
    align_boreholes_to_scene,
    surface_from_sections,
)
from pycsamt.map.geometry import survey_frame, survey_uv


def _survey():
    """Two E-W lines of stations around 22N, 103E."""
    ids, lats, lons, lines = [], [], [], []
    for line_idx, lat in enumerate((22.000, 22.002)):
        for j in range(6):
            ids.append(f"L{line_idx}-{j}")
            lats.append(lat)
            lons.append(103.000 + 0.001 * j)
            lines.append(f"L{line_idx}")
    return ids, lats, lons, lines


def _document(lat, lon, *, crs="EPSG:4326"):
    collar_kwargs = {}
    if crs == "EPSG:4326":
        collar = Collar(lon, lat, 500.0, longitude=lon, latitude=lat)
    else:
        collar = Collar(lon, lat, 500.0)
    return PCBHDocument(
        document_id="test:align",
        created_at="2026-09-03T00:00:00Z",
        created_by="pytest",
        crs=CoordinateReferenceSystem(crs),
        lithologies=[
            VocabularyEntry("A", "Overburden", color="#cccccc"),
            VocabularyEntry("B", "Granite", color="#884422"),
        ],
        boreholes=[
            PCBHBorehole(
                id="BH1",
                name="Borehole 1",
                kind="mining_exploration",
                status="completed",
                collar=collar,
                total_depth_md=120.0,
                interval_logs={
                    "lithology": [
                        LogInterval(0.0, 40.0, code="A"),
                        LogInterval(40.0, 120.0, code="B"),
                    ]
                },
            )
        ],
    )


def test_collar_at_a_known_station_lands_at_that_station_uv():
    ids, lats, lons, lines = _survey()
    frame = survey_frame(lats, lons, lines)
    uv = survey_uv(ids, lats, lons, lines)
    target = uv["L0-3"]

    document = _document(lats[3], lons[3])
    alignment = align_boreholes_to_scene(document, frame, datum="zero")
    hole = alignment.placed[0]

    assert hole.collar_scene[0] == pytest.approx(target[0], abs=1.0)
    assert hole.collar_scene[1] == pytest.approx(target[1], abs=1.0)
    # z=0 datum: collar at surface, trajectory runs downward
    assert hole.collar_scene[2] == pytest.approx(0.0)
    assert hole.centerline[0][2] == pytest.approx(0.0)
    assert hole.centerline[-1][2] == pytest.approx(-120.0, abs=1e-6)


def test_surface_datum_drapes_collar_onto_terrain_line():
    ids, lats, lons, lines = _survey()
    frame = survey_frame(lats, lons, lines)
    u_station = np.array([survey_uv(ids, lats, lons, lines)[i][0] for i in ids])
    elevations = 400.0 + 0.05 * u_station  # a gentle ramp
    surface = surface_from_sections([(u_station, elevations)])

    document = _document(lats[3], lons[3])
    alignment = align_boreholes_to_scene(
        document, frame, datum="surface", surface=surface
    )
    hole = alignment.placed[0]
    expected_top = surface(hole.collar_scene[0])
    assert hole.collar_scene[2] == pytest.approx(expected_top, abs=1e-6)
    assert hole.centerline[-1][2] == pytest.approx(
        expected_top - 120.0, abs=1e-6
    )


def test_segments_carry_vocabulary_colours_and_bounds():
    ids, lats, lons, lines = _survey()
    frame = survey_frame(lats, lons, lines)
    document = _document(lats[0], lons[0])
    alignment = align_boreholes_to_scene(document, frame, datum="zero")
    hole = alignment.placed[0]
    assert [seg.color for seg in hole.segments] == ["#cccccc", "#884422"]
    assert hole.segments[0].from_md == 0.0
    assert hole.segments[1].to_md == 120.0


def test_depth_range_clips_the_trajectory():
    ids, lats, lons, lines = _survey()
    frame = survey_frame(lats, lons, lines)
    document = _document(lats[0], lons[0])
    alignment = align_boreholes_to_scene(
        document, frame, datum="zero", depth_range=(0.0, 50.0)
    )
    hole = alignment.placed[0]
    assert all(point[2] >= -50.0 - 1e-6 for point in hole.centerline)


def test_local_crs_without_lonlat_is_unplaced():
    ids, lats, lons, lines = _survey()
    frame = survey_frame(lats, lons, lines)
    document = _document(0.0, 0.0, crs="LOCAL:mine-grid")
    alignment = align_boreholes_to_scene(document, frame)
    assert alignment.placed == ()
    assert alignment.holes[0].relation == "unplaced"
    assert not alignment.holes[0].placed


def test_offset_shift_and_scale_match_the_volume_line_panels():
    """A hole on a line must land on that line's normalised panel.

    ``pycsamt.map.volume`` draws each line panel at
    ``(median cross-strike v - front-most line's median) * line_spacing``.
    ``align_boreholes_to_scene`` must reproduce that shift + stretch or
    an inserted hole floats in front of / behind its line.
    """
    ids, lats, lons, lines = _survey()
    frame = survey_frame(lats, lons, lines)
    uv = survey_uv(ids, lats, lons, lines)

    line_v = {}
    for sid, ln in zip(ids, lines):
        line_v.setdefault(ln, []).append(uv[sid][1])
    medians = {ln: float(np.median(v)) for ln, v in line_v.items()}
    shift = min(medians.values())
    spacing = 2.5

    # borehole planted at a station on the back line (L1)
    idx = ids.index("L1-3")
    document = _document(lats[idx], lons[idx])
    alignment = align_boreholes_to_scene(
        document,
        frame,
        datum="zero",
        offset_shift=shift,
        offset_scale=spacing,
    )
    hole = alignment.placed[0]
    expected_y = (medians["L1"] - shift) * spacing
    assert hole.collar_scene[1] == pytest.approx(expected_y, abs=1.0)
    # and the whole trajectory sits on that panel (vertical hole)
    assert all(
        point[1] == pytest.approx(expected_y, abs=1.0)
        for point in hole.centerline
    )


def test_lean_tips_a_vertical_hole_and_shortens_its_reach():
    ids, lats, lons, lines = _survey()
    frame = survey_frame(lats, lons, lines)
    document = _document(lats[3], lons[3])

    plumb = align_boreholes_to_scene(
        document, frame, datum="zero"
    ).placed[0]
    leaned = align_boreholes_to_scene(
        document, frame, datum="zero", lean_deg=30.0, lean_azimuth_deg=90.0
    ).placed[0]

    # collar is unmoved; the toe has swung out and risen
    assert leaned.collar_scene == pytest.approx(plumb.collar_scene)
    plumb_toe = plumb.centerline[-1]
    leaned_toe = leaned.centerline[-1]
    horizontal = float(
        np.hypot(
            leaned_toe[0] - plumb_toe[0],
            leaned_toe[1] - plumb_toe[1],
        )
    )
    assert horizontal == pytest.approx(120.0 * np.sin(np.radians(30)), abs=1.0)
    assert leaned_toe[2] == pytest.approx(
        -120.0 * np.cos(np.radians(30)), abs=1.0
    )
    assert "leaning 30" in " ".join(leaned.warnings)


def test_lean_leaves_a_surveyed_deviated_hole_alone():
    ids, lats, lons, lines = _survey()
    frame = survey_frame(lats, lons, lines)
    # a real deviation survey (inclination from vertical) -> skip it
    collar = Collar(
        lons[3], lats[3], 500.0, longitude=lons[3], latitude=lats[3]
    )
    document = PCBHDocument(
        document_id="test:align-dev",
        created_at="2026-09-04T00:00:00Z",
        created_by="pytest",
        crs=CoordinateReferenceSystem("EPSG:4326"),
        lithologies=[VocabularyEntry("A", "Overburden", color="#cccccc")],
        boreholes=[
            PCBHBorehole(
                id="BH1",
                name="Deviated 1",
                kind="mining_exploration",
                status="completed",
                collar=collar,
                total_depth_md=120.0,
                trajectory=Trajectory(
                    method="survey",
                    stations=[
                        SurveyStation(0.0, 40.0, 0.0),
                        SurveyStation(120.0, 40.0, 35.0),
                    ],
                ),
                interval_logs={"lithology": [LogInterval(0.0, 120.0, "A")]},
            )
        ],
    )
    leaned = align_boreholes_to_scene(
        document, frame, datum="zero", lean_deg=30.0, lean_azimuth_deg=90.0
    ).placed[0]
    assert not any("leaning" in w for w in leaned.warnings)


def test_relation_marks_a_far_collar_outside():
    ids, lats, lons, lines = _survey()
    frame = survey_frame(lats, lons, lines)
    document = _document(22.05, 103.05)  # ~5 km away
    alignment = align_boreholes_to_scene(
        document,
        frame,
        datum="zero",
        scene_bounds=(-500.0, 500.0, -500.0, 500.0, -1000.0, 10.0),
    )
    assert alignment.placed[0].relation == "outside"
