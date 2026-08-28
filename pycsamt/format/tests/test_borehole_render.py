"""Geometry-only tests for the shared PCBH render-model builder."""

from __future__ import annotations

import pytest

from pycsamt.format.borehole import (
    BoreholeRenderBuilder,
    Collar,
    CoordinateReferenceSystem,
    DisplayRadiusPolicy,
    LogInterval,
    PCBHBorehole,
    PCBHDocument,
    StructureObservation,
    SurveyStation,
    Trajectory,
    UnitSystem,
    VocabularyEntry,
    build_render_model,
    deterministic_color,
)


def _hole(hole_id="BH-1", *, x=100.0) -> PCBHBorehole:
    return PCBHBorehole(
        id=hole_id,
        name=f"Hole {hole_id}",
        kind="mining_exploration",
        status="completed",
        collar=Collar(x, 200.0, 50.0),
        total_depth_md=100.0,
        diameter=0.2,
        trajectory=Trajectory(
            method="survey",
            north_reference="grid",
            stations=[
                SurveyStation(0.0, 350.0, 0.0),
                SurveyStation(50.0, 0.0, 30.0),
                SurveyStation(100.0, 10.0, 45.0),
            ],
        ),
        interval_logs={
            "lithology": [
                LogInterval(0.0, 37.25, code="CLAY", label="Clay"),
                LogInterval(37.25, 100.0, code="GRAN", label="Granite"),
            ]
        },
        structures=[
            StructureObservation(
                kind="fracture",
                at_md=60.0,
                orientation_representation="global_plane",
                dip_deg=45.0,
                dip_direction_deg=120.0,
            )
        ],
    )


def _document() -> PCBHDocument:
    return PCBHDocument(
        document_id="test:render",
        created_at="2026-08-28T00:00:00Z",
        created_by="pytest",
        crs=CoordinateReferenceSystem("EPSG:32629"),
        boreholes=[_hole(), _hole("BH-2", x=300.0)],
        lithologies=[
            VocabularyEntry("CLAY", "Clay", color="#765432"),
            VocabularyEntry("GRAN", "Granite"),
        ],
    )


def test_render_model_contains_exact_deviated_interval_geometry():
    model = build_render_model(_document())
    hole = model.boreholes[0]

    assert hole.collar.position == (100.0, 200.0, 50.0)
    assert any(point.md == 37.25 for point in hole.centerline.points)
    first, second = hole.interval_segments
    assert first.points[-1].md == 37.25
    assert second.points[0].md == 37.25
    assert first.points[-1] == second.points[0]
    assert first.points[-1].y != 200.0
    assert hole.contacts[0].position == (
        first.points[-1].x,
        first.points[-1].y,
        first.points[-1].z,
    )


def test_colors_use_vocabulary_then_stable_fallback_and_batch():
    model = build_render_model(_document())
    first = model.boreholes[0].interval_segments

    assert first[0].color == "#765432"
    assert first[1].color == deterministic_color("GRAN")
    assert deterministic_color("GRAN") == deterministic_color("gran")
    granite_batch = next(batch for batch in model.batches if batch.color == first[1].color)
    assert len(granite_batch.segment_ids) == 2


def test_selection_and_hover_metadata_are_backend_neutral():
    model = build_render_model(_document(), selected_ids={"BH-1"})
    selected, normal = model.boreholes

    assert selected.collar.selected
    assert selected.centerline.selected
    assert selected.interval_segments[0].selected
    assert not normal.collar.selected
    assert selected.collar.metadata["kind"] == "mining_exploration"
    assert selected.interval_segments[0].metadata["label"] == "Clay"
    assert selected.contacts[0].metadata["tvd"] > 0


def test_structure_glyph_has_exact_position_tangent_and_orientation():
    model = build_render_model(_document())
    glyph = model.boreholes[0].structure_glyphs[0]

    assert glyph.md == 60.0
    assert glyph.representation == "global_plane"
    assert glyph.orientation == {"dip_deg": 45.0, "dip_direction_deg": 120.0}
    assert sum(value * value for value in glyph.borehole_tangent) == pytest.approx(1.0)


def test_radius_modes_do_not_change_physical_diameter():
    document = _document()
    fixed = build_render_model(
        document,
        radius_policy=DisplayRadiusPolicy(mode="fixed", fixed_radius=2.5),
    )
    exaggerated = build_render_model(
        document,
        radius_policy=DisplayRadiusPolicy(mode="exaggeration", exaggeration=10.0),
    )

    assert fixed.boreholes[0].display_radius == 2.5
    assert exaggerated.boreholes[0].display_radius == 1.0
    assert document.boreholes[0].diameter == 0.2


def test_lod_preserves_geological_boundaries_and_reduces_points():
    full = build_render_model(_document(), lod_tolerance=0.0)
    reduced = build_render_model(_document(), lod_tolerance=1000.0)
    full_points = full.boreholes[0].centerline.points
    reduced_points = reduced.boreholes[0].centerline.points

    assert len(reduced_points) < len(full_points)
    assert [point.md for point in reduced_points] == [0.0, 37.25, 100.0]
    assert reduced.boreholes[0].contacts[0].md == 37.25


def test_builder_caches_unsplit_scientific_trajectories():
    builder = BoreholeRenderBuilder()
    document = _document()

    first = builder.build(document)
    cache_size = builder.cache_size
    second = builder.build(document, lod_tolerance=10.0)

    assert cache_size == 2
    assert builder.cache_size == cache_size
    assert first.boreholes[0].source_checksum == second.boreholes[0].source_checksum


def test_sampling_is_bounded_for_large_holes():
    document = _document()
    model = build_render_model(
        document,
        sampling_step_md=0.01,
        max_points_per_hole=20,
    )

    assert len(model.boreholes[0].centerline.points) <= 23
    assert model.max_points_per_hole == 20


def test_rendering_rejects_mismatched_depth_and_coordinate_units():
    document = _document()
    document.units = UnitSystem(depth="ft")
    with pytest.raises(ValueError, match="units to match"):
        build_render_model(document)
