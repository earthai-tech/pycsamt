"""Tests for PCBH and legacy geology compatibility adapters."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.geology import Borehole, Interval, RockDatabase, RockEntry
from pycsamt.interp import ModelCalibrator, ResistivityModel
from pycsamt.format.borehole import (
    Collar,
    CoordinateReferenceSystem,
    LogInterval,
    PCBHBorehole,
    PCBHDocument,
    StructureObservation,
    SurveyStation,
    Trajectory,
    UnitSystem,
    VocabularyEntry,
    from_legacy_borehole,
    legacy_borehole_views,
    legacy_conversion_losses,
    to_legacy_borehole,
)


def _pcbh_hole(**overrides) -> PCBHBorehole:
    values = {
        "id": "BH-01",
        "name": "Calibration hole",
        "kind": "mining_exploration",
        "status": "completed",
        "collar": Collar(421000.0, 900000.0, 250.0),
        "total_depth_md": 100.0,
        "trajectory": Trajectory(method="vertical"),
        "interval_logs": {
            "lithology": [
                LogInterval(0.0, 40.0, code="CLY", resistivity_ohm_m=10.0),
                LogInterval(
                    40.0,
                    100.0,
                    label="Granite",
                    resistivity_ohm_m=1000.0,
                ),
            ]
        },
    }
    values.update(overrides)
    return PCBHBorehole(**values)


def _document(hole: PCBHBorehole | None = None) -> PCBHDocument:
    return PCBHDocument(
        document_id="test:adapters",
        created_at="2026-08-28T00:00:00Z",
        created_by="pytest",
        crs=CoordinateReferenceSystem("EPSG:32629"),
        boreholes=[hole or _pcbh_hole()],
        lithologies=[VocabularyEntry("CLY", "Clay")],
    )


def test_legacy_to_pcbh_preserves_queries_and_explicit_collar():
    legacy = Borehole(
        "BH old",
        125.0,
        [
            Interval(0.0, 20.0, "Clay", 12.0),
            Interval(20.0, 55.0, "Sand", None),
        ],
        collar_elevation=99.0,
    )
    result = from_legacy_borehole(
        legacy,
        collar=Collar(500000.0, 1000000.0, 99.0),
        borehole_id="BH-OLD",
        kind="water",
        status="completed",
    )

    assert result.id == "BH-OLD"
    assert result.collar.x == 500000.0
    assert result.metadata["legacy_profile_x"] == 125.0
    assert result.total_depth_md == 55.0
    assert result.interval_logs["lithology"][0].resistivity_ohm_m == 12.0


def test_pcbh_to_legacy_uses_explicit_profile_mapping_and_vocabulary():
    legacy = to_legacy_borehole(
        _pcbh_hole(),
        profile_x=325.0,
        vocabulary=[VocabularyEntry("CLY", "Clay")],
    )

    assert legacy.x == 325.0
    assert legacy.collar_elevation == 250.0
    assert legacy.lithology_at_depth(0.0) == "Clay"
    assert legacy.lithology_at_depth(39.999) == "Clay"
    assert legacy.lithology_at_depth(40.0) == "Granite"
    assert legacy.tres_at_depth(20.0) == 10.0
    assert legacy.tres_at_depth(75.0) == 1000.0
    assert legacy.interval_at_depth(100.0) is None


def test_rock_database_is_an_opt_in_fallback_for_unresolved_code():
    hole = _pcbh_hole(
        interval_logs={
            "lithology": [
                LogInterval(0.0, 100.0, code="LOCAL", resistivity_ohm_m=25.0)
            ]
        }
    )
    db = RockDatabase([RockEntry("Resolved rock", 20.0, 30.0)])

    unresolved = to_legacy_borehole(hole, profile_x=0.0)
    resolved = to_legacy_borehole(hole, profile_x=0.0, rock_db=db)

    assert unresolved.lithology_at_depth(50.0) == "LOCAL"
    assert resolved.lithology_at_depth(50.0) == "Resolved rock"


def test_calibration_views_accept_mapping_and_callable():
    document = _document()

    mapped = legacy_borehole_views(document, profile_x={"BH-01": 75.0})
    called = legacy_borehole_views(document, profile_x=lambda hole: 125.0)

    assert mapped[0].x == 75.0
    assert called[0].x == 125.0
    assert mapped[0].tres_at_depth(10.0) == 10.0


def test_calibrator_result_is_unchanged_by_pcbh_compatibility_view():
    model = ResistivityModel.from_array(
        np.log10(np.array([[12.0], [900.0]])),
        x_centers=np.array([75.0]),
        z_centers=np.array([20.0, 70.0]),
        station_x=np.array([75.0]),
        station_names=["S1"],
        method="adapter-test",
    )
    direct = Borehole(
        "Calibration hole",
        75.0,
        [
            Interval(0.0, 40.0, "Clay", 10.0),
            Interval(40.0, 100.0, "Granite", 1000.0),
        ],
        collar_elevation=250.0,
    )
    adapted = legacy_borehole_views(
        _document(), profile_x={"BH-01": 75.0}
    )[0]

    direct_result = ModelCalibrator(verbose=False).fit(
        model, [direct]
    ).calibrated_model()
    adapted_result = ModelCalibrator(verbose=False).fit(
        model, [adapted]
    ).calibrated_model()

    np.testing.assert_array_equal(
        adapted_result.rho_2d, direct_result.rho_2d
    )


def test_calibration_views_reject_incompatible_document_units():
    document = _document()
    document.units = UnitSystem(depth="ft")

    with pytest.raises(ValueError, match="convert document units"):
        legacy_borehole_views(document, profile_x={"BH-01": 0.0})


def test_empty_legacy_borehole_requires_explicit_total_depth():
    empty = Borehole("Empty", 0.0, [])
    with pytest.raises(ValueError, match="total_depth_md must be provided"):
        from_legacy_borehole(empty, collar=Collar(0.0, 0.0, 0.0))

    promoted = from_legacy_borehole(
        empty,
        collar=Collar(0.0, 0.0, 0.0),
        total_depth_md=50.0,
    )
    assert promoted.total_depth_md == 50.0
    assert promoted.interval_logs == {}


def test_downgrade_losses_report_source_specific_information():
    hole = _pcbh_hole(
        trajectory=Trajectory(
            method="survey",
            stations=[
                SurveyStation(0.0, 10.0, 2.0),
                SurveyStation(100.0, 20.0, 5.0),
            ],
        ),
        interval_logs={
            **_pcbh_hole().interval_logs,
            "weathering": [LogInterval(0.0, 10.0, label="Weathered")],
        },
        structures=[StructureObservation(kind="fault", at_md=25.0)],
    )
    losses = legacy_conversion_losses(hole)

    assert "deviated trajectory and survey orientations" in losses
    assert "interval log families: weathering" in losses
    assert "structural observations" in losses


def test_pcbh_to_legacy_requires_finite_explicit_profile_x():
    with pytest.raises(ValueError, match="profile_x"):
        to_legacy_borehole(_pcbh_hole(), profile_x=float("nan"))
