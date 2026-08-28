"""Tests for deterministic PCBH trajectory derivation."""

from __future__ import annotations

import math

import pytest

from pycsamt.format.borehole import (
    Collar,
    PCBHBorehole,
    SurveyStation,
    Trajectory,
    desurvey,
    trajectory_checksum,
)


def _borehole(
    stations: list[SurveyStation] | None = None,
    *,
    total_depth: float = 100.0,
) -> PCBHBorehole:
    method = "survey" if stations is not None else "vertical"
    return PCBHBorehole(
        id="BH-1",
        name="Reference hole",
        kind="scientific",
        status="completed",
        collar=Collar(500.0, 1000.0, 80.0),
        total_depth_md=total_depth,
        trajectory=Trajectory(method=method, stations=stations or []),
    )


def test_vertical_trajectory_uses_absolute_xyz_and_positive_tvd():
    result = desurvey(_borehole())

    assert len(result.points) == 2
    assert result.points[0].md == 0.0
    end = result.points[-1]
    assert (end.x, end.y, end.z, end.tvd) == (500.0, 1000.0, -20.0, 100.0)
    assert (end.azimuth_deg, end.inclination_deg) == (0.0, 0.0)


def test_constant_inclined_course_matches_closed_form():
    result = desurvey(
        _borehole(
            [
                SurveyStation(0.0, 90.0, 60.0),
                SurveyStation(100.0, 90.0, 60.0),
            ]
        )
    )
    end = result.points[-1]

    assert end.x == pytest.approx(500.0 + 100.0 * math.sin(math.pi / 3))
    assert end.y == pytest.approx(1000.0)
    assert end.tvd == pytest.approx(50.0)
    assert end.z == pytest.approx(30.0)


def test_quarter_circle_reference_and_interpolation():
    """A 90-degree build has an analytic circular-arc solution."""
    result = desurvey(
        _borehole(
            [
                SurveyStation(0.0, 90.0, 0.0),
                SurveyStation(100.0, 90.0, 90.0),
            ]
        )
    )
    radius = 200.0 / math.pi
    end = result.points[-1]
    midpoint = result.at_md(50.0)

    assert end.x == pytest.approx(500.0 + radius)
    assert end.y == pytest.approx(1000.0)
    assert end.tvd == pytest.approx(radius)
    assert midpoint.x == pytest.approx(
        500.0 + radius * (1.0 - math.cos(math.pi / 4))
    )
    assert midpoint.tvd == pytest.approx(radius * math.sin(math.pi / 4))
    assert midpoint.inclination_deg == pytest.approx(45.0)


def test_azimuth_interpolation_crosses_north_not_south():
    result = desurvey(
        _borehole(
            [
                SurveyStation(0.0, 350.0, 60.0),
                SurveyStation(100.0, 10.0, 60.0),
            ]
        )
    )

    midpoint = result.at_md(50.0)
    assert midpoint.azimuth_deg == pytest.approx(0.0, abs=1e-12)
    assert midpoint.y > 1000.0


def test_split_inserts_exact_boundaries_without_moving_path():
    result = desurvey(
        _borehole(
            [
                SurveyStation(0.0, 30.0, 10.0),
                SurveyStation(100.0, 60.0, 30.0),
            ]
        )
    )
    expected = result.at_md(37.25)
    split = result.split_at([75.0, 37.25, 37.25])

    assert [point.md for point in split.points] == [0.0, 37.25, 75.0, 100.0]
    actual = split.at_md(37.25)
    assert actual == expected
    assert split.at_md(60.0).x == pytest.approx(result.at_md(60.0).x)
    assert split.at_md(60.0).y == pytest.approx(result.at_md(60.0).y)
    assert split.at_md(60.0).z == pytest.approx(result.at_md(60.0).z)


def test_missing_collar_and_terminal_stations_extend_nearest_attitude():
    result = desurvey(
        _borehole(
            [
                SurveyStation(20.0, 90.0, 30.0),
                SurveyStation(80.0, 90.0, 30.0),
            ]
        )
    )

    assert [point.md for point in result.points] == [0.0, 20.0, 80.0, 100.0]
    assert result.points[-1].x == pytest.approx(550.0)
    assert result.points[-1].tvd == pytest.approx(100.0 * math.cos(math.pi / 6))


@pytest.mark.parametrize("md", [-0.01, 100.01, math.nan, math.inf])
def test_interpolation_rejects_out_of_range_or_nonfinite_md(md):
    result = desurvey(_borehole())
    with pytest.raises((TypeError, ValueError)):
        result.at_md(md)


def test_antipodal_dogleg_is_rejected_as_ambiguous():
    borehole = _borehole(
        [
            SurveyStation(0.0, 0.0, 0.0),
            SurveyStation(100.0, 0.0, 180.0),
        ]
    )
    with pytest.raises(ValueError, match="no unique path"):
        desurvey(borehole)


def test_checksum_is_deterministic_and_uses_only_geometry_inputs():
    first = _borehole()
    same_geometry = _borehole()
    same_geometry.name = "Renamed"
    changed_geometry = _borehole(total_depth=101.0)

    assert trajectory_checksum(first) == trajectory_checksum(first)
    assert trajectory_checksum(first) == trajectory_checksum(same_geometry)
    assert trajectory_checksum(first) != trajectory_checksum(changed_geometry)
    assert desurvey(first).source_checksum == trajectory_checksum(first)
