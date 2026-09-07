"""Tests for the PCGL interpretation-legend mini-format."""

from __future__ import annotations

import pytest

from pycsamt.format import (
    GeologyLegend,
    GeologyLegendValidationError,
    legend_from_csv,
    legend_from_dict,
    legend_to_dict,
    read_legend,
    write_legend,
)
from pycsamt.geology import RockDatabase, RockEntry


def _legend() -> GeologyLegend:
    db = RockDatabase(
        [
            RockEntry("Sand/gravel (saturated)", 10.0, 200.0, "#E9C46A"),
            RockEntry(
                "Granodiorite (fresh)", 1000.0, 5000.0, "#8D99AE",
                pattern_id="cross-hatch-01", pattern_source="builtin",
            ),
        ]
    )
    return GeologyLegend.from_rock_database(
        db, document_id="pcgl:test", title="Test legend"
    )


def test_json_round_trip(tmp_path):
    original = _legend()
    path = tmp_path / "legend.pcgl.json"
    write_legend(original, path)
    restored = read_legend(path)
    assert [e.name for e in restored.entries] == [
        "Sand/gravel (saturated)", "Granodiorite (fresh)",
    ]
    assert restored.entries[1].pattern_id == "cross-hatch-01"
    assert restored.entries[1].pattern_source == "builtin"
    assert restored.pcgl_version == "0.1.0"
    assert restored.title == "Test legend"


def test_classify_delegates_to_rock_database():
    legend = _legend()
    assert legend.classify(50.0).name == "Sand/gravel (saturated)"
    assert legend.classify(2000.0).name == "Granodiorite (fresh)"


def test_from_rock_database_default_is_valid():
    legend = GeologyLegend.from_rock_database()
    assert len(legend.entries) == len(RockDatabase.default())
    assert legend.issues() == []


def test_csv_round_trip_carries_pattern_columns(tmp_path):
    original = _legend()
    csv_path = tmp_path / "legend.csv"
    original.to_csv(csv_path)
    restored = legend_from_csv(csv_path)
    assert restored.entries[1].pattern_id == "cross-hatch-01"
    assert restored.entries[1].pattern_source == "builtin"


def test_plain_rock_csv_loads_with_no_pattern(tmp_path):
    csv_path = tmp_path / "plain.csv"
    csv_path.write_text(
        "name,rho_min,rho_max\nClay,1,20\nLimestone,200,5000\n",
        encoding="utf-8",
    )
    legend = legend_from_csv(csv_path)
    assert [e.pattern_id for e in legend.entries] == ["", ""]
    assert legend.issues() == []


def test_validation_rejects_empty_legend():
    bad = GeologyLegend(
        document_id="pcgl:bad",
        created_at="2026-09-04T00:00:00Z",
        created_by="pytest",
        db=RockDatabase([]),
    )
    with pytest.raises(GeologyLegendValidationError):
        bad.validate()


def test_validation_rejects_inverted_range():
    bad = GeologyLegend(
        document_id="pcgl:bad",
        created_at="2026-09-04T00:00:00Z",
        created_by="pytest",
        db=RockDatabase([RockEntry("Bad", 100.0, 10.0)]),
    )
    assert any("rho_min" in issue for issue in bad.issues())


def test_validation_rejects_duplicate_names():
    bad = GeologyLegend(
        document_id="pcgl:bad",
        created_at="2026-09-04T00:00:00Z",
        created_by="pytest",
        db=RockDatabase(
            [RockEntry("Clay", 1.0, 20.0), RockEntry("Clay", 20.0, 40.0)]
        ),
    )
    assert any("duplicate" in issue for issue in bad.issues())


def test_validation_rejects_pattern_id_without_source():
    bad = GeologyLegend(
        document_id="pcgl:bad",
        created_at="2026-09-04T00:00:00Z",
        created_by="pytest",
        db=RockDatabase(
            [RockEntry("Clay", 1.0, 20.0, pattern_id="stipple-01")]
        ),
    )
    assert any("pattern_id" in issue for issue in bad.issues())


def test_legend_from_dict_rejects_non_dict():
    with pytest.raises(GeologyLegendValidationError):
        legend_from_dict([])  # type: ignore[arg-type]


def test_legend_to_dict_is_json_safe_and_stable_keys():
    payload = legend_to_dict(_legend())
    assert payload["pcgl_version"] == "0.1.0"
    assert payload["rho_unit"] == "ohm.m"
    assert payload["entries"][0]["name"] == "Sand/gravel (saturated)"
