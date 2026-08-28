"""Phase 12 PCBH hardening, publication, and benchmark checks."""

from __future__ import annotations

import copy
import json
import random
from importlib import resources
from pathlib import Path

import pytest

from pycsamt.format.borehole import (
    PCBH_SCHEMA_URI,
    benchmark_pcbh,
    pcbh_from_dict,
)


def _example() -> dict:
    path = resources.files("pycsamt.format.borehole").joinpath(
        "examples/minimal-vertical.pcbh.json"
    )
    return json.loads(path.read_text(encoding="utf-8"))


def test_packaged_example_is_canonical_and_readable():
    root = _example()
    assert root["$schema"] == PCBH_SCHEMA_URI
    assert pcbh_from_dict(root).boreholes[0].id == "BH-001"


def test_schema_id_matches_reader_contract():
    path = resources.files("pycsamt.format.borehole").joinpath(
        "schemas/pcbh-0.1.schema.json"
    )
    assert json.loads(path.read_text(encoding="utf-8"))["$id"] == PCBH_SCHEMA_URI


def test_deterministic_mapping_mutations_fail_safely():
    """Exercise arbitrary types at contract boundaries without Hypothesis."""
    rng = random.Random(20260828)
    paths = [
        ("boreholes",),
        ("crs",),
        ("units",),
        ("conventions",),
        ("boreholes", 0, "collar"),
        ("boreholes", 0, "total_depth_md"),
        ("boreholes", 0, "trajectory", "stations"),
    ]
    values = [None, True, -1, float("nan"), "invalid", {}, [], [None]]
    for _ in range(100):
        root = copy.deepcopy(_example())
        path = rng.choice(paths)
        cursor = root
        for key in path[:-1]:
            cursor = cursor[key]
        cursor[path[-1]] = rng.choice(values)
        try:
            pcbh_from_dict(root)
        except (TypeError, ValueError):
            pass


def test_nesting_limit_rejects_adversarial_extension():
    root = _example()
    nested: dict = {}
    cursor = nested
    for _ in range(40):
        cursor["x"] = {}
        cursor = cursor["x"]
    root["extensions"] = {"test:nested": nested}
    with pytest.raises(ValueError, match="nesting"):
        pcbh_from_dict(root)


def test_benchmark_exercises_large_project_contract():
    result = benchmark_pcbh(100)
    assert result.boreholes == 100
    assert result.render_points >= 200
    assert result.json_bytes > 10_000
    assert result.parse_seconds >= 0
    assert result.render_seconds >= 0


def test_ci_matrix_runs_format_tests():
    root = Path(__file__).resolve().parents[3]
    workflow = (root / ".github/workflows/ci.yml").read_text(encoding="utf-8")
    assert 'python-version: ["3.9", "3.10", "3.11", "3.12", "3.13"]' in workflow
    assert "pycsamt/format/tests" in workflow

