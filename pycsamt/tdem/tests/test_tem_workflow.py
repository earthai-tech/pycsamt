"""Tests for high-level TEMAVG workflow helpers."""

from __future__ import annotations

from pathlib import Path

import pytest

from pycsamt.tdem import (
    TEMAVGConversion,
    read_temavg_soundings,
    transform_temavg_survey,
)
from pycsamt.tdem.io import (
    read_temavg_soundings as read_soundings_io,
)

DATA_DIR = Path(__file__).parents[3] / "data" / "TEMAVG" / "JIANGSU"
AVG_FILE = DATA_DIR / "TEM100.AVG"

pytestmark = pytest.mark.skipif(
    not AVG_FILE.exists(),
    reason="TEMAVG sample data not available",
)


def test_read_temavg_soundings_exports_profile_stations():
    """Workflow reader should build soundings from selected profiles."""
    soundings = read_temavg_soundings(DATA_DIR, stems=["TEM100"])

    assert len(soundings) == 51
    assert soundings[0].station_name == "TEM100_100"
    assert soundings[0].n_gates == 25


def test_read_temavg_soundings_io_export():
    """The IO module should expose the workflow sounding reader."""
    soundings = read_soundings_io(DATA_DIR, stems=["TEM100"])

    assert len(soundings) == 51
    assert soundings[0].station_name == "TEM100_100"


def test_transform_temavg_survey_returns_result_bundle():
    """Complete workflow should return survey, soundings, and results."""
    bundle = transform_temavg_survey(DATA_DIR, stems=["TEM100"])

    assert isinstance(bundle, TEMAVGConversion)
    assert bundle.n_soundings == 51
    assert bundle.n_results == 51
    assert bundle.collection is None
    assert bundle.written_paths == []
    assert bundle.results[0]["station_name"] == "TEM100_100"
    assert len(bundle.results[0]["freq"]) == 25
    assert len(bundle.results[0]["rho_a"]) == 25


def test_transform_temavg_survey_can_return_edi_collection():
    """Workflow should create an in-memory EDI collection on request."""
    bundle = transform_temavg_survey(
        DATA_DIR,
        stems=["TEM100"],
        return_collection=True,
    )

    assert bundle.collection is not None
    assert len(bundle.collection) == bundle.n_soundings
