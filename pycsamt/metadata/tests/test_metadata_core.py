# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Unit tests for :mod:`pycsamt.metadata` core value objects:
frequency bands, bounding boxes, and survey metadata.
"""

from __future__ import annotations

import math
from datetime import date

import numpy as np
import pytest

from pycsamt.metadata.frequency import (
    MT_BANDS,
    REGISTRY,
    FrequencyBand,
    band_for_frequency,
    doi_estimate,
    frequency_range,
    register_band,
)
from pycsamt.metadata.survey import BBox, SurveyMeta


# ---------------------------- FrequencyBand ---------------------------

def test_band_period_and_decades():
    b = MT_BANDS["AMT"]
    assert b.f_min == 10.0
    assert b.f_max == 1e5
    assert b.period_min == pytest.approx(1e-5)
    assert b.period_max == pytest.approx(0.1)
    assert b.n_decades == pytest.approx(4.0)
    assert b.frequency_range == (10.0, 1e5)


def test_band_validation_and_swap():
    with pytest.raises(ValueError):
        FrequencyBand("X", "x", -1.0, 10.0, "MT")
    swapped = FrequencyBand("X", "x", 100.0, 1.0, "MT")
    assert swapped.f_min == 1.0
    assert swapped.f_max == 100.0


def test_band_membership_and_clip():
    b = MT_BANDS["AMT"]
    assert 1000.0 in b
    assert not b.contains(1.0)
    assert "spam" not in b
    kept = b.clip_frequencies([1.0, 100.0, 1e6])
    assert kept.tolist() == [100.0]


def test_band_skin_depth_formula():
    b = MT_BANDS["AMT"]
    # delta ~ 503.3 sqrt(rho / f): ratio between two frequencies
    d1 = b.skin_depth_m(100.0, rho=100.0)
    d2 = b.skin_depth_m(1.0, rho=100.0)
    assert d2 / d1 == pytest.approx(10.0)
    assert d1 == pytest.approx(503.3 * math.sqrt(1.0), rel=0.01)

    shallow, deep = b.doi_range_m(rho=100.0)
    assert shallow < deep


def test_band_logspace_and_overlaps():
    b = MT_BANDS["AMT"]
    freqs = b.logspace(11)
    assert freqs.size == 11
    assert freqs[0] == pytest.approx(b.f_min)
    assert freqs[-1] == pytest.approx(b.f_max)
    assert np.all(np.diff(np.log10(freqs)) > 0)

    mt = MT_BANDS["MT"]
    assert b.overlaps(mt)
    inter = b.intersection(mt)
    assert inter == (10.0, 1e3)


def test_band_to_dict_roundtrip_fields():
    d = MT_BANDS["MT"].to_dict()
    assert d["name"] == "MT"
    assert d["f_min"] == 1e-4
    assert "doi_deep_m" in d
    assert "FrequencyBand" in repr(MT_BANDS["MT"])


# ------------------------- registry helpers ---------------------------

def test_band_for_frequency_sorted_by_width():
    matches = band_for_frequency(1000.0)
    assert matches, "1 kHz must belong to at least one band"
    decades = [b.n_decades for b in matches]
    assert decades == sorted(decades)


def test_frequency_range_name_and_method_union():
    assert frequency_range("AMT") == (10.0, 1e5)
    lo, hi = frequency_range("MT")  # direct name lookup
    assert (lo, hi) == (1e-4, 1e3)
    with pytest.raises(KeyError):
        frequency_range("not-a-method")


def test_register_band_adds_to_registry():
    band = FrequencyBand("TESTBAND", "test", 1.0, 2.0, "MT")
    try:
        register_band(band)
        assert "TESTBAND" in REGISTRY
        assert frequency_range("TESTBAND") == (1.0, 2.0)
    finally:
        REGISTRY.pop("TESTBAND", None)


def test_doi_estimate_scalar():
    assert doi_estimate(1.0, rho=100.0) == pytest.approx(5033.0, rel=0.01)


# --------------------------------- BBox --------------------------------

def test_bbox_geometry():
    bbox = BBox(27.8, 28.9, 101.5, 103.2)
    assert bbox.centre == pytest.approx((28.35, 102.35))
    assert bbox.area_deg2 == pytest.approx(1.1 * 1.7)
    assert bbox.span_lat == pytest.approx(1.1)
    assert bbox.span_lon == pytest.approx(1.7)


def test_bbox_membership_forms():
    bbox = BBox(27.8, 28.9, 101.5, 103.2)
    assert 28.0 in bbox                      # latitude scalar
    assert (28.0, 102.0) in bbox             # (lat, lon) tuple
    assert (10.0, 102.0) not in bbox
    assert "junk" not in bbox
    assert bbox.contains(27.8, 101.5)        # boundary inclusive


def test_bbox_validation_and_roundtrip():
    with pytest.raises(ValueError):
        BBox(30.0, 20.0, 0.0, 1.0)
    with pytest.raises(ValueError):
        BBox.from_coords([], [])
    bbox = BBox.from_coords([1.0, 3.0, 2.0], [10.0, 12.0])
    assert bbox.lat_min == 1.0 and bbox.lat_max == 3.0
    assert BBox.from_dict(bbox.to_dict()) == bbox


# ------------------------------ SurveyMeta -----------------------------

def test_survey_meta_validation_and_duration():
    with pytest.raises(ValueError):
        SurveyMeta(name="s", method="SONAR")
    with pytest.raises(ValueError):
        SurveyMeta(
            name="s",
            date_start=date(2024, 5, 2),
            date_end=date(2024, 5, 1),
        )
    meta = SurveyMeta(
        name="s",
        method="amt",  # normalised to upper case
        date_start=date(2024, 5, 1),
        date_end=date(2024, 5, 11),
    )
    assert meta.method == "AMT"
    assert meta.duration_days == 10
    assert SurveyMeta(name="s").duration_days is None


def test_survey_meta_completeness():
    partial = SurveyMeta(name="s")
    assert not partial.is_complete
    full = SurveyMeta(
        name="s",
        bbox=BBox(0.0, 1.0, 0.0, 1.0),
        n_stations=10,
    )
    assert full.is_complete


def test_survey_meta_json_roundtrip(tmp_path):
    meta = SurveyMeta(
        name="WILLY_L18",
        project="Phase-III",
        operator="EarthAI-Tech",
        method="AMT",
        bbox=BBox(27.8, 28.9, 101.5, 103.2),
        date_start=date(2023, 3, 1),
        date_end=date(2023, 3, 20),
        n_stations=128,
        extra={"contractor": "field-co"},
    )
    path = tmp_path / "survey.json"
    meta.to_json(path)
    back = SurveyMeta.from_json(path)
    assert back.name == meta.name
    assert back.method == "AMT"
    assert back.bbox == meta.bbox
    assert back.date_start == meta.date_start
    assert back.duration_days == 19
    assert back.extra.get("contractor") == "field-co"


def test_survey_meta_from_sites_builds_bbox():
    class FakeSite:
        def __init__(self, lat, lon):
            self.coords = (lat, lon, 0.0)

    sites = [FakeSite(28.0, 102.0), FakeSite(28.5, 102.5),
             FakeSite(None, None)]
    meta = SurveyMeta.from_sites(sites, name="from-sites", method="AMT")
    assert meta.n_stations == 3
    assert meta.bbox is not None
    assert meta.bbox.lat_min == 28.0
    assert meta.bbox.lat_max == 28.5
