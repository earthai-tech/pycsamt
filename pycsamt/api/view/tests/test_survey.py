from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from pycsamt.api.view.config import reset_api_view
from pycsamt.api.view.frame import APIFrame
from pycsamt.api.view.io import read_edis
from pycsamt.api.view.survey import APISurvey
from pycsamt.seg.collection import EDICollection

_DATA_ROOT = Path(__file__).resolve().parents[4] / "data"
_THREE_EDIS = _DATA_ROOT / "3edis"


@pytest.fixture(autouse=True)
def _reset_config():
    reset_api_view()
    yield
    reset_api_view()


@pytest.fixture
def real_survey() -> APISurvey:
    return read_edis(_THREE_EDIS, progress=False, strict=False)


# ─────────────────────────────────────────────────────────────────────────
# Construction / coercion
# ─────────────────────────────────────────────────────────────────────────


def test_default_construction_is_empty():
    survey = APISurvey()
    assert survey.n_sites == 0
    assert survey.name == "survey"
    assert isinstance(survey.collection, EDICollection)


def test_construction_from_existing_edicollection_is_not_rewrapped():
    coll = EDICollection()
    survey = APISurvey(coll)
    assert survey.collection is coll


def test_construction_wraps_non_collection_iterable():
    survey = APISurvey([])
    assert isinstance(survey.collection, EDICollection)


def test_construction_sets_name_source_parser_meta():
    survey = APISurvey(name="mine", source="src", parser="p", meta={"x": 1})
    assert survey.name == "mine"
    assert survey.source == "src"
    assert survey.parser == "p"
    assert survey.meta == {"x": 1}


# ─────────────────────────────────────────────────────────────────────────
# Real data: n_sites / stations / paths / iteration / getitem
# ─────────────────────────────────────────────────────────────────────────


@pytest.mark.skipif(not _THREE_EDIS.exists(), reason="sample EDI data not found")
def test_real_survey_basic_properties(real_survey):
    assert real_survey.n_sites == 3
    assert len(real_survey) == 3
    assert len(real_survey.stations) == 3
    assert len(real_survey.paths) == 3
    assert list(iter(real_survey))
    first = real_survey[0]
    assert first is not None


@pytest.mark.skipif(not _THREE_EDIS.exists(), reason="sample EDI data not found")
def test_real_survey_df_property_is_summary(real_survey):
    assert real_survey.df.to_dict() == real_survey.summary().to_dict()


@pytest.mark.skipif(not _THREE_EDIS.exists(), reason="sample EDI data not found")
def test_real_survey_data_property_is_collection(real_survey):
    assert real_survey.data is real_survey.collection


@pytest.mark.skipif(not _THREE_EDIS.exists(), reason="sample EDI data not found")
def test_real_survey_getattr_delegates_to_collection(real_survey):
    assert real_survey.nf_stats() == real_survey.collection.nf_stats()


@pytest.mark.skipif(not _THREE_EDIS.exists(), reason="sample EDI data not found")
def test_real_survey_dir_includes_collection_attrs(real_survey):
    assert "nf_stats" in dir(real_survey)


def test_getattr_raises_for_unknown_attribute():
    survey = APISurvey()
    with pytest.raises(AttributeError):
        survey.totally_unknown_thing


@pytest.mark.skipif(not _THREE_EDIS.exists(), reason="sample EDI data not found")
def test_real_survey_summary_fields_filter(real_survey):
    frame = real_survey.summary(fields=("station", "n_freq"))
    assert isinstance(frame, APIFrame)
    assert list(frame.columns) == ["station", "n_freq"]
    assert frame.kind == "edi.summary"


@pytest.mark.skipif(not _THREE_EDIS.exists(), reason="sample EDI data not found")
def test_real_survey_get_site_found_and_missing(real_survey):
    station = real_survey.stations[0]
    found = real_survey.get_site(station)
    assert found is not None
    missing = real_survey.get_site("does-not-exist-xyz", default="fallback")
    assert missing == "fallback"


@pytest.mark.skipif(not _THREE_EDIS.exists(), reason="sample EDI data not found")
def test_real_survey_to_collection(real_survey):
    assert real_survey.to_collection() is real_survey.collection


@pytest.mark.skipif(not _THREE_EDIS.exists(), reason="sample EDI data not found")
def test_real_survey_update_meta_returns_self(real_survey):
    result = real_survey.update_meta(owner="me")
    assert result is real_survey
    assert real_survey.meta["owner"] == "me"


@pytest.mark.skipif(not _THREE_EDIS.exists(), reason="sample EDI data not found")
def test_real_survey_summary_text_and_str(real_survey):
    text = real_survey.summary_text()
    assert text == str(real_survey)
    assert f"APISurvey: {real_survey.name}" in text
    assert "sites: 3" in text
    assert "stations:" in text


def test_str_without_stations_or_source_or_errors():
    survey = APISurvey(name="empty")
    text = str(survey)
    assert "sites: 0" in text
    assert "stations:" not in text
    assert "source:" not in text
    assert "errors:" not in text


def test_str_truncates_long_station_list():
    class FakeCollection(EDICollection):
        def stations(self):
            return [f"S{i}" for i in range(12)]

    survey = APISurvey(name="many")
    survey.collection = FakeCollection()
    text = str(survey)
    assert "..." in text


def test_str_includes_source_when_given():
    survey = APISurvey(name="n", source="my-source")
    assert "source: my-source" in str(survey)


# ─────────────────────────────────────────────────────────────────────────
# errors()
# ─────────────────────────────────────────────────────────────────────────


def test_errors_returns_empty_list_when_no_parser():
    survey = APISurvey()
    assert survey.errors() == []


def test_errors_returns_empty_list_when_parser_errors_not_callable():
    class FakeParser:
        errors = "not-callable"

    survey = APISurvey(parser=FakeParser())
    assert survey.errors() == []


def test_errors_returns_list_from_callable():
    class FakeParser:
        def errors(self):
            return [("path", ValueError("boom"))]

    survey = APISurvey(parser=FakeParser())
    result = survey.errors()
    assert len(result) == 1
    assert result[0][0] == "path"


def test_str_shows_error_count():
    class FakeParser:
        def errors(self):
            return [("p", ValueError("x"))]

    survey = APISurvey(name="n", parser=FakeParser())
    assert "errors: 1" in str(survey)


# ─────────────────────────────────────────────────────────────────────────
# to_dataframe
# ─────────────────────────────────────────────────────────────────────────


def test_to_dataframe_falls_back_to_summary_when_no_method():
    survey = APISurvey()
    result = survey.to_dataframe()
    assert isinstance(result, APIFrame)


def test_to_dataframe_uses_collection_method_and_wraps_by_default():
    survey = APISurvey()
    survey.collection.to_dataframe = lambda *a, **k: pd.DataFrame({"a": [1]})
    result = survey.to_dataframe()
    assert isinstance(result, APIFrame)
    assert result.kind == "edi.data"
    assert result.name == f"{survey.name}_data"


def test_to_dataframe_api_false_returns_raw():
    survey = APISurvey()
    survey.collection.to_dataframe = lambda *a, **k: pd.DataFrame({"a": [1]})
    result = survey.to_dataframe(api=False)
    assert isinstance(result, pd.DataFrame)
    assert not isinstance(result, APIFrame)


def test_to_dataframe_custom_kind_kwarg():
    survey = APISurvey()
    survey.collection.to_dataframe = lambda *a, **k: pd.DataFrame({"a": [1]})
    result = survey.to_dataframe(kind="custom.kind")
    assert result.kind == "custom.kind"
