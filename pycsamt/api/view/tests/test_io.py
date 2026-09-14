from __future__ import annotations

from pathlib import Path

import pytest

from pycsamt.api.view.config import reset_api_view
from pycsamt.api.view.io import read_edi, read_edis, read_sites
from pycsamt.api.view.survey import APISurvey

_DATA_ROOT = Path(__file__).resolve().parents[4] / "data"
_THREE_EDIS = _DATA_ROOT / "3edis"


@pytest.fixture(autouse=True)
def _reset_config():
    reset_api_view()
    yield
    reset_api_view()


@pytest.mark.skipif(not _THREE_EDIS.exists(), reason="sample EDI data not found")
def test_read_edi_single_file():
    from pycsamt.seg.edi import EDIFile

    one_file = next(_THREE_EDIS.glob("*.edi"))
    result = read_edi(one_file)
    assert isinstance(result, EDIFile)


@pytest.mark.skipif(not _THREE_EDIS.exists(), reason="sample EDI data not found")
def test_read_edis_directory_returns_populated_survey():
    survey = read_edis(_THREE_EDIS, progress=False, strict=False)
    assert isinstance(survey, APISurvey)
    assert survey.n_sites == 3
    assert survey.source == _THREE_EDIS
    assert survey.meta["recursive"] is True
    assert not survey.errors()


@pytest.mark.skipif(not _THREE_EDIS.exists(), reason="sample EDI data not found")
def test_read_edis_list_of_paths():
    files = list(_THREE_EDIS.glob("*.edi"))
    survey = read_edis(files, progress=False, strict=False)
    assert survey.n_sites == len(files)


@pytest.mark.skipif(not _THREE_EDIS.exists(), reason="sample EDI data not found")
def test_read_edis_with_progress_enabled():
    survey = read_edis(_THREE_EDIS, progress=True, strict=False)
    assert survey.n_sites == 3


def test_read_edis_glob_with_no_matches_returns_empty_survey_with_errors(
    tmp_path,
):
    survey = read_edis(tmp_path / "*.edi", progress=False, strict=False)
    assert isinstance(survey, APISurvey)
    assert len(survey) == 0
    assert isinstance(survey.summary(), object)
    assert survey.errors()


@pytest.mark.skipif(not _THREE_EDIS.exists(), reason="sample EDI data not found")
def test_read_sites_is_alias_for_read_edis():
    survey = read_sites(_THREE_EDIS, progress=False, strict=False)
    assert isinstance(survey, APISurvey)
    assert survey.n_sites == 3


@pytest.mark.skipif(not _THREE_EDIS.exists(), reason="sample EDI data not found")
def test_read_edis_on_dup_keep_skips_duplicate_station(tmp_path):
    # Two files sharing the same >HEAD DATAID collapse onto one station;
    # on_dup="keep" must retain the first one read rather than the last.
    source = next(_THREE_EDIS.glob("*.edi"))
    text = source.read_text(encoding="utf-8")
    (tmp_path / "a_first.edi").write_text(text, encoding="utf-8")
    (tmp_path / "b_second.edi").write_text(text, encoding="utf-8")

    survey = read_edis(tmp_path, progress=False, strict=False, on_dup="keep")
    assert survey.n_sites == 1
    kept = survey[0]
    assert kept.path.name == "a_first.edi"
