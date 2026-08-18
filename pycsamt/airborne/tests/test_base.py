from __future__ import annotations

import numpy as np
import pytest

from pycsamt.airborne import (
    AirborneEMDataset,
    AirborneEMLine,
    AirborneEMRecord,
    NavigationTrack,
)
from pycsamt.emtf import EMTF, TransferFunction
from pycsamt.metadata import InstrumentMeta, SurveyMeta


def _admittance_document(scale: float = 1.0) -> EMTF:
    periods = np.array([1.0 / 100.0, 1.0 / 1000.0])
    data = scale * np.arange(12, dtype=float).reshape(2, 3, 2)
    tf = TransferFunction(
        name="admittance",
        data=data.astype(complex),
        input_channels=("Ex", "Ey"),
        output_channels=("Hx", "Hy", "Hz"),
        periods=periods,
        units="S",
    )
    doc = EMTF(periods=periods)
    doc.add_transfer_function(tf)
    return doc


def _line() -> AirborneEMLine:
    nav = NavigationTrack(
        sample_ids=("S1", "S2", "S3"),
        latitude=(5.0, 5.1, 5.2),
        longitude=(-3.0, -3.1, -3.2),
    )
    return AirborneEMLine(line_id="L001", navigation=nav)


def test_record_reuses_general_emtf_payload():
    record = AirborneEMRecord(sample_id="S1", emtf=_admittance_document())
    assert record.transfer_function_names == ("admittance",)
    tf = record.emtf.get_transfer_function("admittance")
    assert tf is not None
    assert tf.shape == (2, 3, 2)
    assert tf.input_channels == ("Ex", "Ey")
    assert tf.output_channels == ("Hx", "Hy", "Hz")


def test_record_rejects_non_emtf_payload():
    with pytest.raises(TypeError):
        AirborneEMRecord(sample_id="S1", emtf=object())


def test_line_supports_sparse_records_and_navigation_order():
    line = _line()
    line.add_emtf("S3", _admittance_document(3.0))
    line.add_emtf("S1", _admittance_document(1.0))
    assert line.n_samples == 3
    assert line.n_records == 2
    assert line.missing_sample_ids == ("S2",)
    assert [record.sample_id for record in line.iter_records()] == ["S1", "S3"]
    assert line.record_at(1) is None
    assert line.transfer_function_names == ("admittance",)


def test_line_rejects_record_not_on_navigation_track():
    line = _line()
    with pytest.raises(KeyError):
        line.add_emtf("S9", _admittance_document())


def test_line_duplicate_requires_replace():
    line = _line()
    first = line.add_emtf("S1", _admittance_document(1.0))
    with pytest.raises(ValueError):
        line.add_record(first)
    replacement = AirborneEMRecord(
        sample_id="S1",
        emtf=_admittance_document(2.0),
    )
    line.add_record(replacement, replace=True)
    assert line.get_record("S1") is replacement


def test_dataset_collection_and_bbox():
    line1 = _line()
    line1.add_emtf("S1", _admittance_document())
    line2 = AirborneEMLine(
        line_id="L002",
        navigation=NavigationTrack(
            sample_ids=("A", "B"),
            latitude=(6.0, 6.2),
            longitude=(-4.0, -4.2),
        ),
    )
    line2.add_emtf("B", _admittance_document(2.0))

    dataset = AirborneEMDataset(name="survey", lines={"L001": line1})
    dataset.add_line(line2)

    assert dataset.line_ids == ("L001", "L002")
    assert dataset.n_lines == 2
    assert dataset.n_samples == 5
    assert dataset.n_records == 2
    assert dataset.transfer_function_names == ("admittance",)
    assert dataset.bbox is not None
    assert dataset.bbox.lat_min == pytest.approx(5.0)
    assert dataset.bbox.lat_max == pytest.approx(6.2)
    assert dataset.bbox.lon_min == pytest.approx(-4.2)
    assert dataset.bbox.lon_max == pytest.approx(-3.0)
    assert set(dataset.emtf_records()) == {("L001", "S1"), ("L002", "B")}


def test_dataset_accepts_reusable_survey_and_instrument_metadata():
    survey = SurveyMeta(name="aem", method="AEM")
    instrument = InstrumentMeta(system="Example airborne EM system")
    dataset = AirborneEMDataset(
        name="aem",
        survey=survey,
        instrument=instrument,
    )
    assert dataset.method == "AEM"
    assert dataset.survey is survey
    assert dataset.instrument is instrument


def test_dataset_duplicate_line_requires_replace():
    line = _line()
    dataset = AirborneEMDataset(name="survey")
    dataset.add_line(line)
    with pytest.raises(ValueError):
        dataset.add_line(line)


def test_dataset_mapping_key_must_match_line_id():
    with pytest.raises(ValueError):
        AirborneEMDataset(name="survey", lines={"wrong": _line()})
