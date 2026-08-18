from __future__ import annotations

from pathlib import Path

import pytest

from pycsamt.airborne import (
    AirborneEMDataset,
    AirborneFormatDefinition,
    AirborneIOError,
    available_airborne_readers,
    available_airborne_writers,
    detect_airborne_format,
    read_airborne,
    register_airborne_format,
    write_airborne,
)


def _dummy_reader(source, **kwargs):
    return AirborneEMDataset(
        name=kwargs.get("name", Path(source).stem),
        attrs={"technology": "MobileMT"},
    )


def _dummy_writer(dataset, target, **kwargs):
    Path(target).write_text(dataset.name + "\n", encoding="utf-8")
    return Path(target)


def _dummy_detector(source):
    return isinstance(source, (str, Path)) and str(source).endswith(".mmtx")


def _register_dummy():
    definition = AirborneFormatDefinition(
        name="test_mobilemt_delivery",
        technology="mobilemt",
        reader=_dummy_reader,
        writer=_dummy_writer,
        detector=_dummy_detector,
        extensions=(".mmtx",),
        aliases=("test_mmt",),
    )
    try:
        register_airborne_format(definition)
    except ValueError:
        pass


def test_no_vendor_reader_is_assumed_for_unknown_file(tmp_path):
    path = tmp_path / "unknown.dat"
    path.write_text("not a verified airborne delivery", encoding="utf-8")
    with pytest.raises(AirborneIOError, match="representative delivery file"):
        read_airborne(path)


def test_custom_format_round_trip_dispatch(tmp_path):
    _register_dummy()
    source = tmp_path / "survey.mmtx"
    source.write_text("dummy", encoding="utf-8")
    assert detect_airborne_format(source) == "test_mobilemt_delivery"

    dataset = read_airborne(source, name="survey_test")
    assert dataset.name == "survey_test"
    assert "test_mobilemt_delivery" in available_airborne_readers(
        technology="mobilemt"
    )
    assert "test_mobilemt_delivery" in available_airborne_writers(
        technology="mobilemt"
    )

    target = tmp_path / "copy.mmtx"
    result = write_airborne(dataset, target)
    assert result == target
    assert target.read_text(encoding="utf-8").strip() == "survey_test"


def test_explicit_format_alias_and_technology_guard(tmp_path):
    _register_dummy()
    source = tmp_path / "anything.bin"
    source.write_text("dummy", encoding="utf-8")
    dataset = read_airborne(
        source,
        format="test_mmt",
        technology="MobileMT",
    )
    assert dataset.attrs["technology"] == "MobileMT"
    with pytest.raises(AirborneIOError, match="not 'ztem'"):
        read_airborne(
            source,
            format="test_mmt",
            technology="ZTEM",
        )


def test_existing_dataset_is_passthrough():
    dataset = AirborneEMDataset(name="already_decoded")
    assert read_airborne(dataset) is dataset
