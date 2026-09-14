from __future__ import annotations

import pytest

from pycsamt.airborne import AirborneEMDataset
from pycsamt.airborne.io import (
    AirborneIOError,
    available_airborne_readers,
    available_airborne_writers,
    read_airborne,
    write_airborne,
)
from pycsamt.airborne.registry import (
    AirborneFormatDefinition,
    register_airborne_format,
)


def _empty_dataset(name: str = "io_test") -> AirborneEMDataset:
    return AirborneEMDataset(name=name)


# ─────────────────────────────────────────────────────────────────────────
# read_airborne / _resolved_format
# ─────────────────────────────────────────────────────────────────────────


def test_read_airborne_unknown_explicit_format_raises():
    with pytest.raises(AirborneIOError):
        read_airborne("whatever", format="not_a_real_format")


def test_read_airborne_ambiguous_detection_wraps_registry_error():
    marker = "__io_test_ambiguous_marker__"

    def _matches(source):
        return source == marker

    register_airborne_format(
        AirborneFormatDefinition(
            name="io_test_ambig_a", technology="mobilemt", detector=_matches,
        )
    )
    register_airborne_format(
        AirborneFormatDefinition(
            name="io_test_ambig_b", technology="mobilemt", detector=_matches,
        )
    )
    with pytest.raises(AirborneIOError):
        read_airborne(marker)


def test_read_airborne_no_match_message_includes_technology():
    with pytest.raises(AirborneIOError, match="ztem"):
        read_airborne("no_extension_source", technology="ztem")


def test_read_airborne_no_match_without_technology():
    with pytest.raises(AirborneIOError):
        read_airborne("no_extension_source")


def test_read_airborne_unknown_technology_raises():
    register_airborne_format(
        AirborneFormatDefinition(
            name="io_test_known_format",
            technology="mobilemt",
            extensions=(".iotestfmt",),
        )
    )
    with pytest.raises(AirborneIOError):
        read_airborne(
            "sample.iotestfmt",
            technology="not_a_real_technology",
        )


def test_read_airborne_no_registered_reader_raises():
    register_airborne_format(
        AirborneFormatDefinition(
            name="io_test_no_reader",
            technology="mobilemt",
            extensions=(".ionoreader",),
        )
    )
    with pytest.raises(AirborneIOError):
        read_airborne("sample.ionoreader")


def test_read_airborne_reader_returns_wrong_type_raises():
    def _bad_reader(source, **kwargs):
        return {"not": "a dataset"}

    register_airborne_format(
        AirborneFormatDefinition(
            name="io_test_bad_reader",
            technology="mobilemt",
            extensions=(".iobadreader",),
            reader=_bad_reader,
        )
    )
    with pytest.raises(AirborneIOError):
        read_airborne("sample.iobadreader")


def test_read_airborne_happy_path_returns_dataset():
    built = _empty_dataset("from_reader")

    def _good_reader(source, **kwargs):
        return built

    register_airborne_format(
        AirborneFormatDefinition(
            name="io_test_good_reader",
            technology="mobilemt",
            extensions=(".iogoodreader",),
            reader=_good_reader,
        )
    )
    result = read_airborne("sample.iogoodreader")
    assert result is built


def test_read_airborne_passthrough_existing_dataset():
    dataset = _empty_dataset("passthrough")
    assert read_airborne(dataset) is dataset


# ─────────────────────────────────────────────────────────────────────────
# write_airborne
# ─────────────────────────────────────────────────────────────────────────


def test_write_airborne_rejects_non_dataset():
    with pytest.raises(TypeError):
        write_airborne("not-a-dataset", "target.xyz")


def test_write_airborne_format_not_inferable_from_non_path_target():
    with pytest.raises(AirborneIOError):
        write_airborne(_empty_dataset(), target=12345)


def test_write_airborne_format_not_inferable_from_unknown_extension():
    with pytest.raises(AirborneIOError):
        write_airborne(_empty_dataset(), target="output.no_such_extension_xyz")


def test_write_airborne_no_registered_writer_raises():
    register_airborne_format(
        AirborneFormatDefinition(
            name="io_test_no_writer",
            technology="mobilemt",
            extensions=(".iowriteronly",),
        )
    )
    with pytest.raises(AirborneIOError):
        write_airborne(_empty_dataset(), "output.iowriteronly")


def test_write_airborne_happy_path_calls_writer():
    calls = {}

    def _writer(dataset, target, **kwargs):
        calls["dataset"] = dataset
        calls["target"] = target
        return "written"

    register_airborne_format(
        AirborneFormatDefinition(
            name="io_test_good_writer",
            technology="mobilemt",
            extensions=(".iogoodwriter",),
            writer=_writer,
        )
    )
    dataset = _empty_dataset()
    result = write_airborne(dataset, "output.iogoodwriter")
    assert result == "written"
    assert calls["dataset"] is dataset
    assert calls["target"] == "output.iogoodwriter"


# ─────────────────────────────────────────────────────────────────────────
# available_airborne_readers / available_airborne_writers
# ─────────────────────────────────────────────────────────────────────────


def test_available_readers_unknown_technology_raises():
    with pytest.raises(AirborneIOError):
        available_airborne_readers(technology="not_a_real_technology")


def test_available_writers_unknown_technology_raises():
    with pytest.raises(AirborneIOError):
        available_airborne_writers(technology="not_a_real_technology")


def test_available_readers_and_writers_list_registered_formats():
    def _reader(source, **kwargs):
        return _empty_dataset()

    def _writer(dataset, target, **kwargs):
        return None

    register_airborne_format(
        AirborneFormatDefinition(
            name="io_test_rw_format",
            technology="mobilemt",
            extensions=(".iorwformat",),
            reader=_reader,
            writer=_writer,
        )
    )
    assert "io_test_rw_format" in available_airborne_readers(
        technology="mobilemt"
    )
    assert "io_test_rw_format" in available_airborne_writers(
        technology="mobilemt"
    )
    assert "io_test_rw_format" in available_airborne_readers()
    assert "io_test_rw_format" in available_airborne_writers()
