# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Phase 2 tests for canonical PCBH JSON I/O."""

from __future__ import annotations

import json
from copy import deepcopy
from pathlib import Path

import pytest

from pycsamt.format import read_pcbh, write_pcbh
from pycsamt.format.borehole import (
    Collar,
    CoordinateReferenceSystem,
    LogInterval,
    PCBHBorehole,
    PCBHDocument,
    PCBHValidationError,
    Trajectory,
    VocabularyEntry,
    pcbh_from_dict,
    pcbh_to_dict,
)

DATA_DIR = Path(__file__).parent / "data" / "borehole"
MINIMAL = DATA_DIR / "minimal_vertical.pcbh.json"
MULTIPLE = DATA_DIR / "multiple_deviated.pcbh.json"


def _raw(path: Path = MINIMAL) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _document(**overrides) -> PCBHDocument:
    hole = PCBHBorehole(
        id="BH-001",
        name="Borehole one",
        kind="water",
        status="completed",
        collar=Collar(x=1.0, y=2.0, z=3.0),
        total_depth_md=20.0,
        trajectory=Trajectory(),
        interval_logs={
            "lithology": [
                LogInterval(
                    0.0,
                    20.0,
                    code="GRAN",
                    data_nature="observed",
                )
            ]
        },
    )
    values = {
        "document_id": "example:phase-2",
        "created_at": "2026-08-28T12:00:00Z",
        "created_by": "pytest",
        "crs": CoordinateReferenceSystem(horizontal="LOCAL:test"),
        "boreholes": [hole],
        "lithologies": [
            VocabularyEntry(code="GRAN", name="Granite")
        ],
    }
    values.update(overrides)
    return PCBHDocument(**values)


class TestMappingConversion:
    def test_fixture_builds_valid_document(self):
        document = pcbh_from_dict(_raw())
        assert document.document_id == "example:minimal-vertical"
        assert document.boreholes[0].interval_logs["lithology"][1].code == (
            "GRAN"
        )

    def test_object_mapping_object_round_trip(self):
        original = _document(
            title="Unicode geology: Côte d’Ivoire",
            extensions={"example:project": {"priority": 2}},
        )
        original.boreholes[0].extensions = {
            "example:owner": "field team"
        }
        mapping = pcbh_to_dict(original)
        restored = pcbh_from_dict(mapping)
        assert pcbh_to_dict(restored) == mapping
        assert restored.extensions == original.extensions
        assert restored.boreholes[0].extensions == {
            "example:owner": "field team"
        }

    def test_unknown_regular_field_is_rejected(self):
        value = _raw()
        value["mystery"] = 1
        with pytest.raises(ValueError, match="unsupported field"):
            pcbh_from_dict(value)

    def test_extension_keys_must_be_namespaced(self):
        document = _document(extensions={"plain": 1})
        with pytest.raises(PCBHValidationError, match="namespaced"):
            pcbh_to_dict(document)

    @pytest.mark.parametrize("bad", [float("nan"), float("inf")])
    def test_non_finite_metadata_is_rejected(self, bad):
        document = _document(metadata={"bad": bad})
        with pytest.raises(ValueError, match="non-finite"):
            pcbh_to_dict(document)

    def test_non_json_metadata_is_rejected(self):
        document = _document(metadata={"bad": object()})
        with pytest.raises(ValueError, match="not JSON serializable"):
            pcbh_to_dict(document)


class TestVersionChecks:
    def test_newer_minor_warns_and_loads(self):
        value = _raw()
        value["pcbh_version"] = "0.2.0"
        with pytest.warns(UserWarning, match="newer"):
            document = pcbh_from_dict(value)
        assert document.pcbh_version == "0.2.0"

    def test_unknown_major_is_rejected(self):
        value = _raw()
        value["pcbh_version"] = "1.0.0"
        with pytest.raises(ValueError, match="unsupported pcbh_version"):
            pcbh_from_dict(value)

    @pytest.mark.parametrize("version", ["0.1", "zero", "0.-1.0"])
    def test_malformed_version_is_rejected(self, version):
        value = _raw()
        value["pcbh_version"] = version
        with pytest.raises(ValueError, match="pcbh_version"):
            pcbh_from_dict(value)


class TestResourceLimits:
    def test_borehole_limit(self):
        value = _raw()
        value["boreholes"].append(deepcopy(value["boreholes"][0]))
        with pytest.raises(ValueError, match="boreholes; limit"):
            pcbh_from_dict(value, max_boreholes=1, validate=False)

    def test_interval_limit(self):
        with pytest.raises(ValueError, match="intervals; limit"):
            pcbh_from_dict(_raw(), max_intervals=1, validate=False)

    def test_nesting_limit(self):
        value = _raw()
        value["metadata"] = {"a": {"b": {"c": 1}}}
        with pytest.raises(ValueError, match="nesting"):
            pcbh_from_dict(value, max_nesting=2, validate=False)

    @pytest.mark.parametrize(
        "keyword",
        ["max_boreholes", "max_intervals", "max_nesting"],
    )
    def test_limits_must_be_positive_integers(self, keyword):
        with pytest.raises(ValueError, match=keyword):
            pcbh_from_dict(_raw(), **{keyword: 0})


class TestFileIO:
    def test_read_write_read_round_trip(self, tmp_path):
        original = read_pcbh(MULTIPLE)
        target = tmp_path / "nested" / "roundtrip.pcbh.json"
        assert write_pcbh(original, target) == target
        restored = read_pcbh(target)
        assert pcbh_to_dict(restored) == pcbh_to_dict(original)
        text = target.read_text(encoding="utf-8")
        assert text.startswith('{\n  "$schema"')
        assert text.endswith("\n")

    def test_utf8_is_not_ascii_escaped(self, tmp_path):
        document = _document(title="Forage à Korhogo")
        target = write_pcbh(document, tmp_path / "unicode.pcbh.json")
        assert "Forage à Korhogo" in target.read_text(encoding="utf-8")

    def test_file_size_limit(self):
        size = MINIMAL.stat().st_size
        with pytest.raises(ValueError, match="max_bytes"):
            read_pcbh(MINIMAL, max_bytes=size - 1)

    def test_malformed_json_has_location(self, tmp_path):
        path = tmp_path / "bad.pcbh.json"
        path.write_text('{"pcbh_version": ', encoding="utf-8")
        with pytest.raises(ValueError, match="line 1, column"):
            read_pcbh(path)

    def test_duplicate_json_key_is_rejected(self, tmp_path):
        path = tmp_path / "duplicate.pcbh.json"
        path.write_text('{"a": 1, "a": 2}', encoding="utf-8")
        with pytest.raises(ValueError, match="duplicate JSON object key"):
            read_pcbh(path)

    def test_non_utf8_is_rejected(self, tmp_path):
        path = tmp_path / "binary.pcbh.json"
        path.write_bytes(b"\xff\xfe")
        with pytest.raises(UnicodeDecodeError):
            read_pcbh(path)

    def test_semantically_invalid_file_can_be_inspected_without_validation(
        self,
    ):
        path = DATA_DIR / "invalid_overlap.pcbh.json"
        with pytest.raises(PCBHValidationError, match="overlaps"):
            read_pcbh(path)
        document = read_pcbh(path, validate=False)
        assert document.boreholes[0].id == "BAD-01"

    def test_atomic_failure_preserves_existing_file(
        self, tmp_path, monkeypatch
    ):
        from pycsamt.format.borehole import jsonio

        target = tmp_path / "existing.pcbh.json"
        target.write_text("original", encoding="utf-8")

        def fail_replace(source, destination):
            raise OSError("simulated replacement failure")

        monkeypatch.setattr(jsonio.os, "replace", fail_replace)
        with pytest.raises(OSError, match="simulated"):
            write_pcbh(_document(), target)
        assert target.read_text(encoding="utf-8") == "original"
        assert list(tmp_path.glob("*.tmp")) == []

    def test_indent_must_be_positive_integer(self, tmp_path):
        with pytest.raises(ValueError, match="indent"):
            write_pcbh(_document(), tmp_path / "x.json", indent=0)

