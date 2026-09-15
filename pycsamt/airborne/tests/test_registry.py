from __future__ import annotations

import numpy as np
import pytest

from pycsamt.airborne import (
    AirborneFormatDefinition,
    AirborneTechnologyAmbiguityError,
    detect_airborne_technology,
    get_airborne_technology,
    identify_airborne_technologies,
    list_airborne_technologies,
    register_airborne_format,
)
from pycsamt.airborne import registry as _registry
from pycsamt.airborne.registry import (
    AirborneFormatDetectionError,
    AirborneTechnologyDefinition,
    detect_airborne_format,
    get_airborne_format,
    list_airborne_formats,
    register_airborne_technology,
)
from pycsamt.airborne.afmag import (
    build_airmt_emtf,
    build_original_afmag_emtf,
)
from pycsamt.airborne.mobilemt import build_mobilemt_emtf
from pycsamt.airborne.ztem import build_ztem_emtf
from pycsamt.emtf import EMTF, TransferFunction


def test_builtin_technology_registry():
    names = {item.name for item in list_airborne_technologies()}
    assert {"mobilemt", "ztem", "afmag", "airmt"} <= names
    assert get_airborne_technology("tensor-afmag").name == "airmt"
    assert get_airborne_technology("original afmag").name == "afmag"


def test_technology_detection_from_adapters():
    mobile = build_mobilemt_emtf(
        np.ones((1, 3, 2), dtype=complex),
        frequency=[100.0],
    )
    ztem = build_ztem_emtf(
        np.ones((1, 2), dtype=complex),
        frequency=[90.0],
    )
    airmt = build_airmt_emtf(
        np.ones((1, 3, 2), dtype=complex),
        frequency=[100.0],
    )
    original = build_original_afmag_emtf(
        [1.0],
        frequency=[150.0],
    )
    assert detect_airborne_technology(mobile) == "mobilemt"
    assert detect_airborne_technology(ztem) == "ztem"
    assert detect_airborne_technology(airmt) == "airmt"
    assert detect_airborne_technology(original) == "afmag"


def test_standard_tipper_alone_is_not_called_ztem():
    doc = EMTF(periods=[1.0])
    doc.add_transfer_function(
        TransferFunction(
            name="T",
            data=np.ones((1, 1, 2), dtype=complex),
            input_channels=("Hx", "Hy"),
            output_channels=("Hz",),
            periods=[1.0],
        )
    )
    assert identify_airborne_technologies(doc) == ()
    assert detect_airborne_technology(doc) is None


def test_mixed_technology_is_explicitly_ambiguous():
    mobile = build_mobilemt_emtf(
        np.ones((1, 3, 2), dtype=complex),
        frequency=[100.0],
    )
    ztem = build_ztem_emtf(
        np.ones((1, 2), dtype=complex),
        frequency=[90.0],
    )

    class Mixed:
        attrs = {}
        records = {"a": type("R", (), {"attrs": {}, "emtf": mobile})(),
                   "b": type("R", (), {"attrs": {}, "emtf": ztem})()}

    assert identify_airborne_technologies(Mixed()) == ("mobilemt", "ztem")
    with pytest.raises(AirborneTechnologyAmbiguityError):
        detect_airborne_technology(Mixed())


def test_register_format_requires_known_technology():
    with pytest.raises(ValueError):
        register_airborne_format(
            AirborneFormatDefinition(
                name="test_unknown_technology_format",
                technology="not_a_real_technology",
            )
        )


def test_identify_technologies_handles_objects_with_no_recognizable_attrs():
    class Bare:
        pass

    assert identify_airborne_technologies(Bare()) == ()


def test_identify_technologies_ignores_non_dict_tf_attrs():
    class FakeTF:
        name = "tipper"
        attrs = None

    class Doc:
        attrs = {}
        subtype = None
        transfer_functions = {"tipper": FakeTF()}

    assert identify_airborne_technologies(Doc()) == ()


def test_detect_airborne_technology_non_strict_returns_none_when_mixed():
    mobile = build_mobilemt_emtf(
        np.ones((1, 3, 2), dtype=complex),
        frequency=[100.0],
    )
    ztem = build_ztem_emtf(
        np.ones((1, 2), dtype=complex),
        frequency=[90.0],
    )

    class Mixed:
        attrs = {}
        records = {"a": type("R", (), {"attrs": {}, "emtf": mobile})(),
                   "b": type("R", (), {"attrs": {}, "emtf": ztem})()}

    assert detect_airborne_technology(Mixed(), strict=False) is None


def test_register_format_rewrites_alias_technology_to_canonical_name():
    definition = register_airborne_format(
        AirborneFormatDefinition(
            name="test_alias_technology_format",
            technology="mobile_mt",
        )
    )
    assert definition.technology == "mobilemt"
    assert get_airborne_format("test_alias_technology_format").technology == (
        "mobilemt"
    )


def test_list_airborne_formats_without_filter_returns_everything():
    all_formats = list_airborne_formats()
    assert isinstance(all_formats, tuple)
    assert any(fmt.name == "test_alias_technology_format" for fmt in all_formats)


def test_list_airborne_formats_unknown_technology_raises():
    with pytest.raises(ValueError):
        list_airborne_formats(technology="not_a_real_technology")


def test_detect_airborne_format_returns_none_for_unrecognized_source():
    assert detect_airborne_format(12345) is None
    assert detect_airborne_format("no_extension_file") is None


def test_register_builtin_technologies_is_idempotent():
    before = {item.name for item in list_airborne_technologies()}
    _registry._register_builtin_technologies()
    after = {item.name for item in list_airborne_technologies()}
    assert before == after


def test_technology_definition_rejects_empty_name_label_or_family():
    with pytest.raises(ValueError):
        AirborneTechnologyDefinition(name="", label="X", family="fam")
    with pytest.raises(ValueError):
        AirborneTechnologyDefinition(name="x", label="", family="fam")
    with pytest.raises(ValueError):
        AirborneTechnologyDefinition(name="x", label="X", family="")


def test_format_definition_rejects_empty_name_or_technology():
    with pytest.raises(ValueError):
        AirborneFormatDefinition(name="", technology="mobilemt")
    with pytest.raises(ValueError):
        AirborneFormatDefinition(name="x", technology="")


def test_register_airborne_technology_rejects_wrong_type():
    with pytest.raises(TypeError):
        register_airborne_technology("not-a-definition")


def test_register_airborne_format_rejects_wrong_type():
    with pytest.raises(TypeError):
        register_airborne_format("not-a-definition")


def test_detect_airborne_format_raises_on_ambiguous_detector_match():
    marker = "__test_detector_marker__"

    def _matches_marker(source):
        return source == marker

    register_airborne_format(
        AirborneFormatDefinition(
            name="test_detector_a",
            technology="mobilemt",
            detector=_matches_marker,
        )
    )
    register_airborne_format(
        AirborneFormatDefinition(
            name="test_detector_b",
            technology="mobilemt",
            detector=_matches_marker,
        )
    )
    with pytest.raises(AirborneFormatDetectionError):
        detect_airborne_format(marker)
    # A source that does not match either detector is unaffected by
    # these two now-permanently-registered detectors.
    assert detect_airborne_format("unrelated_source") is None


def test_detect_airborne_format_raises_on_ambiguous_extension_match():
    register_airborne_format(
        AirborneFormatDefinition(
            name="test_extension_a",
            technology="ztem",
            extensions=(".testext",),
        )
    )
    register_airborne_format(
        AirborneFormatDefinition(
            name="test_extension_b",
            technology="ztem",
            extensions=(".testext",),
        )
    )
    with pytest.raises(AirborneFormatDetectionError):
        detect_airborne_format("file.testext")
