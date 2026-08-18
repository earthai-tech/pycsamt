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
