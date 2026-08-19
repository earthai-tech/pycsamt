from __future__ import annotations

import pytest

from pycsamt.emtf.datatypes import (
    DataTypeDefinition,
    get_emtf_datatype,
    list_emtf_datatypes,
    register_emtf_datatype,
)


def test_primary_registry_matches_fcu_v41_names_and_aliases():
    z = get_emtf_datatype("impedance")
    t = get_emtf_datatype("T")
    zi = get_emtf_datatype("ZI")
    ti = get_emtf_datatype("TI")

    assert z is not None and z.name == "Z" and z.tag == "impedance"
    assert t is not None and t.tag == "tipper"
    assert zi is not None and zi.tag == "interstation_impedance"
    assert ti is not None and ti.tag == "interstation_transfer_functions"

    # aliases preserve compatibility with the 2020 publication table
    assert get_emtf_datatype("Q") is zi
    assert get_emtf_datatype("P") is ti


def test_duplicate_short_name_keeps_impedance_as_default():
    # FCU defines both full and off-diagonal impedance with short code Z.
    assert get_emtf_datatype("Z").tag == "impedance"
    assert (
        get_emtf_datatype("off_diagonal_impedance").name
        == "Z"
    )


def test_registry_contains_primary_and_derived_types():
    primary = list_emtf_datatypes(intention="primary")
    derived = list_emtf_datatypes(intention="derived")
    assert {"impedance", "tipper"}.issubset(primary)
    assert {"apparent_resistivity", "phase_tensor"}.issubset(derived)
    assert get_emtf_datatype("ZEFF").tag == "effective_impedance"
    assert get_emtf_datatype("ZDET").tag == "impedance_determinant"


def test_custom_registration_and_validation():
    definition = DataTypeDefinition(
        name="YTEST",
        tag="test_admittance",
        data_kind="complex",
        input_kind="E",
        output_kind="H",
        units="[test]",
    )
    register_emtf_datatype(definition)
    assert get_emtf_datatype("YTEST") is definition
    assert get_emtf_datatype("test_admittance") is definition

    with pytest.raises(ValueError):
        register_emtf_datatype(definition)

    with pytest.raises(ValueError):
        DataTypeDefinition(
            name="BAD",
            tag="bad",
            data_kind="object",
        )
