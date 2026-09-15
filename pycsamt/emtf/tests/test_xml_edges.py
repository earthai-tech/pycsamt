from __future__ import annotations

from datetime import datetime
import io
import xml.etree.ElementTree as ET

import numpy as np
import pytest

from pycsamt.emtf import (
    EMTFXMLParseError,
    EMTFXMLSerializationError,
    EMTFXMLWarning,
    EMTFXMLWriteWarning,
    EMTFXMLReader,
)
from pycsamt.emtf.xml.parser import (
    XMLDataTypeSpec,
    child,
    children,
    element_to_mapping,
    load_xml_root,
    local_name,
    parse_numeric_text,
    text,
    warn,
)
from pycsamt.emtf.xml.serializer import (
    EMTFXMLSerializer,
    _append_mapping,
    _append_text,
    _value_text,
)


def test_parser_namespace_navigation_and_mapping_roundtrip():
    root = ET.fromstring(
        '<x:Root xmlns:x="urn:x" a="1"> lead '
        '<x:Item k="v">one</x:Item><x:Item>two</x:Item></x:Root>'
    )
    assert local_name(root.tag) == "Root"
    assert local_name("p:Item") == "Item"
    assert children(None, "Item") == []
    assert len(children(root, "Item")) == 2
    assert text(child(root, "Item")) == "one"
    assert text(root, "Missing") is None

    mapping = element_to_mapping(root)
    assert mapping["@attributes"] == {"a": "1"}
    assert mapping["#text"] == "lead"
    assert mapping["Item"] == [
        {"@attributes": {"k": "v"}, "#text": "one"},
        "two",
    ]

    rebuilt = ET.Element("Parent")
    node = _append_mapping(rebuilt, "Root", mapping)
    assert node.attrib == {"a": "1"}
    assert [item.text for item in node] == ["one", "two"]


def test_xml_source_coercion_supports_all_inputs_and_labels(tmp_path):
    payload = b"<EM_TF><Description>x</Description></EM_TF>"
    root, label = load_xml_root(payload)
    assert root.tag == "EM_TF" and label == "<bytes>"
    root, label = load_xml_root(payload.decode())
    assert root.tag == "EM_TF" and label == "<string>"

    stream = io.BytesIO(payload)
    root, label = load_xml_root(stream)
    assert root.tag == "EM_TF" and label == "<stream>"

    path = tmp_path / "input.xml"
    path.write_bytes(payload)
    root, label = load_xml_root(path)
    assert root.tag == "EM_TF" and label == str(path)


@pytest.mark.parametrize(
    ("source", "message"),
    [
        (b"<broken>", "invalid XML bytes"),
        ("<broken>", "invalid XML string"),
        (io.StringIO("<broken>"), "invalid XML stream"),
        ("definitely-missing.xml", "cannot parse"),
    ],
)
def test_xml_source_errors_are_contextual(source, message):
    with pytest.raises(EMTFXMLParseError, match=message):
        load_xml_root(source)


@pytest.mark.parametrize(
    ("value", "complex_", "expected"),
    [
        ("1D+2", False, 100.0),
        ("1, -2", True, 1 - 2j),
        ("3+4i", True, 3 + 4j),
    ],
)
def test_numeric_parser_accepts_fcu_forms(value, complex_, expected):
    assert parse_numeric_text(value, complex_=complex_) == expected


@pytest.mark.parametrize(
    ("value", "complex_"),
    [(None, False), ("1 2", False), ("bad", False), ("1 2 3", True), ("x y", True)],
)
def test_numeric_parser_rejects_ambiguous_values(value, complex_):
    with pytest.raises(EMTFXMLParseError):
        parse_numeric_text(value, complex_=complex_)


def test_parser_warning_and_datatype_scalar_contract():
    with pytest.warns(EMTFXMLWarning, match="recoverable"):
        warn("recoverable")
    assert XMLDataTypeSpec("S", "scalar", "real").is_scalar
    assert not XMLDataTypeSpec(
        "Z", "impedance", "complex", input_kind="H", output_kind="E"
    ).is_scalar


def test_reader_strict_and_permissive_numeric_diagnostics():
    strict = EMTFXMLReader(strict=True)
    with pytest.raises(EMTFXMLParseError, match="missing required"):
        strict._float(None, "x", required=True)
    with pytest.raises(EMTFXMLParseError, match="invalid numeric"):
        strict._float("bad", "x")
    with pytest.raises(EMTFXMLParseError, match="invalid integer"):
        strict._integer("1.2", "count")

    permissive = EMTFXMLReader(strict=False)
    with pytest.warns(EMTFXMLWarning, match="invalid numeric"):
        assert permissive._float("bad", "x") is None
    with pytest.warns(EMTFXMLWarning, match="invalid integer"):
        assert permissive._integer("bad", "count") is None


def test_permissive_reader_drops_bad_periods_and_unknown_datatypes():
    xml = """
    <EM_TF><Data count="3">
      <Period value="1" units="fortnights"><Z><value output="Ex" input="Hx">1 2</value></Z></Period>
      <Period value="0"><Z><value output="Ex" input="Hx">3 4</value></Z></Period>
      <Period value="2"><UNKNOWN><value>5</value></UNKNOWN></Period>
    </Data></EM_TF>
    """
    with pytest.warns(EMTFXMLWarning) as caught:
        doc = EMTFXMLReader(strict=False).read(xml)
    assert len(caught) >= 3
    np.testing.assert_allclose(doc.periods, [1.0, 2.0])
    assert doc.metadata["xml_data_types_inferred"]
    assert "unknown" not in doc.transfer_functions


def test_serializer_low_level_text_mapping_and_formatting():
    root = ET.Element("Root")
    assert _append_text(root, "Skip", None) is None
    assert _append_text(root, "Skip", "") is None
    empty = _append_text(root, "Empty", "", allow_empty=True, attrs={"a": 1, "b": None})
    assert empty is not None and empty.attrib == {"a": "1"}
    assert _value_text(datetime(2024, 1, 2, 3, 4, 5)) == "2024-01-02T03:04:05"

    serializer = EMTFXMLSerializer(precision=6)
    assert serializer._component_name("T", output="Hz", input_="Hy") == "Ty"
    assert serializer._component_name("Z", output="Ex", input_="Hy") == "Zxy"
    assert serializer._format_numeric(1 + 2j, complex_=True) == "1 2"
    assert serializer._format_numeric(1.25, complex_=False) == "1.25"
    assert serializer._format_float(None) is None
    assert serializer._is_missing(np.nan)
    assert serializer._is_missing(1 + np.nan * 1j)
    assert not serializer._is_missing("text")
    assert serializer._is_finite(1 + 2j, complex_=True)
    assert not serializer._is_finite("text", complex_=False)


def test_serializer_problem_policy():
    with pytest.raises(EMTFXMLSerializationError, match="unsafe"):
        EMTFXMLSerializer(strict=True)._problem("unsafe")
    with pytest.warns(EMTFXMLWriteWarning, match="unsafe"):
        EMTFXMLSerializer(strict=False)._problem("unsafe")
