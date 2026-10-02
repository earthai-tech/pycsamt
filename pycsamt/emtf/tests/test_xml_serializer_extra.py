# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Additional pycsamt.emtf.xml.serializer coverage: metadata edge cases
and spec-fallback branches not exercised by test_xml_writer.py's full
round-trip tests."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.emtf import EMTF, EMTFXMLSerializationError, StatisticalEstimate, TransferFunction
from pycsamt.emtf.xml.serializer import EMTFXMLSerializer, _apply_mapping
from pycsamt.metadata import (
    ChannelMeta,
    CopyrightInfo,
    LocationMeta,
    OrientationMeta,
    Person,
    ProcessingMeta,
    ProvenanceMeta,
    QualityComment,
    Reference,
    RemoteReferenceMeta,
    SiteLayout,
    SiteMeta,
    Software,
    TransferFunctionQuality,
)

import xml.etree.ElementTree as ET


def _to_element(doc: EMTF, *, strict: bool = True) -> ET.Element:
    return EMTFXMLSerializer(strict=strict).to_element(doc)


# ─────────────────────────────────────────────────────────────────────────
# _apply_mapping low-level branches
# ─────────────────────────────────────────────────────────────────────────


def test_apply_mapping_skips_none_attribute_values():
    node = ET.Element("x")
    _apply_mapping(node, {"@attributes": {"a": "1", "b": None}})
    assert node.attrib == {"a": "1"}


def test_apply_mapping_none_scalar_leaves_no_text():
    node = ET.Element("x")
    _apply_mapping(node, None)
    assert node.text is None


# ─────────────────────────────────────────────────────────────────────────
# Header
# ─────────────────────────────────────────────────────────────────────────


def test_header_appends_notes_and_repeated_external_url():
    doc = EMTF(
        metadata={
            "notes": {"#text": "free text"},
            "externalurl": [{"@attributes": {"description": "a"}}, {"@attributes": {"description": "b"}}],
        },
    )
    root = _to_element(doc)
    assert root.find("Notes").text == "free text"
    assert len(root.findall("ExternalUrl")) == 2


# ─────────────────────────────────────────────────────────────────────────
# Provenance / Person
# ─────────────────────────────────────────────────────────────────────────


def test_provenance_with_none_creator_omits_creator_node():
    doc = EMTF(provenance=ProvenanceMeta(create_time="2020-01-01"))
    root = _to_element(doc)
    assert root.find("Provenance/Creator") is None


def test_person_with_all_empty_fields_is_omitted():
    doc = EMTF(
        provenance=ProvenanceMeta(creator=Person(name=None, email=None))
    )
    root = _to_element(doc)
    assert root.find("Provenance/Creator") is None


# ─────────────────────────────────────────────────────────────────────────
# Copyright / RotationInfo
# ─────────────────────────────────────────────────────────────────────────


def test_rotation_info_without_copyright_is_written_at_root():
    doc = EMTF(orientation=OrientationMeta(rotation_info="rotated 15deg"))
    root = _to_element(doc)
    assert root.find("Copyright") is None
    assert root.find("RotationInfo").text == "rotated 15deg"


def test_rotation_info_with_copyright_is_written_inside_copyright_node():
    doc = EMTF(
        copyright=CopyrightInfo(
            release_status="open", conditions_of_use="none",
            reference=Reference(title="T", author="A"),
        ),
        orientation=OrientationMeta(rotation_info="rotated 15deg"),
    )
    root = _to_element(doc)
    assert root.find("RotationInfo") is None
    assert root.find("Copyright/RotationInfo").text == "rotated 15deg"


# ─────────────────────────────────────────────────────────────────────────
# Site / Location / Orientation / Quality
# ─────────────────────────────────────────────────────────────────────────


def test_orientation_only_document_writes_site_node_without_site_fields():
    doc = EMTF(
        orientation=OrientationMeta(mode="sitelayout"),
    )
    root = _to_element(doc)
    site_node = root.find("Site")
    assert site_node is not None
    assert site_node.find("Project") is None
    orientation_node = site_node.find("Orientation")
    assert orientation_node.text == "sitelayout"
    assert "angle_to_geographic_north" not in orientation_node.attrib


def test_run_list_as_sequence_is_joined_with_spaces():
    doc = EMTF(
        site=SiteMeta(site_id="S1", extra={"run_list": ["a", "b", "c"]}),
    )
    root = _to_element(doc)
    assert root.find("Site/RunList").text == "a b c"


def test_site_comments_as_single_string_are_written():
    doc = EMTF(site=SiteMeta(site_id="S1", extra={"comments": "one note"}))
    root = _to_element(doc)
    assert [c.text for c in root.findall("Site/Comments")] == ["one note"]


def test_site_without_location_omits_location_node():
    doc = EMTF(site=SiteMeta(site_id="S1"))
    root = _to_element(doc)
    assert root.find("Site/Location") is None


def test_location_without_datum_omits_datum_attribute():
    doc = EMTF(
        site=SiteMeta(
            site_id="S1",
            location=LocationMeta(latitude=1.0, longitude=2.0, datum=None),
        ),
    )
    root = _to_element(doc)
    location = root.find("Site/Location")
    assert "datum" not in location.attrib


def test_location_with_declination_but_no_epoch_omits_epoch_attribute():
    doc = EMTF(
        site=SiteMeta(
            site_id="S1",
            location=LocationMeta(
                latitude=1.0, longitude=2.0, declination=5.0,
            ),
        ),
    )
    root = _to_element(doc)
    node = root.find("Site/Location/Declination")
    assert node.text is not None
    assert "epoch" not in node.attrib


def test_location_without_declination_omits_declination_node():
    doc = EMTF(
        site=SiteMeta(
            site_id="S1",
            location=LocationMeta(latitude=1.0, longitude=2.0),
        ),
    )
    root = _to_element(doc)
    assert root.find("Site/Location/Declination") is None


def test_quality_with_only_warnings_omits_notes_but_writes_warnings():
    quality = TransferFunctionQuality(
        warning_flag=1,
        warnings=[QualityComment(text="careful", author="Reviewer")],
    )
    doc = EMTF(site=SiteMeta(site_id="S1"), quality=quality)
    root = _to_element(doc)
    assert root.find("Site/DataQualityNotes") is None
    warnings_node = root.find("Site/DataQualityWarnings")
    assert warnings_node.find("Flag").text == "1"
    comment = warnings_node.find("Comments")
    assert comment.attrib["author"] == "Reviewer"


def test_quality_with_only_rating_omits_warnings_node():
    quality = TransferFunctionQuality(rating=4)
    doc = EMTF(site=SiteMeta(site_id="S1"), quality=quality)
    root = _to_element(doc)
    assert root.find("Site/DataQualityNotes/Rating").text == "4"
    assert root.find("Site/DataQualityWarnings") is None


# ─────────────────────────────────────────────────────────────────────────
# FieldNotes
# ─────────────────────────────────────────────────────────────────────────


def test_field_notes_key_becomes_run_attribute():
    doc = EMTF(field_notes={"A": {"#text": "note"}})
    root = _to_element(doc)
    node = root.find("FieldNotes")
    assert node.attrib["run"] == "A"


def test_field_notes_default_run_key_is_not_written_as_attribute():
    doc = EMTF(field_notes={"run_1": {"#text": "note"}})
    root = _to_element(doc)
    node = root.find("FieldNotes")
    assert "run" not in node.attrib


# ─────────────────────────────────────────────────────────────────────────
# Processing
# ─────────────────────────────────────────────────────────────────────────


def test_processing_without_remote_reference_type_omits_remote_ref_node():
    doc = EMTF(
        processing=ProcessingMeta(
            remote_reference=RemoteReferenceMeta(site="BASE01"),
        ),
    )
    root = _to_element(doc)
    assert root.find("ProcessingInfo/RemoteRef") is None
    assert root.find("ProcessingInfo/RemoteInfo/Site/Id").text == "BASE01"


def test_processing_without_software_omits_processing_software_node():
    doc = EMTF(processing=ProcessingMeta(processed_by="Analyst"))
    root = _to_element(doc)
    assert root.find("ProcessingInfo/ProcessingSoftware") is None


def test_processing_sign_convention_unknown_value_passes_through():
    doc = EMTF(processing=ProcessingMeta(sign_convention="unknown convention"))
    root = _to_element(doc)
    assert root.find("ProcessingInfo/SignConvention").text == "unknown convention"


def test_processing_sign_convention_negative_form_is_rewritten():
    doc = EMTF(processing=ProcessingMeta(sign_convention="exp(-i ω t)"))
    root = _to_element(doc)
    assert root.find("ProcessingInfo/SignConvention").text == r"exp(- i\omega t)"


def test_processing_software_without_author_name_omits_author_text():
    doc = EMTF(
        processing=ProcessingMeta(
            software=Software(name="EMTF", author=Person(name=None)),
        ),
    )
    root = _to_element(doc)
    assert root.find("ProcessingInfo/ProcessingSoftware/Author") is None


# ─────────────────────────────────────────────────────────────────────────
# Estimate / DataType declaration edge cases
# ─────────────────────────────────────────────────────────────────────────


def test_estimate_specs_adds_default_declaration_for_undeclared_estimate():
    tf = TransferFunction(
        name="impedance",
        data=np.ones((1, 1, 1), dtype=complex),
        input_channels=(),
        output_channels=(),
        periods=[1.0],
    )
    tf.add_estimate(
        StatisticalEstimate(name="VAR", kind="variance", data=np.ones((1, 1, 1)))
    )
    doc = EMTF(periods=[1.0])
    doc.add_transfer_function(tf)
    root = _to_element(doc)
    node = root.find("StatisticalEstimates/Estimate")
    assert node.attrib["name"] == "VAR"


def test_estimate_declaration_spec_without_name_is_skipped():
    doc = EMTF(metadata={"xml_statistical_estimates": [{"tag": "variance"}]})
    root = _to_element(doc)
    assert root.find("StatisticalEstimates/Estimate") is None


def test_estimate_declaration_without_kind_omits_type_attribute():
    doc = EMTF(
        metadata={
            "xml_statistical_estimates": [{"name": "VAR", "tag": "variance"}],
        },
    )
    root = _to_element(doc)
    node = root.find("StatisticalEstimates/Estimate")
    assert "type" not in node.attrib


def test_data_type_declaration_missing_name_or_tag_is_reported():
    doc = EMTF(metadata={"xml_data_types": [{"name": "Z"}]})
    with pytest.warns(UserWarning, match="needs name and tag"):
        root = _to_element(doc, strict=False)
    assert root.find("DataTypes/DataType") is None


def test_data_type_declaration_uses_legacy_output_input_keys():
    doc = EMTF(
        metadata={
            "xml_data_types": [
                {
                    "name": "Z",
                    "tag": "impedance",
                    "output": "E",
                    "input": "H",
                }
            ],
        },
    )
    root = _to_element(doc)
    node = root.find("DataTypes/DataType")
    assert "type" not in node.attrib  # no data_kind or type key present
    assert node.attrib["output"] == "E"
    assert node.attrib["input"] == "H"


def test_data_type_declaration_intention_other_than_primary_or_derived():
    doc = EMTF(
        metadata={
            "xml_data_types": [
                {
                    "name": "Z",
                    "tag": "impedance",
                    "data_kind": "complex",
                    "intention": "custom",
                }
            ],
        },
    )
    root = _to_element(doc)
    assert root.find("DataTypes/DataType/Intention").text == "custom"


def test_channel_group_with_no_channels_writes_nothing():
    doc = EMTF(
        site_layout=SiteLayout(
            input_channels=[
                ChannelMeta(name="Hx", field_type="magnetic"),
            ],
            output_channels=[],
        ),
    )
    root = _to_element(doc)
    assert root.find("SiteLayout/OutputChannels") is None
    assert root.find("SiteLayout/InputChannels/Magnetic") is not None


# ─────────────────────────────────────────────────────────────────────────
# _spec_from_tf fallback paths
# ─────────────────────────────────────────────────────────────────────────


def test_spec_from_tf_uses_tf_name_when_it_is_xml_safe():
    tf = TransferFunction(
        name="custom_tag",
        data=np.ones((1, 1, 1), dtype=complex),
        input_channels=(),
        output_channels=(),
        periods=[1.0],
    )
    doc = EMTF(periods=[1.0])
    doc.add_transfer_function(tf)
    root = _to_element(doc)
    node = root.find("DataTypes/DataType")
    assert node.attrib["name"] == "CUSTOM_TAG"


def test_spec_from_tf_falls_back_to_unknown_for_unsafe_name():
    tf = TransferFunction(
        name="not a safe xml name!",
        data=np.ones((1, 1, 1), dtype=complex),
        input_channels=(),
        output_channels=(),
        periods=[1.0],
    )
    doc = EMTF(periods=[1.0])
    doc.add_transfer_function(tf)
    with pytest.warns(UserWarning, match="no safe EMTF XML code"):
        root = _to_element(doc, strict=False)
    node = root.find("DataTypes/DataType")
    assert node.attrib["name"] == "UNKNOWN"


def test_spec_from_tf_uses_already_set_xml_data_kind():
    tf = TransferFunction(
        name="custom_tag",
        data=np.ones((1, 1, 1), dtype=complex),
        input_channels=(),
        output_channels=(),
        periods=[1.0],
        attrs={"xml_name": "Z", "xml_data_kind": "real"},
    )
    spec = EMTFXMLSerializer()._spec_from_tf(tf)
    assert spec["data_kind"] == "real"


def test_spec_from_tf_infers_real_data_kind_without_definition():
    tf = TransferFunction(
        name="custom_real",
        data=np.ones((1, 1, 1), dtype=float),
        input_channels=(),
        output_channels=(),
        periods=[1.0],
    )
    doc = EMTF(periods=[1.0])
    doc.add_transfer_function(tf)
    root = _to_element(doc)
    node = root.find("DataTypes/DataType")
    assert node.attrib["type"] == "real"


# ─────────────────────────────────────────────────────────────────────────
# Data matrices
# ─────────────────────────────────────────────────────────────────────────


def test_append_data_skipped_for_empty_document():
    doc = EMTF()
    root = _to_element(doc)
    assert root.find("Data") is None


def test_append_data_reports_missing_periods():
    tf = TransferFunction(
        name="Z",
        data=np.ones((1, 1, 1), dtype=complex),
        input_channels=(),
        output_channels=(),
        periods=[1.0],
    )
    doc = EMTF()
    doc.transfer_functions["impedance"] = tf  # bypass periods bookkeeping
    with pytest.warns(UserWarning, match="without periods"):
        root = _to_element(doc, strict=False)
    assert root.find("Data") is None


def test_matrix_with_non_finite_non_missing_value_is_reported():
    tf = TransferFunction(
        name="impedance",
        data=np.array([[[np.inf + 1j]]]),
        input_channels=("Hx",),
        output_channels=("Ex",),
        periods=[1.0],
    )
    doc = EMTF(periods=[1.0])
    doc.add_transfer_function(tf)
    with pytest.warns(UserWarning, match="non-finite value"):
        root = _to_element(doc, strict=False)
    # the non-finite value produced no <value> children
    assert root.find("Data/Period/impedance") is None


def test_data_writing_falls_back_to_spec_from_tf_when_undeclared():
    # metadata declares an unrelated spec; the actual TF's tag has no
    # matching entry, so _append_data must build one on the fly via
    # _spec_from_tf rather than reusing spec_by_tag.
    tf = TransferFunction(
        name="impedance",
        data=np.ones((1, 1, 1), dtype=complex),
        input_channels=(),
        output_channels=(),
        periods=[1.0],
        attrs={"xml_name": "Z"},
    )
    doc = EMTF(
        periods=[1.0],
        metadata={
            "xml_data_types": [
                {"name": "OTHER", "tag": "other_tag", "data_kind": "complex"}
            ],
        },
    )
    doc.add_transfer_function(tf)
    root = _to_element(doc)
    assert root.find("Data/Period/Z") is not None


def test_data_writing_rejects_invalid_declared_name():
    tf = TransferFunction(
        name="impedance",
        data=np.ones((1, 1, 1), dtype=complex),
        input_channels=(),
        output_channels=(),
        periods=[1.0],
    )
    doc = EMTF(
        periods=[1.0],
        metadata={
            "xml_data_types": [
                {
                    "name": "bad name!",
                    "tag": "impedance",
                    "data_kind": "complex",
                }
            ],
        },
    )
    doc.add_transfer_function(tf)
    with pytest.warns(UserWarning, match="invalid XML data-type name"):
        root = _to_element(doc, strict=False)
    assert len(root.find("Data/Period")) == 0


def test_data_writing_skips_unsupported_estimate_code_silently():
    tf = TransferFunction(
        name="impedance",
        data=np.ones((1, 1, 1), dtype=complex),
        input_channels=(),
        output_channels=(),
        periods=[1.0],
    )
    tf.add_estimate(
        StatisticalEstimate(name="BOGUS", kind="other", data=np.ones((1, 1, 1)))
    )
    doc = EMTF(periods=[1.0])
    doc.add_transfer_function(tf)
    with pytest.warns(UserWarning, match="does not yet define a safe mapping"):
        root = _to_element(doc, strict=False)
    period = root.find("Data/Period")
    assert period.find("impedance.BOGUS") is None


def test_estimate_matrix_rejects_wrong_ndim():
    tf = TransferFunction(
        name="impedance",
        data=np.ones((1, 1, 1), dtype=complex),
        input_channels=(),
        output_channels=(),
        periods=[1.0],
    )
    tf.add_estimate(
        StatisticalEstimate(name="VAR", kind="variance", data=np.ones((1, 1)))
    )
    doc = EMTF(periods=[1.0])
    doc.add_transfer_function(tf)
    with pytest.warns(UserWarning, match="must be a 3-D array"):
        root = _to_element(doc, strict=False)
    assert root.find("Data/Period/impedance.VAR") is None


def test_estimate_matrix_rejects_wrong_shape():
    tf = TransferFunction(
        name="impedance",
        data=np.ones((1, 2, 2), dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Ex", "Ey"),
        periods=[1.0],
    )
    tf.add_estimate(
        StatisticalEstimate(
            name="VAR", kind="variance", data=np.ones((1, 3, 3)),
        )
    )
    doc = EMTF(periods=[1.0])
    doc.add_transfer_function(tf)
    with pytest.warns(UserWarning, match="has shape"):
        root = _to_element(doc, strict=False)
    assert root.find("Data/Period/impedance.VAR") is None


def test_estimate_matrix_rejects_complex_variance_with_nonzero_imag():
    tf = TransferFunction(
        name="impedance",
        data=np.ones((1, 1, 1), dtype=complex),
        input_channels=(),
        output_channels=(),
        periods=[1.0],
    )
    tf.add_estimate(
        StatisticalEstimate(
            name="VAR", kind="variance",
            data=np.array([[[1.0 + 1.0j]]]),
        )
    )
    doc = EMTF(periods=[1.0])
    doc.add_transfer_function(tf)
    with pytest.warns(UserWarning, match="non-zero imaginary"):
        root = _to_element(doc, strict=False)
    assert root.find("Data/Period/impedance.VAR") is None


def test_component_name_uses_single_axis_for_vertical_tipper_output():
    name = EMTFXMLSerializer()._component_name("T", output="Hz", input_="Hx")
    assert name == "Tx"


def test_component_name_falls_back_to_full_name_for_non_axis_channel():
    name = EMTFXMLSerializer()._component_name(
        "Z", output="Aux1", input_="Hx",
    )
    assert name == "Zaux1x"


def test_estimate_matrix_converts_zero_imaginary_variance_to_real():
    tf = TransferFunction(
        name="impedance",
        data=np.ones((1, 1, 1), dtype=complex),
        input_channels=(),
        output_channels=(),
        periods=[1.0],
    )
    estimate = StatisticalEstimate(
        name="VAR", kind="variance",
        data=np.array([[[2.0 + 0.0j]]]),
    )
    tf.add_estimate(estimate)
    matrix = EMTFXMLSerializer()._estimate_matrix(tf, estimate, "VAR", 0)
    assert not np.iscomplexobj(matrix)
    assert matrix[0, 0] == pytest.approx(2.0)


def test_period_range_falls_back_to_stored_metadata_when_no_periods():
    doc = EMTF(metadata={"period_range": {"min": "1.0", "max": "10.0"}})
    root = _to_element(doc)
    node = root.find("PeriodRange")
    assert node.attrib == {"min": "1.0", "max": "10.0"}
