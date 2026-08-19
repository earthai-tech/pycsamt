# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Phase 5 EMTF XML writer and scientific round-trip tests."""

from __future__ import annotations

import io
from pathlib import Path

import numpy as np
import pytest

from pycsamt.emtf import (
    EMTF,
    EMTFXMLSerializationError,
    EMTFXMLWriteWarning,
    EMTFXMLWriter,
    StatisticalEstimate,
    TransferFunction,
    write_emtf_xml,
)
from pycsamt.emtf.tests.test_xml_reader import FULL_XML


def _assert_scientific_equal(left: EMTF, right: EMTF) -> None:
    assert left.product_id == right.product_id
    assert left.description == right.description
    assert left.subtype == right.subtype
    assert left.tags == right.tags
    np.testing.assert_allclose(left.periods, right.periods)

    assert left.provenance.create_time == right.provenance.create_time
    assert left.provenance.creating_application == (
        right.provenance.creating_application
    )
    assert left.provenance.creator.name == right.provenance.creator.name
    assert left.provenance.submitter.name == right.provenance.submitter.name
    assert left.copyright.release_status == right.copyright.release_status
    assert (
        left.copyright.conditions_of_use
        == right.copyright.conditions_of_use
    )
    assert left.copyright.reference.title == right.copyright.reference.title
    assert left.copyright.reference.author == right.copyright.reference.author
    assert left.copyright.reference.year == right.copyright.reference.year

    assert left.site.site_id == right.site.site_id
    assert left.site.name == right.site.name
    assert left.site.project == right.site.project
    assert left.site.survey == right.site.survey
    assert left.site.location.latitude == pytest.approx(
        right.site.location.latitude
    )
    assert left.site.location.longitude == pytest.approx(
        right.site.location.longitude
    )
    assert left.site.location.elevation == pytest.approx(
        right.site.location.elevation
    )
    assert left.orientation.mode == right.orientation.mode
    assert left.orientation.angle_to_geographic_north == pytest.approx(
        right.orientation.angle_to_geographic_north
    )
    assert left.site_layout.input_names == right.site_layout.input_names
    assert left.site_layout.output_names == right.site_layout.output_names
    assert left.processing.sign_convention == right.processing.sign_convention
    assert left.processing.processed_by == right.processing.processed_by
    assert left.processing.remote_reference.site == (
        right.processing.remote_reference.site
    )
    assert left.quality.rating == right.quality.rating
    assert left.field_notes == right.field_notes

    assert set(left.transfer_functions) == set(right.transfer_functions)
    for name in left.transfer_functions:
        ltf = left.transfer_functions[name]
        rtf = right.transfer_functions[name]
        assert ltf.input_channels == rtf.input_channels
        assert ltf.output_channels == rtf.output_channels
        assert ltf.units == rtf.units
        np.testing.assert_allclose(ltf.data, rtf.data, equal_nan=True)
        assert set(ltf.estimates) == set(rtf.estimates)
        for key in ltf.estimates:
            lest = ltf.estimates[key]
            rest = rtf.estimates[key]
            assert lest.name == rest.name
            assert lest.kind == rest.kind
            np.testing.assert_allclose(
                lest.data,
                rest.data,
                equal_nan=True,
            )


def test_full_fcu_style_xml_scientific_roundtrip():
    original = EMTF.from_xml(FULL_XML)
    serialized = original.to_xml()
    reread = EMTF.from_xml(serialized)

    _assert_scientific_equal(original, reread)
    assert '<Data count="2">' in serialized
    assert '<Z.INVSIGCOV type="complex" size="2 2">' in serialized
    assert '<T.INVSIGCOV type="complex" size="1 2">' in serialized


def test_missing_component_is_omitted_and_remains_missing():
    original = EMTF.from_xml(FULL_XML)
    serialized = original.to_xml()
    reread = EMTF.from_xml(serialized)

    # Zyy is absent in the second source period. The writer must not invent a
    # zero-valued component while serializing that missing measurement.
    second_period = serialized.split('<Period value="25.6"', 1)[1]
    second_period = second_period.split("</Period>", 1)[0]
    assert 'name="Zyy"' not in second_period
    assert np.isnan(reread.z[1, 1, 1].real)
    assert np.isnan(reread.z[1, 1, 1].imag)


def test_writer_is_deterministic():
    doc = EMTF.from_xml(FULL_XML)
    writer = EMTFXMLWriter()
    assert writer.dumps(doc) == writer.dumps(doc)


def test_metadata_only_roundtrip():
    xml = """
    <EM_TF>
      <Description>metadata only</Description>
      <ProductId>P.S01.2026</ProductId>
      <Tags>impedance</Tags>
      <Site>
        <Project>P</Project>
        <Survey>S</Survey>
        <YearCollected>2026</YearCollected>
        <Id>S01</Id><Name>Site 01</Name>
        <Location datum="WGS84">
          <Latitude>5.1</Latitude><Longitude>-3.2</Longitude>
        </Location>
      </Site>
    </EM_TF>
    """
    doc = EMTF.from_xml(xml)
    out = doc.to_xml()
    reread = EMTF.from_xml(out)
    assert reread.product_id == "P.S01.2026"
    assert reread.is_empty()
    assert reread.site.site_id == "S01"
    assert "<Data " not in out


def test_write_path_and_streams(tmp_path: Path):
    doc = EMTF.from_xml(FULL_XML)
    path = tmp_path / "site.xml"
    returned = write_emtf_xml(doc, path)
    assert returned == path
    assert EMTF.from_xml(path).product_id == doc.product_id

    text_stream = io.StringIO()
    EMTFXMLWriter().write(doc, text_stream)
    assert "<EM_TF>" in text_stream.getvalue()

    binary_stream = io.BytesIO()
    EMTFXMLWriter().write(doc, binary_stream)
    assert b"<EM_TF>" in binary_stream.getvalue()


def test_document_convenience_methods(tmp_path: Path):
    doc = EMTF.from_xml(FULL_XML)
    path = tmp_path / "convenience.xml"
    doc.write(path)
    reread = EMTF.from_xml(path)
    np.testing.assert_allclose(doc.z, reread.z, equal_nan=True)
    edi_path = tmp_path / "convenience.edi"
    with pytest.warns(UserWarning):
        out = doc.write(edi_path, format="edi")
    assert out == edi_path.resolve()
    assert edi_path.exists()


def test_custom_declared_matrix_type_roundtrip():
    xml = """
    <EM_TF>
      <DataTypes>
        <DataType name="Y" type="complex" output="H" input="E" units="[]">
          <Description>Example admittance</Description>
          <Intention>primary data type</Intention>
          <Tag>example_admittance</Tag>
        </DataType>
      </DataTypes>
      <SiteLayout>
        <InputChannels>
          <Electric name="Ex"/><Electric name="Ey"/>
        </InputChannels>
        <OutputChannels>
          <Magnetic name="Hx"/><Magnetic name="Hy"/><Magnetic name="Hz"/>
        </OutputChannels>
      </SiteLayout>
      <Data count="1">
        <Period value="1" units="secs">
          <Y type="complex" size="3 2" units="[]">
            <value output="Hx" input="Ex">1 2</value>
            <value output="Hz" input="Ey">3 4</value>
          </Y>
        </Period>
      </Data>
    </EM_TF>
    """
    doc = EMTF.from_xml(xml)
    out = doc.to_xml()
    reread = EMTF.from_xml(out)
    tf = reread.get_transfer_function("example_admittance")
    assert tf is not None
    assert tf.shape == (1, 3, 2)
    assert tf.data[0, 0, 0] == pytest.approx(1 + 2j)
    assert tf.data[0, 2, 1] == pytest.approx(3 + 4j)
    assert 'name="Y" type="complex" output="H" input="E"' in out


def test_unsupported_estimate_is_rejected_in_strict_mode():
    tf = TransferFunction(
        name="impedance",
        data=np.ones((1, 2, 2), dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Ex", "Ey"),
        periods=[1.0],
    )
    tf.add_estimate(
        StatisticalEstimate(
            name="STD",
            kind="standard_error",
            data=np.ones((1, 2, 2)),
        )
    )
    doc = EMTF(periods=[1.0], transfer_functions={"impedance": tf})
    with pytest.raises(EMTFXMLSerializationError, match="STD"):
        doc.to_xml(strict=True)


def test_unsupported_estimate_can_be_omitted_permissively():
    tf = TransferFunction(
        name="impedance",
        data=np.ones((1, 2, 2), dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Ex", "Ey"),
        periods=[1.0],
    )
    tf.add_estimate(
        StatisticalEstimate(
            name="STD",
            kind="standard_error",
            data=np.ones((1, 2, 2)),
        )
    )
    doc = EMTF(periods=[1.0], transfer_functions={"impedance": tf})
    with pytest.warns(EMTFXMLWriteWarning, match="STD"):
        out = doc.to_xml(strict=False)
    assert ".STD" not in out


def test_writer_recomputes_period_range_from_scientific_periods():
    doc = EMTF.from_xml(FULL_XML)
    out = doc.to_xml()
    assert '<PeriodRange min="7.31429" max="25.6"' in out
