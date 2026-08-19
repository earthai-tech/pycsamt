# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Phase 3 tests for the read-only EMTF XML adapter."""

from __future__ import annotations

import textwrap

import numpy as np
import pytest

from pycsamt.emtf import (
    EMTF,
    EMTFXMLParseError,
    EMTFXMLReader,
    EMTFXMLWarning,
    read_emtf_xml,
)


FULL_XML = """\
<?xml version="1.0" encoding="UTF-8"?>
<EM_TF>
  <Description>Magnetotelluric Transfer Functions</Description>
  <ProductId>USArray.ORL09.2006</ProductId>
  <SubType>MT_TF</SubType>
  <Tags>impedance,tipper</Tags>
  <Provenance>
    <CreateTime>2013-07-19T15:50:21</CreateTime>
    <CreatingApplication>EMTF File Conversion Utilities 3.0</CreatingApplication>
    <Creator>
      <Name>Gary Egbert and Lana Erofeev</Name>
      <Email>egbert@example.edu</Email>
      <Org>Oregon State University</Org>
      <OrgUrl>https://example.edu</OrgUrl>
    </Creator>
    <Submitter><Name>Anna Kelbert</Name></Submitter>
  </Provenance>
  <Copyright>
    <Citation>
      <Title>USArray TA Magnetotelluric Transfer Functions</Title>
      <Authors>Adam Schultz, Gary D. Egbert, Anna Kelbert</Authors>
      <Year>2006</Year>
    </Citation>
    <ReleaseStatus>Unrestricted Release</ReleaseStatus>
    <ConditionsOfUse>Example conditions.</ConditionsOfUse>
  </Copyright>
  <Site>
    <Project>USArray</Project>
    <Survey>TA</Survey>
    <YearCollected>2006</YearCollected>
    <Id>ORL09</Id>
    <Name>Hoppin Springs, OR, USA</Name>
    <Location datum="WGS84">
      <Latitude>42.085064</Latitude>
      <Longitude>-117.552100</Longitude>
      <Elevation units="meters">1978.750</Elevation>
      <Declination epoch="1995.0">15.300</Declination>
    </Location>
    <Orientation angle_to_geographic_north="0.000">orthogonal</Orientation>
    <AcquiredBy>GSY-USA, Inc.</AcquiredBy>
    <Start>2006-10-13T22:50:32</Start>
    <End>2006-10-31T17:47:39</End>
    <RunList>ORL09a ORL09b ORL09c</RunList>
    <DataQualityNotes>
      <Rating>5</Rating>
      <GoodFromPeriod>8.000</GoodFromPeriod>
      <GoodToPeriod>20000.000</GoodToPeriod>
      <Comments author="Gary Egbert">great TF</Comments>
    </DataQualityNotes>
    <DataQualityWarnings>
      <Flag>0</Flag>
      <Comments author="Gary Egbert"></Comments>
    </DataQualityWarnings>
  </Site>
  <FieldNotes run="ORL09b">
    <Instrument><Manufacturer>Barry Narod</Manufacturer><Name>NIMS</Name></Instrument>
    <SamplingRate units="Hz">1.000</SamplingRate>
  </FieldNotes>
  <ProcessingInfo>
    <SignConvention>exp(+ i\\omega t)</SignConvention>
    <RemoteRef type="Robust Remote Reference"/>
    <RemoteInfo><Site><Id>ORJ09</Id></Site></RemoteInfo>
    <ProcessedBy>Gary Egbert and Prasanta Patro</ProcessedBy>
    <ProcessingSoftware>
      <Name>EMTF</Name>
      <LastMod>1998-03-24</LastMod>
      <Author>Gary Egbert</Author>
    </ProcessingSoftware>
    <ProcessingTag>ORL09bc_J9</ProcessingTag>
  </ProcessingInfo>
  <StatisticalEstimates>
    <Estimate name="VAR" type="real"><Tag>variance</Tag></Estimate>
    <Estimate name="INVSIGCOV" type="complex"><Tag>inverse_signal_covariance</Tag></Estimate>
    <Estimate name="RESIDCOV" type="complex"><Tag>residual_covariance</Tag></Estimate>
  </StatisticalEstimates>
  <DataTypes>
    <DataType name="Z" type="complex" output="E" input="H" units="[mV/km]/[nT]">
      <Description>MT impedance</Description>
      <Intention>primary data type</Intention>
      <Tag>impedance</Tag>
    </DataType>
    <DataType name="T" type="complex" output="H" input="H" units="[]">
      <Description>Vertical Field Transfer Functions (Tipper)</Description>
      <Intention>primary data type</Intention>
      <Tag>tipper</Tag>
    </DataType>
  </DataTypes>
  <SiteLayout>
    <InputChannels ref="site" units="m">
      <Magnetic name="Hx" orientation="15.8" x="0.0" y="0.0" z="0.0"/>
      <Magnetic name="Hy" orientation="105.8" x="0.0" y="0.0" z="0.0"/>
    </InputChannels>
    <OutputChannels ref="site" units="m">
      <Magnetic name="Hz" orientation="0.0" x="0.0" y="0.0" z="0.0"/>
      <Electric name="Ex" orientation="15.8" x="-50.0" y="0.0" z="0.0" x2="50.0" y2="0.0" z2="0.0"/>
      <Electric name="Ey" orientation="105.8" x="0.0" y="-50.0" z="0.0" x2="0.0" y2="50.0" z2="0.0"/>
    </OutputChannels>
  </SiteLayout>
  <Data count="2">
    <Period value="7.31429" units="secs">
      <Z type="complex" size="2 2" units="[mV/km]/[nT]">
        <value name="Zxx" output="Ex" input="Hx">-3.092013e0 -2.721066e0</value>
        <value name="Zxy" output="Ex" input="Hy">1.321225e1 9.352477e0</value>
        <value name="Zyx" output="Ey" input="Hx">-4.331750e0 -6.302522e0</value>
        <value name="Zyy" output="Ey" input="Hy">-8.166870e-1 1.109706e0</value>
      </Z>
      <Z.VAR type="real" size="2 2">
        <value name="Zxx" output="Ex" input="Hx">1.128847e-1</value>
        <value name="Zxy" output="Ex" input="Hy">1.580089e-1</value>
        <value name="Zyx" output="Ey" input="Hx">6.515949e-2</value>
        <value name="Zyy" output="Ey" input="Hy">9.120612e-2</value>
      </Z.VAR>
      <Z.INVSIGCOV type="complex" size="2 2">
        <value output="Hx" input="Hx">3.444128e1 -2.222731e-7</value>
        <value output="Hx" input="Hy">1.809302e0 -2.094002e-1</value>
        <value output="Hy" input="Hx">1.809302e0 2.093998e-1</value>
        <value output="Hy" input="Hy">4.820872e1 4.466731e-7</value>
      </Z.INVSIGCOV>
      <Z.RESIDCOV type="complex" size="2 2">
        <value output="Ex" input="Ex">6.555198e-3 0.0</value>
        <value output="Ex" input="Ey">-3.697359e-4 -1.465000e-4</value>
        <value output="Ey" input="Ex">-3.697359e-4 1.465000e-4</value>
        <value output="Ey" input="Ey">3.783802e-3 0.0</value>
      </Z.RESIDCOV>
      <T type="complex" size="1 2" units="[]">
        <value name="Tx" output="Hz" input="Hx">7.031946e-2 -6.107178e-2</value>
        <value name="Ty" output="Hz" input="Hy">3.806275e-1 1.202571e-2</value>
      </T>
      <T.VAR type="real" size="1 2">
        <value name="Tx" output="Hz" input="Hx">5.365952e-4</value>
        <value name="Ty" output="Hz" input="Hy">7.510919e-4</value>
      </T.VAR>
      <T.INVSIGCOV type="complex" size="1 2">
        <value output="Hx" input="Hx">3.444128e1 -2.222731e-7</value>
        <value output="Hx" input="Hy">1.809302e0 -2.094002e-1</value>
        <value output="Hy" input="Hx">1.809302e0 2.093998e-1</value>
        <value output="Hy" input="Hy">4.820872e1 4.466731e-7</value>
      </T.INVSIGCOV>
      <T.RESIDCOV type="complex" size="1 2">
        <value output="Hz" input="Hz">3.116000e-5 0.0</value>
      </T.RESIDCOV>
    </Period>
    <Period value="25.60000" units="secs">
      <Z type="complex" size="2 2" units="[mV/km]/[nT]">
        <value name="Zxx" output="Ex" input="Hx">-2.658900e-1 -1.369750e0</value>
        <value name="Zxy" output="Ex" input="Hy">7.020470e0 7.808172e0</value>
        <value name="Zyx" output="Ey" input="Hx">-1.869531e0 -3.936828e0</value>
      </Z>
      <Z.VAR type="real" size="2 2">
        <value name="Zxx" output="Ex" input="Hx">2.372779e-4</value>
      </Z.VAR>
      <T type="complex" size="1 2" units="[]">
        <value name="Tx" output="Hz" input="Hx">1.0e-1 2.0e-2</value>
        <value name="Ty" output="Hz" input="Hy">2.0e-1 3.0e-2</value>
      </T>
    </Period>
  </Data>
  <PeriodRange min="7.31429" max="25.60000"/>
</EM_TF>
"""


def test_read_full_fcu_style_xml():
    doc = read_emtf_xml(FULL_XML)

    assert doc.product_id == "USArray.ORL09.2006"
    assert doc.site.site_id == "ORL09"
    assert doc.site.location.latitude == pytest.approx(42.085064)
    assert doc.orientation.mode == "orthogonal"
    assert doc.orientation.angle_to_geographic_north == pytest.approx(0.0)
    assert doc.site_layout.input_names == ("Hx", "Hy")
    assert doc.site_layout.output_names == ("Hz", "Ex", "Ey")
    assert doc.processing.sign_convention == "exp(+i ω t)"
    assert doc.processing.remote_reference.site == "ORJ09"
    assert doc.quality.rating == 5
    assert "ORL09b" in doc.field_notes

    np.testing.assert_allclose(doc.periods, [7.31429, 25.6])
    assert doc.impedance.shape == (2, 2, 2)
    assert doc.tipper_tf.shape == (2, 1, 2)
    assert doc.impedance.input_channels == ("Hx", "Hy")
    assert doc.impedance.output_channels == ("Ex", "Ey")

    assert doc.z[0, 0, 1] == pytest.approx(13.21225 + 9.352477j)
    assert doc.z[1, 1, 1].real != doc.z[1, 1, 1].real  # missing -> NaN
    assert doc.tipper[0, 0, 0] == pytest.approx(0.07031946 - 0.06107178j)

    var = doc.impedance.get_estimate("VAR")
    inv = doc.impedance.get_estimate("INVSIGCOV")
    resid = doc.impedance.get_estimate("RESIDCOV")
    assert var.shape == (2, 2, 2)
    assert inv.shape == (2, 2, 2)
    assert resid.shape == (2, 2, 2)
    assert var.data[0, 0, 1] == pytest.approx(1.580089e-1)
    assert inv.data[0, 0, 1] == pytest.approx(1.809302 - 0.2094002j)

    t_inv = doc.tipper_tf.get_estimate("INVSIGCOV")
    t_resid = doc.tipper_tf.get_estimate("RESIDCOV")
    assert t_inv.shape == (2, 2, 2)
    assert t_resid.shape == (2, 1, 1)


def test_existing_z_and_tipper_compatibility_objects_are_available():
    doc = EMTF.from_xml(FULL_XML)
    assert doc.Z is not None
    assert doc.Tip is not None
    np.testing.assert_allclose(doc.Z.z[:, 0, 1], doc.z[:, 0, 1])


def test_frequency_units_are_converted_to_periods():
    xml = """
    <EM_TF>
      <DataTypes>
        <DataType name="Z" type="complex" output="E" input="H">
          <Tag>impedance</Tag>
        </DataType>
      </DataTypes>
      <SiteLayout>
        <InputChannels><Magnetic name="Hx"/><Magnetic name="Hy"/></InputChannels>
        <OutputChannels><Electric name="Ex"/><Electric name="Ey"/></OutputChannels>
      </SiteLayout>
      <Data count="1">
        <Period value="10" units="Hz">
          <Z><value output="Ex" input="Hy">1 2</value></Z>
        </Period>
      </Data>
    </EM_TF>
    """
    doc = read_emtf_xml(textwrap.dedent(xml))
    assert doc.periods[0] == pytest.approx(0.1)
    assert doc.frequency[0] == pytest.approx(10.0)


def test_reader_is_namespace_safe():
    xml = FULL_XML.replace(
        "<EM_TF>", '<EM_TF xmlns="http://example.org/emtf">', 1
    )
    doc = read_emtf_xml(xml)
    assert doc.product_id == "USArray.ORL09.2006"
    assert doc.z.shape == (2, 2, 2)


def test_missing_component_is_nan_not_zero():
    doc = read_emtf_xml(FULL_XML)
    assert np.isnan(doc.z[1, 1, 1].real)
    assert np.isnan(doc.z[1, 1, 1].imag)


def test_reader_can_infer_registered_types_without_datatypes_block():
    xml = """
    <EM_TF>
      <SiteLayout>
        <InputChannels><Magnetic name="Hx"/><Magnetic name="Hy"/></InputChannels>
        <OutputChannels><Electric name="Ex"/><Electric name="Ey"/></OutputChannels>
      </SiteLayout>
      <Data count="1">
        <Period value="1" units="secs">
          <Z><value output="Ex" input="Hy">1 2</value></Z>
        </Period>
      </Data>
    </EM_TF>
    """
    doc = read_emtf_xml(textwrap.dedent(xml))
    assert doc.impedance is not None
    assert doc.z[0, 0, 1] == pytest.approx(1 + 2j)


def test_reader_infers_channels_when_layout_is_absent():
    xml = """
    <EM_TF>
      <DataTypes>
        <DataType name="Z" type="complex" output="E" input="H">
          <Tag>impedance</Tag>
        </DataType>
      </DataTypes>
      <Data count="1">
        <Period value="1" units="secs">
          <Z>
            <value output="Ex" input="Hx">1 0</value>
            <value output="Ex" input="Hy">2 0</value>
            <value output="Ey" input="Hx">3 0</value>
            <value output="Ey" input="Hy">4 0</value>
          </Z>
        </Period>
      </Data>
    </EM_TF>
    """
    doc = read_emtf_xml(textwrap.dedent(xml))
    assert doc.impedance.input_channels == ("Hx", "Hy")
    assert doc.impedance.output_channels == ("Ex", "Ey")


def test_strict_reader_rejects_unknown_channel_reference():
    bad = FULL_XML.replace(
        'output="Ex" input="Hx">-3.092013e0 -2.721066e0',
        'output="BAD" input="Hx">-3.092013e0 -2.721066e0',
        1,
    )
    with pytest.raises(EMTFXMLParseError, match="unknown channel"):
        EMTFXMLReader(strict=True).read(bad)


def test_permissive_reader_skips_unknown_channel_reference():
    bad = FULL_XML.replace(
        'output="Ex" input="Hx">-3.092013e0 -2.721066e0',
        'output="BAD" input="Hx">-3.092013e0 -2.721066e0',
        1,
    )
    with pytest.warns(EMTFXMLWarning, match="unknown channel"):
        doc = EMTFXMLReader(strict=False).read(bad)
    assert np.isnan(doc.z[0, 0, 0].real)


def test_metadata_only_document_is_valid():
    xml = """
    <EM_TF>
      <Description>Magnetotelluric Transfer Functions</Description>
      <ProductId>GrandProject.ID001.2006</ProductId>
      <SubType>MT_TF</SubType>
      <Tags>impedance,tipper</Tags>
      <Site>
        <Project>GrandProject</Project>
        <Survey>Example</Survey>
        <YearCollected>2006</YearCollected>
        <Id>ID001</Id>
        <Name>Example Site</Name>
        <Location datum="WGS84">
          <Latitude>42.0</Latitude><Longitude>-117.0</Longitude>
          <Elevation units="meters">1000</Elevation>
        </Location>
      </Site>
    </EM_TF>
    """
    doc = read_emtf_xml(textwrap.dedent(xml))
    assert doc.product_id == "GrandProject.ID001.2006"
    assert doc.is_empty()
    assert doc.periods is None


def test_wrong_root_is_rejected():
    with pytest.raises(EMTFXMLParseError, match="root element"):
        read_emtf_xml("<Configuration/>")
