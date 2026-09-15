# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Additional pycsamt.emtf.xml.reader coverage: metadata edge cases,
permissive-mode recoverable problems, and the scalar/channel-inference
branches not exercised by test_xml_reader.py's full-document round
trip."""

from __future__ import annotations

import textwrap

import numpy as np
import pytest

from pycsamt.emtf import EMTFXMLParseError, EMTFXMLReader, EMTFXMLWarning


def _read(xml: str, *, strict: bool = True):
    return EMTFXMLReader(strict=strict).read(textwrap.dedent(xml))


# ─────────────────────────────────────────────────────────────────────────
# Provenance / Copyright
# ─────────────────────────────────────────────────────────────────────────


def test_creator_with_all_empty_fields_is_none():
    doc = _read(
        """
        <EM_TF>
          <Provenance><Creator></Creator></Provenance>
        </EM_TF>
        """
    )
    assert doc.provenance.creator is None


def test_copyright_handles_invalid_year_and_doi_forms():
    doc = _read(
        """
        <EM_TF>
          <Copyright>
            <Citation>
              <Title>T</Title>
              <Authors>A</Authors>
              <Year>not-a-year</Year>
              <DOI>doi:10.1234/abc</DOI>
              <SurveyDOI>10.9999/survey</SurveyDOI>
            </Citation>
            <SelectedPublications>Pub1</SelectedPublications>
            <Acknowledgement>Thanks</Acknowledgement>
            <AdditionalInfo>Extra</AdditionalInfo>
          </Copyright>
        </EM_TF>
        """
    )
    ref = doc.copyright.reference
    assert ref.year is None
    assert ref.extra["year_text"] == "not-a-year"
    assert ref.doi == "10.1234/abc"  # "doi:" prefix stripped
    assert ref.extra["survey_doi"] == "10.9999/survey"
    assert ref.extra["selectedpublications"] == "Pub1"
    assert ref.extra["acknowledgement"] == "Thanks"
    assert ref.extra["additionalinfo"] == "Extra"


def test_copyright_rejects_malformed_doi_into_extra():
    doc = _read(
        """
        <EM_TF>
          <Copyright><Citation><DOI>not-a-doi</DOI></Citation></Copyright>
        </EM_TF>
        """
    )
    ref = doc.copyright.reference
    assert ref.doi is None
    assert ref.extra["doi"] == "not-a-doi"


# ─────────────────────────────────────────────────────────────────────────
# Site / Location
# ─────────────────────────────────────────────────────────────────────────


def test_site_comments_are_collected_into_extra():
    doc = _read(
        """
        <EM_TF>
          <Site>
            <Id>S1</Id>
            <Comments>first note</Comments>
            <Comments>second note</Comments>
          </Site>
        </EM_TF>
        """
    )
    assert doc.site.extra["comments"] == ["first note", "second note"]


def test_site_without_location_has_none_location():
    doc = _read("<EM_TF><Site><Id>S1</Id></Site></EM_TF>")
    assert doc.site.location is None


def test_invalid_site_metadata_is_reported(strict=True):
    with pytest.raises(EMTFXMLParseError, match="invalid <Site>"):
        _read(
            "<EM_TF><Site><Id>S1</Id>"
            "<YearCollected>-5</YearCollected></Site></EM_TF>",
            strict=True,
        )


def test_permissive_reader_drops_invalid_site_metadata():
    with pytest.warns(EMTFXMLWarning, match="invalid <Site>"):
        doc = _read(
            "<EM_TF><Site><Id>S1</Id>"
            "<YearCollected>-5</YearCollected></Site></EM_TF>",
            strict=False,
        )
    assert doc.site is None


def test_invalid_location_metadata_is_reported():
    with pytest.raises(EMTFXMLParseError, match="invalid <Location>"):
        _read(
            """
            <EM_TF><Site><Id>S1</Id>
              <Location><Latitude>999</Latitude></Location>
            </Site></EM_TF>
            """,
            strict=True,
        )


# ─────────────────────────────────────────────────────────────────────────
# Orientation / Quality
# ─────────────────────────────────────────────────────────────────────────


def test_orientation_falls_back_to_rotation_info_without_orientation_node():
    doc = _read(
        """
        <EM_TF>
          <Copyright><RotationInfo>rotated to strike</RotationInfo></Copyright>
          <Site><Id>S1</Id></Site>
        </EM_TF>
        """
    )
    assert doc.orientation.rotation_info == "rotated to strike"
    assert doc.orientation.mode is None


def test_invalid_orientation_mode_is_reported():
    with pytest.raises(EMTFXMLParseError, match="invalid orientation"):
        _read(
            "<EM_TF><Site><Id>S1</Id>"
            "<Orientation>bogus_mode</Orientation></Site></EM_TF>",
            strict=True,
        )


def test_invalid_quality_rating_is_reported():
    with pytest.raises(EMTFXMLParseError, match="invalid data-quality"):
        _read(
            """
            <EM_TF><Site><Id>S1</Id>
              <DataQualityNotes><Rating>9</Rating></DataQualityNotes>
            </Site></EM_TF>
            """,
            strict=True,
        )


# ─────────────────────────────────────────────────────────────────────────
# Processing / RemoteRef
# ─────────────────────────────────────────────────────────────────────────


def test_processing_without_remote_ref_or_info_has_no_remote_reference():
    doc = _read(
        "<EM_TF><ProcessingInfo><ProcessedBy>X</ProcessedBy>"
        "</ProcessingInfo></EM_TF>"
    )
    assert doc.processing.remote_reference is None


def test_processing_records_process_date_in_extra():
    doc = _read(
        "<EM_TF><ProcessingInfo><ProcessDate>2020-01-01</ProcessDate>"
        "</ProcessingInfo></EM_TF>"
    )
    assert doc.processing.extra["process_date"] == "2020-01-01"


# ─────────────────────────────────────────────────────────────────────────
# SiteLayout / channels
# ─────────────────────────────────────────────────────────────────────────


def test_no_sitelayout_and_no_legacy_channels_gives_none_layout():
    doc = _read("<EM_TF><Site><Id>S1</Id></Site></EM_TF>")
    assert doc.site_layout is None


def test_legacy_xml_without_sitelayout_wrapper_still_reads_channels():
    with pytest.warns(EMTFXMLWarning, match="legacy EMTF XML"):
        doc = _read(
            """
            <EM_TF>
              <InputChannels><Magnetic name="Hx"/></InputChannels>
              <OutputChannels><Electric name="Ex"/></OutputChannels>
            </EM_TF>
            """,
            strict=False,
        )
    assert doc.site_layout.input_names == ("Hx",)
    assert doc.site_layout.extra["legacy_without_sitelayout"] is True


def test_sitelayout_rejects_duplicate_channel_names_case_insensitively():
    with pytest.raises(EMTFXMLParseError, match="invalid <SiteLayout>"):
        _read(
            """
            <EM_TF><SiteLayout>
              <InputChannels>
                <Magnetic name="Hx"/><Magnetic name="hx"/>
              </InputChannels>
            </SiteLayout></EM_TF>
            """,
            strict=True,
        )


def test_sitelayout_with_missing_channel_groups_yields_empty_lists():
    doc = _read("<EM_TF><SiteLayout></SiteLayout></EM_TF>")
    assert doc.site_layout.input_channels == []
    assert doc.site_layout.output_channels == []


def test_unknown_channel_element_kind_is_skipped():
    doc = _read(
        """
        <EM_TF><SiteLayout>
          <InputChannels>
            <Other name="X"/><Magnetic name="Hx"/>
          </InputChannels>
        </SiteLayout></EM_TF>
        """
    )
    assert doc.site_layout.input_names == ("Hx",)


def test_channel_missing_name_is_reported_and_skipped():
    with pytest.warns(EMTFXMLWarning, match="missing required name"):
        doc = _read(
            """
            <EM_TF><SiteLayout>
              <InputChannels>
                <Magnetic orientation="0"/><Magnetic name="Hx"/>
              </InputChannels>
            </SiteLayout></EM_TF>
            """,
            strict=False,
        )
    assert doc.site_layout.input_names == ("Hx",)


def test_channel_with_non_finite_orientation_is_reported():
    with pytest.raises(EMTFXMLParseError, match="invalid channel"):
        _read(
            """
            <EM_TF><SiteLayout>
              <InputChannels><Magnetic name="Hx" orientation="nan"/></InputChannels>
            </SiteLayout></EM_TF>
            """,
            strict=True,
        )


# ─────────────────────────────────────────────────────────────────────────
# FieldNotes
# ─────────────────────────────────────────────────────────────────────────


def test_duplicate_field_notes_run_key_becomes_a_list():
    doc = _read(
        """
        <EM_TF>
          <FieldNotes run="A"><SamplingRate units="Hz">1</SamplingRate></FieldNotes>
          <FieldNotes run="A"><SamplingRate units="Hz">2</SamplingRate></FieldNotes>
          <FieldNotes run="A"><SamplingRate units="Hz">3</SamplingRate></FieldNotes>
        </EM_TF>
        """
    )
    assert isinstance(doc.field_notes["A"], list)
    assert len(doc.field_notes["A"]) == 3


def test_auxiliary_metadata_collects_notes_and_repeated_urls():
    doc = _read(
        """
        <EM_TF>
          <Notes>free text</Notes>
          <ExternalUrl>http://a</ExternalUrl>
          <ExternalUrl>http://b</ExternalUrl>
        </EM_TF>
        """
    )
    assert doc.metadata["notes"]["#text"] == "free text"
    assert len(doc.metadata["externalurl"]) == 2


# ─────────────────────────────────────────────────────────────────────────
# DataType declarations
# ─────────────────────────────────────────────────────────────────────────


def test_datatype_name_and_tag_fall_back_to_registry():
    # No name= attribute, but the Tag matches a registered type -> name filled.
    doc = _read(
        """
        <EM_TF><DataTypes>
          <DataType type="complex"><Tag>impedance</Tag></DataType>
        </DataTypes></EM_TF>
        """
    )
    assert doc.metadata["xml_data_types"][0]["name"] == "Z"


def test_datatype_tag_falls_back_to_registry_from_name():
    # No <Tag>, but name="Z" matches a registered type -> tag filled.
    doc = _read(
        '<EM_TF><DataTypes><DataType name="Z" type="complex">'
        "</DataType></DataTypes></EM_TF>"
    )
    assert doc.metadata["xml_data_types"][0]["tag"] == "impedance"


def test_datatype_with_neither_name_nor_resolvable_tag_is_dropped():
    with pytest.warns(EMTFXMLWarning, match="requires a name"):
        doc = _read(
            "<EM_TF><DataTypes><DataType type='complex'></DataType>"
            "</DataTypes></EM_TF>",
            strict=False,
        )
    assert doc.metadata["xml_data_types"] == []


def test_datatype_invalid_type_falls_back_to_registered_kind():
    with pytest.warns(EMTFXMLWarning, match="invalid/missing type"):
        doc = _read(
            '<EM_TF><DataTypes><DataType name="Z" type="bogus">'
            "<Tag>impedance</Tag></DataType></DataTypes></EM_TF>",
            strict=False,
        )
    assert doc.metadata["xml_data_types"][0]["data_kind"] == "complex"


def test_duplicate_datatype_declaration_is_reported():
    with pytest.warns(EMTFXMLWarning, match="duplicate DataType"):
        doc = _read(
            """
            <EM_TF><DataTypes>
              <DataType name="Z" type="complex"><Tag>impedance</Tag></DataType>
              <DataType name="Z" type="complex"><Tag>impedance</Tag></DataType>
            </DataTypes></EM_TF>
            """,
            strict=False,
        )
    assert len(doc.metadata["xml_data_types"]) == 1


# ─────────────────────────────────────────────────────────────────────────
# Periods
# ─────────────────────────────────────────────────────────────────────────


def test_data_count_mismatch_is_reported_permissively():
    with pytest.warns(EMTFXMLWarning, match="does not match"):
        doc = _read(
            """
            <EM_TF><Data count="5">
              <Period value="1"><Z><value output="Ex" input="Hx">1 2</value></Z></Period>
            </Data></EM_TF>
            """,
            strict=False,
        )
    assert doc.periods is not None


def test_all_periods_invalid_yields_no_usable_periods():
    with pytest.warns(EMTFXMLWarning, match="no usable"):
        doc = _read(
            '<EM_TF><Data count="1"><Period value="0"/></Data></EM_TF>',
            strict=False,
        )
    assert doc.periods is None


def test_infer_data_type_specs_skips_dotted_children_and_dedupes():
    # Two periods both carrying <Z>; the second must not re-add "Z", and
    # the dotted <Z.VAR> child must never become its own inferred type.
    doc = _read(
        """
        <EM_TF><Data count="2">
          <Period value="1">
            <Z><value output="Ex" input="Hx">1 2</value></Z>
            <Z.VAR><value output="Ex" input="Hx">0.1</value></Z.VAR>
          </Period>
          <Period value="2">
            <Z><value output="Ex" input="Hx">3 4</value></Z>
          </Period>
        </Data></EM_TF>
        """,
        strict=False,
    )
    names = [d["name"] for d in doc.metadata["xml_data_types"]]
    assert names.count("Z") == 1
    assert "VAR" not in names


# ─────────────────────────────────────────────────────────────────────────
# Transfer function / channel resolution
# ─────────────────────────────────────────────────────────────────────────


def test_cannot_determine_channels_without_layout_or_value_attrs():
    with pytest.warns(EMTFXMLWarning, match="cannot determine input/output"):
        doc = _read(
            """
            <EM_TF><DataTypes>
              <DataType name="Z" type="complex" output="E" input="H">
                <Tag>impedance</Tag>
              </DataType>
            </DataTypes>
            <Data count="1">
              <Period value="1"><Z><value>1 2</value></Z></Period>
            </Data></EM_TF>
            """,
            strict=False,
        )
    assert doc.get_transfer_function("impedance") is None


def test_declared_datatype_never_present_in_any_period_yields_no_tf():
    doc = _read(
        """
        <EM_TF><DataTypes>
          <DataType name="Z" type="complex" output="E" input="H">
            <Tag>impedance</Tag>
          </DataType>
        </DataTypes>
        <SiteLayout>
          <InputChannels><Magnetic name="Hx"/></InputChannels>
          <OutputChannels><Electric name="Ex"/></OutputChannels>
        </SiteLayout>
        <Data count="1">
          <Period value="1"><T><value output="Hx" input="Hx">1 2</value></T></Period>
        </Data></EM_TF>
        """,
        strict=False,
    )
    assert doc.get_transfer_function("impedance") is None


def test_channel_names_for_spec_with_no_output_kind_returns_empty_tuple():
    # input="H" resolves to Hx via layout; output kind is unset entirely,
    # so family() short-circuits on the falsy-kind branch for the output
    # side while the type remains non-scalar (input_kind is set).
    doc = _read(
        """
        <EM_TF><DataTypes>
          <DataType name="X" type="complex" input="H">
            <Tag>custom_x</Tag>
          </DataType>
        </DataTypes>
        <SiteLayout>
          <InputChannels><Magnetic name="Hx"/></InputChannels>
          <OutputChannels><Electric name="Ex"/></OutputChannels>
        </SiteLayout>
        <Data count="1">
          <Period value="1">
            <X><value output="Ex" input="Hx">1 2</value></X>
          </Period>
        </Data></EM_TF>
        """,
        strict=False,
    )
    tf = doc.get_transfer_function("custom_x")
    # output channels come from _infer_component_channels since the
    # declared spec's own output_kind resolved to no channels.
    assert tf.input_channels == ("Hx",)


def test_channel_kind_other_than_h_or_e_returns_every_channel():
    doc = _read(
        """
        <EM_TF><DataTypes>
          <DataType name="X" type="complex" output="ALL" input="H">
            <Tag>custom_x</Tag>
          </DataType>
        </DataTypes>
        <SiteLayout>
          <InputChannels><Magnetic name="Hx"/></InputChannels>
          <OutputChannels>
            <Magnetic name="Hz"/><Electric name="Ex"/>
          </OutputChannels>
        </SiteLayout>
        <Data count="1">
          <Period value="1">
            <X>
              <value output="Hz" input="Hx">1 2</value>
              <value output="Ex" input="Hx">3 4</value>
            </X>
          </Period>
        </Data></EM_TF>
        """,
        strict=False,
    )
    tf = doc.get_transfer_function("custom_x")
    assert tf.output_channels == ("Hz", "Ex")


def test_matrix_value_missing_channel_attrs_is_reported():
    with pytest.warns(EMTFXMLWarning, match="requires input and output"):
        _read(
            """
            <EM_TF><DataTypes>
              <DataType name="Z" type="complex" output="E" input="H">
                <Tag>impedance</Tag>
              </DataType>
            </DataTypes>
            <SiteLayout>
              <InputChannels><Magnetic name="Hx"/></InputChannels>
              <OutputChannels><Electric name="Ex"/></OutputChannels>
            </SiteLayout>
            <Data count="1">
              <Period value="1"><Z><value>1 2</value></Z></Period>
            </Data></EM_TF>
            """,
            strict=False,
        )


def test_matrix_value_with_duplicate_component_is_reported():
    with pytest.warns(EMTFXMLWarning, match="duplicate component"):
        doc = _read(
            """
            <EM_TF><DataTypes>
              <DataType name="Z" type="complex" output="E" input="H">
                <Tag>impedance</Tag>
              </DataType>
            </DataTypes>
            <SiteLayout>
              <InputChannels><Magnetic name="Hx"/></InputChannels>
              <OutputChannels><Electric name="Ex"/></OutputChannels>
            </SiteLayout>
            <Data count="1">
              <Period value="1">
                <Z>
                  <value output="Ex" input="Hx">1 2</value>
                  <value output="Ex" input="Hx">3 4</value>
                </Z>
              </Period>
            </Data></EM_TF>
            """,
            strict=False,
        )
    tf = doc.get_transfer_function("impedance")
    # last value wins
    assert tf.data[0, 0, 0] == pytest.approx(3 + 4j)


def test_matrix_value_with_bad_numeric_text_is_reported():
    with pytest.warns(EMTFXMLWarning):
        doc = _read(
            """
            <EM_TF><DataTypes>
              <DataType name="Z" type="complex" output="E" input="H">
                <Tag>impedance</Tag>
              </DataType>
            </DataTypes>
            <SiteLayout>
              <InputChannels><Magnetic name="Hx"/></InputChannels>
              <OutputChannels><Electric name="Ex"/></OutputChannels>
            </SiteLayout>
            <Data count="1">
              <Period value="1">
                <Z><value output="Ex" input="Hx">not numeric</value></Z>
              </Period>
            </Data></EM_TF>
            """,
            strict=False,
        )
    tf = doc.get_transfer_function("impedance")
    assert np.isnan(tf.data[0, 0, 0].real)


# ─────────────────────────────────────────────────────────────────────────
# Scalar data types
# ─────────────────────────────────────────────────────────────────────────


def test_scalar_datatype_reads_a_single_value():
    doc = _read(
        """
        <EM_TF><DataTypes>
          <DataType name="ZSTRIKE" type="real">
            <Tag>impedance_strike</Tag>
          </DataType>
        </DataTypes>
        <Data count="1">
          <Period value="1"><ZSTRIKE><value>12.5</value></ZSTRIKE></Period>
        </Data></EM_TF>
        """
    )
    tf = doc.get_transfer_function("impedance_strike")
    assert tf.data[0, 0, 0] == pytest.approx(12.5)


def test_scalar_datatype_with_no_values_stays_nan():
    doc = _read(
        """
        <EM_TF><DataTypes>
          <DataType name="ZSTRIKE" type="real">
            <Tag>impedance_strike</Tag>
          </DataType>
        </DataTypes>
        <Data count="1">
          <Period value="1"><ZSTRIKE></ZSTRIKE></Period>
        </Data></EM_TF>
        """
    )
    tf = doc.get_transfer_function("impedance_strike")
    assert np.isnan(tf.data[0, 0, 0])


def test_scalar_datatype_with_bad_numeric_text_is_reported():
    with pytest.warns(EMTFXMLWarning):
        doc = _read(
            """
            <EM_TF><DataTypes>
              <DataType name="ZSTRIKE" type="real">
                <Tag>impedance_strike</Tag>
              </DataType>
            </DataTypes>
            <Data count="1">
              <Period value="1"><ZSTRIKE><value>abc</value></ZSTRIKE></Period>
            </Data></EM_TF>
            """,
            strict=False,
        )
    tf = doc.get_transfer_function("impedance_strike")
    assert np.isnan(tf.data[0, 0, 0])


def test_scalar_datatype_with_multiple_values_is_reported():
    with pytest.warns(EMTFXMLWarning, match="multiple scalar values"):
        doc = _read(
            """
            <EM_TF><DataTypes>
              <DataType name="ZSTRIKE" type="real">
                <Tag>impedance_strike</Tag>
              </DataType>
            </DataTypes>
            <Data count="1">
              <Period value="1">
                <ZSTRIKE><value>1.0</value><value>2.0</value></ZSTRIKE>
              </Period>
            </Data></EM_TF>
            """,
            strict=False,
        )
    tf = doc.get_transfer_function("impedance_strike")
    assert tf.data[0, 0, 0] == pytest.approx(2.0)  # last value wins


def test_scalar_datatype_has_no_invsigcov_or_residcov_estimates():
    doc = _read(
        """
        <EM_TF><DataTypes>
          <DataType name="ZSTRIKE" type="real">
            <Tag>impedance_strike</Tag>
          </DataType>
        </DataTypes>
        <StatisticalEstimates>
          <Estimate name="VAR" type="real"><Tag>variance</Tag></Estimate>
          <Estimate name="INVSIGCOV" type="complex">
            <Tag>inverse_signal_covariance</Tag>
          </Estimate>
        </StatisticalEstimates>
        <Data count="1">
          <Period value="1">
            <ZSTRIKE><value>1.0</value></ZSTRIKE>
            <ZSTRIKE.VAR><value>0.1</value></ZSTRIKE.VAR>
            <ZSTRIKE.INVSIGCOV><value output="a" input="a">1 2</value></ZSTRIKE.INVSIGCOV>
          </Period>
        </Data></EM_TF>
        """
    )
    tf = doc.get_transfer_function("impedance_strike")
    assert tf.get_estimate("VAR") is not None
    assert tf.get_estimate("INVSIGCOV") is None
