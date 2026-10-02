from __future__ import annotations

import re

import pytest

from pycsamt.zonge.heads import Header
from pycsamt.zonge.property import (
    Hardware,
    Receiver,
    SkipFlag,
    SurveyAnnotation,
    SurveyConfiguration,
    Transmitter,
)


# -------------------------------------------------------------------- #
# SkipFlag                                                             #
# -------------------------------------------------------------------- #
def test_skipflag_codes_and_labels():
    s = SkipFlag()  # default "2"
    assert s.code == "2"
    assert s.get() == "good"

    s.set(1)
    assert s.code == "1"
    assert s.get() == "skip"

    s.set(0)
    assert s.code == "0"
    assert s.get() == "reject"

    s.set("*")
    assert s.code == "*"
    assert s.get() == "nodata"

    with pytest.raises(ValueError):
        s.set("3")  # invalid code


# -------------------------------------------------------------------- #
# Hardware                                                             #
# -------------------------------------------------------------------- #


def test_hardware_keywords_roundtrip():
    hw = Hardware()
    hw.update_from_keywords(
        {
            "version": "7.76",
            "source_file": "LCS01.fld",
            "dated": "99-01-01",
            "processed": "22 Jul 16",
            "astatic_ver": "v3.60d",
            "updated": "22/07/16",
            "tma_points": "5",
            "tma_freq": "1024",
        }
    )

    kw = hw.to_keywords()
    # types/values normalized
    assert kw["version"] == "7.76"
    assert str(kw["source_file"]).endswith("LCS01.fld")
    assert kw["dated"] == "99-01-01"
    assert kw["processed"] == "22 Jul 16"
    assert kw["astatic_ver"] == "v3.60d"
    assert kw["updated"] == "22/07/16"
    assert kw["tma_points"] == 5
    assert kw["tma_freq"] == 1024.0

    # JSON serialization should not explode
    js = hw.to_json(indent=2)
    assert isinstance(js, str)
    assert "LCS01.fld" in js


# -------------------------------------------------------------------- #
# SurveyAnnotation / SurveyConfiguration                               #
# -------------------------------------------------------------------- #


def test_annotation_update_and_export():
    ann = SurveyAnnotation()
    ann.update_from_keywords(
        {
            "$Job.Name": "North Silverbell",
            "$Job.Area": "Tucson, AZ",
            "$Job.For": "Zonge Engineering",
            "$Job.By": "Zonge",
            "$Job.Number": "9309",
            "$Job.Date": "Nov 93",
        }
    )
    kw = ann.to_keywords()
    assert kw["Job.Name"] == "North Silverbell"
    assert kw["Job.Area"] == "Tucson, AZ"
    assert kw["Job.For"] == "Zonge Engineering"
    assert kw["Job.By"] == "Zonge"
    assert kw["Job.Number"] == "9309"
    assert kw["Job.Date"] == "Nov 93"


def test_config_update_and_export():
    cfg = SurveyConfiguration()
    cfg.update_from_keywords(
        {
            "$Survey.Type": "CSAMT",
            "$Survey.Array": "Scalar",
            "$Line.Name": "LCS01",
            "$Line.Number": "0",
            "$Line.Azimuth": "0",
            "$Stn.GdpBeg": "0",
            "$Stn.GdpInc": "50",
            "$Stn.Beg": "0",
            "$Stn.Inc": "50",
            "$Stn.Left": "25",
            "$Stn.Right": "1375",
            "$Unit.Length": "m",
            "$Unit.E": "nV/Am",
            "$Unit.B": "pT/A",
            "$Unit.Phase": "mrad",
        }
    )
    kw = cfg.to_keywords()
    assert kw["Survey.Type"] == "CSAMT"
    assert kw["Survey.Array"] == "Scalar"
    assert kw["Line.Name"] == "LCS01"
    assert kw["Line.Number"] == 0.0
    assert kw["Line.Azimuth"] == 0.0
    assert kw["Stn.Left"] == 25.0
    assert kw["Stn.Right"] == 1375.0
    assert kw["Stn.GdpInc"] == 50.0
    assert kw["Unit.Length"] == "m"
    assert kw["Unit.E"] == "nV/Am"
    assert kw["Unit.B"] == "pT/A"
    assert kw["Unit.Phase"] == "mrad"


# Receiver / Transmitter
def test_rx_update_and_export_minimal():
    rx = Receiver()
    rx.update_from_keywords(
        {
            "$Rx.GdpStn": "25",
            "$Rx.Stn": "25",
            "$Rx.Length": "50 m",
            "$Rx.Cmp": "ExHy",
        }
    )
    assert rx.station == 25
    assert rx.gdp_station == 25
    assert rx.length_m == 50.0
    assert rx.unit == "m"
    assert rx.comps == "ExHy"

    kw = rx.to_keywords()
    assert kw["Rx.Stn"] == 25
    # length should be rendered with unit
    assert isinstance(kw["Rx.Length"], str)
    assert "50" in kw["Rx.Length"]
    assert "m" in kw["Rx.Length"]
    assert kw["Rx.Cmp"] == "ExHy"
    assert kw["Rx.GdpStn"] == 25


def test_tx_update_and_export_minimal():
    tx = Transmitter()
    tx.update_from_keywords(
        {
            "$Tx.Type": "Natural",
            "$Tx.GdpStn": "20",
            "$Tx.Stn": "20",
            "$Tx.Length": "5000",
        }
    )
    print(tx)
    assert tx.tx_type == "Natural"
    print(type(tx.gdp_station))
    assert tx.gdp_station == 20
    assert tx.station == 20
    assert tx.length_m == 5000.0

    kw = tx.to_keywords()
    assert kw["Tx.Type"] == "Natural"
    assert kw["Tx.GdpStn"] == 20
    assert kw["Tx.Stn"] == 20
    assert kw["Tx.Length"] == "5000 m"


# Header facade

LEGACY_HEADER_LINES = [
    r'\ AMTAVG 7.76: "LCS01.fld", Dated 99-01-01, Processed 22 Jul 16',
    r"\ ASTATIC v3.60d updated data on 22/07/16",
    r"\ 5-point TMA Filter at 1024 hertz",
    "",
    "$Survey.Type=CSAMT",
    "$Survey.Array=Scalar",
    "$Line.Name=LCS01",
    "$Line.Number=0",
    "$Line.Azimuth=0",
    "$Stn.GdpBeg=0",
    "$Stn.GdpInc=50",
    "$Stn.Beg=0",
    "$Stn.Inc=50",
    "$Stn.Left=25",
    "$Stn.Right=1375",
    "$Unit.Length=m",
    "$Unit.E=nV/Am",
    "$Unit.B=pT/A",
    "$Unit.Phase=mrad",
    "$Tx.GdpStn=20",
    "$Tx.Stn=20",
    "$Tx.Type=Natural",
    "$Rx.GdpStn=25",
    "$Rx.Stn=25",
    "$Rx.Length=50 m",
    "$Rx.Cmp=ExHy",
]


def test_header_from_lines_parses_banner_and_keywords():
    hdr = Header.from_lines(LEGACY_HEADER_LINES)

    # Hardware from banner
    assert hdr.hardware.version == "7.76"
    assert str(hdr.hardware.source_file).endswith("LCS01.fld")
    assert hdr.hardware.tma_points == 5
    assert hdr.hardware.tma_freq == 1024.0
    assert hdr.hardware.astatic_ver.startswith("v3.60")

    # Config & units
    assert hdr.config.survey_type == "CSAMT"
    assert hdr.config.array_type == "Scalar"
    assert hdr.config.unit_length == "m"
    assert hdr.config.unit_emag == "nV/Am"
    assert hdr.config.unit_hfield == "pT/A"
    assert hdr.config.unit_phase == "mrad"

    # Rx / Tx blocks
    assert hdr.rx.station == 25
    assert hdr.rx.length_m == 50.0
    assert hdr.rx.comps == "ExHy"
    assert hdr.tx.tx_type == "Natural"
    assert hdr.tx.station == 20

    # Writing the header should emit $Written and expected keys
    out = hdr.write()
    text = "\n".join(out)
    assert "$Survey.Type=CSAMT" in text
    assert "$Rx.Stn=25" in text
    assert "$Tx.Stn=20" in text
    assert re.search(r"^\$Written=\d{4}-\d{2}-\d{2}T", text, re.M)


def test_header_read_from_meta_mapping_directly():
    meta = {
        "$Job.Name": "North Silverbell",
        "$Survey.Type": "CSAMT",
        "$Rx.Stn": "75",
        "$Rx.Length": "200 ft",
        "$Rx.Cmp": "ExHy",
    }
    hdr = Header()
    hdr.read(meta=meta)

    assert hdr.annotation.project_name == "North Silverbell"
    assert hdr.config.survey_type == "CSAMT"
    assert hdr.rx.station == 75
    # unit conversion policy is class-specific; here we only
    # check that a number was parsed for length.
    assert pytest.approx(hdr.rx.length_m, rel=0, abs=1e-6) == 200.0


# -------------------------------------------------------------------- #
# smoke: __str__                                                       #
# -------------------------------------------------------------------- #
def test_strs_are_informative():
    assert "SkipFlag" in str(SkipFlag())
    assert "Hardware" in Hardware().__str__()
    assert "Survey(" in SurveyConfiguration().__str__()
    assert "Receiver" in Receiver().__str__()
    assert "Transmitter" in Transmitter().__str__()
    assert "Header(" in Header().__str__()


# -------------------------------------------------------------------- #
# Coverage: set()/get() direct APIs + unknown-key -> _extra fallback   #
# -------------------------------------------------------------------- #
def test_skipflag_set_none_is_noop():
    s = SkipFlag()
    s.set(None)
    assert s.code == "2"


def test_hardware_set_get_unknown_key_lands_in_extra():
    hw = Hardware()
    hw.set(version="9.0", custom_field="abc")
    assert hw.version == "9.0"
    assert hw.get("custom_field") == "abc"
    assert hw.get("missing", "fallback") == "fallback"


def test_hardware_to_keywords_defaults_omit_unset_fields():
    hw = Hardware()
    kw = hw.to_keywords()
    # only version/astatic_ver have non-None defaults
    assert kw["version"] == "7.76"
    assert kw["astatic_ver"] == "v3.60"
    assert "dated" not in kw
    assert "processed" not in kw
    assert "updated" not in kw
    assert "tma_points" not in kw
    assert "tma_freq" not in kw
    assert kw["source_file"] is None


def test_hardware_update_from_keywords_empty_dict_is_noop():
    hw = Hardware()
    hw.update_from_keywords({})
    assert hw.version == "7.76"
    assert hw.source_file is None
    assert hw.tma_points is None


def test_receiver_set_get_direct():
    rx = Receiver()
    rx.set(station=7, notes="edge site")
    assert rx.get("station") == 7
    assert rx.get("notes") == "edge site"
    assert rx.get("nope", "d") == "d"


def test_receiver_update_from_keywords_hpr_tuple_and_string_variants():
    rx1 = Receiver()
    rx1.update_from_keywords({"$Rx.HPR": (10.0, 2.0, 0.5)})
    assert rx1.hpr == (10.0, 2.0, 0.5)
    assert rx1.azimuth_deg == 10.0

    rx2 = Receiver()
    rx2.update_from_keywords({"$Rx.HPR": "15,3,1"})
    assert rx2.hpr == (15.0, 3.0, 1.0)

    rx3 = Receiver()
    rx3.update_from_keywords({"$Rx.HPR": "15;3;1"})
    assert rx3.hpr == (15.0, 3.0, 1.0)

    # too few parts -> the parsed-tuple branch is skipped, but the
    # generic keymap pass has already stored the raw string verbatim.
    rx4 = Receiver()
    rx4.update_from_keywords({"$Rx.HPR": "15,3"})
    assert rx4.hpr == "15,3"


def test_receiver_update_from_keywords_gps_and_missing_cmp():
    rx = Receiver()
    rx.update_from_keywords({"$GPS.Lat": "38.5", "$GPS.Lon": "-119.1"})
    assert rx.latitude == 38.5
    assert rx.longitude == -119.1
    # comps left at default since Rx.Cmp absent
    assert rx.comps == "ExHy"
    # station/gdp_station left None since absent
    assert rx.station is None
    assert rx.gdp_station is None


def test_receiver_to_keywords_empty_instance_omits_all():
    rx = Receiver(comps=None)
    kw = rx.to_keywords()
    assert kw == {}


def test_transmitter_set_get_direct():
    tx = Transmitter()
    tx.set(tx_type="Natural", notes="loop A")
    assert tx.get("tx_type") == "Natural"
    assert tx.get("notes") == "loop A"
    assert tx.get("nope", "d") == "d"


def test_transmitter_update_from_keywords_xmtr_alias_nonnumeric():
    tx = Transmitter()
    tx.update_from_keywords({"XMTR": "GDP-9"})
    # non-numeric XMTR value stored verbatim
    assert tx.gdp_station == "GDP-9"


def test_transmitter_update_from_keywords_xmtr_numeric_fallback():
    # Tx.GdpStn absent entirely (not merely non-numeric) so the XMTR
    # legacy-alias fallback path is exercised on its numeric-success leg.
    tx = Transmitter()
    tx.update_from_keywords({"XMTR": "42"})
    assert tx.gdp_station == 42


def test_transmitter_update_from_keywords_center_and_hpr_variants():
    tx1 = Transmitter()
    tx1.update_from_keywords({"$Tx.Center": (100.0, 200.0, 0.0)})
    assert tx1.center == (100.0, 200.0, 0.0)

    tx2 = Transmitter()
    tx2.update_from_keywords({"$Tx.Center": "100,200,0"})
    assert tx2.center == (100.0, 200.0, 0.0)

    tx3 = Transmitter()
    tx3.update_from_keywords({"$Tx.HPR": (0.0, 0.0, 0.0)})
    assert tx3.hpr == (0.0, 0.0, 0.0)

    tx4 = Transmitter()
    tx4.update_from_keywords({"$Tx.HPR": "1;2;3"})
    assert tx4.hpr == (1.0, 2.0, 3.0)

    # no Tx.Stn -> station stays None
    assert tx4.station is None


def test_transmitter_to_keywords_with_center_and_hpr():
    tx = Transmitter()
    tx.set(center=(1.0, 2.0, 3.0), hpr=(4.0, 5.0, 6.0))
    kw = tx.to_keywords()
    assert kw["Tx.Center"] == "1,2,3"
    assert kw["Tx.HPR"] == "4,5,6"


def test_transmitter_to_keywords_empty_instance_omits_all():
    tx = Transmitter()
    kw = tx.to_keywords()
    assert kw == {}


def test_survey_configuration_set_get_unknown_key():
    cfg = SurveyConfiguration()
    cfg.set(survey_type="TDEM", odd_field=1.5)
    assert cfg.survey_type == "TDEM"
    assert cfg.get("odd_field") == 1.5
    assert cfg.get("nope", "d") == "d"


def test_survey_configuration_to_json_roundtrips():
    cfg = SurveyConfiguration()
    cfg.set(extra_note="hello")
    js = cfg.to_json(indent=2)
    assert isinstance(js, str)
    assert "CSAMT" in js


def test_survey_annotation_set_get_unknown_key():
    ann = SurveyAnnotation()
    ann.set(project_name="X-Line", odd_field=3)
    assert ann.project_name == "X-Line"
    assert ann.get("odd_field") == 3
    assert ann.get("nope", "d") == "d"


def test_survey_annotation_str_and_to_json():
    ann = SurveyAnnotation(project_area="Area 51")
    s = str(ann)
    assert "Annotation(project=" in s
    assert "Area 51" in s

    js = ann.to_json(indent=2)
    assert isinstance(js, str)
    assert "pyCSAMT" in js


def test_survey_annotation_str_without_area():
    ann = SurveyAnnotation()
    assert "area=-" in str(ann)


# -------------------------------------------------------------------- #
# Coverage: module-private helper functions                            #
# -------------------------------------------------------------------- #
def test_to_number_sentinels_and_non_numeric():
    from pycsamt.zonge.property import _to_number

    assert _to_number(None) is None
    assert _to_number("") is None
    assert _to_number("*") is None
    assert _to_number("nan") is None
    assert _to_number("null") is None
    assert _to_number("not-a-number") == "not-a-number"
    assert _to_number("5") == 5
    assert _to_number("5.5") == 5.5


def test_parse_length_empty_string():
    from pycsamt.zonge.property import _parse_length

    assert _parse_length("") == (None, "m")
    assert _parse_length(None) == (None, "m")
    assert _parse_length("200") == (200.0, "m")
    assert _parse_length("200 ft") == (200.0, "ft")


def test_kv_roundtrip_skips_missing_and_none():
    from pycsamt.zonge.property import _kv_roundtrip

    keymap = {"a": "A.Key", "b": "B.Key", "c": "C.Key"}
    data = {"a": 1, "b": None}  # "c" intentionally absent from data
    out = _kv_roundtrip(data, keymap)
    assert out == {"A.Key": 1}


def test_apply_keywords_with_explicit_aliases():
    from pycsamt.zonge.property import _apply_keywords

    class _Obj:
        foo = None

    obj = _Obj()
    _apply_keywords(
        obj,
        keymap={"foo": "Foo.Key"},
        meta={"LegacyFoo": "99"},
        aliases={"LegacyFoo": "foo"},
    )
    assert obj.foo == "99"


if __name__ == "__main__":  # pragma: no-cover
    pytest.main([__file__])
