
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
    s = SkipFlag()                      # default "2"
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
        s.set("3")                      # invalid code

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
    assert kw["Tx.Length"] == '5000 m'



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

if __name__=='__main__': # pragma: no-cover
   pytest.main( [__file__])
