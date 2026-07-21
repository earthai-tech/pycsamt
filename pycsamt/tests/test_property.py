from __future__ import annotations

import pytest

from pycsamt.exceptions import FileHandlingError
from pycsamt.property import FieldAliases, FileRecognizer, TermDefinitions


def test_term_definitions_attributes_are_strings():
    assert isinstance(TermDefinitions.reference_frequency, str)
    assert "frequency" in TermDefinitions.reference_frequency
    assert isinstance(TermDefinitions.j_format, str)
    assert "Jones" in TermDefinitions.j_format


def test_field_aliases_lists():
    assert "lon" in FieldAliases.longitude
    assert "lat" in FieldAliases.latitude
    assert "e" in FieldAliases.easting
    assert "n" in FieldAliases.northing
    assert "sta" in FieldAliases.station
    assert "elev" in FieldAliases.elevation
    assert "azim" in FieldAliases.azimuth
    assert None in FieldAliases.missing_values


def test_recognize_invalid_path_raises():
    with pytest.raises(FileHandlingError):
        FileRecognizer.recognize("does/not/exist.edi")

    with pytest.raises(FileHandlingError):
        FileRecognizer.recognize("")


def test_recognize_shallow_uses_extension_only(tmp_path):
    p = tmp_path / "site.edi"
    p.write_text("garbage content that matches nothing", encoding="utf-8")
    assert FileRecognizer.recognize(str(p), deep=False) == "edi"


def test_recognize_shallow_unknown_extension_falls_back_to_deep(tmp_path):
    p = tmp_path / "site.unknownext"
    p.write_text(">HEAD\nsome\n>END\n", encoding="utf-8")
    # extension not in signatures -> falls through to deep content scan
    assert FileRecognizer.recognize(str(p), deep=False) == "edi"


def test_recognize_deep_edi(tmp_path):
    p = tmp_path / "site.dat"
    p.write_text(">HEAD\nfoo\n>END\n", encoding="utf-8")
    assert FileRecognizer.recognize(str(p)) == "edi"


def test_recognize_deep_j(tmp_path):
    p = tmp_path / "site.dat"
    p.write_text(">AZIMUTH 0.0\nRXX\nRXY\n", encoding="utf-8")
    assert FileRecognizer.recognize(str(p)) == "j"


def test_recognize_deep_avg(tmp_path):
    p = tmp_path / "site.dat"
    p.write_text("AMTAVG\nskp\n%Rho\n", encoding="utf-8")
    assert FileRecognizer.recognize(str(p)) == "avg"


def test_recognize_deep_stn(tmp_path):
    p = tmp_path / "site.dat"
    p.write_text("Station Freq\n1 2\n", encoding="utf-8")
    assert FileRecognizer.recognize(str(p)) == "stn"


def test_recognize_deep_mesh(tmp_path):
    p = tmp_path / "site.dat"
    p.write_text("ZZZZZZZZZZZZ\n????????????\n", encoding="utf-8")
    assert FileRecognizer.recognize(str(p)) == "mesh"


def test_recognize_deep_model(tmp_path):
    p = tmp_path / "site.dat"
    p.write_text("MODEL NAME\nNUM LAYERS\n", encoding="utf-8")
    assert FileRecognizer.recognize(str(p)) == "model"


def test_recognize_deep_startup(tmp_path):
    p = tmp_path / "site.dat"
    p.write_text("STARTUP\nITERATION\n", encoding="utf-8")
    assert FileRecognizer.recognize(str(p)) == "startup"


def test_recognize_deep_iter(tmp_path):
    p = tmp_path / "site.dat"
    p.write_text("Iteration 1\nMisfit 0.5\n", encoding="utf-8")
    assert FileRecognizer.recognize(str(p)) == "iter"


def test_recognize_deep_logfile(tmp_path):
    p = tmp_path / "site.dat"
    p.write_text("MISFIT\nROUGHNESS\n", encoding="utf-8")
    assert FileRecognizer.recognize(str(p)) == "logfile"


def test_recognize_deep_numeric_resp(tmp_path):
    p = tmp_path / "site.dat"
    p.write_text("1.0 2.0 3.0\n4.0 5.0 6.0\n", encoding="utf-8")
    assert FileRecognizer.recognize(str(p)) == "resp"


def test_recognize_deep_unrecognized(tmp_path):
    p = tmp_path / "site.dat"
    p.write_text("just some random text\nnot numeric either\n", encoding="utf-8")
    with pytest.raises(FileHandlingError):
        FileRecognizer.recognize(str(p))


def test_recognize_deep_unreadable_binary(tmp_path):
    p = tmp_path / "site.dat"
    p.write_bytes(b"\xff\xfe\x00\x01\x80\x81")
    with pytest.raises(FileHandlingError):
        FileRecognizer.recognize(str(p))


def test_is_numeric_file_empty_line_is_false(tmp_path):
    p = tmp_path / "site.dat"
    # first line empty (no tokens) -> not numeric -> falls through to
    # unrecognized since no other signature matches either.
    p.write_text("\nnot numeric\n", encoding="utf-8")
    with pytest.raises(FileHandlingError):
        FileRecognizer.recognize(str(p))
