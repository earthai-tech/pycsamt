from __future__ import annotations

from datetime import datetime
from types import SimpleNamespace

import numpy as np
import pytest

from pycsamt.emtf import EMTF, StatisticalEstimate, TransferFunction
from pycsamt.emtf.converters.edi import (
    DataLossWarning,
    EMTFEDIConversionError,
    _array_or_none,
    _as_edi_date,
    _canonical_channel_name,
    _elevation_units_to_edi,
    _elevation_units_to_emtf,
    _explicit_keys,
    _extract_edi_info_lines,
    _first,
    _legacy_error,
    _loss,
    _measurement_angle,
    _preferred_variance,
    _raw_file_has_block,
    _rotation_vector,
    _tf_to_edi_error,
    _year_from_date,
    write_edi,
)


def _tf(*, estimate: StatisticalEstimate | None = None) -> TransferFunction:
    tf = TransferFunction(
        name="impedance",
        data=np.ones((1, 2, 2), dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Ex", "Ey"),
        periods=np.array([1.0]),
    )
    if estimate is not None:
        tf.add_estimate(estimate)
    return tf


def test_loss_policies_are_explicit():
    with pytest.warns(DataLossWarning, match="lost"):
        _loss(" WARN ", "lost")
    with pytest.raises(EMTFEDIConversionError, match="lost"):
        _loss("raise", "lost")
    _loss("ignore", "lost")
    with pytest.raises(ValueError, match="on_loss"):
        _loss("invalid", "lost")


def test_scalar_normalization_helpers_cover_empty_and_fallback_values():
    assert _first(None, "", "None", 0) == 0
    assert _first(None, "") is None
    assert _explicit_keys([" LAT = 1", "comment", "LONG=-3=ignored"]) == {
        "lat",
        "long",
    }

    assert _year_from_date(None) is None
    assert _year_from_date(datetime(2024, 1, 2)) == 2024
    assert _year_from_date("recorded/1998/archive") == 1998
    assert _year_from_date("not-a-date") is None
    assert _as_edi_date(datetime(2024, 1, 2)) == "01/02/24"
    assert _as_edi_date("2024-03-04") == "03/04/24"
    assert _as_edi_date("unknown") == "unknown"
    assert _as_edi_date("  ") is None


@pytest.mark.parametrize(
    ("value", "emtf", "edi"),
    [
        (None, "meters", "m"),
        ("metre", "meters", "m"),
        ("FEET", "feet", "ft"),
        ("furlong", "furlong", "furlong"),
    ],
)
def test_elevation_unit_mappings(value, emtf, edi):
    assert _elevation_units_to_emtf(value) == emtf
    assert _elevation_units_to_edi(value) == edi


def test_array_and_rotation_helpers():
    assert _array_or_none(None) is None
    assert _array_or_none([]) is None
    np.testing.assert_array_equal(_array_or_none([1]), [1])
    assert _rotation_vector(None, 2) is None
    np.testing.assert_allclose(_rotation_vector(12.0, 2), [12.0, 12.0])
    np.testing.assert_allclose(_rotation_vector([1.0, 2.0], 2), [1.0, 2.0])
    assert _rotation_vector([1.0, 2.0, 3.0], 2) is None


def test_raw_block_detection_handles_paths_and_read_failures(tmp_path):
    path = tmp_path / "sample.edi"
    path.write_text("  >ZROT // 1\n0\n", encoding="utf-8")
    assert _raw_file_has_block(SimpleNamespace(path=path), "zrot")
    assert not _raw_file_has_block(SimpleNamespace(path=path), "trot")
    assert not _raw_file_has_block(SimpleNamespace(path=None), "zrot")
    assert not _raw_file_has_block(
        SimpleNamespace(path=tmp_path / "missing.edi"), "zrot"
    )


def test_measurement_angles_and_channel_names():
    assert _measurement_angle(SimpleNamespace(azm="90")) == 90.0
    assert _measurement_angle(SimpleNamespace(azm="bad")) is None
    assert _measurement_angle(
        SimpleNamespace(azm=None, x=0, y=0, x2=0, y2=1)
    ) == 90.0
    assert _measurement_angle(
        SimpleNamespace(azm=None, x=0, y=0, x2=0, y2=0)
    ) is None
    assert _canonical_channel_name(None) is None
    assert _canonical_channel_name("rhx") == "RHx"
    assert _canonical_channel_name("aux") == "aux"


def test_info_lines_accept_strings_and_xml_text_dicts():
    doc = EMTF(field_notes={"edi_info": ["plain", {"#text": 3}, {"x": 1}, 4]})
    assert _extract_edi_info_lines(doc) == ["plain", "3"]
    doc.field_notes["edi_info"] = "single"
    assert _extract_edi_info_lines(doc) == ["single"]


def test_variance_is_preferred_and_converted_to_writer_error():
    variance = StatisticalEstimate(
        name="VAR", kind="variance", data=np.full((1, 2, 2), 9.0)
    )
    tf = _tf(estimate=variance)
    np.testing.assert_allclose(_preferred_variance(tf), 9.0)
    np.testing.assert_allclose(_tf_to_edi_error(tf), 3.0)
    assert _legacy_error(tf) is None

    bad = _tf(
        estimate=StatisticalEstimate(
            name="VAR", kind="variance", data=np.ones((1, 1, 1))
        )
    )
    with pytest.raises(EMTFEDIConversionError, match="VAR shape"):
        _preferred_variance(bad)


def test_legacy_error_requires_shape_and_explicit_semantics():
    generic = _tf(
        estimate=StatisticalEstimate(
            name="legacy",
            kind="standard_error",
            data=np.ones((1, 2, 2)),
        )
    )
    with pytest.raises(EMTFEDIConversionError, match="statistical convention"):
        _legacy_error(generic)

    explicit = _tf(
        estimate=StatisticalEstimate(
            name="legacy",
            kind="standard_error",
            data=np.full((1, 2, 2), 0.25),
            attrs={"semantics": "pycsamt_legacy_z_err"},
        )
    )
    np.testing.assert_allclose(_tf_to_edi_error(explicit), 0.25)

    explicit.get_estimate("standard_error").data = np.ones((1, 1, 1))
    with pytest.raises(EMTFEDIConversionError, match="shape"):
        _legacy_error(explicit)


def test_write_edi_rejects_invalid_target_and_object(tmp_path):
    with pytest.raises(TypeError, match="filesystem target"):
        write_edi(EMTF(), object())
    with pytest.raises(TypeError, match="accepts EMTF or EDIFile"):
        write_edi(object(), tmp_path / "out.edi")
