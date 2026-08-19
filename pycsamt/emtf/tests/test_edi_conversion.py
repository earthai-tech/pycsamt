# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Phase 6 regression tests for explicit EDI <-> EMTF interoperability."""

from __future__ import annotations

from pathlib import Path
import warnings

import numpy as np
import pytest

from pycsamt.emtf import (
    DataLossWarning,
    EMTF,
    EMTFEDIConversionError,
    StatisticalEstimate,
    TransferFunction,
)
from pycsamt.io import list_tf_formats, write_transfer_function
from pycsamt.metadata import LocationMeta, SiteMeta
from pycsamt.seg import EDIFile


MINIMAL_EDI = """\
>HEAD
 DATAID="S01"
 ACQBY="Demo"
 ACQDATE=01/02/26
 ENDDATE=01/03/26
 LAT=05:06:00.000
 LONG=-003:12:00.000
 ELEV=100
 UNITS=M
 STDVERS="SEG 1.0"
 COORDSYS="Geographic North"
 DATUM=WGS84
 EMPTY=1.0E+32

>INFO
 PROJECT=DemoProject
 SURVEY=DemoSurvey
 PROCESSEDBY="Tester"
 PROCESSINGSOFTWARE="DemoProc"

>=DEFINEMEAS
 MAXCHAN=7
 MAXRUN=1
 MAXMEAS=7
 UNITS=M
 REFTYPE=CART
 REFLAT=05:06:00.000
 REFLONG=-003:12:00.000
 REFELEV=100
>HMEAS ID=1 CHTYPE=HX X=0 Y=0 Z=0 AZM=0 ACQCHAN=H1
>HMEAS ID=2 CHTYPE=HY X=0 Y=0 Z=0 AZM=90 ACQCHAN=H2
>HMEAS ID=3 CHTYPE=HZ X=0 Y=0 Z=0 AZM=0 ACQCHAN=H3
>EMEAS ID=4 CHTYPE=EX X=-50 Y=0 X2=50 Y2=0 ACQCHAN=E1
>EMEAS ID=5 CHTYPE=EY X=0 Y=-50 X2=0 Y2=50 ACQCHAN=E2
>HMEAS ID=6 CHTYPE=HX X=0 Y=1000 AZM=0 ACQCHAN=RH1
>HMEAS ID=7 CHTYPE=HY X=0 Y=1000 AZM=90 ACQCHAN=RH2

>=MTSECT
 SECTID="S01"
 NFREQ=2
>! dynamic defaults
 HX=1
 HY=2
 HZ=3
 EX=4
 EY=5
 RX=6
 RY=7

>FREQ //2
 100 10
>ZROT //2
 0 0
>ZXXR ROT=ZROT //2
 1 2
>ZXXI ROT=ZROT //2
 0.1 0.2
>ZXX.VAR ROT=ZROT //2
 4 9
>ZXYR ROT=ZROT //2
 3 4
>ZXYI ROT=ZROT //2
 0.3 0.4
>ZXY.VAR ROT=ZROT //2
 16 25
>ZYXR ROT=ZROT //2
 5 6
>ZYXI ROT=ZROT //2
 0.5 0.6
>ZYX.VAR ROT=ZROT //2
 36 49
>ZYYR ROT=ZROT //2
 7 8
>ZYYI ROT=ZROT //2
 0.7 0.8
>ZYY.VAR ROT=ZROT //2
 64 81
>TROT //2
 0 0
>TXR.EXP ROT=TROT //2
 0.1 0.2
>TXI.EXP ROT=TROT //2
 0.01 0.02
>TXVAR.EXP ROT=TROT //2
 0.0004 0.0009
>TYR.EXP ROT=TROT //2
 -0.1 -0.2
>TYI.EXP ROT=TROT //2
 -0.01 -0.02
>TYVAR.EXP ROT=TROT //2
 0.0016 0.0025
>END
"""


def _edi_path(tmp_path: Path) -> Path:
    path = tmp_path / "minimal.edi"
    path.write_text(MINIMAL_EDI, encoding="utf-8")
    return path


def _simple_document(*, include_tipper: bool = True) -> EMTF:
    periods = np.array([0.01, 0.1])
    z = np.array(
        [
            [[1 + 0.1j, 3 + 0.3j], [5 + 0.5j, 7 + 0.7j]],
            [[2 + 0.2j, 4 + 0.4j], [6 + 0.6j, 8 + 0.8j]],
        ]
    )
    ztf = TransferFunction(
        name="impedance",
        data=z,
        input_channels=("Hx", "Hy"),
        output_channels=("Ex", "Ey"),
        periods=periods,
    )
    ztf.add_estimate(
        StatisticalEstimate(
            name="VAR",
            kind="variance",
            data=np.array(
                [
                    [[4.0, 16.0], [36.0, 64.0]],
                    [[9.0, 25.0], [49.0, 81.0]],
                ]
            ),
        )
    )
    tfs = {"impedance": ztf}
    if include_tipper:
        tip = np.array(
            [
                [[0.1 + 0.01j, -0.1 - 0.01j]],
                [[0.2 + 0.02j, -0.2 - 0.02j]],
            ]
        )
        ttf = TransferFunction(
            name="tipper",
            data=tip,
            input_channels=("Hx", "Hy"),
            output_channels=("Hz",),
            periods=periods,
        )
        ttf.add_estimate(
            StatisticalEstimate(
                name="VAR",
                kind="variance",
                data=np.array(
                    [
                        [[0.0004, 0.0016]],
                        [[0.0009, 0.0025]],
                    ]
                ),
            )
        )
        tfs["tipper"] = ttf
    return EMTF(
        periods=periods,
        transfer_functions=tfs,
        site=SiteMeta(
            site_id="S01",
            project="DemoProject",
            survey="DemoSurvey",
            location=LocationMeta(
                latitude=5.1,
                longitude=-3.2,
                elevation=100.0,
            ),
        ),
    )


def test_edi_to_emtf_preserves_tf_variance_and_remote_layout(tmp_path: Path):
    edi = EDIFile(_edi_path(tmp_path))
    doc = EMTF.from_edi(edi)

    assert doc.site.site_id == "S01"
    assert doc.site.project == "DemoProject"
    assert doc.site.survey == "DemoSurvey"
    assert doc.site.location.latitude == pytest.approx(5.1)
    assert doc.site.location.longitude == pytest.approx(-3.2)
    assert doc.site_layout.input_names == ("Hx", "Hy")
    assert doc.site_layout.output_names == ("Hz", "Ex", "Ey")

    np.testing.assert_allclose(doc.frequency, edi.Z.freq)
    np.testing.assert_allclose(doc.z, edi.Z.z)
    var = doc.impedance.get_estimate("VAR")
    np.testing.assert_allclose(var.data, np.square(edi.Z.z_err))

    tvar = doc.tipper_tf.get_estimate("VAR")
    np.testing.assert_allclose(tvar.data, np.square(edi.Tip.tipper_err))

    remote = doc.processing.remote_reference
    assert remote.extra["edi_rx"] == "6"
    assert remote.extra["edi_ry"] == "7"
    assert len(doc.metadata["edi_measurements"]) == 7


def test_edi_to_emtf_does_not_invent_measurement_dates(tmp_path: Path):
    doc = EMTF.from_edi(_edi_path(tmp_path))
    records = doc.metadata["edi_measurements"]
    assert all(record["measdate"] is None for record in records)


def test_emtf_to_edi_restores_legacy_error_carrier_and_remote(tmp_path: Path):
    original = EDIFile(_edi_path(tmp_path))
    doc = EMTF.from_edi(original)
    converted = doc.to_edi()

    np.testing.assert_allclose(converted.Z.z, original.Z.z)
    np.testing.assert_allclose(converted.Z.z_err, original.Z.z_err)
    np.testing.assert_allclose(converted.Tip.tipper, original.Tip.tipper)
    np.testing.assert_allclose(
        converted.Tip.tipper_err,
        original.Tip.tipper_err,
    )
    mtsect = converted.get_section("mtsect")
    assert mtsect.rx == "6"
    assert mtsect.ry == "7"
    measurements = converted.get_section("definemeas").all_meas()
    assert {str(item.id) for item in measurements} == {
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
    }


def test_edi_file_roundtrip_preserves_numeric_tf(tmp_path: Path):
    original = EDIFile(_edi_path(tmp_path))
    doc = EMTF.from_edi(original)
    target = tmp_path / "roundtrip.edi"
    doc.write(target, format="edi")
    reread = EDIFile(target)

    np.testing.assert_allclose(reread.Z.freq, original.Z.freq)
    np.testing.assert_allclose(reread.Z.z, original.Z.z, rtol=2e-6)
    np.testing.assert_allclose(
        reread.Z.z_err,
        original.Z.z_err,
        rtol=2e-6,
    )
    np.testing.assert_allclose(
        reread.Tip.tipper,
        original.Tip.tipper,
        rtol=2e-6,
    )
    np.testing.assert_allclose(
        reread.Tip.tipper_err,
        original.Tip.tipper_err,
        rtol=2e-6,
    )


def test_edi_to_xml_roundtrip_preserves_standard_scientific_subset(
    tmp_path: Path,
):
    doc = EMTF.from_edi(_edi_path(tmp_path))
    reread = EMTF.from_xml(doc.to_xml())

    np.testing.assert_allclose(reread.frequency, doc.frequency)
    np.testing.assert_allclose(reread.z, doc.z)
    np.testing.assert_allclose(reread.tipper, doc.tipper)
    np.testing.assert_allclose(
        reread.impedance.get_estimate("VAR").data,
        doc.impedance.get_estimate("VAR").data,
    )
    assert reread.site.site_id == "S01"
    assert (
        reread.processing.remote_reference.reference_type
        == "Remote Reference"
    )


def test_full_covariance_and_custom_tf_emit_data_loss_warning():
    doc = _simple_document()
    invsig = StatisticalEstimate(
        name="INVSIGCOV",
        kind="inverse_signal_covariance",
        data=np.ones((2, 2, 2), dtype=complex),
    )
    doc.impedance.add_estimate(invsig)
    custom = TransferFunction(
        name="custom_response",
        data=np.ones((2, 1, 1), dtype=complex),
        input_channels=("H",),
        output_channels=("E",),
        periods=doc.periods,
    )
    doc.add_transfer_function(custom)

    with pytest.warns(DataLossWarning):
        edi = doc.to_edi(on_loss="warn")
    assert edi.Z.n_freq == 2

    with pytest.raises(EMTFEDIConversionError):
        doc.to_edi(on_loss="raise")


def test_generic_standard_error_is_not_silently_interpreted_as_edi_var():
    periods = np.array([1.0])
    tf = TransferFunction(
        name="impedance",
        data=np.ones((1, 2, 2), dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Ex", "Ey"),
        periods=periods,
    )
    tf.add_estimate(
        StatisticalEstimate(
            name="STD",
            kind="standard_error",
            data=np.ones((1, 2, 2)),
        )
    )
    doc = EMTF(
        periods=periods,
        transfer_functions={"impedance": tf},
        site=SiteMeta(site_id="S01"),
    )
    with pytest.raises(EMTFEDIConversionError, match="statistical convention"):
        doc.to_edi(on_loss="ignore")


def test_neutral_edi_writer_preserves_zero_and_marks_nan_missing(
    tmp_path: Path,
):
    doc = _simple_document(include_tipper=False)
    data = np.array(doc.impedance.data, copy=True)
    data[0, 0, 0] = 0.0 + 0.0j
    data[0, 1, 1] = np.nan + 1j * np.nan
    doc.impedance.data = data

    target = tmp_path / "zero_missing.edi"
    write_transfer_function(doc, target, format="edi")
    text = target.read_text(encoding="utf-8")

    # Physical zero must survive the neutral converter.
    zxx_block = text.split(">ZXXR", 1)[1].split(">ZXXI", 1)[0]
    assert "0.000000E+00" in zxx_block
    # Missing data remain the EDI EMPTY sentinel, not a physical zero.
    zyy_block = text.split(">ZYYR", 1)[1].split(">ZYYI", 1)[0]
    assert "1.000000E+32" in zyy_block


def test_all_zero_tipper_is_not_dropped_by_neutral_writer(tmp_path: Path):
    doc = _simple_document(include_tipper=True)
    doc.tipper_tf.data[:] = 0.0
    target = tmp_path / "zero_tipper.edi"
    write_transfer_function(doc, target, format="edi")
    text = target.read_text(encoding="utf-8")
    assert ">TXR.EXP" in text
    assert ">TYR.EXP" in text


def test_phase6_registry_advertises_edi_write_support():
    formats = list_tf_formats()
    assert formats["edi"]["writable"] is True
    assert formats["emtf_xml"]["writable"] is True


def test_direct_emtf_to_edi_keeps_explicit_all_zero_tipper(tmp_path: Path):
    doc = _simple_document(include_tipper=True)
    doc.tipper_tf.data[:] = 0.0

    edi = doc.to_edi(on_loss="ignore")
    assert getattr(edi, "_force_tipper_on_write", False) is True

    out_dir = tmp_path / "direct"
    written = edi.write(
        new_edifn="zero_tipper_direct.edi",
        savepath=out_dir,
        preserve_zero=True,
        stamp_headers=False,
    )
    text = Path(written).read_text(encoding="utf-8")
    assert ">TXR.EXP" in text
    assert ">TYR.EXP" in text
