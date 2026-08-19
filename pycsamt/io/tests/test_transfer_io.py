# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

import io
from pathlib import Path

import numpy as np
import pytest

from pycsamt.emtf import EMTF
from pycsamt.io import (
    TransferFunctionFormatError,
    detect_tf_format,
    get_tf_format,
    list_tf_formats,
    read_transfer_function,
    register_tf_format,
)
from pycsamt.seg.edi import EDIFile


XML = """\
<?xml version="1.0"?>
<EM_TF>
  <ProductId>TEST.S01.2026</ProductId>
  <DataTypes>
    <DataType name="Z" type="complex" output="E" input="H">
      <Tag>impedance</Tag>
    </DataType>
  </DataTypes>
  <SiteLayout>
    <InputChannels>
      <Magnetic name="Hx"/><Magnetic name="Hy"/>
    </InputChannels>
    <OutputChannels>
      <Electric name="Ex"/><Electric name="Ey"/>
    </OutputChannels>
  </SiteLayout>
  <Data count="1">
    <Period value="1" units="secs">
      <Z>
        <value output="Ex" input="Hy">1 2</value>
        <value output="Ey" input="Hx">-3 -4</value>
      </Z>
    </Period>
  </Data>
</EM_TF>
"""


def _edi_text(station: str = "S01") -> str:
    return "\n".join(
        [
            ">HEAD",
            f"  DATAID={station}",
            "  LAT=26:00:00N",
            "  LONG=010:00:00E",
            "  ELEV=1000",
            "",
            ">INFO",
            "  PROJECT=SIM",
            "",
            ">=MTSECT",
            f"  SECTID={station}",
            "  NFREQ=2",
            "",
            ">FREQ  //2",
            "  1.000000E+02  1.000000E+01",
            "",
            ">ZXYR ROT=ZROT  //2",
            "  1.000000E+00  2.000000E+00",
            ">ZXYI ROT=ZROT  //2",
            "  3.000000E+00  4.000000E+00",
            "",
            ">END",
        ]
    )


def test_builtin_registry_and_aliases():
    formats = list_tf_formats()
    assert set(formats) >= {"edi", "emtf_xml"}
    assert get_tf_format("xml").name == "emtf_xml"
    assert get_tf_format("emtf").name == "emtf_xml"
    assert get_tf_format("seg-edi").name == "edi"


def test_detect_emtf_xml_from_inline_and_stream():
    assert detect_tf_format(XML) == "emtf_xml"
    assert detect_tf_format(io.StringIO(XML)) == "emtf_xml"


def test_detect_namespaced_emtf_xml():
    xml = XML.replace(
        "<EM_TF>", '<EM_TF xmlns="http://example.org/emtf">', 1
    )
    assert detect_tf_format(xml) == "emtf_xml"


def test_detect_edi_from_content_not_extension(tmp_path: Path):
    path = tmp_path / "site.dat"
    path.write_text(_edi_text(), encoding="utf-8")
    assert detect_tf_format(path) == "edi"


def test_bad_xml_extension_is_not_accepted_as_format(tmp_path: Path):
    path = tmp_path / "not_emtf.xml"
    path.write_text("<something_else/>", encoding="utf-8")
    with pytest.raises(TransferFunctionFormatError, match="content"):
        detect_tf_format(path)


def test_unknown_format_raises_clear_error():
    with pytest.raises(TransferFunctionFormatError):
        detect_tf_format("plain text with no transfer function")


def test_read_transfer_function_xml_returns_emtf():
    obj = read_transfer_function(XML)
    assert isinstance(obj, EMTF)
    assert obj.product_id == "TEST.S01.2026"
    assert obj.Z is not None
    assert obj.Z.z[0, 0, 1] == pytest.approx(1 + 2j)


def test_read_transfer_function_xml_forced_alias():
    obj = read_transfer_function(XML, format="xml", strict=True)
    assert isinstance(obj, EMTF)


def test_read_transfer_function_edi_preserves_established_backend(
    tmp_path: Path,
):
    path = tmp_path / "S01.edi"
    path.write_text(_edi_text(), encoding="utf-8")
    obj = read_transfer_function(path)
    assert isinstance(obj, EDIFile)
    assert obj.station == "S01"
    np.testing.assert_allclose(obj.Z.freq, [100.0, 10.0])
    np.testing.assert_allclose(obj.Z.z[:, 0, 1], [1 + 3j, 2 + 4j])


def test_explicit_format_can_read_xml_bytes():
    obj = read_transfer_function(XML.encode(), format="emtf_xml")
    assert isinstance(obj, EMTF)


def test_edi_in_memory_source_fails_with_boundary_message():
    with pytest.raises(TypeError, match="filesystem path"):
        read_transfer_function(_edi_text(), format="edi")


def test_custom_format_registration_is_additive():
    token = "phase4_dummy"
    register_tf_format(
        token,
        reader=lambda source, **kwargs: (source, kwargs),
        detector=lambda source: source == "PHASE4-DUMMY",
        aliases=("p4dummy",),
    )
    assert detect_tf_format("PHASE4-DUMMY") == token
    out = read_transfer_function("x", format="p4dummy", answer=42)
    assert out == ("x", {"answer": 42})


def test_registry_reports_phase6_write_capability():
    formats = list_tf_formats()
    assert formats["emtf_xml"]["writable"] is True
    assert formats["edi"]["writable"] is True


def test_write_transfer_function_xml_by_extension(tmp_path: Path):
    from pycsamt.io import write_transfer_function

    doc = read_transfer_function(XML)
    path = tmp_path / "written.xml"
    out = write_transfer_function(doc, path)
    assert out == path
    assert detect_tf_format(path) == "emtf_xml"
    reread = read_transfer_function(path)
    assert isinstance(reread, EMTF)
    np.testing.assert_allclose(doc.z, reread.z, equal_nan=True)


def test_write_transfer_function_stream_requires_explicit_format():
    from pycsamt.io import write_transfer_function

    doc = read_transfer_function(XML)
    stream = io.StringIO()
    with pytest.raises(TransferFunctionFormatError, match="explicit"):
        write_transfer_function(doc, stream)
    write_transfer_function(doc, stream, format="xml")
    assert "<EM_TF>" in stream.getvalue()


def test_write_transfer_function_edi_is_registered(tmp_path: Path):
    from pycsamt.io import write_transfer_function

    doc = read_transfer_function(XML)
    path = tmp_path / "site.edi"
    out = write_transfer_function(doc, path)
    assert out == path.resolve()
    assert path.exists()
    assert detect_tf_format(path) == "edi"
