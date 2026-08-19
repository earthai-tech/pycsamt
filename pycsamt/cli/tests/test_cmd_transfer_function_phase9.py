# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Phase 9 CLI tests for generic EDI/EMTF XML workflows."""

from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

from pycsamt.cli import main
from pycsamt.emtf import EMTF


_MINIMAL_XML = """\
<EM_TF>
  <Description>Phase 9 CLI fixture</Description>
  <ProductId>P9.S01.2026</ProductId>
  <SubType>MT_TF</SubType>
  <Tags>impedance</Tags>
  <Site>
    <Project>P9</Project><Survey>CLI</Survey>
    <YearCollected>2026</YearCollected>
    <Id>S01</Id><Name>Phase 9 Site</Name>
    <Location datum="WGS84">
      <Latitude>5.1</Latitude><Longitude>-3.2</Longitude>
      <Elevation units="meters">100</Elevation>
    </Location>
    <Orientation angle_to_geographic_north="0">orthogonal</Orientation>
  </Site>
  <SiteLayout>
    <InputChannels ref="site" units="m">
      <Magnetic name="Hx" orientation="0"/>
      <Magnetic name="Hy" orientation="90"/>
    </InputChannels>
    <OutputChannels ref="site" units="m">
      <Electric name="Ex" orientation="0"/>
      <Electric name="Ey" orientation="90"/>
    </OutputChannels>
  </SiteLayout>
  <Data count="2">
    <Period value="1" units="secs">
      <Z type="complex" size="2 2" units="[mV/km]/[nT]">
        <value output="Ex" input="Hx">1 2</value>
        <value output="Ex" input="Hy">3 4</value>
        <value output="Ey" input="Hx">5 6</value>
        <value output="Ey" input="Hy">7 8</value>
      </Z>
    </Period>
    <Period value="10" units="secs">
      <Z type="complex" size="2 2" units="[mV/km]/[nT]">
        <value output="Ex" input="Hx">2 1</value>
        <value output="Ex" input="Hy">4 3</value>
        <value output="Ey" input="Hx">6 5</value>
        <value output="Ey" input="Hy">8 7</value>
      </Z>
    </Period>
  </Data>
</EM_TF>
"""


def _write_xml(path: Path) -> Path:
    path.write_text(_MINIMAL_XML, encoding="utf-8")
    return path


def _write_edi(path: Path) -> Path:
    doc = EMTF.from_xml(_MINIMAL_XML)
    doc.write(path, format="edi", on_loss="ignore")
    return path


def test_convert_edi_to_xml_positional(runner: CliRunner, tmp_path: Path):
    src = _write_edi(tmp_path / "source.edi")
    dst = tmp_path / "converted.xml"

    result = runner.invoke(main, ["convert", str(src), str(dst)])

    assert result.exit_code == 0, result.output
    assert dst.exists()
    doc = EMTF.from_xml(dst)
    assert doc.station == "S01"
    assert doc.z.shape == (2, 2, 2)


def test_convert_xml_to_edi_positional(runner: CliRunner, tmp_path: Path):
    src = _write_xml(tmp_path / "source.xml")
    dst = tmp_path / "converted.edi"

    result = runner.invoke(
        main,
        ["convert", str(src), str(dst), "--on-loss", "ignore"],
    )

    assert result.exit_code == 0, result.output
    assert dst.exists()
    reread = EMTF.from_edi(dst, prefer_spectra=False)
    assert reread.station == "S01"
    assert reread.z.shape == (2, 2, 2)


def test_convert_explicit_from_to_and_derived_name(
    runner: CliRunner, tmp_path: Path
):
    src = _write_edi(tmp_path / "source.dat")
    out = tmp_path / "out"

    result = runner.invoke(
        main,
        [
            "convert",
            str(src),
            "--from",
            "edi",
            "--to",
            "emtf-xml",
            "--output-dir",
            str(out),
        ],
    )

    assert result.exit_code == 0, result.output
    assert (out / "source.xml").exists()


def test_convert_generic_dry_run_writes_nothing(
    runner: CliRunner, tmp_path: Path
):
    src = _write_xml(tmp_path / "source.xml")
    dst = tmp_path / "converted.edi"

    result = runner.invoke(
        main,
        ["convert", str(src), str(dst), "--dry-run"],
    )

    assert result.exit_code == 0
    assert "Dry run" in result.output
    assert not dst.exists()


def test_info_accepts_emtf_xml(runner: CliRunner, tmp_path: Path):
    src = _write_xml(tmp_path / "source.xml")

    result = runner.invoke(main, ["info", str(src), "--format", "json"])

    assert result.exit_code == 0, result.output
    data = json.loads(result.output)
    assert data[0]["format"] == "emtf_xml"
    assert data[0]["station"] == "S01"
    assert "impedance" in data[0]["components"]


def test_convert_help_documents_modern_syntax(runner: CliRunner):
    result = runner.invoke(main, ["convert", "--help"])

    assert result.exit_code == 0
    assert "SOURCE [TARGET]" in result.output
    assert "station.edi station.xml" in result.output
    assert "--from" in result.output
    assert "--on-loss" in result.output
