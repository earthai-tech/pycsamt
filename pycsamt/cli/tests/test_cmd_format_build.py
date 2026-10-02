# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for the PCBH/PCGL/PCGS/PCPT builder and EDI/EMTF-XML sub-commands
of ``pycsamt format`` (:mod:`pycsamt.cli.commands.format`).

Test strategy
-------------
PCGL/PCGS/PCPT/PCBH-CSV fixtures are small, self-contained CSVs written
via ``tmp_path`` (they don't need any bundled sample data). The EDI ⇄
EMTF-XML tests use the real, bundled Broken Hill EDI survey under
``data/MT/broken-hill/edis`` and skip automatically when that directory
is absent.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from pycsamt.cli import main

_ROOT = Path(__file__).resolve().parents[3]
_EDI_DIR = _ROOT / "data" / "MT" / "broken-hill" / "edis"


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


def _edi_dir() -> Path:
    if not _EDI_DIR.exists() or not any(_EDI_DIR.glob("*.edi")):
        pytest.skip(f"No Broken Hill EDI sample data found at {_EDI_DIR}")
    return _EDI_DIR


# ---------------------------------------------------------------------------
# group help lists the new sub-commands
# ---------------------------------------------------------------------------


def test_group_help_lists_build_and_edi_commands(runner):
    result = runner.invoke(main, ["format", "--help"])
    assert result.exit_code == 0
    for sub in (
        "build-pcbh",
        "build-pcgl",
        "build-pcgs",
        "build-pcpt",
        "edi-to-xml",
        "xml-to-edi",
    ):
        assert sub in result.output


# ---------------------------------------------------------------------------
# build-pcgl
# ---------------------------------------------------------------------------


def test_build_pcgl_from_csv(runner, tmp_path):
    src = tmp_path / "units.csv"
    src.write_text(
        "name,rho_min,rho_max,color\n"
        "Topsoil,1,50,#8B5A2B\n"
        "Fresh granite,300,10000,#7D7D7D\n"
    )
    dst = tmp_path / "units.pcgl.json"
    result = runner.invoke(
        main,
        ["format", "build-pcgl", str(src), "-o", str(dst), "--title", "T"],
    )
    assert result.exit_code == 0, result.output
    assert dst.exists()
    payload = json.loads(dst.read_text())
    assert payload["title"] == "T"


def test_build_pcgl_json_output(runner, tmp_path):
    src = tmp_path / "units.csv"
    src.write_text("name,rho_min,rho_max\nA,1,10\n")
    dst = tmp_path / "units.pcgl.json"
    result = runner.invoke(
        main,
        ["format", "build-pcgl", str(src), "-o", str(dst), "-f", "json"],
    )
    assert result.exit_code == 0, result.output
    data = json.loads(result.output)
    assert data["n_entries"] == 1


def test_build_pcgl_existing_output_requires_overwrite(runner, tmp_path):
    src = tmp_path / "units.csv"
    src.write_text("name,rho_min,rho_max\nA,1,10\n")
    dst = tmp_path / "units.pcgl.json"
    dst.write_text("{}")
    result = runner.invoke(main, ["format", "build-pcgl", str(src), "-o", str(dst)])
    assert result.exit_code != 0
    assert "overwrite" in result.output.lower()


# ---------------------------------------------------------------------------
# build-pcgs
# ---------------------------------------------------------------------------


def test_build_pcgs_requires_one_source(runner, tmp_path):
    dst = tmp_path / "s.pcgs.json"
    result = runner.invoke(main, ["format", "build-pcgs", "-o", str(dst)])
    assert result.exit_code != 0


def test_build_pcgs_from_planar_csv(runner, tmp_path):
    planar = tmp_path / "planar.csv"
    planar.write_text(
        "x,kind,strike_deg,dip_deg,dip_direction_deg\n"
        "100,bedding,45,60,135\n"
    )
    dst = tmp_path / "s.pcgs.json"
    result = runner.invoke(
        main,
        ["format", "build-pcgs", "--planar", str(planar), "-o", str(dst)],
    )
    assert result.exit_code == 0, result.output
    assert dst.exists()


# ---------------------------------------------------------------------------
# build-pcpt
# ---------------------------------------------------------------------------


def test_build_pcpt_from_csv(runner, tmp_path):
    src = tmp_path / "targets.csv"
    src.write_text("name,x,y,z\nA,100,200,10\nB,150,250,12\n")
    dst = tmp_path / "targets.pcpt.json"
    result = runner.invoke(main, ["format", "build-pcpt", str(src), "-o", str(dst)])
    assert result.exit_code == 0, result.output
    payload = json.loads(dst.read_text())
    assert len(payload["points"]) == 2


# ---------------------------------------------------------------------------
# build-pcbh
# ---------------------------------------------------------------------------


def test_build_pcbh_from_combined_csv(runner, tmp_path):
    src = tmp_path / "boreholes.csv"
    src.write_text(
        "borehole_id,x,y,z,crs,kind,status,total_depth_md,from_md,to_md,"
        "lithology\n"
        "BH-1,100,200,10,EPSG:4326,water,completed,20,0,10,Sand\n"
    )
    dst = tmp_path / "boreholes.pcbh.json"
    result = runner.invoke(main, ["format", "build-pcbh", str(src), "-o", str(dst)])
    assert result.exit_code == 0, result.output
    payload = json.loads(dst.read_text())
    assert len(payload["boreholes"]) == 1


def test_build_pcbh_las_requires_collar_options(runner, tmp_path):
    fake_las = tmp_path / "hole.las"
    fake_las.write_text("~V\nVERS. 2.0 :\n")
    dst = tmp_path / "hole.pcbh.json"
    result = runner.invoke(main, ["format", "build-pcbh", str(fake_las), "-o", str(dst)])
    assert result.exit_code != 0
    assert "--collar-id" in result.output


# ---------------------------------------------------------------------------
# edi-to-xml / xml-to-edi
# ---------------------------------------------------------------------------


def test_edi_to_xml_single_file(runner, tmp_path):
    edi_dir = _edi_dir()
    first_edi = sorted(edi_dir.glob("*.edi"))[0]
    out_dir = tmp_path / "xml"
    result = runner.invoke(
        main, ["format", "edi-to-xml", str(first_edi), "-o", str(out_dir)]
    )
    assert result.exit_code == 0, result.output
    xml_files = list(out_dir.glob("*.xml"))
    assert len(xml_files) == 1


def test_edi_to_xml_directory_matches_file_count(runner, tmp_path):
    edi_dir = _edi_dir()
    n_edi = len(list(edi_dir.glob("*.edi")))
    out_dir = tmp_path / "xml"
    result = runner.invoke(
        main, ["format", "edi-to-xml", str(edi_dir), "-o", str(out_dir), "-f", "json"]
    )
    assert result.exit_code == 0, result.output
    written = json.loads(result.output)
    assert len(written) == n_edi
    assert len(list(out_dir.glob("*.xml"))) == n_edi


def test_edi_xml_edi_round_trip_preserves_station_names(runner, tmp_path):
    edi_dir = _edi_dir()
    xml_dir = tmp_path / "xml"
    back_dir = tmp_path / "edi_back"

    r1 = runner.invoke(
        main, ["format", "edi-to-xml", str(edi_dir), "-o", str(xml_dir)]
    )
    assert r1.exit_code == 0, r1.output

    r2 = runner.invoke(
        main, ["format", "xml-to-edi", str(xml_dir), "-o", str(back_dir)]
    )
    assert r2.exit_code == 0, r2.output

    original_names = {p.stem for p in edi_dir.glob("*.edi")}
    round_tripped_names = {p.stem for p in back_dir.glob("*.edi")}
    assert original_names == round_tripped_names


def test_xml_to_edi_nonexistent_source_fails(runner, tmp_path):
    result = runner.invoke(
        main, ["format", "xml-to-edi", str(tmp_path / "nope")]
    )
    assert result.exit_code != 0
