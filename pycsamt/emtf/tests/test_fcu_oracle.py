# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Numerical validation of pycsamt's EDI->EMTF-XML conversion against a
real, locally-compiled EMTF-FCU v4.1 (Fortran) reference implementation.

EMTF-FCU is vendored, gitignored, and never shipped -- see
``pycsamt/emtf/fcu-v4.1/`` and its ``.gitignore`` entry. It is not
guaranteed to exist (or be built) in any given environment, so this
whole module is gated behind ``requires_real_fcu`` and skips, rather
than fails, when no built ``edi2xml`` executable is found.

Build recipe (Linux/WSL2; needs gfortran):

    cd pycsamt/emtf/fcu-v4.1
    git clone --depth 1 https://github.com/andreww/fox.git FoX
    (cd FoX && ./configure && make)
    (cd f90 && sed -i 's#^BIN_DIR = .*#BIN_DIR = /tmp/fcu-bin#' Makefile \\
        && make -k && cp z2edi xml2z z2xml xml2edi edi2xml z2z edi2edi \\
           xml2xml *.pl /tmp/fcu-bin)
    mkdir -p bin && cp /tmp/fcu-bin/* bin/ \\
        && cp -r f90/DATATYPES f90/COPYRIGHT bin/

On Windows, the built (Linux ELF) executables are invoked through
``wsl.exe``; on Linux/macOS they are invoked directly. See
:func:`_run_fcu` and :func:`_fcu_available`.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from pycsamt.emtf.converters.edi import edi_to_emtf
from pycsamt.emtf.document import EMTF
from pycsamt.seg.edi import EDIFile

_FCU_BIN_DIR = Path(__file__).resolve().parents[1] / "fcu-v4.1" / "bin"

_CONFIG_XML = """\
<Configuration>
  <TimeSeriesArchived>0</TimeSeriesArchived>
  <Network>EM</Network>
  <Project>PYCSAMT</Project>
  <Survey>FCU-Oracle</Survey>
  <YearCollected>2026</YearCollected>
  <Country>USA</Country>
  <Tags>
  impedance,tipper
  </Tags>
  <ReleaseStatus>Unrestricted Release</ReleaseStatus>
  <ProcessedBy>pyCSAMT</ProcessedBy>
  <DateFormat>MM/DD/YY</DateFormat>
  <ChannelsOnTwoLines>0</ChannelsOnTwoLines>
  <ParseEDIInfo>0</ParseEDIInfo>
  <WriteEDIInfo>0</WriteEDIInfo>
  <MetadataOnly>0</MetadataOnly>
</Configuration>
"""

_WILLY_EDI = Path("data") / "AMT" / "WILLY_DATA" / "L18PLT" / "18-001A.edi"
_GV100_EDI = Path("data") / "gv_data" / "gv_final_edi" / "gv100.edi"


def _win_to_wsl_path(p: Path) -> str:
    """``C:\\foo\\bar`` -> ``/mnt/c/foo/bar`` (standard WSL2 drive mount)."""

    s = str(p.resolve())
    drive, rest = s.split(":", 1)
    return f"/mnt/{drive.lower()}{rest.replace(chr(92), '/')}"


def _fcu_available() -> bool:
    if not (_FCU_BIN_DIR / "edi2xml").exists():
        return False
    if sys.platform == "win32":
        return shutil.which("wsl.exe") is not None
    return True


requires_real_fcu = pytest.mark.skipif(
    not _fcu_available(),
    reason=(
        "no built EMTF-FCU v4.1 reference found at "
        f"{_FCU_BIN_DIR} (and/or wsl.exe unavailable on Windows); "
        "see this module's docstring for the local build recipe"
    ),
)


def _run_fcu(tool: str, *args: str, cwd: Path, timeout: int = 30) -> str:
    """Run a built FCU executable, returning combined stdout+stderr."""

    if sys.platform == "win32":
        wsl_cwd = _win_to_wsl_path(cwd)
        cmd = str(cwd / tool).replace("\\", "/")
        # `cwd` already holds a copy of the tool + DATATYPES/COPYRIGHT
        # (see `_stage_fcu_workdir`), so a plain relative "./tool" run
        # from inside that directory is what FCU expects for its
        # homedir-relative COPYRIGHT/DATATYPES lookup.
        shell_cmd = f"cd '{wsl_cwd}' && ./{tool} " + " ".join(args)
        result = subprocess.run(
            ["wsl.exe", "-e", "bash", "-lc", shell_cmd],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    else:
        result = subprocess.run(
            [f"./{tool}", *args],
            cwd=str(cwd),
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    return result.stdout + result.stderr


def _stage_fcu_workdir(tmp_path: Path, edi_source: Path) -> Path:
    """Copy the FCU binaries + DATATYPES/COPYRIGHT + a config/EDI into a
    throwaway working directory (FCU resolves DATATYPES/COPYRIGHT and
    config.xml relative to the current directory)."""

    work = tmp_path / "fcu_work"
    work.mkdir()
    for name in ("edi2xml", "DATATYPES", "COPYRIGHT"):
        src = _FCU_BIN_DIR / name
        dst = work / name
        if src.is_dir():
            shutil.copytree(src, dst)
        else:
            shutil.copy2(src, dst)
    (work / "config.xml").write_text(_CONFIG_XML, encoding="utf-8")
    edi_dst = work / edi_source.name
    edi_dst.write_text(edi_source.read_text(encoding="utf-8"), encoding="utf-8")
    return work


def _fcu_edi2xml(tmp_path: Path, edi_source: Path) -> EMTF:
    """Run the real FCU ``edi2xml`` and read its output back through
    pycsamt's own (separately tested) XML reader."""

    work = _stage_fcu_workdir(tmp_path, edi_source)
    xml_name = edi_source.stem + "_fcu.xml"
    log = _run_fcu(
        "edi2xml", edi_source.name, xml_name, "silent", "0.0", cwd=work
    )
    xml_path = work / xml_name
    assert xml_path.exists(), f"FCU edi2xml produced no output:\n{log}"
    return EMTF.from_xml(xml_path, strict=False)


def _sorted_by_period(tf: EMTF) -> tuple[np.ndarray, np.ndarray]:
    order = np.argsort(tf.periods)
    return order, tf.periods[order]


# ---------------------------------------------------------------------------
# Axis-aligned real field EDI: direct Z/variance/site comparison, no
# rotation ambiguity (channels are AZM=0/90 already).
# ---------------------------------------------------------------------------


@requires_real_fcu
def test_fcu_edi2xml_periods_match(tmp_path: Path) -> None:
    edi = EDIFile(_WILLY_EDI)
    tf_py = edi_to_emtf(edi)
    tf_fcu = _fcu_edi2xml(tmp_path, _WILLY_EDI)

    _, p_py = _sorted_by_period(tf_py)
    _, p_fcu = _sorted_by_period(tf_fcu)
    assert p_py.size == p_fcu.size
    np.testing.assert_allclose(p_py, p_fcu, rtol=1e-4)


@requires_real_fcu
def test_fcu_edi2xml_impedance_matches(tmp_path: Path) -> None:
    edi = EDIFile(_WILLY_EDI)
    tf_py = edi_to_emtf(edi)
    tf_fcu = _fcu_edi2xml(tmp_path, _WILLY_EDI)

    order_py, _ = _sorted_by_period(tf_py)
    order_fcu, _ = _sorted_by_period(tf_fcu)
    z_py = tf_py.z[order_py]
    z_fcu = tf_fcu.z[order_fcu]
    # FCU's own text serialization rounds to a handful of significant
    # digits, so this is a numerical-precision check, not bit-exact.
    np.testing.assert_allclose(z_py, z_fcu, rtol=1e-4, atol=1e-2)


@requires_real_fcu
def test_fcu_edi2xml_variance_matches(tmp_path: Path) -> None:
    edi = EDIFile(_WILLY_EDI)
    tf_py = edi_to_emtf(edi)
    tf_fcu = _fcu_edi2xml(tmp_path, _WILLY_EDI)

    order_py, _ = _sorted_by_period(tf_py)
    order_fcu, _ = _sorted_by_period(tf_fcu)
    v_py = tf_py.get_transfer_function("impedance").estimates["variance"].data
    v_fcu = (
        tf_fcu.get_transfer_function("impedance").estimates["variance"].data
    )
    np.testing.assert_allclose(
        v_py[order_py], v_fcu[order_fcu], rtol=1e-4, atol=1e-2
    )


@requires_real_fcu
def test_fcu_edi2xml_site_coordinates_match(tmp_path: Path) -> None:
    edi = EDIFile(_WILLY_EDI)
    tf_py = edi_to_emtf(edi)
    tf_fcu = _fcu_edi2xml(tmp_path, _WILLY_EDI)

    assert tf_py.lat == pytest.approx(tf_fcu.lat, abs=1e-4)
    assert tf_py.lon == pytest.approx(tf_fcu.lon, abs=1e-4)
    assert tf_py.elev == pytest.approx(tf_fcu.elev, abs=1e-2)


@requires_real_fcu
def test_fcu_edi2xml_channel_geometry_matches(tmp_path: Path) -> None:
    edi = EDIFile(_WILLY_EDI)
    tf_py = edi_to_emtf(edi)
    tf_fcu = _fcu_edi2xml(tmp_path, _WILLY_EDI)

    py_out = {c.name.upper(): c for c in tf_py.site_layout.output_channels}
    fcu_out = {c.name.upper(): c for c in tf_fcu.site_layout.output_channels}
    assert set(py_out) == set(fcu_out)
    for name, c_py in py_out.items():
        c_fcu = fcu_out[name]
        assert c_py.orientation == pytest.approx(c_fcu.orientation, abs=1e-3)
        assert c_py.x == pytest.approx(c_fcu.x, abs=1e-3)
        assert c_py.y == pytest.approx(c_fcu.y, abs=1e-3)


# ---------------------------------------------------------------------------
# Non-orthogonal real field EDI (AZM=-12.5/77.5) with tipper: validates
# pycsamt's own rotation math (EMTF.rotate) against FCU's rotation-to-
# geographic-north, matching the implementation plan's Phase 7 goal.
# ---------------------------------------------------------------------------


@requires_real_fcu
def test_fcu_edi2xml_rotation_matches(tmp_path: Path) -> None:
    edi = EDIFile(_GV100_EDI)
    tf_py = edi_to_emtf(edi)
    with pytest.warns(UserWarning):
        tf_py_rot = tf_py.rotate(
            0.0, target="orthogonal", use_legacy_edi_rotation=True
        )

    work = _stage_fcu_workdir(tmp_path, _GV100_EDI)
    log = _run_fcu(
        "edi2xml", _GV100_EDI.name, "gv100_fcu.xml", "silent", "0.0", cwd=work
    )
    xml_path = work / "gv100_fcu.xml"
    assert xml_path.exists(), f"FCU edi2xml produced no output:\n{log}"
    tf_fcu = EMTF.from_xml(xml_path, strict=False)

    order_py, _ = _sorted_by_period(tf_py_rot)
    order_fcu, _ = _sorted_by_period(tf_fcu)
    z_py = tf_py_rot.z[order_py]
    z_fcu = tf_fcu.z[order_fcu]
    np.testing.assert_allclose(z_py, z_fcu, rtol=1e-3, atol=1e-1)
