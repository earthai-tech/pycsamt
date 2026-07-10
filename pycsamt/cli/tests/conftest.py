# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Shared fixtures for pycsamt.cli tests.

Fixtures
--------
edi_dir         Path to edi_out/ with real .edi files (skips if absent).
single_edi      First .edi file in edi_dir.
site_edi_dir    Small real-EDI dir: prefers data/3edis/, falls back to
                data/AMT/TIPPER/, then edi_out/ (skips if all absent).
site_edi_dir_writable  Copy of site_edi_dir in tmp_path (for edit/export).
runner          click.testing.CliRunner instance.
isolated_home   Patch pycsamt.cli.survey path constants to a tmp dir.
occam_workdir   Tmp dir with minimal Occam2D signature files.
modem_workdir   Tmp dir with minimal ModEM signature files.
occam_workdir_with_iters  occam_workdir + iter/resp/log files.
fake_sites      Picklable _FakeSites object (no EDI I/O required).
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

# ---------------------------------------------------------------------------
# Project paths
# ---------------------------------------------------------------------------
# Path layout:  <repo>/pycsamt/cli/tests/conftest.py
#   parents[0] = <repo>/pycsamt/cli/tests/
#   parents[1] = <repo>/pycsamt/cli/
#   parents[2] = <repo>/pycsamt/
#   parents[3] = <repo>/                ← repo root

_PROJECT_ROOT  = Path(__file__).resolve().parents[3]   # pycsamt/ repo root
_EDI_OUT       = _PROJECT_ROOT / "edi_out"
_DATA_3EDIS    = _PROJECT_ROOT / "data" / "3edis"
_DATA_AMT_TIP  = _PROJECT_ROOT / "data" / "AMT" / "TIPPER"
_DATA_WILLY    = _PROJECT_ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _first_edi_dir(*candidates: Path) -> Path | None:
    """Return the first candidate that contains at least one .edi file."""
    for d in candidates:
        if d.exists() and list(d.glob("*.edi")):
            return d
    return None


# ---------------------------------------------------------------------------
# EDI directory fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def edi_dir() -> Path:
    """Directory containing real EDI files (edi_out/ only)."""
    if not _EDI_OUT.exists() or not list(_EDI_OUT.glob("*.edi")):
        pytest.skip("edi_out/ not found — skipping live EDI tests")
    return _EDI_OUT


@pytest.fixture(scope="session")
def single_edi(edi_dir: Path) -> Path:
    """First .edi file in edi_dir."""
    return sorted(edi_dir.glob("*.edi"))[0]


@pytest.fixture(scope="session")
def site_edi_dir() -> Path:
    """Small (2–5 station) real-EDI directory for site command tests.

    Resolution order:
      1. data/3edis/          — 3 EDIs, no tipper (good for basic tests)
      2. data/AMT/TIPPER/     — 2 EDIs with tipper
      3. edi_out/             — 3 mixed EDIs from repo root
    Skips when none are available.
    """
    found = _first_edi_dir(_DATA_3EDIS, _DATA_AMT_TIP, _EDI_OUT)
    if found is None:
        pytest.skip("No small EDI directory found — skipping site live tests")
    return found


@pytest.fixture
def site_edi_dir_writable(site_edi_dir: Path, tmp_path: Path) -> Path:
    """Copy of site_edi_dir in a writable tmp directory.

    Used by edit / export tests that write modified EDIs to disk.
    """
    dst = tmp_path / "edis"
    shutil.copytree(site_edi_dir, dst)
    return dst


# ---------------------------------------------------------------------------
# CLI runner
# ---------------------------------------------------------------------------

@pytest.fixture
def runner() -> CliRunner:
    """Click test runner."""
    return CliRunner()


# ---------------------------------------------------------------------------
# Isolated home directory (avoids polluting ~/.pycsamt during tests)
# ---------------------------------------------------------------------------

@pytest.fixture
def isolated_home(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Redirect pycsamt.cli.survey path constants to a tmp directory.

    ``_PYCSAMT_DIR``, ``_CONTEXT_FILE``, and ``_CACHE_ROOT`` are module-level
    constants computed once at import time, so we patch them directly rather
    than patching ``Path.home()``.
    """
    import pycsamt.cli.commands.survey as _cmd_survey
    import pycsamt.cli.survey as _survey

    fake_home   = tmp_path / "home"
    pycsamt_dir = fake_home / ".pycsamt"
    cache_root  = pycsamt_dir / "cache"
    fake_home.mkdir()
    pycsamt_dir.mkdir()

    monkeypatch.setattr(_survey, "_PYCSAMT_DIR",  pycsamt_dir)
    monkeypatch.setattr(_survey, "_CONTEXT_FILE", pycsamt_dir / "context.json")
    monkeypatch.setattr(_survey, "_CACHE_ROOT",   cache_root)
    monkeypatch.setattr(_cmd_survey, "_CACHE_ROOT",   cache_root)
    monkeypatch.setattr(_cmd_survey, "_PYCSAMT_DIR",  pycsamt_dir)

    return fake_home


# ---------------------------------------------------------------------------
# Fake Sites (picklable, no EDI I/O needed)
# ---------------------------------------------------------------------------

class _FakeSite:
    """Minimal picklable Site-like object for unit tests."""
    def __init__(self, name: str, lat: float = 0.0, lon: float = 0.0,
                 elev: float = 0.0, nfreq: int = 10) -> None:
        self.name  = name
        self._lat  = lat
        self._lon  = lon
        self._elev = elev
        self._nfreq= nfreq

    @property
    def coords(self):
        return (self._lat, self._lon, self._elev)


class _FakeSites:
    """Minimal picklable Sites-like container for unit tests."""
    def __init__(self, names: list[str]) -> None:
        self._items = [_FakeSite(n, lat=float(i), lon=float(i + 100))
                       for i, n in enumerate(names)]

    def __len__(self) -> int:
        return len(self._items)

    def __iter__(self):
        return iter(self._items)

    def __getitem__(self, key):
        if isinstance(key, int):
            return self._items[key]
        for s in self._items:
            if s.name == key:
                return s
        raise KeyError(key)


def make_fake_sites(n: int = 4, prefix: str = "S") -> _FakeSites:
    """Return a _FakeSites with *n* stations named S01, S02, …"""
    return _FakeSites([f"{prefix}{i:02d}" for i in range(1, n + 1)])


# ---------------------------------------------------------------------------
# CLI runner
# ---------------------------------------------------------------------------

@pytest.fixture
def runner() -> CliRunner:
    """Click test runner."""
    return CliRunner()


# ---------------------------------------------------------------------------
# Isolated home directory (avoids polluting ~/.pycsamt during tests)
# ---------------------------------------------------------------------------

@pytest.fixture
def isolated_home(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Redirect pycsamt.cli.survey path constants to a tmp directory.

    ``_PYCSAMT_DIR``, ``_CONTEXT_FILE``, and ``_CACHE_ROOT`` are module-level
    constants computed once at import time, so we patch them directly rather
    than patching ``Path.home()``.
    """
    import pycsamt.cli.commands.survey as _cmd_survey
    import pycsamt.cli.survey as _survey

    fake_home   = tmp_path / "home"
    pycsamt_dir = fake_home / ".pycsamt"
    cache_root  = pycsamt_dir / "cache"
    fake_home.mkdir()
    pycsamt_dir.mkdir()

    # patch the core module
    monkeypatch.setattr(_survey, "_PYCSAMT_DIR",  pycsamt_dir)
    monkeypatch.setattr(_survey, "_CONTEXT_FILE", pycsamt_dir / "context.json")
    monkeypatch.setattr(_survey, "_CACHE_ROOT",   cache_root)

    # patch the command module (it imports the same names)
    monkeypatch.setattr(_cmd_survey, "_CACHE_ROOT",   cache_root)
    monkeypatch.setattr(_cmd_survey, "_PYCSAMT_DIR",  pycsamt_dir)

    return fake_home


# ---------------------------------------------------------------------------
# Minimal inversion workdir fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def occam_workdir(tmp_path: Path) -> Path:
    """Create a workdir with Occam2D signature files."""
    wd = tmp_path / "occam_run"
    wd.mkdir()
    (wd / "Occam2DMesh").touch()
    (wd / "Occam2DModel").touch()
    (wd / "OccamData.dat").touch()
    (wd / "OccamStartup").write_text("Occam2D startup file\n")
    return wd


@pytest.fixture
def modem_workdir(tmp_path: Path) -> Path:
    """Create a workdir with ModEM signature files."""
    wd = tmp_path / "modem_run"
    wd.mkdir()
    (wd / "ModEMData.dat").touch()
    (wd / "ModEM.inv").touch()
    (wd / "ModEM.cov").touch()
    return wd


@pytest.fixture
def occam_workdir_with_iters(occam_workdir: Path) -> Path:
    """Occam2D workdir that also has iteration output files."""
    for i in range(1, 6):
        (occam_workdir / f"OCCAM2MT_{i:03d}.iter").touch()
        (occam_workdir / f"OCCAM2MT_{i:03d}.resp").touch()
    log = occam_workdir / "logFile.logFile"
    log.write_text(
        "Iteration 1: RMS = 3.421\n"
        "Iteration 2: RMS = 2.189\n"
        "Iteration 3: RMS = 1.542\n"
        "Iteration 4: RMS = 1.231\n"
        "Iteration 5: RMS = 1.087\n"
    )
    return occam_workdir
