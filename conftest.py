# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
Pytest configuration and shared fixtures.

Locates project data folders and exposes paths to legacy AVG,
STN, and EDI files. Also provides a simulated minimal EDI file
for tests that do not rely on repository data.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest


def _terminate_process(code: int) -> None:
    """Terminate without running CPython extension finalizers."""

    sys.stdout.flush()
    sys.stderr.flush()
    if sys.platform == "win32":
        import ctypes

        kernel32 = ctypes.windll.kernel32
        kernel32.TerminateProcess(kernel32.GetCurrentProcess(), code)
    else:
        os._exit(code)


def _is_qt_interface_run(config) -> bool:
    """Return whether this pytest invocation owns the desktop Qt tests."""

    args = {str(arg).replace("\\", "/").rstrip("/") for arg in config.args}
    return "PySide6" in sys.modules and any(
        arg == "pycsamt/app/tests" or "/pycsamt/app/tests" in arg
        for arg in args
    )


@pytest.hookimpl(trylast=True)
def pytest_terminal_summary(terminalreporter, exitstatus, config):  # noqa: ARG001
    """Exit after reports are persisted but before Qt plugin finalization."""

    if _is_qt_interface_run(config):
        _terminate_process(int(exitstatus))


@pytest.hookimpl(trylast=True)
def pytest_sessionfinish(session, exitstatus):
    """Bypass unsafe Qt/Shiboken finalization after interface test runs."""

    if not _is_qt_interface_run(session.config):
        return

    # Root conftest is loaded during pytest's initial configuration, so this
    # session hook cannot be unregistered with a nested test directory.
    # ``trylast`` lets coverage and terminal reporters persist results first.
    _terminate_process(int(exitstatus))


# Root conftest.py is imported during pytest's initial-conftest phase,
# before pytest-cov's coverage tracer starts (that happens in
# pytest_configure). Some test modules (e.g. test_web_callbacks_inversion.py)
# import torch at module scope for the first time mid-session; initializing
# torch's C extension while coverage.py is already tracing corrupts memory
# and causes unrelated-looking segfaults later on (see
# pycsamt/app/tests/_cov_runner_scratch.py for the original diagnosis).
# Pre-importing here, before tracing begins, avoids it.
try:
    import torch  # noqa: F401
except ImportError:
    pass


@pytest.fixture(autouse=True)
def _close_matplotlib_figures():
    """
    Close every pyplot figure after each test to bound figure accumulation.

    Nothing in the suite closes figures globally: the interfaces shard alone
    makes ~160 ``plt.subplots``/``plt.figure`` calls, most of which never get
    a matching ``plt.close``. Left open, each figure keeps its Axes, transform
    graphs, and Agg render buffers alive in pyplot's ``Gcf`` registry, so the
    live graph grows monotonically across the whole session.

    On Python 3.9 -- the only interpreter that resolves the *terminal*
    matplotlib 3.9.x line (3.10+ pull matplotlib 3.10+, which is unaffected) --
    that unbounded pile eventually tips matplotlib's C extensions into a
    segmentation fault deep in ``matplotlib.transforms`` while a brand-new
    Axes is being created. The tell is the crash point: it lands ~15% into the
    session, not at the first ``plt.subplots`` call, which is the signature of
    resource accumulation rather than an intrinsic per-call bug (CI:
    "Python 3.9 / interfaces" shard, test_batch_export_tool.py).

    Bounding the live figure set to what a single test creates removes the
    trigger, and it is good hygiene on every interpreter. Runs in teardown
    (after the test body), and only touches pyplot if the test actually
    imported it -- pure parsing/CLI tests never pull matplotlib in.
    """
    yield
    plt = sys.modules.get("matplotlib.pyplot")
    if plt is not None:
        try:
            plt.close("all")
        except Exception:  # pragma: no cover - defensive teardown
            pass


def get_project_root() -> Path:
    """Find repo root by locating the 'pycsamt' dir."""
    cur = Path(__file__).resolve()
    for parent in cur.parents:
        if (parent / "pycsamt").is_dir():
            return parent
    raise FileNotFoundError("Project root not found.")


@pytest.fixture(scope="session")
def project_root() -> Path:
    """Project root directory."""
    return get_project_root()


# --------------------- AVG / legacy data ----------------------


@pytest.fixture(scope="session")
def data_path(project_root: Path) -> Path:
    """Base path to bundled AVG data."""
    return project_root / "data" / "avg"


@pytest.fixture(scope="session")
def legacy_data_file(data_path: Path) -> Path:
    """Path to legacy K1.AVG; skip if missing."""
    p = data_path / "K1.AVG"
    if not p.exists():
        pytest.skip(f"Missing legacy data: {p}")
    return p


@pytest.fixture(scope="session")
def modern_data_file(data_path: Path) -> Path:
    """Path to modern K2.AVG; skip if missing."""
    p = data_path / "K2.AVG"
    if not p.exists():
        pytest.skip(f"Missing modern data: {p}")
    return p


@pytest.fixture(scope="session")
def stn_file_k1(data_path: Path) -> Path:
    """Path to K1.stn topography; skip if missing."""
    p = data_path / "K1.stn"
    if not p.exists():
        pytest.skip(f"Missing STN file: {p}")
    return p


@pytest.fixture(scope="session")
def stn_file_k2(data_path: Path) -> Path:
    """Path to K2.stn topography; skip if missing."""
    p = data_path / "K2.stn"
    if not p.exists():
        pytest.skip(f"Missing STN file: {p}")
    return p


# --------------------------- EDI data -------------------------


@pytest.fixture(scope="session")
def edi_path(project_root: Path) -> Path:
    """Base path to bundled EDI data."""
    return project_root / "data"


@pytest.fixture(scope="session")
def edi_imp_file(edi_path: Path) -> Path:
    """Path to a bundled MT impedance EDI; skip if missing."""
    p = edi_path / "MT" / "kap03lmt_edis" / "kap103.edi"
    if not p.exists():
        pytest.skip(f"Missing EDI: {p}")
    return p


@pytest.fixture(scope="session")
def edi_spe_file(edi_path: Path) -> Path:
    """Path to a bundled spectra EDI; skip if missing."""
    p = edi_path / "MT" / "SPECTRA" / "spectra01.edi"
    if not p.exists():
        pytest.skip(f"Missing EDI: {p}")
    return p


@pytest.fixture(scope="session")
def edi_csamt_file(edi_path: Path) -> Path:
    """Path to a bundled CSAMT EDI; skip if missing."""
    p = edi_path / "CSAMT" / "csa000.edi"
    if not p.exists():
        pytest.skip(f"Missing EDI: {p}")
    return p


@pytest.fixture(
    scope="session",
    params=[
        Path("MT") / "kap03lmt_edis" / "kap103.edi",
        Path("MT") / "SPECTRA" / "spectra01.edi",
        Path("CSAMT") / "csa000.edi",
    ],
)
def any_edi_file(request: pytest.FixtureRequest, edi_path: Path) -> Path:
    """
    Parametrized EDI path. Skips the case if a file is absent.
    """
    p = edi_path / request.param
    if not p.exists():
        pytest.skip(f"Missing EDI: {p}")
    return p


# --------------- Synthetic/minimal EDI for tests --------------


@pytest.fixture()
def simulated_edi(tmp_path: Path) -> Path:
    """
    Create a minimal valid EDI on the fly for isolated tests.
    """
    lines = [
        ">HEAD",
        "  DATAID=SIM01",
        "  LAT=26:00:00N",
        "  LONG=010:00:00E",
        "  ELEV=1000",
        "  STDVERS=SEG 1.0",
        "",
        ">INFO",
        "  PROJECT=SIM",
        "  PROCESSEDBY=pyCSAMT",
        "  PROCESSINGSOFTWARE=pyCSAMT",
        "",
        ">=MTSECT",
        "  SECTID=SIM01",
        "  NFREQ=2",
        "",
        ">!****FREQUENCIES****!",
        ">FREQ  //2",
        "  1.000000E+02  2.000000E+02",
        "",
        ">ZXXR ROT=ZROT  //2",
        "  1.000000E+00  1.000000E+00",
        "",
        ">END",
    ]
    p = tmp_path / "simulated.edi"
    p.write_text("\n".join(lines), encoding="utf-8")
    return p


# ----------------------- EDIS collection -----------------------


@pytest.fixture(scope="session")
def edis_path(project_root: Path) -> Path:
    """
    Base path to a folder holding multiple EDI files.
    """
    p = project_root / "data" / "edis"
    if not p.exists():
        pytest.skip(f"Missing EDIS folder: {p}")
    return p


@pytest.fixture(scope="session")
def edis_files(edis_path: Path) -> list[Path]:
    """
    All *.edi files discovered under the EDIS folder.
    """
    files = sorted(edis_path.rglob("*.edi"))
    if not files:
        pytest.skip(f"No EDI files found in: {edis_path}")
    return files


@pytest.fixture(scope="session")
def any_edis_file(edis_files: list[Path]) -> Path:
    """
    A single EDI from the collection (first after sorting).
    """
    return edis_files[0]


@pytest.fixture(scope="session")
def edi_collection(edis_path: Path):
    """
    Parsed EDICollection built from the EDIS folder.
    """
    try:
        from pycsamt.seg.cbase import CoreParser
        from pycsamt.seg.collection import EDICollection
    except Exception as exc:  # pragma: no cover
        pytest.skip(f"Cannot import EDICollection: {exc}")

    cor = CoreParser(recursive=True, verbose=0)
    edis = cor.parse([edis_path])
    if not edis:
        pytest.skip("EDICollection parsed 0 items.")
    col = EDICollection(items=edis)
    return col


# ------------------------ J (Jones) data ----------------------
@pytest.fixture(scope="session")
def j_path(project_root: Path) -> Path:
    """Base path to J single-file samples."""
    return project_root / "data" / "j"


@pytest.fixture(scope="session")
def j_single_file(j_path: Path) -> Path:
    """Single J file: kb0-s001.txt; skip if missing."""
    p = j_path / "kb0-s001.txt"
    if not p.exists():
        pytest.skip(f"Missing J file: {p}")
    return p


@pytest.fixture(scope="session")
def jc_path(project_root: Path) -> Path:
    """Base path to a J collection (S00.dat, ...)."""
    p = project_root / "data" / "j" / "nia"
    if not p.exists():
        pytest.skip(f"Missing J collection: {p}")
    return p


@pytest.fixture(scope="session")
def jc_files(jc_path: Path) -> list[Path]:
    """All S*.dat files discovered under jc_path."""
    files = sorted(jc_path.glob("S*.dat"))
    if not files:
        pytest.skip(f"No J files found in: {jc_path}")
    return files


@pytest.fixture(scope="session")
def any_jc_file(jc_files: list[Path]) -> Path:
    """A single J file from the collection (first sorted)."""
    return jc_files[0]
