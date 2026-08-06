# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Download and build the MARE2DEM Fortran source.

This module provides :class:`SourceManager`, which handles the
full lifecycle of the MARE2DEM source tree:

1. **Locate** – resolve where the sources should live, using a
   four-level fallback that works for both development installs
   and PyPI installs on Linux, macOS, and Windows (WSL).
2. **Download** – clone from Bitbucket with ``git`` or fetch a
   ``tar.gz`` archive via ``requests`` with a ``tqdm`` progress
   bar drawn from :func:`~pycsamt.api.view.progress.iter_progress`.
3. **Build** – auto-generate a platform-specific Make include
   file and run ``make INCLUDE=<inc>``, or accept a user-supplied
   include file for cluster environments.
4. **Locate binary** – resolve the compiled ``MARE2DEM``
   executable so :class:`~pycsamt.models.mare2dem.runner.Mare2DEMRunner`
   can use it immediately.

Source-directory resolution order
----------------------------------
1. ``source_dir`` argument to :class:`SourceManager`.
2. ``config.source_dir`` field of :class:`Mare2DEMConfig`.
3. ``PYCSAMT_MARE2DEM_SOURCE`` environment variable.
4. Bundled ``_source/`` directory inside the installed package —
   only used when it is writable (development/editable installs).
5. Platform user-data directory — the robust fallback for PyPI
   installs where site-packages are read-only:

   * Linux / WSL:   ``~/.local/share/pycsamt/mare2dem/``
   * macOS:         ``~/Library/Application Support/pycsamt/mare2dem/``
   * Windows:       ``%LOCALAPPDATA%\\pycsamt\\mare2dem\\``
"""

from __future__ import annotations

import os
import platform
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

from .base import Mare2DEMBase
from .doc import _mare2dem_param_docs as _params

if TYPE_CHECKING:
    from .config import Mare2DEMConfig

__all__ = ["SourceManager"]

# ---------------------------------------------------------------------------
# Bitbucket URLs
# ---------------------------------------------------------------------------

_REPO_URL = "https://bitbucket.org/mare2dem/mare2dem_source"
_ARCHIVE_URL = (
    "https://bitbucket.org/mare2dem/mare2dem_source/get/master.tar.gz"
)

# Marker file that proves source tree is populated (not just the .gitkeep)
_SOURCE_MARKER = "Makefile"

# Expected binary name (case-sensitive on Linux/macOS)
_BINARY_NAME = "MARE2DEM"

# Triangle mesh generator is compiled alongside MARE2DEM itself (see
# ``TRICOPTS`` in :func:`_generate_inc`), so its binary name is not
# fixed by pycsamt the way ``_BINARY_NAME`` is -- try both the
# lowercase (upstream Shewchuk convention) and MARE2DEM's own
# capitalized build product name.
_TRIANGLE_BINARY_NAMES = ("triangle", "Triangle")


# ---------------------------------------------------------------------------
# Platform user-data directory (no external dependencies)
# ---------------------------------------------------------------------------


def _user_data_dir() -> Path:
    """Return the platform-appropriate user-data directory for pycsamt."""
    system = platform.system()
    if system == "Darwin":
        base = Path.home() / "Library" / "Application Support"
    elif system == "Windows":
        base = Path(
            os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local")
        )
    else:
        # Linux, FreeBSD, WSL
        base = Path(
            os.environ.get("XDG_DATA_HOME", Path.home() / ".local" / "share")
        )
    return base / "pycsamt" / "mare2dem"


# ---------------------------------------------------------------------------
# Compiler auto-detection
# ---------------------------------------------------------------------------


def _detect_fc() -> str:
    """Return the best available MPI-Fortran compiler name.

    ``mpiifx`` (wrapping ``ifx``) is tried before the classic ``mpiifort``
    (wrapping ``ifort``): current Intel oneAPI releases (2024.0+) either do
    not ship ``ifort``/``mpiifort`` at all, or ship a ``mpiifort`` that
    internally tries to exec the now-removed classic ``ifort`` and fails
    with ``eval: ifort: not found`` -- confirmed by a real build against a
    2026.1 oneAPI install. ``mpiifort`` is still tried as a fallback for
    older oneAPI installations where it genuinely works.
    """
    for candidate in ("mpiifx", "mpiifort", "mpifort", "mpif90"):
        if shutil.which(candidate):
            return candidate
    return "mpifort"


def _detect_cc() -> str:
    """Return the best available MPI-C compiler name.

    Same reasoning as :func:`_detect_fc`: ``mpiicx`` (wrapping ``icx``)
    before the classic ``mpiicc`` (wrapping the now-removed ``icc``).
    """
    for candidate in ("mpiicx", "mpiicc", "mpicc"):
        if shutil.which(candidate):
            return candidate
    return "mpicc"


def _detect_mkl() -> str | None:
    """Return the MKLROOT path from the environment or None."""
    mklroot = os.environ.get("MKLROOT")
    if mklroot and Path(mklroot).is_dir():
        return mklroot
    for candidate in (
        "/opt/intel/oneapi/mkl/latest",
        "/opt/intel/mkl",
        "/usr/local/intel/mkl",
    ):
        if Path(candidate).is_dir():
            return candidate
    return None


def _is_intel_compiler(name: str) -> bool:
    """Return True when *name* wraps the Intel Fortran compiler."""
    try:
        out = subprocess.check_output(
            [name, "--version"], stderr=subprocess.STDOUT, text=True
        )
        return "Intel" in out or "ifort" in out.lower()
    except Exception:
        return False


# ---------------------------------------------------------------------------
# Include-file generator
# ---------------------------------------------------------------------------


def _generate_inc(
    fc: str,
    cc: str,
    mklroot: str | None,
    dest: Path,
) -> Path:
    """Write a Make include file and return its path.

    Parameters
    ----------
    fc :
        MPI-Fortran compiler (e.g. ``"mpiifx"``, ``"mpiifort"``, or
        ``"mpifort"``).
    cc :
        MPI-C compiler (e.g. ``"mpiicx"``, ``"mpiicc"``, or ``"mpicc"``).
    mklroot :
        Path to the Intel MKL root directory or ``None``.
    dest :
        Destination directory for the generated file.

    Returns
    -------
    pathlib.Path
        Path of the written include file.
    """
    dest.mkdir(parents=True, exist_ok=True)
    inc_path = dest / "pycsamt_auto.inc"

    intel = _is_intel_compiler(fc)

    # Real, confirmed-by-build flag differences from the vendored example
    # include files (habanero.inc/macos.inc), all found by actually
    # compiling against a current oneAPI release, not inferred:
    #
    # * Fortran needs ``-fpp`` (not GNU ``-cpp``) for Intel, since
    #   call_triangle.f90's ``#IF DEFINED(...)``-style directives are
    #   uppercase and only ``-fpp`` accepts them case-insensitively.
    # * ``ifx`` specifically (not the classic ``ifort``) needs
    #   ``-cxxlib`` or the final link fails; untested for the classic
    #   compiler since it no longer exists in a current oneAPI install to
    #   test against, so it is only added for the modern path.
    # * The vendored BLACS/ScaLAPACK C sources (e.g. ``igsum2d_.c``) rely
    #   on 1990s-style implicit function declarations that both ``icx``
    #   and modern GCC reject by default under C99+ -- ``-std=gnu89``
    #   permits them, matching how this code was originally written.
    # * ``TRICOPTS``'s vendored example value (``-fp-model precise
    #   -fp-model source``) is classic ``icc`` syntax; ``icx`` rejects it
    #   outright (``unsupported argument 'source' to option
    #   '-ffp-model='``). Not required for correctness, so dropped rather
    #   than translated.
    if intel:
        fflags = "-cxxlib -O2 -fpp -fPIC" if "ifx" in fc else "-O2 -fpp -fPIC"
        cflags = "-O2 -fPIC -std=gnu89"
    else:
        fflags = "-O2 -cpp -fPIC"
        cflags = "-O2 -fPIC"
    tricopts = "-O2 -fPIC"

    # ``xiar`` (the classic Intel archiver) no longer exists in current
    # oneAPI releases; plain GNU ``ar`` links ifx-compiled .o files on
    # Linux/macOS fine -- confirmed by a real build, so this is no longer
    # conditional on the detected compiler.
    arch_tool = "ar"
    ranlib = "ranlib"

    if mklroot:
        mkllib = (
            f"-L{mklroot}/lib -Wl,-rpath,{mklroot}/lib -I{mklroot}/include "
            "-lmkl_intel_lp64 -lmkl_sequential -lmkl_core "
            "-lpthread -lm -ldl"
        )
    else:
        # Fallback: OpenBLAS or system LAPACK — may not link fully
        mkllib = "-lopenblas -lpthread -lm -ldl"

    system = platform.system()
    macos_sign = (
        "codesign --force --deep --sign - ./MARE2DEM"
        if system == "Darwin"
        else ""
    )

    lines = [
        "# Auto-generated by pycsamt SourceManager — edit if needed",
        f"   FC        = {fc}",
        f"   FFLAGS    = {fflags}",
        f"   CC        = {cc}",
        f"   CFLAGS    = {cflags}",
        f"    TRICOPTS = {tricopts}",
        f"   ARCH      = {arch_tool}",
        "   ARCHFLAGS = ruv",
        f"   RANLIB    = {ranlib}",
        f"    MKLLIB   = {mkllib}",
    ]
    if macos_sign:
        lines.append(f"    MACOSCODESIGN={macos_sign}")

    inc_path.write_text("\n".join(lines) + "\n")
    return inc_path


# ---------------------------------------------------------------------------
# SourceManager
# ---------------------------------------------------------------------------


class SourceManager(Mare2DEMBase):
    """Download, build, and locate the MARE2DEM binary.

    See module docstring for the full source-directory resolution
    order and usage examples.
    """

    REPO_URL: str = _REPO_URL
    ARCHIVE_URL: str = _ARCHIVE_URL

    def __init__(
        self,
        config: Mare2DEMConfig | None = None,
        source_dir: str | Path | None = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        # Lazy import to avoid circular dependency
        from .config import Mare2DEMConfig

        self.config: Mare2DEMConfig = config or Mare2DEMConfig()
        self._explicit_source_dir: Path | None = (
            Path(source_dir) if source_dir is not None else None
        )

    # ------------------------------------------------------------------
    # Directory resolution
    # ------------------------------------------------------------------

    def resolve_source_dir(self) -> Path:
        """Return the directory where MARE2DEM sources should live.

        The resolution order is:

        1. Explicit ``source_dir`` argument passed to the constructor.
        2. ``config.source_dir`` field.
        3. ``PYCSAMT_MARE2DEM_SOURCE`` environment variable.
        4. Bundled ``_source/`` inside the package (only when writable
           — i.e. editable / development installs).
        5. Platform user-data directory (always writable).

        Returns
        -------
        pathlib.Path
            Resolved directory (created if it does not yet exist).
        """
        # 1. explicit constructor arg
        if self._explicit_source_dir is not None:
            d = self._explicit_source_dir
            d.mkdir(parents=True, exist_ok=True)
            return d

        # 2. config field
        if self.config.source_dir is not None:
            d = Path(self.config.source_dir)
            d.mkdir(parents=True, exist_ok=True)
            return d

        # 3. environment variable
        env = os.environ.get("PYCSAMT_MARE2DEM_SOURCE")
        if env:
            d = Path(env)
            d.mkdir(parents=True, exist_ok=True)
            return d

        # 4. bundled _source/ — only if writable (dev install)
        pkg_source = Path(__file__).parent / "_source"
        if pkg_source.is_dir():
            try:
                test_file = pkg_source / ".write_test"
                test_file.touch()
                test_file.unlink()
                return pkg_source
            except OSError:
                pass

        # 5. platform user-data fallback
        d = _user_data_dir()
        d.mkdir(parents=True, exist_ok=True)
        return d

    def resolve_binary(self) -> Path | None:
        """Return the path to the compiled MARE2DEM binary or ``None``.

        Resolution order:

        1. Binary name found on ``PATH``.
        2. ``<source_dir>/MARE2DEM``.
        3. Platform user-data directory ``MARE2DEM``.
        """
        name = self.config.binary
        if shutil.which(name):
            return Path(name)
        src = self.resolve_source_dir()
        for candidate in (src / name, _user_data_dir() / name):
            if candidate.is_file() and os.access(candidate, os.X_OK):
                return candidate
        return None

    def resolve_triangle_binary(self, name: str | None = None) -> Path | None:
        """Return the path to a Triangle mesh-generator executable or ``None``.

        Triangle meshing is a lighter-weight, independent concern from the
        MARE2DEM solver itself (no MPI/Fortran/LAPACK build required), so
        this resolves separately from :meth:`resolve_binary` rather than
        assuming a full MARE2DEM build has happened.

        Parameters
        ----------
        name : str, optional
            Explicit executable name to look for. Tries both
            ``"triangle"`` and ``"Triangle"`` when omitted.

        Returns
        -------
        pathlib.Path or None
            Resolved executable, or ``None`` if not found anywhere.

        Resolution order (per candidate name):

        1. ``PATH`` lookup via :func:`shutil.which`.
        2. ``<source_dir>/<name>`` -- Triangle is typically a byproduct of
           :meth:`build`, since MARE2DEM's own Makefile compiles it.
        3. Platform user-data directory ``<name>``.

        Each directory candidate is checked with ``shutil.which(name,
        path=directory)`` before a literal-file fallback, so a bare name
        like ``"triangle"`` also resolves to ``"triangle.exe"`` via
        ``PATHEXT`` on Windows -- the same fix already required for
        :func:`pycsamt.forward.maxwell.external.resolve_executable`'s
        equivalent search-path loop.
        """
        names = (name,) if name else _TRIANGLE_BINARY_NAMES
        src = self.resolve_source_dir()
        search_dirs = (src, _user_data_dir())
        for candidate_name in names:
            found = shutil.which(candidate_name)
            if found:
                return Path(found)
            for directory in search_dirs:
                found = shutil.which(candidate_name, path=str(directory))
                if found:
                    return Path(found)
                candidate = directory / candidate_name
                if candidate.is_file() and os.access(candidate, os.X_OK):
                    return candidate
        return None

    # ------------------------------------------------------------------
    # Status checks
    # ------------------------------------------------------------------

    def is_downloaded(self) -> bool:
        """Return ``True`` when the source tree appears populated."""
        src = self.resolve_source_dir()
        return (src / _SOURCE_MARKER).exists()

    def is_built(self) -> bool:
        """Return ``True`` when the MARE2DEM binary exists."""
        return self.resolve_binary() is not None

    # ------------------------------------------------------------------
    # Download
    # ------------------------------------------------------------------

    def download(
        self,
        *,
        method: str = "auto",
        force: bool = False,
    ) -> Path:
        """Download the MARE2DEM source tree.

        Parameters
        ----------
        method : {"auto", "git", "archive"}, default "auto"
            Download strategy. ``"auto"`` tries ``"git"`` first and
            falls back to ``"archive"``.
        force : bool, default False
            Re-download even if the source tree is already present.

        Returns
        -------
        pathlib.Path
            Source directory containing the downloaded tree.

        Raises
        ------
        RuntimeError
            When neither ``git`` nor ``requests`` is available, or
            when the download fails.
        """
        src = self.resolve_source_dir()

        if not force and self.is_downloaded():
            if self.verbose:
                self.logger.info(
                    "SourceManager: sources already at %s — skipping download "
                    "(use force=True to re-download).",
                    src,
                )
            return src

        method = method.lower()
        if method == "auto":
            method = "git" if shutil.which("git") else "archive"

        if method == "git":
            self._download_git(src)
        elif method == "archive":
            self._download_archive(src)
        else:
            raise ValueError(
                f"Unknown download method '{method}'. "
                "Choose 'auto', 'git', or 'archive'."
            )

        if self.verbose:
            self.logger.info(
                "SourceManager: MARE2DEM source ready at %s.", src
            )
        return src

    def _download_git(self, dest: Path) -> None:
        """Clone the MARE2DEM repository into *dest*."""
        if dest.exists() and any(dest.iterdir()):
            shutil.rmtree(dest)
        dest.mkdir(parents=True, exist_ok=True)

        cmd = ["git", "clone", "--depth", "1", self.REPO_URL, str(dest)]
        if self.verbose:
            self.logger.info("SourceManager: %s", " ".join(cmd))

        print(f"Cloning MARE2DEM source from {self.REPO_URL} …")
        proc = subprocess.run(cmd, check=False)
        if proc.returncode != 0:
            raise RuntimeError(
                f"git clone failed (exit code {proc.returncode}). "
                "Check your network connection or try method='archive'."
            )

        # Strip the clone's own .git/ so *dest* holds plain files, matching
        # what the "archive" method produces. Left in place, a nested .git
        # turns *dest* into an embedded repository from the parent repo's
        # point of view -- `git add` on a dev checkout would either skip it
        # silently (it also matches the _source/* ignore pattern) or, for
        # anyone forcing it, record a dangling gitlink with no
        # corresponding .gitmodules entry.
        shutil.rmtree(dest / ".git", ignore_errors=True)

    def _download_archive(self, dest: Path) -> None:
        """Download and extract the MARE2DEM tar.gz archive into *dest*."""
        try:
            import requests
        except ImportError as exc:
            raise RuntimeError(
                "The 'requests' package is required for archive download. "
                "Install it with:  pip install requests"
            ) from exc

        dest.mkdir(parents=True, exist_ok=True)
        url = self.ARCHIVE_URL

        print(f"Downloading MARE2DEM archive from {url} …")

        with requests.get(url, stream=True, timeout=120) as resp:
            resp.raise_for_status()
            total = int(resp.headers.get("content-length", 0)) or None
            chunk_size = 1 << 14  # 16 KiB

            try:
                from tqdm import tqdm as _tqdm

                progress = _tqdm(
                    total=total,
                    unit="B",
                    unit_scale=True,
                    unit_divisor=1024,
                    desc="MARE2DEM source",
                    leave=False,
                )
                ctx = progress
            except ImportError:
                ctx = None

            with tempfile.NamedTemporaryFile(
                suffix=".tar.gz", delete=False
            ) as tmp:
                tmp_path = Path(tmp.name)
                for chunk in resp.iter_content(chunk_size=chunk_size):
                    if chunk:
                        tmp.write(chunk)
                        if ctx is not None:
                            ctx.update(len(chunk))
            if ctx is not None:
                ctx.close()

        print("Extracting archive …")
        with tarfile.open(tmp_path, "r:gz") as tar:
            members = tar.getmembers()
            # Strip the top-level directory that Bitbucket adds
            prefix = members[0].name.split("/")[0] if members else ""
            for member in members:
                if prefix and member.name.startswith(prefix + "/"):
                    member.name = member.name[len(prefix) + 1 :]
                if member.name:
                    tar.extract(member, dest)  # nosec — controlled source

        tmp_path.unlink(missing_ok=True)

    # ------------------------------------------------------------------
    # Build
    # ------------------------------------------------------------------

    def build(
        self,
        *,
        inc_file: str | Path | None = None,
        fc: str | None = None,
        cc: str | None = None,
        clean_first: bool = False,
    ) -> Path:
        """Compile MARE2DEM from the downloaded source.

        Parameters
        ----------
        inc_file : path-like or None, default None
            Explicit Make include file. When ``None``, one is
            auto-generated from the detected compilers and MKL.
        fc : str or None, default None
            MPI-Fortran compiler override. Falls back to
            ``config.fc_compiler`` then auto-detection.
        cc : str or None, default None
            MPI-C compiler override. Falls back to
            ``config.cc_compiler`` then auto-detection.
        clean_first : bool, default False
            Run ``make clean_all`` before compiling to start fresh.

        Returns
        -------
        pathlib.Path
            Path to the compiled ``MARE2DEM`` binary.

        Raises
        ------
        FileNotFoundError
            When the source tree is not present. Call
            :meth:`download` first.
        RuntimeError
            When the build fails or the binary is not found after
            the build.
        """
        if sys.platform == "win32":
            raise RuntimeError(
                "MARE2DEM cannot be compiled natively on Windows. "
                "Use Windows Subsystem for Linux (WSL) and run "
                "SourceManager from within the WSL terminal."
            )

        src = self.resolve_source_dir()
        if not self.is_downloaded():
            raise FileNotFoundError(
                f"MARE2DEM source not found at {src}. "
                "Call SourceManager.download() first."
            )

        # ---- resolve compilers ----
        _fc = fc or self.config.fc_compiler or _detect_fc()
        _cc = cc or self.config.cc_compiler or _detect_cc()
        mklroot = _detect_mkl()

        if mklroot is None:
            print(
                "\n[pycsamt] WARNING: Intel MKL not found. MARE2DEM requires "
                "MKL for ScaLAPACK/BLACS. Set the MKLROOT environment variable "
                "or run the Intel oneAPI setvars.sh script. The build may fail "
                "or produce an unusable binary without MKL.\n"
            )

        # ---- resolve include file ----
        if inc_file is not None:
            _inc = Path(inc_file)
        else:
            _inc = _generate_inc(_fc, _cc, mklroot, src / "_pycsamt_build")
            if self.verbose:
                self.logger.info(
                    "SourceManager: generated include file at %s", _inc
                )

        # ---- optional clean ----
        if clean_first:
            print("Running make clean_all …")
            subprocess.run(
                ["make", "clean_all"],
                cwd=str(src),
                check=False,
            )

        # ---- build ----
        cmd = ["make", f"INCLUDE={_inc}"]
        if self.verbose:
            self.logger.info("SourceManager: %s (cwd=%s)", " ".join(cmd), src)

        print("Building MARE2DEM (this may take several minutes) …")
        print(f"  FC = {_fc}")
        print(f"  CC = {_cc}")
        print(f"  MKLROOT = {mklroot or '(not found)'}")
        print(f"  Include file = {_inc}")

        proc = subprocess.run(cmd, cwd=str(src), check=False)
        if proc.returncode != 0:
            raise RuntimeError(
                f"MARE2DEM build failed (make exit code {proc.returncode}). "
                "Check the compiler output above. Common fixes:\n"
                "  • Ensure Intel oneAPI compilers are loaded "
                "(source /opt/intel/oneapi/setvars.sh)\n"
                "  • Set MKLROOT to your MKL installation directory\n"
                "  • Pass inc_file= with a custom make include file for "
                "your cluster"
            )

        binary = src / _BINARY_NAME
        if not binary.exists():
            raise RuntimeError(
                f"Build succeeded (exit code 0) but '{_BINARY_NAME}' was not "
                f"found in {src}. Check the Makefile output for the actual "
                "binary location."
            )

        print(f"\nMARE2DEM binary ready: {binary}\n")
        if self.verbose:
            self.logger.info("SourceManager: binary at %s", binary)
        return binary

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------

    def status(self) -> dict[str, object]:
        """Return a summary dictionary of source and build status.

        Returns
        -------
        dict
            Keys: ``source_dir``, ``downloaded``, ``binary_path``,
            ``built``, ``fc``, ``cc``, ``mklroot``.
        """
        binary = self.resolve_binary()
        return {
            "source_dir": self.resolve_source_dir(),
            "downloaded": self.is_downloaded(),
            "built": self.is_built(),
            "binary_path": binary,
            "fc": self.config.fc_compiler or _detect_fc(),
            "cc": self.config.cc_compiler or _detect_cc(),
            "mklroot": _detect_mkl(),
        }

    def print_status(self) -> None:
        """Print a human-readable source and build status report."""
        s = self.status()
        lines = [
            "MARE2DEM SourceManager status",
            "─" * 40,
            f"  source_dir  : {s['source_dir']}",
            f"  downloaded  : {s['downloaded']}",
            f"  built       : {s['built']}",
            f"  binary      : {s['binary_path'] or '(not found)'}",
            f"  FC compiler : {s['fc']}",
            f"  CC compiler : {s['cc']}",
            f"  MKLROOT     : {s['mklroot'] or '(not found — required)'}",
        ]
        print("\n".join(lines))


SourceManager.__doc__ = rf"""
Manage the MARE2DEM Fortran source: download, build, and locate.

``SourceManager`` is the entry point for obtaining a working
MARE2DEM binary. It separates source management from inversion
execution so that users who already have a compiled binary can
skip directly to :class:`~pycsamt.models.mare2dem.runner.Mare2DEMRunner`.

Source-directory resolution
---------------------------
The directory where sources are stored follows this priority:

1. ``source_dir`` constructor argument.
2. ``config.source_dir`` field.
3. ``PYCSAMT_MARE2DEM_SOURCE`` environment variable.
4. Bundled ``_source/`` inside the installed package (writable
   dev installs only).
5. Platform user-data directory — the safe PyPI fallback:

   * Linux / WSL:   ``~/.local/share/pycsamt/mare2dem/``
   * macOS:         ``~/Library/Application Support/pycsamt/mare2dem/``
   * Windows:       ``%LOCALAPPDATA%\\pycsamt\\mare2dem\\``

Parameters
----------
config : Mare2DEMConfig, optional
    Configuration object. Supplies ``source_dir``,
    ``fc_compiler``, ``cc_compiler``, and ``binary``.
source_dir : path-like, optional
    Explicit path that overrides ``config.source_dir`` and the
    environment variable.
{_params.common.verbose}
{_params.common.logger}

Notes
-----
MARE2DEM requires Intel compilers (``mpiifx``/``mpiicx`` on current
oneAPI releases, or the classic ``mpiifort``/``mpiicc`` on older
ones) and the Intel MKL. When Intel oneAPI is installed, source the
``setvars.sh`` script before calling :meth:`build`:

.. code-block:: bash

    source /opt/intel/oneapi/setvars.sh
    python -c "
    from pycsamt.models.mare2dem import SourceManager
    sm = SourceManager(verbose=1)
    sm.download()
    sm.build()
    "

Windows native builds are not supported. Use WSL2.

Examples
--------
Quick-start — download and build:

>>> from pycsamt.models.mare2dem import SourceManager
>>> sm = SourceManager(verbose=1)
>>> sm.download()          # clones from Bitbucket
>>> sm.build()             # compiles with auto-detected compilers

Check status without downloading:

>>> sm.print_status()

Point to sources already on disk:

>>> sm = SourceManager(source_dir="/data/mare2dem_source")
>>> sm.build(inc_file="/data/mare2dem_source/include/habanero.inc")

See Also
--------
Mare2DEMRunner
    Launch the compiled MARE2DEM executable for inversion.
Mare2DEMConfig
    Configuration controlling compiler selection and binary name.

References
----------
.. [SourceManager-1] Key, K. (2016). MARE2DEM: A 2-D inversion code for
   controlled-source electromagnetic and magnetotelluric data.
   *Geophysical Journal International*, 207(1), 571-588.
   doi:10.1093/gji/ggw290.
"""
