# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
solver_build — check, install and build pycsamt's external solvers.

Qt-free engine behind the desktop **Solver Builder** (and reusable from
scripts): Occam2D, ModEM 2-D/3-D and MARE2DEM.

What each route was verified with (2026-09-25):

* **Occam2D / ModEM** — vendored Fortran sources, GNU ``make`` driven from
  Python (no bash needed: the only recipe command outside ``clean`` is
  ``mkdir -p objs``, pre-created here).  Toolchain, in order: a
  *pycsamt-managed* one (created with micromamba -- no conda or admin
  rights needed), the ``pycsamt-fortran`` conda env the bash scripts
  create, then ``gfortran`` + ``make`` on PATH.  Windows: ModEM 2-D built
  in 19 s with the managed MinGW toolchain and ran standalone with its
  five runtime DLLs copied alongside.
* **MARE2DEM** — not vendored (downloaded) and needs MPI + Intel oneMKL
  (its forward solvers call MKL's sparse direct solver).  oneMKL is free,
  so no Intel *compiler* is needed: gfortran + OpenMPI + MKL from a
  managed conda-forge toolchain builds it after the source fixes in
  :mod:`pycsamt.models.mare2dem._gnu_compat`, and the result passes the
  half-space benchmark.  Built natively on Linux/macOS, and inside WSL2
  from Windows (MARE2DEM cannot be built as a native Windows program);
  the resulting binary is registered as ``wsl:<linux path>``.

Built binaries are recorded in ``~/.pycsamt/solvers.json`` so the
Inversion window can pick them up automatically.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import time
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

_IS_WIN = sys.platform == "win32"
_MODELS = Path(__file__).resolve().parent

MICROMAMBA_URL = (
    "https://github.com/mamba-org/micromamba-releases/releases/latest/"
    "download/micromamba-{platform}"
)
MARE2DEM_ARCHIVE_URL = (
    "https://bitbucket.org/mare2dem/mare2dem_source/get/master.zip"
)
# conda-forge packages of the managed toolchains
WIN_FORTRAN_PKGS = ["m2w64-gcc-fortran", "m2w64-openblas", "make"]
POSIX_FORTRAN_PKGS = ["gfortran", "make", "openblas"]
MARE2DEM_PKGS = ["gfortran", "gcc", "openmpi", "make", "mkl-devel",
                 "mkl-include", "python"]
WIN_RUNTIME_DLLS = ["libgcc_s_seh-1.dll", "libgfortran-3.dll",
                    "libopenblas.dll", "libquadmath-0.dll",
                    "libwinpthread-1.dll"]


# ── Solver catalogue ──────────────────────────────────────────────────────


@dataclass(frozen=True)
class SolverSpec:
    key: str
    label: str
    summary: str
    binary_stem: str
    source_rel: str  # relative to pycsamt/models
    vendored: bool
    engine: str  # "fortran" | "mare2dem"
    make_style: str = ""  # "modem" | "occam"
    license_note: str = ""

    def default_source_dir(self) -> Path:
        return _MODELS / self.source_rel


SOLVERS: dict[str, SolverSpec] = {
    "occam2d": SolverSpec(
        "occam2d", "Occam2D", "2-D MT smooth-model inversion (Occam).",
        "Occam2D", "occam2d/_source", True, "fortran", "occam",
    ),
    "modem2d": SolverSpec(
        "modem2d", "ModEM 2-D", "2-D MT NLCG inversion (Mod2DMT).",
        "Mod2DMT", "modem/_source/2D", True, "fortran", "modem",
    ),
    "modem3d": SolverSpec(
        "modem3d", "ModEM 3-D", "3-D MT NLCG inversion (Mod3DMT).",
        "Mod3DMT", "modem/_source/3D", True, "fortran", "modem",
    ),
    "mare2dem": SolverSpec(
        "mare2dem", "MARE2DEM",
        "2-D adaptive finite-element MT/CSEM inversion (Key, 2016).",
        "MARE2DEM", "mare2dem/_source", False, "mare2dem",
        license_note="Source downloaded from bitbucket.org/mare2dem "
                     "under its own license (not bundled with pycsamt).",
    ),
}


# ── Locations & registry ──────────────────────────────────────────────────


def toolchain_root() -> Path:
    """Per-user folder for pycsamt-managed toolchains (no admin rights)."""
    if _IS_WIN:
        base = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData"
                                   / "Local"))
        return base / "pycsamt" / "toolchain"
    base = Path(os.environ.get("XDG_DATA_HOME",
                               Path.home() / ".local" / "share"))
    return base / "pycsamt" / "toolchain"


def managed_fortran_prefix() -> Path:
    return toolchain_root() / "fortran"


def micromamba_path() -> Path:
    return toolchain_root() / "bin" / (
        "micromamba.exe" if _IS_WIN else "micromamba")


def registry_path() -> Path:
    return Path.home() / ".pycsamt" / "solvers.json"


def load_registry(path: Path | None = None) -> dict:
    p = path or registry_path()
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return {}


def register_binary(key: str, binary: str, *, path: Path | None = None,
                    **extra) -> None:
    """Record a built binary so other windows (Inversion) can use it."""
    p = path or registry_path()
    reg = load_registry(p)
    reg[key] = {"binary": str(binary), "built": time.strftime(
        "%Y-%m-%d %H:%M"), **extra}
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(reg, indent=2), encoding="utf-8")


def is_wsl_binary(binary: str | None) -> bool:
    return bool(binary) and str(binary).startswith("wsl:")


def registered_binary(key: str, *, path: Path | None = None) -> str | None:
    """Registered binary for *key* if it still exists (WSL ones trusted)."""
    entry = load_registry(path).get(key) or {}
    b = entry.get("binary")
    if not b:
        return None
    if is_wsl_binary(b) or Path(b).is_file():
        return b
    return None


def expected_binary(key: str, source_dir: Path | None = None) -> Path:
    spec = SOLVERS[key]
    src = Path(source_dir) if source_dir else spec.default_source_dir()
    name = spec.binary_stem + (".exe" if _IS_WIN and
                               spec.engine == "fortran" else "")
    return src / name


def find_binary(key: str, source_dir: Path | None = None) -> str | None:
    """Registry first, then a binary already present in the source tree."""
    reg = registered_binary(key)
    if reg:
        return reg
    b = expected_binary(key, source_dir)
    return str(b) if b.is_file() else None


@dataclass(frozen=True)
class BinaryCandidate:
    """A solver executable found on this machine and where it came from."""

    path: str
    origin: str  # "Solver Builder" | "source tree" | "PATH" | "WSL build"


# Where the Solver Builder's MARE2DEM script leaves its binary (Linux/WSL).
_M2D_BUILD_REL = "pycsamt/mare2dem/build/MARE2DEM"


def _wsl_mare2dem_build() -> str | None:
    """``wsl:<path>`` of a MARE2DEM built in WSL, found by probing WSL."""
    script = ('p="${XDG_DATA_HOME:-$HOME/.local/share}/'
              + _M2D_BUILD_REL + '"; [ -x "$p" ] && echo "$p"')
    try:
        out = subprocess.run(["wsl", "-e", "bash", "-lc", script],
                             capture_output=True, timeout=60)
    except Exception:
        return None
    path = out.stdout.decode("utf-8", "replace").strip().splitlines()
    return f"wsl:{path[-1].strip()}" if out.returncode == 0 and path else None


def discover_binaries(key: str, source_dir: Path | None = None, *,
                      probe_wsl: bool = True) -> list[BinaryCandidate]:
    """Every usable *key* executable, best first, without duplicates.

    Looks in the Solver Builder registry, the source tree, ``PATH`` and --
    for MARE2DEM -- the standard Linux/WSL build folder (the build may have
    been made outside the desktop app, so it is not always registered).
    ``probe_wsl=False`` skips the WSL call (about a second).
    """
    spec = SOLVERS[key]
    found: list[BinaryCandidate] = []

    def add(path, origin):
        if path and all(c.path != str(path) for c in found):
            found.append(BinaryCandidate(str(path), origin))

    add(registered_binary(key), "Solver Builder")
    b = expected_binary(key, source_dir)
    if b.is_file():
        add(b, "source tree")
    add(shutil.which(spec.binary_stem), "PATH")
    if key == "mare2dem":
        if _IS_WIN:
            if probe_wsl and shutil.which("wsl"):
                add(_wsl_mare2dem_build(), "WSL build")
        else:
            data = Path(os.environ.get("XDG_DATA_HOME",
                                       Path.home() / ".local" / "share"))
            local = data / _M2D_BUILD_REL
            if local.is_file():
                add(local, "local build")
    return found


def stage_binary(binary: str, dest_dir: str | Path, *,
                 move: bool = False) -> str:
    """Copy (or move) *binary* into *dest_dir* and return the new path.

    Makes an inversion folder self-contained and portable.  On Windows the
    MinGW runtime DLLs sitting next to the binary travel with it, otherwise
    the copy would not start on a machine without that toolchain.  Binaries
    living inside WSL (``wsl:``) cannot be staged into a Windows folder.
    """
    if is_wsl_binary(binary):
        raise ValueError("A WSL binary runs inside Linux and cannot be "
                         "copied into a Windows folder.")
    src = Path(binary)
    if not src.is_file():
        raise FileNotFoundError(f"Binary not found: {src}")
    dest = Path(dest_dir)
    dest.mkdir(parents=True, exist_ok=True)
    target = dest / src.name
    if target.resolve() == src.resolve():
        return str(target)
    op = shutil.move if move else shutil.copy2
    op(str(src), str(target))
    if _IS_WIN:
        for dll in WIN_RUNTIME_DLLS:
            if (src.parent / dll).is_file() and not (dest / dll).exists():
                shutil.copy2(str(src.parent / dll), str(dest / dll))
    if move:
        # Keep the Solver Builder registry pointing at the binary.
        for key, entry in load_registry().items():
            if isinstance(entry, dict) and entry.get("binary") == str(src):
                register_binary(key, str(target))
    return str(target)


# ── Toolchain detection ───────────────────────────────────────────────────


@dataclass
class FortranToolchain:
    f90: str
    make: str
    bin_dirs: list[str]
    libs: str | None  # LIBS= override, or None for the Makefile default
    origin: str  # "managed" | "conda env" | "PATH"


def _win_prefix_toolchain(prefix: Path, origin: str) -> FortranToolchain | None:
    mingw = prefix / "Library" / "mingw-w64"
    f90 = mingw / "bin" / "gfortran.exe"
    make = prefix / "Library" / "bin" / "make.exe"
    blas = mingw / "lib" / "libopenblas.a"
    if f90.is_file() and make.is_file() and blas.is_file():
        return FortranToolchain(
            str(f90), str(make), [str(mingw / "bin"),
                                  str(prefix / "Library" / "bin")],
            f"-L{(mingw / 'lib').as_posix()} -lopenblas", origin)
    return None


def _conda_env_dirs(name: str) -> list[Path]:
    roots = [Path.home() / d for d in ("anaconda3", "miniconda3",
                                        "miniforge3", "mambaforge")]
    roots.append(Path.home() / ".conda")
    for var in ("CONDA_PREFIX", "CONDA_EXE"):
        v = os.environ.get(var)
        if v:
            p = Path(v)
            roots += [p, p.parent, p.parent.parent]
    return [r / "envs" / name for r in roots]


def fortran_toolchain() -> FortranToolchain | None:
    """The toolchain a Fortran build would use, or ``None``."""
    managed = managed_fortran_prefix()
    if _IS_WIN:
        tc = _win_prefix_toolchain(managed, "managed")
        if tc:
            return tc
        env = os.environ.get("PYCSAMT_FORTRAN_ENV", "pycsamt-fortran")
        for d in _conda_env_dirs(env):
            tc = _win_prefix_toolchain(d, "conda env")
            if tc:
                return tc
    else:
        f90, make = managed / "bin" / "gfortran", managed / "bin" / "make"
        if f90.is_file() and make.is_file():
            lib = managed / "lib"
            return FortranToolchain(
                str(f90), str(make), [str(managed / "bin")],
                f"-L{lib} -Wl,-rpath,{lib} -lopenblas", "managed")
    make = shutil.which("make")
    f90 = shutil.which("gfortran") or next(
        (shutil.which(f"gfortran-{v}") for v in range(16, 8, -1)
         if shutil.which(f"gfortran-{v}")), None)
    if make and f90:
        return FortranToolchain(f90, make, [], None, "PATH")
    return None


def wsl_available() -> bool:
    """A usable WSL distribution exists (Windows only)."""
    if not _IS_WIN or not shutil.which("wsl"):
        return False
    try:
        out = subprocess.run(["wsl", "-e", "true"], capture_output=True,
                             timeout=60)
        return out.returncode == 0
    except Exception:
        return False


def to_wsl_path(path: str | Path) -> str:
    """``C:\\a\\b`` -> ``/mnt/c/a/b`` (paths already POSIX pass through)."""
    s = str(path)
    if len(s) > 1 and s[1] == ":":
        return f"/mnt/{s[0].lower()}{s[2:].replace(chr(92), '/')}"
    return s.replace("\\", "/")


# ── Readiness (dependency checklist) ──────────────────────────────────────


@dataclass
class Check:
    key: str
    label: str
    ok: bool
    detail: str = ""
    hint: str = ""  # manual instructions when missing
    auto: bool = True  # can "Install automatically" provide it?


@dataclass
class Readiness:
    checks: list[Check] = field(default_factory=list)

    @property
    def ready(self) -> bool:
        return all(c.ok for c in self.checks)

    @property
    def can_auto_install(self) -> bool:
        return all(c.ok or c.auto for c in self.checks)

    @property
    def missing(self) -> list[Check]:
        return [c for c in self.checks if not c.ok]


def _fortran_hint() -> str:
    if _IS_WIN:
        return ("Install automatically (a private MinGW gfortran toolchain, "
                "no admin rights), or install MSYS2 and add gfortran + make "
                "to PATH.")
    if sys.platform == "darwin":
        return "Install automatically, or: brew install gcc openblas make"
    return ("Install automatically (no sudo), or: sudo apt install gfortran "
            "make libopenblas-dev")


def _check_fortran(spec: SolverSpec, source_dir: Path | None) -> Readiness:
    r = Readiness()
    tc = fortran_toolchain()
    r.checks.append(Check(
        "fortran", "Fortran compiler (gfortran)", tc is not None,
        f"{tc.f90}  ({tc.origin})" if tc else "not found", _fortran_hint()))
    r.checks.append(Check(
        "make", "GNU make", tc is not None,
        tc.make if tc else "not found", _fortran_hint()))
    r.checks.append(Check(
        "blas", "LAPACK / BLAS", tc is not None and (
            tc.libs is not None or not _IS_WIN),
        "OpenBLAS (toolchain)" if tc and tc.libs else
        ("system -llapack -lblas" if tc else "not found"), _fortran_hint()))
    src = Path(source_dir) if source_dir else spec.default_source_dir()
    has_src = (src / "Makefile").is_file()
    r.checks.append(Check(
        "source", "Source code", has_src,
        str(src) if has_src else f"no Makefile in {src}",
        "Point “Source folder” at a directory containing the solver's "
        "Makefile and Fortran sources.", auto=False))
    return r


_M2D_PROBE = r"""
TC="${XDG_DATA_HOME:-$HOME/.local/share}/pycsamt/toolchain/mare2dem"
[ -d "$TC/bin" ] && export PATH="$TC/bin:$PATH"
echo "managed=$([ -x "$TC/bin/mpifort" ] && echo 1 || echo 0)"
for t in mpifort mpicc mpiifx mpiicx make python3 curl; do
  echo "$t=$(command -v $t 2>/dev/null || true)"; done
mkl=""
for r in "$TC" "${MKLROOT:-}" /opt/intel/oneapi/mkl/latest /usr; do
  [ -n "$r" ] || continue
  for inc in "$r/include" "$r/include/mkl"; do
    [ -f "$inc/mkl_dss.f90" ] && { mkl="$r"; break 2; }; done; done
echo "mkl=$mkl"
echo "arch=$(uname -m)"
echo "home=$HOME"
"""


def probe_mare2dem_env() -> dict:
    """key=value facts about the build environment (native or WSL)."""
    argv = (["wsl", "-e", "bash", "-s"] if _IS_WIN else ["bash", "-s"])
    try:
        # Bytes, not text: a text-mode pipe on Windows turns LF into
        # CRLF, and bash (in WSL) then sees every line ending in CR.
        out = subprocess.run(argv, input=_M2D_PROBE.encode("utf-8"),
                             capture_output=True, timeout=120)
    except Exception as exc:
        return {"error": str(exc)}
    facts = {}
    for line in out.stdout.decode("utf-8", "replace").splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            facts[k.strip()] = v.strip()
    return facts


def _check_mare2dem(source_dir: Path | None) -> Readiness:
    r = Readiness()
    if _IS_WIN:
        ok = wsl_available()
        r.checks.append(Check(
            "wsl", "Linux environment (WSL2)", ok,
            "available" if ok else "not installed",
            "MARE2DEM builds only as a Linux program. In an administrator "
            "terminal run:  wsl --install   then restart Windows.",
            auto=False))
        if not ok:
            return r
    facts = probe_mare2dem_env()
    if "error" in facts:
        r.checks.append(Check("env", "Build environment", False,
                              facts["error"], auto=False))
        return r
    where = " (WSL)" if _IS_WIN else ""
    mpi = bool(facts.get("mpifort") and facts.get("mpicc")) or bool(
        facts.get("mpiifx") and facts.get("mpiicx"))
    hint = ("Install automatically: a private toolchain (gfortran, OpenMPI, "
            "make, free oneMKL) is created with micromamba, no sudo needed.")
    r.checks.append(Check(
        "mpi", f"MPI compilers{where}", mpi,
        facts.get("mpifort") or facts.get("mpiifx") or "not found", hint))
    r.checks.append(Check("make", f"GNU make{where}",
                          bool(facts.get("make")),
                          facts.get("make") or "not found", hint))
    r.checks.append(Check("mkl", f"Intel oneMKL (free){where}",
                          bool(facts.get("mkl")),
                          facts.get("mkl") or "not found", hint))
    r.checks.append(Check(
        "python", f"Python 3{where} (applies source fixes)",
        bool(facts.get("python3")), facts.get("python3") or "not found",
        hint))
    src = Path(source_dir) if source_dir else SOLVERS[
        "mare2dem"].default_source_dir()
    local = (src / "Makefile").is_file()
    r.checks.append(Check(
        "source", "Source code", local or bool(facts.get("curl") or
                                               facts.get("python3")),
        str(src) if local else "will be downloaded from bitbucket.org",
        "Needs internet access (curl) to download the MARE2DEM source, or "
        "set “Source folder” to an existing copy.",
        auto=False))
    return r


def check(key: str, source_dir: Path | None = None) -> Readiness:
    spec = SOLVERS[key]
    if spec.engine == "mare2dem":
        return _check_mare2dem(source_dir)
    return _check_fortran(spec, source_dir)


# ── Execution plumbing ────────────────────────────────────────────────────


class Cancelled(Exception):
    """The user stopped the operation."""


@dataclass
class Context:
    """Callbacks a runner (GUI worker, CLI) supplies."""

    log: Callable[[str], None] = print
    progress: Callable[[float | None], None] = lambda f: None  # 0..1 / None
    stage: Callable[[str], None] = lambda s: None
    cancelled: Callable[[], bool] = lambda: False


@dataclass
class Step:
    label: str
    run: Callable[[Context], None]


def _stream(argv: list[str], ctx: Context, *, cwd: str | None = None,
            env: dict | None = None, stdin_text: str | None = None,
            on_line: Callable[[str], None] | None = None) -> int:
    """Run *argv*, forwarding output lines; honours cancellation."""
    ctx.log("$ " + " ".join(argv[:6]) + (" …" if len(argv) > 6 else ""))
    flags = subprocess.CREATE_NO_WINDOW if _IS_WIN else 0
    # Binary pipes: text mode on Windows would rewrite the script's LF
    # as CRLF on its way into bash (WSL), breaking every line.
    proc = subprocess.Popen(
        argv, cwd=cwd, env=env, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, stdin=subprocess.PIPE if stdin_text
        is not None else subprocess.DEVNULL, creationflags=flags)
    if stdin_text is not None:
        proc.stdin.write(stdin_text.encode("utf-8"))
        proc.stdin.close()
    for raw in iter(proc.stdout.readline, b""):
        line = raw.decode("utf-8", "replace").rstrip()
        if on_line:
            on_line(line)
        ctx.log(line)
        if ctx.cancelled():
            proc.terminate()
            proc.wait(timeout=30)
            raise Cancelled()
    return proc.wait()


def _download(url: str, dest: Path, ctx: Context) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    ctx.log(f"Downloading {url}")
    with urllib.request.urlopen(url, timeout=60) as resp, open(tmp, "wb") as f:
        total = int(resp.headers.get("Content-Length") or 0)
        got = 0
        while True:
            if ctx.cancelled():
                raise Cancelled()
            chunk = resp.read(1 << 16)
            if not chunk:
                break
            f.write(chunk)
            got += len(chunk)
            ctx.progress(got / total if total else None)
    tmp.replace(dest)
    if not _IS_WIN:
        dest.chmod(0o755)
    ctx.log(f"Saved {dest} ({got / 1e6:.1f} MB)")


def _micromamba_platform() -> str:
    if _IS_WIN:
        return "win-64"
    import platform

    mach = platform.machine().lower()
    if sys.platform == "darwin":
        return "osx-arm64" if mach in ("arm64", "aarch64") else "osx-64"
    return "linux-aarch64" if mach in ("arm64", "aarch64") else "linux-64"


# ── Install plans (root-free, via micromamba) ─────────────────────────────


def _ensure_micromamba(ctx: Context) -> Path:
    exe = micromamba_path()
    if not exe.is_file():
        ctx.stage("Downloading package manager (micromamba)")
        _download(MICROMAMBA_URL.format(platform=_micromamba_platform()),
                  exe, ctx)
    return exe


def _create_env(prefix: Path, pkgs: list[str], ctx: Context) -> None:
    exe = _ensure_micromamba(ctx)
    ctx.stage("Installing compilers and libraries")
    ctx.progress(None)
    env = dict(os.environ, MAMBA_ROOT_PREFIX=str(toolchain_root() / "root"))
    rc = _stream([str(exe), "create", "-y", "-p", str(prefix), "-c",
                  "conda-forge", *pkgs], ctx, env=env)
    if rc != 0:
        raise RuntimeError(f"Toolchain installation failed (exit {rc}).")


def install_plan(key: str) -> list[Step]:
    """Steps that provide every auto-installable missing dependency."""
    spec = SOLVERS[key]
    if spec.engine == "fortran":
        pkgs = WIN_FORTRAN_PKGS if _IS_WIN else POSIX_FORTRAN_PKGS
        return [Step("Install Fortran toolchain",
                     lambda ctx: _create_env(managed_fortran_prefix(), pkgs,
                                             ctx))]
    return [Step("Install MARE2DEM toolchain",
                 lambda ctx: _run_mare2dem_script(ctx, install=True,
                                                  build=False))]


# ── Build plans ───────────────────────────────────────────────────────────


def _clean_tree(src: Path, spec: SolverSpec, ctx: Context) -> None:
    n = 0
    for pat in ("*.o", "*.mod", "objs/*.o", "objs/*.mod"):
        for f in src.glob(pat):
            f.unlink(missing_ok=True)
            n += 1
    for f in (src / spec.binary_stem, src / f"{spec.binary_stem}.exe"):
        if f.exists():
            f.unlink()
    ctx.log(f"Cleaned {n} object/module files.")


def _build_fortran(spec: SolverSpec, src: Path, clean: bool,
                   ctx: Context) -> str:
    tc = fortran_toolchain()
    if tc is None:
        raise RuntimeError("No Fortran toolchain — install it first.")
    if not (src / "Makefile").is_file():
        raise RuntimeError(f"No Makefile in {src}")
    if clean:
        _clean_tree(src, spec, ctx)
    env = dict(os.environ)
    env["PATH"] = os.pathsep.join(tc.bin_dirs + [env.get("PATH", "")])
    if spec.make_style == "modem":
        (src / "objs").mkdir(exist_ok=True)  # recipe is `mkdir -p` (no sh)
        args = [tc.make, "-k", f"F90={tc.f90}"]
        if tc.libs:
            args.append(f"LIBS={tc.libs}")
    else:  # Occam2D: default target runs `clean` (rm) first -- skip it
        args = [tc.make, "-k", f"FC90={tc.f90}",
                "FCFLAGS=-O2 -ffree-line-length-none"]
        if tc.libs:
            args.append(f"LIBS={tc.libs}")
        args.append("Occam2D")
    total = max(len(list(src.glob("*.f90"))), 1)
    done = [0]

    def on_line(line: str) -> None:
        if " -c " in line and (".f90" in line or ".f" in line):
            done[0] += 1
            ctx.progress(min(done[0] / total, 0.97))

    binary = expected_binary(spec.key, src)
    ctx.stage("Compiling")
    # The vendored ModEM Makefiles are not in strict module-dependency
    # order; a few `make -k` passes converge (see _solver_build/README.md).
    for attempt in range(1, 6):
        if ctx.cancelled():
            raise Cancelled()
        ctx.log(f"── make pass {attempt}")
        _stream(args, ctx, cwd=str(src), env=env, on_line=on_line)
        if binary.is_file():
            break
    if not binary.is_file():
        raise RuntimeError("Compilation did not produce "
                           f"{binary.name}; see the log above.")
    if _IS_WIN and tc.bin_dirs:
        for dll in WIN_RUNTIME_DLLS:
            s = Path(tc.bin_dirs[0]) / dll
            if s.is_file():
                shutil.copy2(s, src / dll)
        ctx.log("Copied MinGW runtime DLLs next to the binary.")
    ctx.progress(1.0)
    return str(binary)


# The MARE2DEM script runs natively (Linux/macOS) or inside WSL (fed on
# stdin), so the Windows app needs no pycsamt install on the Linux side.
_M2D_SCRIPT = r"""
set -eu
DATA="${XDG_DATA_HOME:-$HOME/.local/share}/pycsamt"
TC="$DATA/toolchain/mare2dem"
BUILD="$DATA/mare2dem/build"
DO_INSTALL=@@INSTALL@@
DO_BUILD=@@BUILD@@
CLEAN=@@CLEAN@@
SRC_IN='@@SOURCE@@'
stage() { echo "@@STAGE $*"; }

if [ "$DO_INSTALL" = 1 ] && [ ! -x "$TC/bin/mpifort" ]; then
  MM="$DATA/toolchain/bin/micromamba"
  if [ ! -x "$MM" ]; then
    stage "Downloading package manager (micromamba)"
    mkdir -p "$(dirname "$MM")"
    curl -fsSL -o "$MM" "@@MM_URL@@"; chmod +x "$MM"
  fi
  stage "Installing compilers, OpenMPI and oneMKL"
  export MAMBA_ROOT_PREFIX="$DATA/toolchain/root"
  "$MM" create -y -p "$TC" -c conda-forge @@PKGS@@
fi
[ "$DO_BUILD" = 1 ] || { echo "@@DONE install"; exit 0; }

[ -d "$TC/bin" ] && export PATH="$TC/bin:$PATH"
FC=mpifort; CC=mpicc; INTEL=""
if ! command -v mpifort >/dev/null && command -v mpiifx >/dev/null; then
  FC=mpiifx; CC=mpiicx; INTEL="--intel"; fi
command -v "$FC" >/dev/null || { echo "MPI Fortran compiler not found"; exit 3; }
PY=$(command -v python3 || command -v python || true)
[ -n "$PY" ] || { echo "python3 not found"; exit 3; }

stage "Preparing source"
if [ "$CLEAN" = 1 ]; then rm -rf "$BUILD"; fi
if [ ! -f "$BUILD/Makefile" ]; then
  mkdir -p "$BUILD"
  if [ -n "$SRC_IN" ] && [ -f "$SRC_IN/Makefile" ]; then
    cp -r "$SRC_IN"/. "$BUILD"/
  else
    stage "Downloading MARE2DEM source"
    "$PY" - "$BUILD" <<'PYDL'
import io, shutil, sys, urllib.request, zipfile, pathlib
dest = pathlib.Path(sys.argv[1])
data = urllib.request.urlopen("@@ARCHIVE_URL@@", timeout=120).read()
with zipfile.ZipFile(io.BytesIO(data)) as z:
    root = z.namelist()[0].split("/")[0]
    z.extractall(dest.parent / "_dl")
for p in (dest.parent / "_dl" / root).iterdir():
    shutil.move(str(p), dest / p.name)
shutil.rmtree(dest.parent / "_dl")
print("source downloaded")
PYDL
  fi
fi

stage "Applying source fixes"
mkdir -p "$BUILD/_pycsamt_build"
cat > "$BUILD/_pycsamt_build/gnu_compat.py" <<'PYCOMPAT'
@@COMPAT@@
PYCOMPAT
"$PY" "$BUILD/_pycsamt_build/gnu_compat.py" "$BUILD" $INTEL

MKL=""
for r in "$TC" "${MKLROOT:-}" /opt/intel/oneapi/mkl/latest /usr; do
  [ -n "$r" ] || continue
  for inc in "$r/include" "$r/include/mkl"; do
    [ -f "$inc/mkl_dss.f90" ] && { MKL="$r"; MKLINC="$inc"; break 2; }; done
done
[ -n "$MKL" ] || { echo "oneMKL not found"; exit 3; }
MKLLIB=""
for d in "$MKL/lib" "$MKL/lib/intel64" "$MKL/lib/x86_64-linux-gnu"; do
  ls "$d"/libmkl_core.* >/dev/null 2>&1 && { MKLLIB="$d"; break; }; done
lk() { if [ -e "$MKLLIB/lib$1.so" ] || [ -e "$MKLLIB/lib$1.a" ]; then echo "-l$1";
       else echo "-l:$(basename "$(ls "$MKLLIB"/lib$1.so.* | head -1)")"; fi; }
if [ -n "$INTEL" ]; then IFACE=mkl_intel_lp64
  FFLAGS="-O2 -fpp -fPIC"; [ "$FC" = mpiifx ] && FFLAGS="-cxxlib $FFLAGS"
  TRIC="-O2 -fPIC"
else IFACE=mkl_gf_lp64
  FFLAGS="-O2 -cpp -fPIC -fallow-argument-mismatch -fdec-format-defaults -I$MKLINC"
  TRIC="-O2 -fPIC -std=gnu89"
fi
cat > "$BUILD/_pycsamt_build/auto.inc" <<INC
   FC        = $FC
   FFLAGS    = $FFLAGS
   CC        = $CC
   CFLAGS    = -O2 -fPIC -std=gnu89
    TRICOPTS = $TRIC
   ARCH      = ar
   ARCHFLAGS = ruv
   RANLIB    = ranlib
    MKLLIB   = -L$MKLLIB -Wl,-rpath,$MKLLIB -I$MKLINC $(lk $IFACE) $(lk mkl_sequential) $(lk mkl_core) -lpthread -lm -ldl
INC

cd "$BUILD"
echo "@@TOTAL $(find . -name '*.f' -o -name '*.f90' -o -name '*.c' | wc -l)"
stage "Compiling (ScaLAPACK + MARE2DEM)"
make INCLUDE=_pycsamt_build/auto.inc
[ -x "$BUILD/MARE2DEM" ] || { echo "MARE2DEM binary not produced"; exit 4; }
echo "@@BINARY $BUILD/MARE2DEM"
"""


def _mare2dem_script(*, install: bool, build: bool, clean: bool = False,
                     source_dir: Path | None = None) -> str:
    compat = (_MODELS / "mare2dem" / "_gnu_compat.py").read_text(
        encoding="utf-8")
    src = ""
    if source_dir and (Path(source_dir) / "Makefile").is_file():
        src = to_wsl_path(source_dir) if _IS_WIN else str(source_dir)
    elif not source_dir:
        vend = SOLVERS["mare2dem"].default_source_dir()
        if (vend / "Makefile").is_file():
            src = to_wsl_path(vend) if _IS_WIN else str(vend)
    repl = {
        "@@INSTALL@@": "1" if install else "0",
        "@@BUILD@@": "1" if build else "0",
        "@@CLEAN@@": "1" if clean else "0",
        "@@SOURCE@@": src.replace("'", ""),
        "@@MM_URL@@": MICROMAMBA_URL.format(platform=(
            "linux-64" if _IS_WIN else _micromamba_platform())),
        "@@PKGS@@": " ".join(MARE2DEM_PKGS),
        "@@ARCHIVE_URL@@": MARE2DEM_ARCHIVE_URL,
        "@@COMPAT@@": compat,
    }
    text = _M2D_SCRIPT
    for k, v in repl.items():
        text = text.replace(k, v)
    return text


def _run_mare2dem_script(ctx: Context, *, install: bool, build: bool,
                         clean: bool = False,
                         source_dir: Path | None = None) -> str | None:
    script = _mare2dem_script(install=install, build=build, clean=clean,
                              source_dir=source_dir)
    argv = ["wsl", "-e", "bash", "-s"] if _IS_WIN else ["bash", "-s"]
    total = [0]
    done = [0]
    binary = [None]

    def on_line(line: str) -> None:
        if line.startswith("@@STAGE "):
            ctx.stage(line[8:])
            ctx.progress(None)
        elif line.startswith("@@TOTAL "):
            total[0] = int(line.split()[1] or 0)
        elif line.startswith("@@BINARY "):
            binary[0] = line[9:].strip()
        elif total[0] and " -c " in line:
            done[0] += 1
            ctx.progress(min(done[0] / total[0], 0.98))

    rc = _stream(argv, ctx, stdin_text=script, on_line=on_line)
    if rc != 0:
        raise RuntimeError(f"MARE2DEM {'build' if build else 'setup'} "
                           f"failed (exit {rc}); see the log above.")
    if build and not binary[0]:
        raise RuntimeError("Build finished but no MARE2DEM binary was "
                           "reported.")
    return binary[0]


def build_plan(key: str, source_dir: Path | None = None, *,
               clean: bool = False, auto_install: bool = False,
               result: dict | None = None) -> list[Step]:
    """Steps to build *key*; ``result["binary"]`` is set on success."""
    spec = SOLVERS[key]
    out = result if result is not None else {}
    steps: list[Step] = []
    if spec.engine == "fortran":
        if auto_install and fortran_toolchain() is None:
            steps += install_plan(key)
        src = Path(source_dir) if source_dir else spec.default_source_dir()

        def build(ctx: Context) -> None:
            b = _build_fortran(spec, src, clean, ctx)
            out["binary"] = b
            register_binary(key, b, source=str(src))

        steps.append(Step(f"Compile {spec.label}", build))
    else:
        def build_m2d(ctx: Context) -> None:
            b = _run_mare2dem_script(ctx, install=auto_install, build=True,
                                     clean=clean, source_dir=source_dir)
            binary = f"wsl:{b}" if _IS_WIN else b
            out["binary"] = binary
            register_binary(key, binary, runtime="wsl" if _IS_WIN else
                            "native")

        steps.append(Step("Build MARE2DEM", build_m2d))
    return steps


def run_steps(steps: list[Step], ctx: Context) -> None:
    for i, step in enumerate(steps, 1):
        if ctx.cancelled():
            raise Cancelled()
        ctx.log(f"━━ Step {i}/{len(steps)}: {step.label}")
        step.run(ctx)


__all__ = [
    "Cancelled", "Check", "Context", "FortranToolchain", "Readiness",
    "SOLVERS", "SolverSpec", "Step", "build_plan", "check",
    "expected_binary", "find_binary", "fortran_toolchain", "install_plan",
    "is_wsl_binary", "load_registry", "managed_fortran_prefix",
    "micromamba_path", "probe_mare2dem_env", "register_binary",
    "registered_binary", "registry_path", "run_steps", "to_wsl_path",
    "toolchain_root", "wsl_available",
]
