"""Tests for pycsamt.models.solver_build (Qt-free Solver Builder engine).

Real compilation is verified separately (see the module docstring); these
tests cover detection, readiness, the registry, script generation and the
run/cancel plumbing without installing or compiling anything.  Every test
redirects the registry and toolchain locations to ``tmp_path``.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

import pycsamt.models.solver_build as sb


@pytest.fixture(autouse=True)
def _isolated(monkeypatch, tmp_path):
    monkeypatch.setattr(sb, "registry_path", lambda: tmp_path / "solvers.json")
    monkeypatch.setattr(sb, "managed_fortran_prefix",
                        lambda: tmp_path / "toolchain" / "fortran")
    yield


# ── catalogue & registry ──────────────────────────────────────────────────


def test_catalogue_has_four_solvers_with_bundled_fortran_sources():
    assert list(sb.SOLVERS) == ["occam2d", "modem2d", "modem3d", "mare2dem"]
    for key in ("occam2d", "modem2d", "modem3d"):
        spec = sb.SOLVERS[key]
        assert spec.vendored and spec.engine == "fortran"
        assert (spec.default_source_dir() / "Makefile").is_file(), key
    assert not sb.SOLVERS["mare2dem"].vendored


def test_registry_round_trip_and_stale_entries(tmp_path):
    exe = tmp_path / "Mod2DMT.exe"
    exe.write_bytes(b"x")
    sb.register_binary("modem2d", str(exe), source="bundled")
    assert sb.registered_binary("modem2d") == str(exe)
    assert sb.load_registry()["modem2d"]["source"] == "bundled"
    exe.unlink()
    assert sb.registered_binary("modem2d") is None  # deleted -> not offered
    sb.register_binary("mare2dem", "wsl:/root/x/MARE2DEM")
    assert sb.registered_binary("mare2dem") == "wsl:/root/x/MARE2DEM"
    assert sb.is_wsl_binary("wsl:/a") and not sb.is_wsl_binary("C:/a.exe")


def test_find_binary_prefers_registry_then_source_tree(tmp_path):
    src = tmp_path / "src"
    src.mkdir()
    built = sb.expected_binary("occam2d", src)
    built.write_bytes(b"x")
    assert sb.find_binary("occam2d", src) == str(built)
    other = tmp_path / "custom.exe"
    other.write_bytes(b"x")
    sb.register_binary("occam2d", str(other))
    assert sb.find_binary("occam2d", src) == str(other)


def test_to_wsl_path():
    assert sb.to_wsl_path(r"C:\Users\me\run") == "/mnt/c/Users/me/run"
    assert sb.to_wsl_path("/home/me") == "/home/me"


# ── toolchain detection & readiness ───────────────────────────────────────


def _fake_managed_toolchain(prefix: Path) -> None:
    if sys.platform == "win32":
        mingw = prefix / "Library" / "mingw-w64"
        (mingw / "bin").mkdir(parents=True)
        (mingw / "lib").mkdir(parents=True)
        (mingw / "bin" / "gfortran.exe").write_bytes(b"")
        (prefix / "Library" / "bin").mkdir(parents=True, exist_ok=True)
        (prefix / "Library" / "bin" / "make.exe").write_bytes(b"")
        (mingw / "lib" / "libopenblas.a").write_bytes(b"")
    else:
        (prefix / "bin").mkdir(parents=True)
        (prefix / "bin" / "gfortran").write_bytes(b"")
        (prefix / "bin" / "make").write_bytes(b"")


def test_managed_toolchain_is_detected_first(tmp_path):
    prefix = sb.managed_fortran_prefix()
    _fake_managed_toolchain(prefix)
    tc = sb.fortran_toolchain()
    assert tc is not None and tc.origin == "managed"
    assert "openblas" in tc.libs
    r = sb.check("modem2d")
    assert r.ready and not r.missing


def test_missing_toolchain_is_auto_installable(monkeypatch):
    monkeypatch.setattr(sb, "fortran_toolchain", lambda: None)
    r = sb.check("occam2d")
    labels = {c.key: c for c in r.checks}
    assert not r.ready and r.can_auto_install
    assert not labels["fortran"].ok and labels["fortran"].hint
    assert labels["source"].ok  # the bundled source is present


def test_custom_source_without_makefile_is_not_auto_fixable(tmp_path,
                                                             monkeypatch):
    monkeypatch.setattr(sb, "fortran_toolchain", lambda: None)
    r = sb.check("modem3d", tmp_path / "empty")
    assert not r.can_auto_install
    assert any(c.key == "source" and not c.ok and not c.auto
               for c in r.checks)


def test_install_plan_uses_micromamba_env(monkeypatch):
    calls = []
    monkeypatch.setattr(sb, "_create_env",
                        lambda prefix, pkgs, ctx: calls.append((prefix, pkgs)))
    sb.run_steps(sb.install_plan("modem2d"), sb.Context(log=lambda s: None))
    prefix, pkgs = calls[0]
    assert prefix == sb.managed_fortran_prefix()
    assert "make" in pkgs and any("fortran" in p for p in pkgs)


# ── MARE2DEM build script ─────────────────────────────────────────────────


def test_mare2dem_script_is_complete_and_valid_bash():
    script = sb._mare2dem_script(install=True, build=True, clean=True)
    leftovers = [t for t in ("@@INSTALL@@", "@@BUILD@@", "@@CLEAN@@",
                             "@@SOURCE@@", "@@PKGS@@", "@@COMPAT@@",
                             "@@MM_URL@@", "@@ARCHIVE_URL@@")
                 if t in script]
    assert not leftovers
    assert "def patch_source_tree" in script  # fixes embedded, stdlib-only
    assert "mkl-devel" in script and "openmpi" in script
    assert "-fallow-argument-mismatch" in script
    bash = shutil.which("bash")
    if bash:
        # bash -n: syntax only, nothing runs
        out = subprocess.run([bash, "-n"], input=script.encode(),
                             capture_output=True)
        assert out.returncode == 0, out.stderr.decode()


def test_mare2dem_readiness_on_windows_without_wsl(monkeypatch):
    monkeypatch.setattr(sb, "_IS_WIN", True)
    monkeypatch.setattr(sb, "wsl_available", lambda: False)
    r = sb.check("mare2dem")
    (c,) = r.checks
    assert c.key == "wsl" and not c.ok and not c.auto
    assert "wsl --install" in c.hint


def test_mare2dem_readiness_from_probe(monkeypatch):
    monkeypatch.setattr(sb, "_IS_WIN", False)
    monkeypatch.setattr(sb, "probe_mare2dem_env", lambda: {
        "mpifort": "/x/mpifort", "mpicc": "/x/mpicc", "make": "/x/make",
        "mkl": "", "python3": "/x/python3", "curl": "/x/curl"})
    r = sb.check("mare2dem")
    missing = [c.key for c in r.missing]
    assert missing == ["mkl"] and r.can_auto_install


# ── run plumbing ──────────────────────────────────────────────────────────


def test_run_steps_logs_and_honours_cancel():
    log = []
    ran = []
    steps = [sb.Step("a", lambda ctx: ran.append("a")),
             sb.Step("b", lambda ctx: ran.append("b"))]
    sb.run_steps(steps, sb.Context(log=log.append))
    assert ran == ["a", "b"] and any("Step 2/2" in m for m in log)
    with pytest.raises(sb.Cancelled):
        sb.run_steps(steps, sb.Context(log=log.append,
                                       cancelled=lambda: True))


def test_stream_sends_lf_only_stdin():
    """Regression: text-mode stdin on Windows turned LF into CRLF, which
    broke every line of the script bash received inside WSL."""
    py = sys.executable
    seen = []
    rc = sb._stream(
        [py, "-c", "import sys; d=sys.stdin.buffer.read(); "
                   "print('CR' if b'\\r' in d else 'LF-only')"],
        sb.Context(log=lambda s: None), stdin_text="a\nb\n",
        on_line=seen.append)
    assert rc == 0 and seen == ["LF-only"]


def test_build_without_toolchain_explains(monkeypatch):
    monkeypatch.setattr(sb, "fortran_toolchain", lambda: None)
    steps = sb.build_plan("modem2d")
    with pytest.raises(RuntimeError, match="install it first"):
        sb.run_steps(steps, sb.Context(log=lambda s: None))


# ── discovery & staging (Inversion Studio binary row) ──────────────────────


def test_discover_binaries_orders_registry_source_path(monkeypatch, tmp_path):
    reg = tmp_path / "reg.json"
    monkeypatch.setattr(sb, "registry_path", lambda: reg)
    built = tmp_path / "built" / "Occam2D.exe"
    built.parent.mkdir()
    built.write_bytes(b"x")
    sb.register_binary("occam2d", str(built))
    src = tmp_path / "src"
    src.mkdir()
    in_tree = sb.expected_binary("occam2d", src)
    in_tree.write_bytes(b"x")
    monkeypatch.setattr(sb.shutil, "which", lambda name: str(built))
    found = sb.discover_binaries("occam2d", src, probe_wsl=False)
    assert [c.origin for c in found] == ["Solver Builder", "source tree"]
    assert found[0].path == str(built)  # PATH duplicate not repeated


def test_discover_mare2dem_probes_wsl_build(monkeypatch, tmp_path):
    monkeypatch.setattr(sb, "registry_path", lambda: tmp_path / "r.json")
    monkeypatch.setattr(sb, "_IS_WIN", True)
    monkeypatch.setattr(sb.shutil, "which",
                        lambda name: "wsl.exe" if name == "wsl" else None)
    monkeypatch.setattr(sb, "_wsl_mare2dem_build",
                        lambda: "wsl:/root/.local/share/pycsamt/mare2dem/"
                                "build/MARE2DEM")
    found = sb.discover_binaries("mare2dem", tmp_path)
    assert found and found[0].origin == "WSL build"
    assert found[0].path.startswith("wsl:")
    assert sb.discover_binaries("mare2dem", tmp_path, probe_wsl=False) == []


def test_stage_binary_copies_runtime_dlls(monkeypatch, tmp_path):
    monkeypatch.setattr(sb, "_IS_WIN", True)
    src = tmp_path / "bin"
    src.mkdir()
    exe = src / "Mod2DMT.exe"
    exe.write_bytes(b"exe")
    (src / "libopenblas.dll").write_bytes(b"dll")
    out = sb.stage_binary(str(exe), tmp_path / "run")
    assert Path(out) == tmp_path / "run" / "Mod2DMT.exe" and exe.is_file()
    assert (tmp_path / "run" / "libopenblas.dll").is_file()
    assert sb.stage_binary(out, tmp_path / "run") == out  # already there


def test_stage_binary_move_updates_registry(monkeypatch, tmp_path):
    reg = tmp_path / "reg.json"
    monkeypatch.setattr(sb, "registry_path", lambda: reg)
    exe = tmp_path / "Occam2D.exe"
    exe.write_bytes(b"exe")
    sb.register_binary("occam2d", str(exe))
    out = sb.stage_binary(str(exe), tmp_path / "run", move=True)
    assert not exe.exists() and Path(out).is_file()
    assert sb.registered_binary("occam2d") == out


def test_stage_binary_rejects_wsl_and_missing(tmp_path):
    with pytest.raises(ValueError, match="WSL"):
        sb.stage_binary("wsl:/opt/MARE2DEM", tmp_path)
    with pytest.raises(FileNotFoundError):
        sb.stage_binary(str(tmp_path / "nope.exe"), tmp_path / "run")
