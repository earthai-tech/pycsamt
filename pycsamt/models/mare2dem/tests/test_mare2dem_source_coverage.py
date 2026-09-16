# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.models.mare2dem.source`.

Everything here is self-contained: no network access, no compiler
toolchain, no bundled ``data/mare2dem/`` example data. All subprocess,
network, and filesystem-permission behaviour is monkeypatched or driven
through real temporary directories, so these contribute real coverage
in CI (unlike anything that would depend on ``data/mare2dem/``, which is
gitignored and absent there).
"""

from __future__ import annotations

import subprocess
import tarfile
from pathlib import Path

import pytest

from pycsamt.models.mare2dem.config import Mare2DEMConfig
from pycsamt.models.mare2dem.source import (
    SourceManager,
    _detect_mkl,
    _is_intel_compiler,
    _user_data_dir,
)

pytestmark = pytest.mark.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# _user_data_dir / _detect_mkl / _is_intel_compiler
# ---------------------------------------------------------------------------


def test_user_data_dir_windows(monkeypatch):
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.platform.system", lambda: "Windows"
    )
    monkeypatch.setenv("LOCALAPPDATA", r"C:\fake\local")
    d = _user_data_dir()
    assert str(d).endswith("pycsamt\\mare2dem") or str(d).endswith(
        "pycsamt/mare2dem"
    )
    assert "fake" in str(d)


def test_user_data_dir_darwin(monkeypatch):
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.platform.system", lambda: "Darwin"
    )
    d = _user_data_dir()
    assert "Application Support" in str(d)


def test_user_data_dir_linux(monkeypatch):
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.platform.system", lambda: "Linux"
    )
    monkeypatch.delenv("XDG_DATA_HOME", raising=False)
    d = _user_data_dir()
    assert ".local" in str(d) and "share" in str(d)


def test_detect_mkl_from_env(tmp_path, monkeypatch):
    monkeypatch.setenv("MKLROOT", str(tmp_path))
    assert _detect_mkl() == str(tmp_path)


def test_detect_mkl_env_set_but_not_dir_falls_back(monkeypatch, tmp_path):
    missing = tmp_path / "nope"
    monkeypatch.setenv("MKLROOT", str(missing))
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.Path.is_dir",
        lambda self: self.as_posix() == "/opt/intel/oneapi/mkl/latest",
    )
    assert _detect_mkl() == "/opt/intel/oneapi/mkl/latest"


def test_detect_mkl_none_found(monkeypatch):
    monkeypatch.delenv("MKLROOT", raising=False)
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.Path.is_dir", lambda self: False
    )
    assert _detect_mkl() is None


def test_is_intel_compiler_true(monkeypatch):
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.subprocess.check_output",
        lambda *a, **k: "Intel(R) Fortran Compiler",
    )
    assert _is_intel_compiler("mpiifx") is True


def test_is_intel_compiler_false_on_gfortran(monkeypatch):
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.subprocess.check_output",
        lambda *a, **k: "GNU Fortran (GCC) 13.0.0",
    )
    assert _is_intel_compiler("mpifort") is False


def test_is_intel_compiler_exception_returns_false(monkeypatch):
    def _raise(*a, **k):
        raise FileNotFoundError()

    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.subprocess.check_output", _raise
    )
    assert _is_intel_compiler("nope") is False


# ---------------------------------------------------------------------------
# resolve_source_dir
# ---------------------------------------------------------------------------


def test_resolve_source_dir_explicit_arg(tmp_path):
    d = tmp_path / "explicit"
    sm = SourceManager(source_dir=d)
    resolved = sm.resolve_source_dir()
    assert resolved == d
    assert d.is_dir()


def test_resolve_source_dir_config_field(tmp_path):
    d = tmp_path / "cfg_dir"
    cfg = Mare2DEMConfig(source_dir=str(d))
    sm = SourceManager(config=cfg)
    resolved = sm.resolve_source_dir()
    assert resolved == Path(d)
    assert d.is_dir()


def test_resolve_source_dir_env_var(tmp_path, monkeypatch):
    d = tmp_path / "env_dir"
    monkeypatch.setenv("PYCSAMT_MARE2DEM_SOURCE", str(d))
    sm = SourceManager()
    resolved = sm.resolve_source_dir()
    assert resolved == Path(d)
    assert d.is_dir()


def test_resolve_source_dir_bundled_writable(tmp_path, monkeypatch):
    monkeypatch.delenv("PYCSAMT_MARE2DEM_SOURCE", raising=False)
    bundled = tmp_path / "_source"
    bundled.mkdir()
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.Path.__file__"
        if False
        else "pycsamt.models.mare2dem.source.__file__",
        str(tmp_path / "source.py"),
    )
    sm = SourceManager()
    resolved = sm.resolve_source_dir()
    assert resolved == bundled


def test_resolve_source_dir_bundled_unwritable_falls_back(
    tmp_path, monkeypatch
):
    monkeypatch.delenv("PYCSAMT_MARE2DEM_SOURCE", raising=False)
    bundled = tmp_path / "_source"
    bundled.mkdir()
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.__file__", str(tmp_path / "source.py")
    )

    orig_touch = Path.touch

    def _fail_touch(self, *a, **k):
        if self.name == ".write_test":
            raise OSError("read-only filesystem")
        return orig_touch(self, *a, **k)

    monkeypatch.setattr("pycsamt.models.mare2dem.source.Path.touch", _fail_touch)

    fallback = tmp_path / "userdata"
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source._user_data_dir", lambda: fallback
    )
    sm = SourceManager()
    resolved = sm.resolve_source_dir()
    assert resolved == fallback
    assert fallback.is_dir()


def test_resolve_source_dir_platform_fallback_when_no_bundled(
    tmp_path, monkeypatch
):
    monkeypatch.delenv("PYCSAMT_MARE2DEM_SOURCE", raising=False)
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.__file__",
        str(tmp_path / "nonexistent_pkg" / "source.py"),
    )
    fallback = tmp_path / "userdata2"
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source._user_data_dir", lambda: fallback
    )
    sm = SourceManager()
    resolved = sm.resolve_source_dir()
    assert resolved == fallback


# ---------------------------------------------------------------------------
# resolve_binary
# ---------------------------------------------------------------------------


def test_resolve_binary_found_on_path(tmp_path, monkeypatch):
    # resolve_binary() returns the bare name (relying on PATH resolution
    # at subprocess-launch time), matching runner.py's own _resolve_binary
    # helper -- not the absolute path which() found it at.
    sm = SourceManager(source_dir=tmp_path)
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.shutil.which",
        lambda name: "/usr/bin/MARE2DEM" if name == "MARE2DEM" else None,
    )
    assert sm.resolve_binary() == Path("MARE2DEM")


def test_resolve_binary_found_in_source_dir(tmp_path, monkeypatch):
    sm = SourceManager(source_dir=tmp_path)
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.shutil.which", lambda name: None
    )
    binfile = tmp_path / "MARE2DEM"
    binfile.write_text("#!/bin/sh\n")
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.os.access", lambda p, m: True
    )
    assert sm.resolve_binary() == binfile


def test_resolve_binary_not_found(tmp_path, monkeypatch):
    sm = SourceManager(source_dir=tmp_path)
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.shutil.which", lambda name: None
    )
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source._user_data_dir", lambda: tmp_path / "ud"
    )
    assert sm.resolve_binary() is None


# ---------------------------------------------------------------------------
# resolve_triangle_binary
# ---------------------------------------------------------------------------


def test_resolve_triangle_binary_found_on_path(tmp_path, monkeypatch):
    sm = SourceManager(source_dir=tmp_path)
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.shutil.which",
        lambda name, path=None: "/usr/bin/triangle"
        if name == "triangle" and path is None
        else None,
    )
    assert sm.resolve_triangle_binary() == Path("/usr/bin/triangle")


def test_resolve_triangle_binary_found_in_dir_via_which(tmp_path, monkeypatch):
    sm = SourceManager(source_dir=tmp_path)

    def _which(name, path=None):
        if path is not None and name == "triangle":
            return str(Path(path) / "triangle")
        return None

    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.shutil.which", _which
    )
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source._user_data_dir", lambda: tmp_path / "ud"
    )
    result = sm.resolve_triangle_binary()
    assert result == Path(str(tmp_path / "triangle"))


def test_resolve_triangle_binary_found_as_literal_file(tmp_path, monkeypatch):
    sm = SourceManager(source_dir=tmp_path)
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.shutil.which", lambda name, path=None: None
    )
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source._user_data_dir", lambda: tmp_path / "ud"
    )
    triangle_bin = tmp_path / "Triangle"
    triangle_bin.write_text("bin")
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.os.access", lambda p, m: True
    )
    result = sm.resolve_triangle_binary()
    assert result == triangle_bin


def test_resolve_triangle_binary_explicit_name_not_found(tmp_path, monkeypatch):
    sm = SourceManager(source_dir=tmp_path)
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.shutil.which", lambda name, path=None: None
    )
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source._user_data_dir", lambda: tmp_path / "ud"
    )
    assert sm.resolve_triangle_binary(name="mytriangle") is None


# ---------------------------------------------------------------------------
# is_downloaded / is_built
# ---------------------------------------------------------------------------


def test_is_downloaded_true_and_false(tmp_path):
    sm = SourceManager(source_dir=tmp_path)
    assert sm.is_downloaded() is False
    (tmp_path / "Makefile").write_text("all:\n")
    assert sm.is_downloaded() is True


def test_is_built_delegates_to_resolve_binary(tmp_path, monkeypatch):
    sm = SourceManager(source_dir=tmp_path)
    monkeypatch.setattr(sm, "resolve_binary", lambda: None)
    assert sm.is_built() is False
    monkeypatch.setattr(sm, "resolve_binary", lambda: tmp_path / "MARE2DEM")
    assert sm.is_built() is True


# ---------------------------------------------------------------------------
# download() / _download_git / _download_archive
# ---------------------------------------------------------------------------


def test_download_skips_when_already_downloaded(tmp_path, capsys):
    (tmp_path / "Makefile").write_text("all:\n")
    sm = SourceManager(source_dir=tmp_path, verbose=1)
    result = sm.download()
    assert result == tmp_path


def test_download_invalid_method_raises(tmp_path):
    sm = SourceManager(source_dir=tmp_path)
    with pytest.raises(ValueError, match="Unknown download method"):
        sm.download(method="ftp")


def test_download_auto_prefers_git_when_available(tmp_path, monkeypatch):
    sm = SourceManager(source_dir=tmp_path)
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.shutil.which",
        lambda name: "/usr/bin/git" if name == "git" else None,
    )
    called = {}

    def _fake_git(self, dest):
        called["git"] = dest

    monkeypatch.setattr(SourceManager, "_download_git", _fake_git)
    sm.download(method="auto")
    assert called["git"] == tmp_path


def test_download_auto_falls_back_to_archive_without_git(tmp_path, monkeypatch):
    sm = SourceManager(source_dir=tmp_path)
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.shutil.which", lambda name: None
    )
    called = {}
    monkeypatch.setattr(
        SourceManager, "_download_archive", lambda self, dest: called.setdefault("a", dest)
    )
    sm.download(method="auto")
    assert called["a"] == tmp_path


def test_download_git_method_explicit(tmp_path, monkeypatch):
    sm = SourceManager(source_dir=tmp_path)
    called = {}
    monkeypatch.setattr(
        SourceManager, "_download_git", lambda self, dest: called.setdefault("d", dest)
    )
    sm.download(method="GIT")
    assert called["d"] == tmp_path


def test_download_git_removes_existing_nonempty_dest(tmp_path, monkeypatch):
    sm = SourceManager(source_dir=tmp_path)
    (tmp_path / "stale.txt").write_text("old")

    def _fake_run(cmd, check=False):
        Path(cmd[-1]).mkdir(parents=True, exist_ok=True)
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.subprocess.run", _fake_run
    )
    sm._download_git(tmp_path)
    assert not (tmp_path / "stale.txt").exists()


def test_download_git_failure_raises_runtime_error(tmp_path, monkeypatch):
    sm = SourceManager(source_dir=tmp_path / "dest")
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.subprocess.run",
        lambda cmd, check=False: subprocess.CompletedProcess(cmd, 1),
    )
    with pytest.raises(RuntimeError, match="git clone failed"):
        sm._download_git(tmp_path / "dest")


def test_download_archive_without_requests_raises(tmp_path, monkeypatch):
    sm = SourceManager(source_dir=tmp_path)
    import builtins

    real_import = builtins.__import__

    def _fake_import(name, *a, **k):
        if name == "requests":
            raise ImportError("no requests")
        return real_import(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", _fake_import)
    with pytest.raises(RuntimeError, match="requests"):
        sm._download_archive(tmp_path)


def test_download_archive_success(tmp_path, monkeypatch):
    sm = SourceManager(source_dir=tmp_path)

    src_dir = tmp_path / "srcbuild"
    src_dir.mkdir()
    (src_dir / "Makefile").write_text("all:\n")
    tar_path = tmp_path / "archive.tar.gz"
    with tarfile.open(tar_path, "w:gz") as tar:
        tar.add(src_dir, arcname="mare2dem-abcdef")

    class _FakeResp:
        headers = {"content-length": str(tar_path.stat().st_size)}

        def raise_for_status(self):
            pass

        def iter_content(self, chunk_size=None):
            with open(tar_path, "rb") as fh:
                while True:
                    chunk = fh.read(chunk_size)
                    if not chunk:
                        break
                    yield chunk

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    class _FakeRequests:
        @staticmethod
        def get(url, stream=True, timeout=120):
            return _FakeResp()

    import sys as _sys

    monkeypatch.setitem(_sys.modules, "requests", _FakeRequests())

    dest = tmp_path / "extracted"
    sm._download_archive(dest)
    assert (dest / "Makefile").exists()


# ---------------------------------------------------------------------------
# build()
# ---------------------------------------------------------------------------


def test_build_raises_on_windows(tmp_path, monkeypatch):
    monkeypatch.setattr("pycsamt.models.mare2dem.source.sys.platform", "win32")
    sm = SourceManager(source_dir=tmp_path)
    with pytest.raises(RuntimeError, match="Windows"):
        sm.build()


def test_build_raises_when_not_downloaded(tmp_path, monkeypatch):
    monkeypatch.setattr("pycsamt.models.mare2dem.source.sys.platform", "linux")
    sm = SourceManager(source_dir=tmp_path)
    with pytest.raises(FileNotFoundError, match="Call SourceManager.download"):
        sm.build()


def test_build_success_with_generated_inc(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr("pycsamt.models.mare2dem.source.sys.platform", "linux")
    (tmp_path / "Makefile").write_text("all:\n")
    sm = SourceManager(source_dir=tmp_path, verbose=1)

    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source._detect_fc", lambda: "mpifort"
    )
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source._detect_cc", lambda: "mpicc"
    )
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source._detect_mkl", lambda: None
    )
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source._is_intel_compiler", lambda fc: False
    )

    def _fake_run(cmd, cwd=None, check=False):
        if cmd[:2] == ["make", "clean_all"]:
            return subprocess.CompletedProcess(cmd, 0)
        (Path(cwd) / "MARE2DEM").write_text("bin")
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.subprocess.run", _fake_run
    )

    binary = sm.build(clean_first=True)
    assert binary == tmp_path / "MARE2DEM"
    out = capsys.readouterr().out
    assert "WARNING: Intel MKL not found" in out


def test_build_with_explicit_inc_file_and_compilers(tmp_path, monkeypatch):
    monkeypatch.setattr("pycsamt.models.mare2dem.source.sys.platform", "linux")
    (tmp_path / "Makefile").write_text("all:\n")
    sm = SourceManager(source_dir=tmp_path)
    inc = tmp_path / "custom.inc"
    inc.write_text("FC=mpiifx\n")

    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source._detect_mkl",
        lambda: str(tmp_path),
    )

    def _fake_run(cmd, cwd=None, check=False):
        (Path(cwd) / "MARE2DEM").write_text("bin")
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.subprocess.run", _fake_run
    )
    binary = sm.build(inc_file=inc, fc="mpiifx", cc="mpiicx")
    assert binary.exists()


def test_build_failure_raises_runtime_error(tmp_path, monkeypatch):
    monkeypatch.setattr("pycsamt.models.mare2dem.source.sys.platform", "linux")
    (tmp_path / "Makefile").write_text("all:\n")
    sm = SourceManager(source_dir=tmp_path)
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source._detect_fc", lambda: "mpifort"
    )
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source._detect_cc", lambda: "mpicc"
    )
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source._detect_mkl", lambda: str(tmp_path)
    )
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.subprocess.run",
        lambda cmd, cwd=None, check=False: subprocess.CompletedProcess(cmd, 1),
    )
    with pytest.raises(RuntimeError, match="build failed"):
        sm.build()


def test_build_binary_missing_after_success_raises(tmp_path, monkeypatch):
    monkeypatch.setattr("pycsamt.models.mare2dem.source.sys.platform", "linux")
    (tmp_path / "Makefile").write_text("all:\n")
    sm = SourceManager(source_dir=tmp_path)
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source._detect_fc", lambda: "mpifort"
    )
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source._detect_cc", lambda: "mpicc"
    )
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source._detect_mkl", lambda: str(tmp_path)
    )
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.subprocess.run",
        lambda cmd, cwd=None, check=False: subprocess.CompletedProcess(cmd, 0),
    )
    with pytest.raises(RuntimeError, match="was not"):
        sm.build()


# ---------------------------------------------------------------------------
# status() / print_status()
# ---------------------------------------------------------------------------


def test_status_and_print_status(tmp_path, monkeypatch, capsys):
    sm = SourceManager(source_dir=tmp_path)
    monkeypatch.setattr(sm, "resolve_binary", lambda: None)
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source._detect_fc", lambda: "mpifort"
    )
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source._detect_cc", lambda: "mpicc"
    )
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source._detect_mkl", lambda: None
    )
    s = sm.status()
    assert s["source_dir"] == tmp_path
    assert s["downloaded"] is False
    assert s["built"] is False
    assert s["binary_path"] is None
    assert s["fc"] == "mpifort"
    assert s["mklroot"] is None

    sm.print_status()
    out = capsys.readouterr().out
    assert "MARE2DEM SourceManager status" in out
    assert "not found — required" in out
