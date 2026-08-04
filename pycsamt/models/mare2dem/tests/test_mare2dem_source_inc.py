# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for :mod:`pycsamt.models.mare2dem.source`'s compiler-detection and
Make-include generation (``_detect_fc``/``_detect_cc``/``_generate_inc``).

These were rewritten 2026-07-31 after a real WSL2 + Intel oneAPI build
surfaced multiple genuine bugs in the previous toolchain choices (classic
``mpiifort``/``mpiicc``/``xiar`` either missing or broken on current oneAPI
releases; missing ``-std=gnu89``/``-cxxlib`` flags; invalid ``TRICOPTS``
syntax for ``icx``) -- see ``pycsamt/forward/maxwell/mare2dem.py``'s module
docstring and ``docs/source/user_guide/models/compilation.rst`` for the
full build story these fixes came from. No compiled MARE2DEM binary is
required to run these: they only check the generated include-file text.
"""

from __future__ import annotations

from pycsamt.models.mare2dem.source import _detect_cc, _detect_fc, _generate_inc


def test_detect_fc_prefers_mpiifx_over_classic_mpiifort(monkeypatch):
    available = {"mpiifx", "mpiifort", "mpifort"}
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.shutil.which",
        lambda name: name if name in available else None,
    )
    assert _detect_fc() == "mpiifx"


def test_detect_fc_falls_back_to_classic_mpiifort_when_mpiifx_absent(
    monkeypatch,
):
    available = {"mpiifort", "mpifort"}
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.shutil.which",
        lambda name: name if name in available else None,
    )
    assert _detect_fc() == "mpiifort"


def test_detect_cc_prefers_mpiicx_over_classic_mpiicc(monkeypatch):
    available = {"mpiicx", "mpiicc"}
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source.shutil.which",
        lambda name: name if name in available else None,
    )
    assert _detect_cc() == "mpiicx"


def test_generate_inc_for_modern_ifx_toolchain(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source._is_intel_compiler", lambda fc: True
    )
    inc = _generate_inc("mpiifx", "mpiicx", "/opt/intel/oneapi/mkl/latest", tmp_path)
    text = inc.read_text()
    assert "FC        = mpiifx" in text
    assert "-cxxlib" in text
    assert "-fpp" in text
    # BLACS/ScaLAPACK's legacy implicit-declaration C sources need this,
    # confirmed by a real build failing without it.
    assert "-std=gnu89" in text
    # xiar no longer exists in current oneAPI releases; always plain ar.
    assert "ARCH      = ar" in text
    assert "xiar" not in text
    # invalid icx syntax, dropped rather than translated.
    assert "-fp-model" not in text
    assert "-Wl,-rpath,/opt/intel/oneapi/mkl/latest/lib" in text


def test_generate_inc_for_classic_ifort_omits_cxxlib(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source._is_intel_compiler", lambda fc: True
    )
    inc = _generate_inc("mpiifort", "mpiicc", None, tmp_path)
    text = inc.read_text()
    assert "FC        = mpiifort" in text
    # untested/unconfirmed for the classic compiler -- not added.
    assert "-cxxlib" not in text
    assert "-std=gnu89" in text
    assert "ARCH      = ar" in text


def test_generate_inc_for_non_intel_toolchain(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "pycsamt.models.mare2dem.source._is_intel_compiler", lambda fc: False
    )
    inc = _generate_inc("mpifort", "mpicc", None, tmp_path)
    text = inc.read_text()
    assert "FC        = mpifort" in text
    assert "-cpp" in text
    assert "-fpp" not in text
    assert "ARCH      = ar" in text
