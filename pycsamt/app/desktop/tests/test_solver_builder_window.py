# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for SolverBuilderWindow (Tools ▸ Solver Builder).

The engine (:mod:`pycsamt.models.solver_build`) is faked: readiness comes
from a stub ``check`` and "building" runs scripted steps, so nothing is
installed or compiled.  ``SolverBuildWorker.start`` runs ``run()``
synchronously (no QThread under the offscreen platform); the registry is
redirected to ``tmp_path``.
"""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

import pycsamt.app.desktop.workers.solver_build_worker as wmod
import pycsamt.models.solver_build as sb
from pycsamt.app.desktop.windows.solver_builder_window import (
    SolverBuilderWindow,
)


def _ready(ok=True, auto=True):
    return sb.Readiness([
        sb.Check("fortran", "Fortran compiler (gfortran)", ok, "x", "hint",
                 auto=auto),
        sb.Check("source", "Source code", True, "src", auto=False),
    ])


@pytest.fixture
def env(monkeypatch, tmp_path):
    monkeypatch.setattr(wmod.SolverBuildWorker, "start",
                        lambda self: self.run())
    monkeypatch.setattr(sb, "registry_path", lambda: tmp_path / "reg.json")
    state = {"ready": _ready(), "built": {}}
    monkeypatch.setattr(sb, "check",
                        lambda key, source_dir=None: state["ready"])
    monkeypatch.setattr(sb, "find_binary",
                        lambda key, source_dir=None: state["built"].get(key))

    def build_plan(key, source_dir=None, *, clean=False, auto_install=False,
                   result=None):
        def step(ctx):
            ctx.stage("Compiling")
            ctx.progress(0.5)
            ctx.log("gfortran -c a.f90")
            if state.get("fail"):
                raise RuntimeError("make failed")
            binary = str(tmp_path / f"{key}.exe")
            result["binary"] = binary
            state["built"][key] = binary

        return [sb.Step("Compile", step)]

    monkeypatch.setattr(sb, "build_plan", build_plan)
    return state


@pytest.fixture
def win(qapp, env):
    w = SolverBuilderWindow()
    w.show()
    yield w
    w.close()


def test_lists_four_solvers_and_checks_the_first(win):
    assert win._list.count() == 4
    assert win._key == "occam2d"
    assert win._deps.topLevelItemCount() == 2
    assert "ready to build" in win._deps_summary.text()
    assert win._btn_build.isEnabled() and not win._btn_install.isEnabled()


def test_missing_auto_installable_dependency(win, env):
    env["ready"] = _ready(ok=False)
    win._recheck()
    assert win._deps.topLevelItem(0).text(0) == "Missing"
    assert "installed automatically" in win._deps_summary.text()
    assert win._btn_build.isEnabled()  # Build installs first (auto mode)
    assert win._btn_install.isEnabled()
    win._rb_manual.setChecked(True)  # "I will install them myself"
    assert not win._btn_build.isEnabled()
    assert win._manual.isVisible() and "hint" in win._manual.text()


def test_manual_only_dependency_blocks_build(win, env):
    env["ready"] = _ready(ok=False, auto=False)
    win._recheck()
    assert not win._btn_build.isEnabled()
    assert "manually" in win._deps_summary.text()


def test_successful_build_registers_and_signals(win, env):
    got = []
    win.binary_built.connect(lambda k, b: got.append((k, b)))
    win.open_solver("modem3d")
    win._btn_build.click()
    assert got and got[0][0] == "modem3d"
    assert win._stage.text() == "✓ Build complete"
    assert win._result.isVisible() and "modem3d.exe" in win._result_msg.text()
    assert win._rows["modem3d"].pill.text() == "Built"
    assert win._btn_build.isEnabled() and not win._btn_stop.isEnabled()
    assert "gfortran -c a.f90" in win._log.toPlainText()


def test_failed_build_shows_reason_and_log(win, env):
    env["fail"] = True
    win.open_solver("modem2d")
    win._btn_build.click()
    assert win._stage.text().startswith("✕") and "make failed" in \
        win._stage.text()
    assert win._rows["modem2d"].pill.text() == "Error"
    assert win._btn_log.isChecked() and win._log.isVisible()
    assert not win._result.isVisible()


def test_busy_state_disables_controls(win, monkeypatch):
    # A started-but-unfinished worker must lock the UI (regression: the
    # buttons were refreshed before the thread started, leaving Build
    # clickable and Stop disabled).
    monkeypatch.setattr(wmod.SolverBuildWorker, "start", lambda self: None)
    win._btn_build.click()
    assert not win._btn_build.isEnabled() and win._btn_stop.isEnabled()
    assert not win._list.isEnabled()


def test_custom_source_folder_is_passed(win, env, monkeypatch, tmp_path):
    seen = []
    monkeypatch.setattr(sb, "check",
                        lambda key, source_dir=None: (seen.append(source_dir),
                                                      env["ready"])[1])
    win._rb_custom.setChecked(True)
    win._src_edit.setText(str(tmp_path))
    win._recheck()
    assert seen[-1] == tmp_path


def test_open_solver_and_session_round_trip(win):
    win.open_solver("mare2dem")
    assert win._key == "mare2dem" and "MARE2DEM" in win._btn_build.text()
    store = {}
    win.save_geometry_to(store)
    assert "solver_builder" in store
    win.restore_geometry_from(store)
