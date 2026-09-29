# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for WorkflowController — the Qt-free model behind Pipeline Studio.

Runs the library pipeline engine for real on the small WILLY L18 survey;
history/export paths always point at pytest's tmp_path.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import pytest

from pycsamt.app.desktop.controllers.workflow_controller import (
    STATUS_STYLE,
    RunStatus,
    WorkflowController,
    parse_literal,
)

_ROOT = Path(__file__).parents[4]
_WILLY = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"


@pytest.fixture(scope="module")
def sites():
    if not (_WILLY.exists() and any(_WILLY.glob("*.edi"))):
        pytest.skip("WILLY L18PLT data not available")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(_WILLY))


@pytest.fixture
def ctrl():
    return WorkflowController()


def _run_all(c):
    for k, i in enumerate(c.plan("all"), 1):
        _, _ok, stop = c.execute_step(i, c.input_for(i), k)
        if stop:
            break


class TestCatalogue:
    def test_registry_and_presets_exposed(self, ctrl):
        assert len(ctrl.catalogue()) >= 55
        assert "static_shift" in ctrl.categories()
        assert "full_processing" in {p.name for p in ctrl.presets()}

    def test_status_colours_have_white_text_contrast(self):
        """Every pill colour must reach WCAG AA (4.5:1) against white."""

        def lum(hex_):
            rgb = [int(hex_[i:i + 2], 16) / 255 for i in (1, 3, 5)]
            lin = [c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
                   for c in rgb]
            return 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]

        for status, (_text, colour) in STATUS_STYLE.items():
            ratio = 1.05 / (lum(colour) + 0.05)
            assert ratio >= 4.5, (status, colour, round(ratio, 2))


class TestEditing:
    def test_load_preset_builds_library_steps(self, ctrl):
        ctrl.load_preset("basic_qc")
        assert [s.code for s in ctrl.steps][:2] == ["NR001", "FREQ002"]

    def test_add_move_remove_and_unique_labels(self, ctrl):
        i = ctrl.add_step("SS001")
        j = ctrl.add_step("SS001")
        assert ctrl.steps[j].label != ctrl.steps[i].label
        ctrl.add_step("NR001", index=0)
        assert ctrl.steps[0].code == "NR001"
        assert ctrl.move_step(0, +1) == 1
        ctrl.remove_step(1)
        assert [s.code for s in ctrl.steps] == ["SS001", "SS001"]

    def test_param_fields_split_defaults_and_signature(self, ctrl):
        ctrl.add_step("NR001")
        fields = ctrl.param_fields(0)
        basic = {f.name for f in fields if not f.advanced}
        assert {"mains_hz", "n_harm", "tol_hz"} <= basic
        assert all(f.name not in {"sites", "inplace", "verbose"}
                   for f in fields)

    def test_set_param_and_reset(self, ctrl):
        ctrl.add_step("NR001")
        ctrl.set_param(0, "mains_hz", 60)
        assert ctrl.steps[0].step.params["mains_hz"] == 60
        ctrl.reset_params(0)
        assert ctrl.steps[0].step.params["mains_hz"] == 50

    def test_parse_literal(self):
        assert parse_literal("(1e-3, 100.0)") == (0.001, 100.0)
        assert parse_literal("") is None
        assert parse_literal("tri") == "tri"

    def test_yaml_and_json_round_trip(self, ctrl, tmp_path):
        ctrl.load_preset("full_processing")
        ctrl.set_enabled(4, False)  # disabled steps are not saved
        for ext in ("yaml", "json"):
            path = tmp_path / f"wf.{ext}"
            ctrl.save(path)
            other = WorkflowController()
            other.load(path)
            assert [s.code for s in other.steps] == [
                s.code for s in ctrl.steps if s.enabled
            ]


class TestRunning:
    def test_run_all_and_outputs(self, ctrl, sites):
        ctrl.set_input(sites, "test")
        ctrl.load_preset("full_processing")
        ctrl.reset_run()
        _run_all(ctrl)
        assert ctrl.is_complete
        assert all(s.status is RunStatus.DONE for s in ctrl.steps)
        assert len(ctrl.output_sites) == len(sites)
        result = ctrl.build_result()
        assert len(result.step_results) == len(ctrl.steps)
        assert result.ok

    def test_disabled_step_is_skipped_and_chain_continues(self, ctrl, sites):
        ctrl.set_input(sites)
        ctrl.load_preset("full_processing")
        ctrl.reset_run()
        ctrl.set_enabled(2, False)
        _run_all(ctrl)
        assert ctrl.steps[2].status is RunStatus.DISABLED
        assert ctrl.steps[3].status is RunStatus.DONE

    def test_editing_marks_downstream_outdated(self, ctrl, sites):
        ctrl.set_input(sites)
        ctrl.load_preset("basic_qc")
        ctrl.reset_run()
        _run_all(ctrl)
        ctrl.reset_params(1)  # any edit of step 1 invalidates 1..end
        assert ctrl.steps[0].status is RunStatus.DONE
        assert all(s.status is RunStatus.OUTDATED for s in ctrl.steps[1:])
        assert ctrl.input_for(2) is None  # step 1 must be re-run first

    def test_error_policy_raise_stops(self, ctrl, sites, monkeypatch):
        from pycsamt.api.pipe import PYCSAMT_PIPE

        ctrl.set_input(sites)
        ctrl.add_step("NR001")

        def boom(_sites):
            raise RuntimeError("synthetic")

        monkeypatch.setattr(ctrl.steps[0].step, "transform", boom)
        monkeypatch.setattr(PYCSAMT_PIPE, "on_step_error", "raise")
        out, ok, stop = ctrl.execute_step(0, sites)
        assert (ok, stop) == (False, True)
        assert out is sites and "synthetic" in ctrl.steps[0].error
        monkeypatch.setattr(PYCSAMT_PIPE, "on_step_error", "warn")
        assert ctrl.execute_step(0, sites)[2] is False

    def test_export_and_history(self, ctrl, sites, tmp_path, monkeypatch):
        from pycsamt.api.pipe import PYCSAMT_PIPE

        monkeypatch.setattr(PYCSAMT_PIPE, "report_formats",
                            ("html", "txt", "dashboard"))
        ctrl.set_input(sites)
        ctrl.load_preset("basic_qc")
        ctrl.reset_run()
        _run_all(ctrl)
        written = ctrl.export(tmp_path / "out", figures=False)
        assert len(written["edis"]) == len(sites)
        assert {p.name for p in written["reports"]} == {
            "summary.txt", "report.html", "dashboard.html"}
        hist = tmp_path / "history.jsonl"
        ctrl.record_history(hist)
        rows = ctrl.load_history(hist)
        assert len(rows) == 1 and rows[0]["n_sites_out"] == len(sites)
