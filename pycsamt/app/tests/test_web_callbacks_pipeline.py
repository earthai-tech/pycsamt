# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.web.callbacks.pipeline (Processing Pipeline page).

Real data
---------
data/AMT/WILLY_DATA/L18PLT/ -- a small (first 6 stations) real-EDI subset
is used to drive ``PipelineController.execute_step()`` for real through the
Dash callbacks, mirroring the existing desktop-side
``test_pipeline_controller.py`` convention. The one thing intentionally
*not* exercised for real is the Export step's ``et.export_edis()`` call --
it is already documented elsewhere (``test_pipeline_controller.py``,
"KNOWN BUG #3") as always raising ``TypeError: ... missing ... 'new_z'``;
here we only confirm that the *callback* surfaces that failure as an error
toast instead of crashing.

Bugs found while writing these tests (documented via test, not fixed)
-----------------------------------------------------------------------
1. ``_generate_preview()``'s diag-less fallback calls
   ``sites.plot_pseudosection(ax=ax, component="res")``, but ``Sites``
   has no such method (nor any ``plot_*`` method at all) -- the call
   always raises ``AttributeError``, silently swallowed by
   ``except Exception: plt.close("all")``, so the preview always falls
   through to the 1x1 transparent placeholder. This affects step 0
   (Load Data), the only pipeline step with ``diag_fn=None``.
2. ``_filter_by_active_lines()`` fails to detect the "every line muted"
   case when ``active`` is an explicit empty list. It calls
   ``filter_sites_by_lines(sites, records, list(_active))`` with an
   empty list, and ``filter_sites_by_lines`` treats an empty
   ``active_lines`` argument as "no filter requested" (its own
   ``if not active_lines: return sites`` short-circuit) rather than
   "filter down to nothing" -- so it hands back the *unfiltered* sites
   and ``all_muted`` comes back ``False``. The "No active lines" guard
   only works when ``active`` names lines that genuinely don't match
   any station (see ``test_genuinely_nonexistent_active_line_...``).
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest
from dash import no_update
from dash.exceptions import PreventUpdate

from pycsamt.app.desktop.controllers.pipeline_controller import StepStatus
from pycsamt.app.web.callbacks.pipeline import (
    _filter_by_active_lines,
    _generate_preview,
    _result_to_figure,
    _stepper_classes,
)
from pycsamt.app.web.layout import IDs
from pycsamt.app.web.utils import empty_src

_ROOT = Path(__file__).parents[3]
_WILLY_L18 = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
_HAS_WILLY = _WILLY_L18.exists() and any(_WILLY_L18.glob("*.edi"))


# ── Dash callback-lookup helpers (shared pattern across web-callback tests) ──


def _unwrap(entry):
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def _cb(web_app, output_id_prop):
    return _unwrap(web_app.callback_map[output_id_prop])


def _cb_multi(web_app, *substrings):
    key = next(
        k for k in web_app.callback_map if all(s in k for s in substrings)
    )
    return _unwrap(web_app.callback_map[key])


def _cb_by_input(web_app, output_substr, input_id):
    for k, v in web_app.callback_map.items():
        if output_substr not in k:
            continue
        if any(input_id in str(i.get("id")) for i in v.get("inputs", [])):
            return _unwrap(v)
    raise AssertionError(
        f"no callback found for output~={output_substr!r} input={input_id!r}"
    )


# ── Fixtures ─────────────────────────────────────────────────────────────────


@pytest.fixture(scope="session")
def willy_sites():
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_WILLY:
        pytest.skip("WILLY L18PLT data not available")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(_WILLY_L18))


@pytest.fixture(scope="session")
def willy_subset(willy_sites):
    """Small (6-station) real Sites subset -- keeps per-step execution fast."""
    from pycsamt.emtools import ensure_sites

    return ensure_sites(willy_sites.as_list()[:6])


@pytest.fixture
def willy_subset_records(willy_subset):
    return [
        {"ID": ed.station, "Line": "L1" if i % 2 == 0 else "L2"}
        for i, ed in enumerate(willy_subset.as_list())
    ]


@pytest.fixture
def cached_session(willy_subset):
    from pycsamt.app.web.cache import cache_set

    session_id = "test-pipeline-session"
    cache_set(session_id, willy_subset)
    return session_id


@pytest.fixture(autouse=True)
def _reset_pipeline_ctrl():
    """Isolate tests from each other: pipeline.py's ``_CTRL`` is a module-
    level singleton mutated by every run_step/run_all/skip/reset call."""
    import pycsamt.app.web.callbacks.pipeline as pipeline_mod
    from pycsamt.app.desktop.controllers.pipeline_controller import (
        PipelineController,
    )

    pipeline_mod._CTRL = PipelineController()
    yield


# ── 1. _stepper_classes ───────────────────────────────────────────────────────


class TestStepperClasses:
    def test_fresh_controller_first_step_active(self):
        classes = _stepper_classes()
        assert classes[0] == "pipe-step-node step-active"
        assert all(c == "pipe-step-node step-pending" for c in classes[1:])
        assert len(classes) == 8

    def test_done_step_advances_active_marker(self):
        import pycsamt.app.web.callbacks.pipeline as pipeline_mod

        pipeline_mod._CTRL.steps[0].status = StepStatus.DONE
        classes = _stepper_classes()
        assert classes[0] == "pipe-step-node step-done"
        assert classes[1] == "pipe-step-node step-active"
        assert classes[2] == "pipe-step-node step-pending"

    def test_error_step_blocks_active_marker(self):
        import pycsamt.app.web.callbacks.pipeline as pipeline_mod

        pipeline_mod._CTRL.steps[0].status = StepStatus.ERROR
        classes = _stepper_classes()
        assert classes[0] == "pipe-step-node step-failed"
        assert "step-active" not in " ".join(classes)

    def test_skipped_step_counts_as_done_for_active_marker(self):
        import pycsamt.app.web.callbacks.pipeline as pipeline_mod

        pipeline_mod._CTRL.steps[0].status = StepStatus.SKIPPED
        classes = _stepper_classes()
        assert classes[0] == "pipe-step-node step-skipped"
        assert classes[1] == "pipe-step-node step-active"


# ── 2. _result_to_figure ─────────────────────────────────────────────────────


class TestResultToFigure:
    def test_figure_passthrough(self):
        import matplotlib.pyplot as plt

        fig = plt.figure()
        try:
            assert _result_to_figure(fig) is fig
        finally:
            plt.close(fig)

    def test_axes_resolves_to_its_figure(self):
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        try:
            assert _result_to_figure(ax) is fig
        finally:
            plt.close(fig)

    def test_object_without_get_figure_returns_none(self):
        assert _result_to_figure(object()) is None

    def test_get_figure_raising_is_swallowed(self):
        class Bad:
            def get_figure(self):
                raise RuntimeError("boom")

        assert _result_to_figure(Bad()) is None

    def test_get_figure_returning_non_figure_returns_none(self):
        class Weird:
            def get_figure(self):
                return "not-a-figure"

        assert _result_to_figure(Weird()) is None


# ── 3. _generate_preview ─────────────────────────────────────────────────────


class TestGeneratePreview:
    def test_sites_none_returns_placeholder(self):
        import pycsamt.app.web.callbacks.pipeline as pipeline_mod

        step = pipeline_mod._CTRL.steps[1]
        assert _generate_preview(step, None, True) == empty_src()

    def test_step_without_diag_fn_hits_missing_pseudosection_bug(
        self, willy_subset
    ):
        """See bug #1 in the module docstring above."""
        import pycsamt.app.web.callbacks.pipeline as pipeline_mod

        step = pipeline_mod._CTRL.steps[0]  # Load Data: diag_fn is None
        assert step.diag_fn is None
        out = _generate_preview(step, willy_subset, True)
        assert out == empty_src()

    def test_step_with_real_diag_fn_renders_a_real_image(self, willy_subset):
        import pycsamt.app.web.callbacks.pipeline as pipeline_mod

        step = pipeline_mod._CTRL.steps[1]  # QC: plot_confidence_band_summary
        out = _generate_preview(step, willy_subset, False)
        assert out.startswith("data:image/png;base64,")
        assert len(out) > 2000  # much larger than the placeholder

    def test_unknown_diag_fn_name_falls_back_to_broken_pseudosection(
        self, willy_subset
    ):
        import pycsamt.app.web.callbacks.pipeline as pipeline_mod

        step = pipeline_mod._CTRL.steps[1]
        step.diag_fn = "does_not_exist_on_emtools"
        out = _generate_preview(step, willy_subset, True)
        assert out == empty_src()

    def test_diag_fn_raising_falls_back_to_broken_pseudosection(
        self, willy_subset, monkeypatch
    ):
        import pycsamt.emtools as et
        import pycsamt.app.web.callbacks.pipeline as pipeline_mod

        step = pipeline_mod._CTRL.steps[1]
        monkeypatch.setattr(
            et,
            step.diag_fn,
            lambda sites: (_ for _ in ()).throw(RuntimeError("x")),
        )
        out = _generate_preview(step, willy_subset, True)
        assert out == empty_src()

    def test_diag_fn_returning_non_figure_falls_back(
        self, willy_subset, monkeypatch
    ):
        import pycsamt.emtools as et
        import pycsamt.app.web.callbacks.pipeline as pipeline_mod

        step = pipeline_mod._CTRL.steps[1]
        monkeypatch.setattr(et, step.diag_fn, lambda sites: "not-a-figure")
        out = _generate_preview(step, willy_subset, True)
        assert out == empty_src()


# ── 4. _filter_by_active_lines ───────────────────────────────────────────────


class TestFilterByActiveLines:
    def test_no_active_lines_data_returns_sites_unchanged(self):
        sites = object()
        out_sites, all_muted = _filter_by_active_lines(sites, {}, None)
        assert out_sites is sites
        assert all_muted is False

    def test_active_equals_all_returns_sites_unchanged(self):
        sites = object()
        als = {"active": ["L1", "L2"], "all": ["L1", "L2"]}
        out_sites, all_muted = _filter_by_active_lines(sites, {}, als)
        assert out_sites is sites
        assert all_muted is False

    def test_empty_active_list_bug_does_not_report_all_muted(
        self, willy_subset, willy_subset_records
    ):
        """See bug #2 in the module docstring above."""
        store_data = {"station_records": willy_subset_records}
        als = {"active": [], "all": ["L1", "L2"]}
        out_sites, all_muted = _filter_by_active_lines(
            willy_subset, store_data, als
        )
        assert all_muted is False
        assert out_sites is willy_subset

    def test_real_subset_of_lines_filters_correctly(
        self, willy_subset, willy_subset_records
    ):
        store_data = {"station_records": willy_subset_records}
        als = {"active": ["L1"], "all": ["L1", "L2"]}
        out_sites, all_muted = _filter_by_active_lines(
            willy_subset, store_data, als
        )
        assert all_muted is False
        assert out_sites is not None
        assert len(out_sites.as_list()) < len(willy_subset.as_list())

    def test_genuinely_nonexistent_active_line_reports_all_muted(
        self, willy_subset, willy_subset_records
    ):
        store_data = {"station_records": willy_subset_records}
        als = {"active": ["GhostLine"], "all": ["GhostLine", "L1"]}
        out_sites, all_muted = _filter_by_active_lines(
            willy_subset, store_data, als
        )
        assert out_sites is None
        assert all_muted is True


# ── 5. sync_step_info callback ───────────────────────────────────────────────


class TestSyncStepInfo:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.PIPE_STEP_INFO}.children")

    def test_none_or_empty_step_id_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None)
        with pytest.raises(PreventUpdate):
            self._fn(web_app)("")

    def test_normal_step_hides_export_card(self, web_app):
        desc, opts, default_method, export_style = self._fn(web_app)("1")
        assert isinstance(desc, str) and desc
        assert {
            "label": "Drop low-confidence frequencies",
            "value": "drop_low_conf",
        } in opts
        assert default_method == "drop_low_conf"
        assert export_style == {"display": "none"}

    def test_export_step_shows_export_card(self, web_app):
        _, opts, default_method, export_style = self._fn(web_app)("7")
        assert default_method == "edi"
        assert export_style == {"display": "block"}
        assert opts == [{"label": "Export EDI files", "value": "edi"}]


# ── 6. run_step callback ─────────────────────────────────────────────────────


class TestRunStep:
    def _fn(self, web_app):
        return _cb_by_input(web_app, f"{IDs.IMG_PIPE}.src", IDs.BTN_PIPE_RUN)

    def test_no_clicks_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(
                None, "0", "current", "", "sid", None, None, None, "dark"
            )

    def test_all_muted_returns_warning(self, web_app, willy_subset_records):
        store_data = {"station_records": willy_subset_records}
        als = {"active": ["GhostLine"], "all": ["GhostLine", "L1", "L2"]}
        out = self._fn(web_app)(
            1,
            "1",
            "drop_low_conf",
            "",
            "no-such-session",
            None,
            als,
            store_data,
            "dark",
        )
        log, src, store, view_val, step_val, spinner, toast, body = out
        assert "No active lines" in log
        assert store == {}
        assert view_val is no_update
        assert step_val is no_update
        assert toast is False

    def test_step0_current_method_hits_missing_pseudosection_bug(
        self, web_app, willy_subset, cached_session
    ):
        out = self._fn(web_app)(
            1, "0", "current", "", cached_session, None, None, None, "dark"
        )
        log, src, store, view_val, step_val, spinner, toast, body = out
        assert "Using currently loaded survey data." in log
        assert toast is False
        # Bug #1 (see module docstring): step 0 has no diag_fn and the
        # pseudosection fallback always fails, so even a successfully
        # processed real survey only ever gets the placeholder preview.
        assert src == empty_src()
        assert store["0"] == src
        assert view_val == "0"
        assert step_val == "1"

    def test_step1_qc_produces_a_real_diagnostic_image(
        self, web_app, willy_subset, cached_session
    ):
        out = self._fn(web_app)(
            1,
            "1",
            "drop_low_conf",
            "",
            cached_session,
            None,
            None,
            None,
            "dark",
        )
        log, src, store, view_val, step_val, spinner, toast, body = out
        assert toast is False
        assert src.startswith("data:image/png;base64,")
        assert len(src) > 2000
        assert store["1"] == src
        assert view_val == "1"
        assert step_val == "2"

    def test_step7_export_folder_injected_and_hits_new_z_bug(
        self, web_app, willy_subset, cached_session, tmp_path
    ):
        # KNOWN BUG (documented in test_pipeline_controller.py): export_edis()
        # requires a positional `new_z` arg the controller never supplies, so
        # the Export step always raises TypeError. At the callback level this
        # must surface as an error toast, not crash the app.
        out_dir = tmp_path / "export_out"
        out = self._fn(web_app)(
            1,
            "7",
            "edi",
            str(out_dir),
            cached_session,
            None,
            None,
            None,
            "dark",
        )
        log, src, store, view_val, step_val, spinner, toast, body = out
        assert "ERROR" in log
        assert toast is True
        assert "new_z" in body
        assert src == empty_src()

    def test_invalid_step_id_hits_outer_exception_handler(self, web_app):
        out = self._fn(web_app)(1, "99", "x", "", None, None, None, None, "dark")
        log, src, store, view_val, step_val, spinner, toast, body = out
        assert toast is True
        assert "ERROR" in log
        assert src == empty_src()
        assert view_val is no_update
        assert step_val is no_update

    def test_cache_get_failure_hits_outer_exception_handler(
        self, web_app, monkeypatch
    ):
        import pycsamt.app.web.callbacks.pipeline as pipeline_mod

        monkeypatch.setattr(
            pipeline_mod,
            "cache_get",
            lambda sid: (_ for _ in ()).throw(RuntimeError("boom")),
        )
        out = self._fn(web_app)(
            1, "0", "current", "", "sid", None, None, None, "light"
        )
        *_, toast, body = out
        assert toast is True
        assert "boom" in body


# ── 7. run_all callback ──────────────────────────────────────────────────────


class TestRunAll:
    def _fn(self, web_app):
        return _cb_by_input(
            web_app, f"{IDs.PIPE_SPINNER}.children", IDs.BTN_PIPE_ALL
        )

    def test_no_clicks_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None, "", "sid", None, None, None, "dark")

    def test_all_muted_returns_warning(self, web_app, willy_subset_records):
        store_data = {"station_records": willy_subset_records}
        als = {"active": ["Ghost"], "all": ["Ghost", "L1", "L2"]}
        out = self._fn(web_app)(
            1, "", "no-session", None, als, store_data, "dark"
        )
        log, store, spinner, toast, body = out
        assert "No active lines" in log
        assert toast is False

    def test_runs_all_steps_and_stops_at_export_without_a_folder(
        self, web_app, willy_subset, cached_session, monkeypatch
    ):
        import pycsamt.app.web.callbacks.pipeline as pipeline_mod

        def execute_step(index, log_cb):
            if index == 7:
                raise ValueError("No output folder selected")
            pipeline_mod._CTRL._sites_chain[index] = willy_subset
            log_cb(f"step {index + 1}")

        monkeypatch.setattr(pipeline_mod._CTRL, "execute_step", execute_step)
        monkeypatch.setattr(
            pipeline_mod, "_generate_preview", lambda *_a, **_k: empty_src()
        )
        out = self._fn(web_app)(
            1, "", cached_session, None, None, None, "dark"
        )
        log, store, spinner, toast, body = out
        assert toast is False
        for i in range(1, 8):
            assert f"Step {i} done" in log
        assert "Step 8 ERROR" in log
        assert "No output folder selected" in log
        assert store["0"] == empty_src()
        assert set(store.keys()) == {"0", "1", "2", "3", "4", "5", "6"}

    def test_runs_all_steps_with_export_folder_hits_new_z_bug(
        self, web_app, willy_subset, cached_session, tmp_path, monkeypatch
    ):
        import pycsamt.app.web.callbacks.pipeline as pipeline_mod

        def execute_step(index, log_cb):
            if index == 7:
                raise TypeError("missing required keyword-only argument: new_z")
            pipeline_mod._CTRL._sites_chain[index] = willy_subset
            log_cb(f"step {index + 1}")

        monkeypatch.setattr(pipeline_mod._CTRL, "execute_step", execute_step)
        monkeypatch.setattr(
            pipeline_mod, "_generate_preview", lambda *_a, **_k: empty_src()
        )
        out_dir = tmp_path / "run_all_export"
        out = self._fn(web_app)(
            1, str(out_dir), cached_session, None, None, None, "dark"
        )
        log, store, spinner, toast, body = out
        assert toast is False  # per-step errors are swallowed, loop just stops
        assert "Step 8 ERROR" in log
        assert "new_z" in log

    def test_outer_exception_handler(self, web_app, monkeypatch):
        import pycsamt.app.web.callbacks.pipeline as pipeline_mod

        monkeypatch.setattr(
            pipeline_mod,
            "cache_get",
            lambda sid: (_ for _ in ()).throw(RuntimeError("kaboom")),
        )
        out = self._fn(web_app)(1, "", "sid", None, None, None, "dark")
        log, store, spinner, toast, body = out
        assert toast is True
        assert "kaboom" in body


# ── 8. view_step_preview callback ────────────────────────────────────────────


class TestViewStepPreview:
    def _fn(self, web_app):
        return _cb_by_input(
            web_app, f"{IDs.IMG_PIPE}.src", IDs.PIPE_VIEW_STEP
        )

    def test_none_or_empty_step_id_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None, {}, "dark")
        with pytest.raises(PreventUpdate):
            self._fn(web_app)("", {}, "dark")

    def test_missing_cached_step_returns_placeholder(self, web_app):
        out = self._fn(web_app)("3", {}, "dark")
        assert out == empty_src()

    def test_returns_cached_preview(self, web_app):
        out = self._fn(web_app)(
            "2", {"2": "data:image/png;base64,XYZ"}, "light"
        )
        assert out == "data:image/png;base64,XYZ"


# ── 9. skip_step callback ────────────────────────────────────────────────────


class TestSkipStep:
    def _fn(self, web_app):
        return _cb_by_input(web_app, f"{IDs.PIPE_LOG}.children", IDs.BTN_PIPE_SKIP)

    def test_no_clicks_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None, "2")

    def test_marks_step_skipped(self, web_app):
        import pycsamt.app.web.callbacks.pipeline as pipeline_mod

        out = self._fn(web_app)(1, "2")
        assert out == "Step 3 skipped."
        assert pipeline_mod._CTRL.steps[2].status == StepStatus.SKIPPED


# ── 10. reset_pipeline callback ──────────────────────────────────────────────


class TestResetPipeline:
    def _fn(self, web_app):
        return _cb_by_input(
            web_app, f"{IDs.PIPE_STATUS}.children", IDs.BTN_PIPE_RESET
        )

    def test_no_clicks_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None)

    def test_resets_state_and_node_classes(self, web_app):
        import pycsamt.app.web.callbacks.pipeline as pipeline_mod

        pipeline_mod._CTRL.steps[0].status = StepStatus.DONE
        out = self._fn(web_app)(1)
        log, status, store, *node_classes = out
        assert log == "Pipeline reset."
        assert status == ""
        assert store == {}
        assert len(node_classes) == 8
        assert node_classes[0] == "pipe-step-node step-active"
        assert all(c == "pipe-step-node step-pending" for c in node_classes[1:])
        assert pipeline_mod._CTRL.steps[0].status == StepStatus.PENDING


# ── 11. update_status callback ───────────────────────────────────────────────


class TestUpdateStatus:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "pipe-node-0.className", IDs.PIPE_LOG)

    def test_builds_rows_and_classes(self, web_app):
        import pycsamt.app.web.callbacks.pipeline as pipeline_mod

        pipeline_mod._CTRL.steps[0].status = StepStatus.DONE
        pipeline_mod._CTRL.steps[0].elapsed_s = 1.234
        out = self._fn(web_app)("whatever log text")
        rows, *node_classes = out
        assert len(rows) == 8
        assert len(node_classes) == 8
        assert rows[0].children[0].children == "1. Load Data"
        assert "1.23s" in rows[0].children[2].children


# ── 12. Folder-browser modal callbacks ───────────────────────────────────────


class TestToggleBrowseModal:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.PIPE_BROWSE_MODAL}.is_open")

    def _ctx(self, monkeypatch, triggered_id):
        import pycsamt.app.web.callbacks.pipeline as pipeline_mod

        monkeypatch.setattr(
            pipeline_mod, "ctx", type("C", (), {"triggered_id": triggered_id})()
        )

    def test_cancel_closes(self, web_app, monkeypatch):
        self._ctx(monkeypatch, "pipe-browse-cancel")
        is_open, path = self._fn(web_app)(None, 1, "", False)
        assert is_open is False
        assert path is no_update

    def test_open_with_valid_current_folder(self, web_app, monkeypatch, tmp_path):
        self._ctx(monkeypatch, IDs.BTN_PIPE_BROWSE)
        is_open, path = self._fn(web_app)(1, None, str(tmp_path), False)
        assert is_open is True
        assert path == str(tmp_path)

    def test_open_with_file_path_uses_parent_dirname(
        self, web_app, monkeypatch, tmp_path
    ):
        self._ctx(monkeypatch, IDs.BTN_PIPE_BROWSE)
        sub = tmp_path / "sub"
        sub.mkdir()
        fake_file = sub / "output.edi"
        is_open, path = self._fn(web_app)(1, None, str(fake_file), False)
        assert is_open is True
        assert path == str(sub)

    def test_open_with_missing_folder_falls_back_to_home(
        self, web_app, monkeypatch
    ):
        self._ctx(monkeypatch, IDs.BTN_PIPE_BROWSE)
        is_open, path = self._fn(web_app)(
            1, None, "Z:/definitely/not/a/real/path", False
        )
        assert is_open is True
        assert path == os.path.expanduser("~")

    def test_open_with_no_current_folder_defaults_to_home(
        self, web_app, monkeypatch
    ):
        self._ctx(monkeypatch, IDs.BTN_PIPE_BROWSE)
        is_open, path = self._fn(web_app)(1, None, "", False)
        assert is_open is True
        assert path == os.path.expanduser("~")


class TestNavigateInto:
    def _fn(self, web_app):
        return _cb_by_input(
            web_app, f"{IDs.PIPE_BROWSE_PATH}.data", "pipe-dir-item"
        )

    def test_no_clicks_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)([0, 0])
        with pytest.raises(PreventUpdate):
            self._fn(web_app)([])

    def test_navigates_into_valid_dir(self, web_app, monkeypatch, tmp_path):
        import pycsamt.app.web.callbacks.pipeline as pipeline_mod

        sub = tmp_path / "child"
        sub.mkdir()
        monkeypatch.setattr(
            pipeline_mod,
            "ctx",
            type(
                "C",
                (),
                {"triggered_id": {"type": "pipe-dir-item", "index": str(sub)}},
            )(),
        )
        out = self._fn(web_app)([1])
        assert out == str(sub)

    def test_nonexistent_dir_prevents_update(self, web_app, monkeypatch):
        import pycsamt.app.web.callbacks.pipeline as pipeline_mod

        monkeypatch.setattr(
            pipeline_mod,
            "ctx",
            type(
                "C",
                (),
                {
                    "triggered_id": {
                        "type": "pipe-dir-item",
                        "index": "Z:/nope/nope",
                    }
                },
            )(),
        )
        with pytest.raises(PreventUpdate):
            self._fn(web_app)([1])

    def test_non_dict_triggered_id_prevents_update(self, web_app, monkeypatch):
        import pycsamt.app.web.callbacks.pipeline as pipeline_mod

        monkeypatch.setattr(
            pipeline_mod, "ctx", type("C", (), {"triggered_id": "something"})()
        )
        with pytest.raises(PreventUpdate):
            self._fn(web_app)([1])


class TestNavigateUp:
    def _fn(self, web_app):
        return _cb_by_input(
            web_app, f"{IDs.PIPE_BROWSE_PATH}.data", IDs.BTN_PIPE_BROWSE_UP
        )

    def test_no_clicks_or_no_path_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None, "/a/b")
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(1, None)

    def test_goes_to_parent(self, web_app, tmp_path):
        sub = tmp_path / "child"
        sub.mkdir()
        out = self._fn(web_app)(1, str(sub))
        assert out == str(tmp_path)

    def test_root_prevents_update(self, web_app):
        root = os.path.abspath(os.sep)
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(1, root)


class TestRefreshListing:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.PIPE_BROWSE_LIST}.children")

    def test_none_path_defaults_to_home(self, web_app):
        children, cwd = self._fn(web_app)(None)
        assert cwd == os.path.expanduser("~")

    def test_lists_subfolders_only(self, web_app, tmp_path):
        (tmp_path / "alpha").mkdir()
        (tmp_path / "beta").mkdir()
        (tmp_path / "afile.txt").write_text("x")
        children, cwd = self._fn(web_app)(str(tmp_path))
        assert cwd == str(tmp_path)
        listgroup = children[0]
        assert len(listgroup.children) == 2

    def test_empty_folder_shows_hint(self, web_app, tmp_path):
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()
        children, cwd = self._fn(web_app)(str(empty_dir))
        assert "No sub-folders" in children[0].children


class TestMakeDir:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.PIPE_BROWSE_NEW}.value")

    def test_missing_args_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None, "x", "/tmp")
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(1, "", "/tmp")
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(1, "x", "")

    def test_creates_dir_and_navigates_into_it(self, web_app, tmp_path):
        new_path, cleared = self._fn(web_app)(1, "My New/Folder", str(tmp_path))
        assert cleared == ""
        assert os.path.isdir(new_path)
        assert "My New_Folder" in new_path

    def test_makedirs_failure_prevents_update(self, web_app, monkeypatch):
        import pycsamt.app.web.callbacks.pipeline as pipeline_mod

        monkeypatch.setattr(
            pipeline_mod.os,
            "makedirs",
            lambda *a, **k: (_ for _ in ()).throw(OSError("nope")),
        )
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(1, "x", "/tmp")


class TestSelectFolder:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.PIPE_EXPORT_FOLDER}.value")

    def test_missing_args_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None, "/x")
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(1, None)

    def test_selects_and_closes(self, web_app):
        value, is_open = self._fn(web_app)(1, "/some/path")
        assert value == "/some/path"
        assert is_open is False
