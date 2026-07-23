# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for PipelineController — step catalogue, state machine, and the
8-step ``execute_step()`` dispatch logic in
``pycsamt.app.desktop.controllers.pipeline_controller``.

Real data
---------
data/MT/kap03lmt_edis/  — 26 KP TIPPER EDIs (flat folder — used for the
                           Load step's "folder" method).
data/AMT/WILLY_DATA/    — 128 WILLY AMT EDIs across 5 line subfolders;
                           a small 5-station subset of L18PLT is used for
                           the (fast) processing-step sweep.

Strategy
--------
* Cheap dataclass/enum surface (``StepStatus``, ``ParamSpec``,
  ``MethodSpec``, ``PipelineStep``) is covered directly with synthetic
  data — no I/O needed.
* ``PipelineController``'s state machine (``set_input_sites``,
  ``is_ready``, ``reset``, ``to_config``/``from_config``,
  ``get_output_sites``) and the ``_n``/``_coerce_sites`` staticmethods
  are covered directly.
* ``_build_steps()`` is called to obtain the real step catalogue, and
  every step/method combination is exercised through
  ``PipelineController.execute_step()`` against a small real-EDI subset.

Known bugs found while writing these tests (NOT fixed — tests document
the current, real behaviour instead)
--------------------------------------------------------------------------
1. Step 2 "Frequency Edit" / method "regrid" calls
   ``et.regrid_logspace(sites_in, n_freq=n)``, but ``regrid_logspace``
   has no ``n_freq`` parameter (its real knobs are ``per_decade`` /
   ``n_per_decade`` / ``fmin`` / ``fmax``). Always raises ``TypeError``.
2. Step 4 "Noise Removal" / method "emap": the ParamSpec offers options
   ``["median", "mean", "gaussian"]`` for ``apply_emap_filter``'s
   ``method`` kwarg, but that function only accepts
   ``'ama'``/``'flma'``/``'tma'`` (aliases ``'adaptive'``/``'fixed'``/
   ``'trimmed'``). Every offered UI option raises
   ``ValueError: method must be one of 'ama', 'flma', or 'tma'.``
3. Step 7 "Export" / method "edi" calls
   ``et.export_edis(edis, savepath=folder)``, but ``export_edis``'s real
   signature is ``export_edis(edi_objs, new_z, savepath=None, **kws)`` —
   the required ``new_z`` positional argument is never supplied, so the
   Export step always raises ``TypeError`` once a folder is set.
4. Step 1 "QC Screening" / method "drop_low_conf": the "Score method"
   ParamSpec (coverage/snr/coherence) is defined but never forwarded to
   ``et.drop_low_confidence_frequencies()`` — a dead UI parameter; the
   function always runs with its own default (``"composite"``).
5. Step 2 "Frequency Edit" / method "close_gaps": the "max_gap"
   ParamSpec is defined but never forwarded to ``et.close_skew_gaps()``
   — dead UI parameter, the function's own default (``max_gap=1``) is
   always used regardless of the value shown/set in the UI.
6. Step 5 "Strike Analysis": the ``result_info`` angle-formatting logic
   looks for a DataFrame column named "strike"/"angle"/"theta"/
   "azimuth", but ``estimate_strike_consensus`` / ``_phase_tensor`` /
   ``_sweep`` all return the angle in a column named "ang". The
   estimated angle is therefore never shown — the generic
   "Strike table computed." message is always used instead of
   "Strike ~ X deg".
7. ``PipelineController.from_config()`` sets ``step.active_method``
   from the saved config but never calls ``step.reset_params()``
   first. It then only copies a saved param over if its key already
   exists in the *current* (stale, previous-method) ``params`` dict.
   So restoring a config whose step used a non-default method drops
   every one of that method's params that don't happen to share a key
   name with the step's default method, leaving the step with a mix
   of the correct method but the wrong (previous/default method's)
   param values.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from pycsamt.app.desktop.controllers.pipeline_controller import (
    STATUS_ICON,
    MethodSpec,
    ParamSpec,
    PipelineController,
    PipelineStep,
    StepStatus,
    _build_steps,
)

# ── Paths ─────────────────────────────────────────────────────────────────────

_ROOT = Path(__file__).parents[3]  # pycsamt/
_TIPPER = _ROOT / "data" / "MT" / "kap03lmt_edis"
_WILLY_L18 = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"

_HAS_TIPPER = _TIPPER.exists() and any(_TIPPER.glob("*.edi"))
_HAS_WILLY = _WILLY_L18.exists() and any(_WILLY_L18.glob("*.edi"))


# ── Fixtures ──────────────────────────────────────────────────────────────────


@pytest.fixture(scope="session")
def willy_subset():
    """Small (5-station) real Sites subset — first 5 EDIs of WILLY L18PLT."""
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_WILLY:
        pytest.skip("WILLY data not available")
    from pycsamt.emtools import ensure_sites
    from pycsamt.site.base import Sites

    full = ensure_sites(str(_WILLY_L18))
    return Sites(full.as_list()[:5])


@pytest.fixture
def controller():
    return PipelineController()


@pytest.fixture
def ready_controller(controller, willy_subset):
    controller.set_input_sites(willy_subset)
    return controller


# ── StepStatus / STATUS_ICON ─────────────────────────────────────────────────


class TestStepStatus:
    def test_values(self):
        assert StepStatus.PENDING.value == "pending"
        assert StepStatus.QUEUED.value == "queued"
        assert StepStatus.RUNNING.value == "running"
        assert StepStatus.DONE.value == "done"
        assert StepStatus.ERROR.value == "error"
        assert StepStatus.SKIPPED.value == "skipped"

    def test_six_members(self):
        assert len(StepStatus) == 6

    def test_status_icon_covers_every_member(self):
        for s in StepStatus:
            assert s in STATUS_ICON
            assert isinstance(STATUS_ICON[s], str)
            assert STATUS_ICON[s]


# ── ParamSpec ─────────────────────────────────────────────────────────────────


class TestParamSpec:
    def test_required_fields(self):
        p = ParamSpec("float", "Label", 1.0)
        assert p.kind == "float"
        assert p.label == "Label"
        assert p.default == 1.0

    def test_defaults(self):
        p = ParamSpec("float", "Label", 1.0)
        assert p.options == []
        assert p.min_val == 0.0
        assert p.max_val == 1e9
        assert p.decimals == 3

    def test_combo_options_stored(self):
        p = ParamSpec("combo", "Mode", "a", options=["a", "b", "c"])
        assert p.options == ["a", "b", "c"]

    def test_independent_default_option_lists(self):
        """default_factory=list must give each instance its own list."""
        p1 = ParamSpec("combo", "A", "x")
        p2 = ParamSpec("combo", "B", "y")
        p1.options.append("mutated")
        assert p2.options == []

    def test_bounds_and_decimals(self):
        p = ParamSpec(
            "float", "L", 0.5, min_val=0.0, max_val=1.0, decimals=2
        )
        assert p.min_val == 0.0
        assert p.max_val == 1.0
        assert p.decimals == 2


# ── MethodSpec ────────────────────────────────────────────────────────────────


class TestMethodSpec:
    def test_defaults_empty_params(self):
        m = MethodSpec("foo", "Foo")
        assert m.name == "foo"
        assert m.label == "Foo"
        assert m.params == {}

    def test_with_params(self):
        m = MethodSpec(
            "foo", "Foo", params={"x": ParamSpec("int", "X", 1)}
        )
        assert "x" in m.params
        assert m.params["x"].default == 1

    def test_independent_default_params_dicts(self):
        m1 = MethodSpec("a", "A")
        m2 = MethodSpec("b", "B")
        m1.params["x"] = ParamSpec("int", "X", 1)
        assert m2.params == {}


# ── PipelineStep ──────────────────────────────────────────────────────────────


def _make_step(active_method: str = "") -> PipelineStep:
    methods = [
        MethodSpec(
            "a", "Method A", params={"x": ParamSpec("int", "X", 5)}
        ),
        MethodSpec(
            "b", "Method B", params={"y": ParamSpec("float", "Y", 1.5)}
        ),
    ]
    return PipelineStep(
        id=0,
        name="Test Step",
        description="d",
        methods=methods,
        default_method="a",
        active_method=active_method,
    )


class TestPipelineStepPostInit:
    def test_active_method_filled_from_default_when_blank(self):
        step = _make_step()
        assert step.active_method == "a"

    def test_explicit_active_method_not_overridden(self):
        step = _make_step(active_method="b")
        assert step.active_method == "b"

    def test_params_seeded_from_default_method(self):
        step = _make_step()
        assert step.params == {"x": 5}

    def test_params_seeded_from_explicit_active_method(self):
        step = _make_step(active_method="b")
        assert step.params == {"y": 1.5}

    def test_status_starts_pending(self):
        step = _make_step()
        assert step.status == StepStatus.PENDING

    def test_output_sites_starts_none(self):
        step = _make_step()
        assert step.output_sites is None

    def test_can_skip_default_true(self):
        step = _make_step()
        assert step.can_skip is True

    def test_abort_on_error_default_false(self):
        step = _make_step()
        assert step.abort_on_error is False


class TestPipelineStepResetParams:
    def test_reset_after_switching_method(self):
        step = _make_step()
        step.active_method = "b"
        step.reset_params()
        assert step.params == {"y": 1.5}

    def test_reset_no_match_is_a_noop(self):
        """get_method() returns None for an unknown active_method, and
        reset_params() only reseeds when a match is found — so stale
        params are left untouched."""
        step = _make_step()
        step.active_method = "does-not-exist"
        step.params = {"stale": True}
        step.reset_params()
        assert step.params == {"stale": True}

    def test_reset_does_not_mutate_other_instances(self):
        step1 = _make_step()
        step2 = _make_step()
        step1.params["x"] = 999
        assert step2.params["x"] == 5


class TestPipelineStepGetMethod:
    def test_match_returns_spec(self):
        step = _make_step()
        m = step.get_method("b")
        assert m is not None
        assert m.name == "b"
        assert m.label == "Method B"

    def test_no_match_returns_none(self):
        step = _make_step()
        assert step.get_method("ghost") is None

    def test_empty_name_no_match(self):
        step = _make_step()
        assert step.get_method("") is None


class TestPipelineStepActiveMethodSpec:
    def test_matches_active_method(self):
        step = _make_step()
        spec = step.active_method_spec
        assert spec is not None
        assert spec.name == "a"

    def test_reflects_method_switch(self):
        step = _make_step()
        step.active_method = "b"
        assert step.active_method_spec.name == "b"

    def test_none_when_active_method_invalid(self):
        step = _make_step()
        step.active_method = "ghost"
        assert step.active_method_spec is None


# ── _build_steps() catalogue ──────────────────────────────────────────────────


class TestBuildSteps:
    def test_returns_8_steps(self):
        assert len(_build_steps()) == 8

    def test_step_ids_sequential_0_to_7(self):
        assert [s.id for s in _build_steps()] == list(range(8))

    def test_every_default_method_exists_in_its_methods_list(self):
        for step in _build_steps():
            names = [m.name for m in step.methods]
            assert step.default_method in names, step.name

    def test_load_and_export_cannot_be_skipped(self):
        steps = _build_steps()
        assert steps[0].can_skip is False
        assert steps[7].can_skip is False

    def test_middle_steps_can_be_skipped(self):
        steps = _build_steps()
        for step in steps[1:7]:
            assert step.can_skip is True, step.name

    def test_every_step_has_at_least_one_method(self):
        for step in _build_steps():
            assert len(step.methods) >= 1, step.name

    def test_fresh_independent_instances_each_call(self):
        a = _build_steps()
        b = _build_steps()
        assert a[0] is not b[0]
        a[0].params["mutated"] = True
        assert "mutated" not in b[0].params

    def test_expected_step_names(self):
        names = [s.name for s in _build_steps()]
        assert names == [
            "Load Data",
            "QC Screening",
            "Frequency Edit",
            "Static Shift Correction",
            "Noise Removal",
            "Strike Analysis",
            "Strike Rotation",
            "Export",
        ]


# ── PipelineController construction ──────────────────────────────────────────


class TestControllerConstruction:
    def test_has_8_steps(self, controller):
        assert len(controller.steps) == 8

    def test_sites_chain_starts_as_8_nones(self, controller):
        assert controller._sites_chain == [None] * 8

    def test_sites_input_starts_none(self, controller):
        assert controller._sites_input is None

    def test_is_ready_false_initially(self, controller):
        assert controller.is_ready() is False

    def test_two_instances_have_independent_steps(self):
        c1 = PipelineController()
        c2 = PipelineController()
        c1.steps[0].params["x"] = "only-on-c1"
        assert "x" not in c2.steps[0].params


# ── set_input_sites / is_ready ────────────────────────────────────────────────


class TestSetInputSites:
    def test_sets_sites_input(self, controller, willy_subset):
        controller.set_input_sites(willy_subset)
        assert controller._sites_input is willy_subset

    def test_is_ready_true_after_set(self, controller, willy_subset):
        controller.set_input_sites(willy_subset)
        assert controller.is_ready() is True

    def test_none_sites_keeps_not_ready(self, controller):
        controller.set_input_sites(None)
        assert controller.is_ready() is False

    def test_prefills_step0_result_info_for_current_method(
        self, controller, willy_subset
    ):
        assert controller.steps[0].active_method == "current"
        controller.set_input_sites(willy_subset)
        assert "Ready" in controller.steps[0].result_info
        assert "stations" in controller.steps[0].result_info

    def test_no_prefill_when_step0_method_not_current(
        self, controller, willy_subset
    ):
        controller.steps[0].active_method = "folder"
        controller.set_input_sites(willy_subset)
        assert controller.steps[0].result_info == ""

    def test_no_prefill_when_sites_is_none(self, controller):
        controller.set_input_sites(None)
        assert controller.steps[0].result_info == ""


# ── reset() ───────────────────────────────────────────────────────────────────


class TestReset:
    def test_reset_clears_status_for_all_steps(self, ready_controller):
        ready_controller.execute_step(0)
        ready_controller.reset()
        assert all(
            s.status == StepStatus.PENDING for s in ready_controller.steps
        )

    def test_reset_clears_sites_chain(self, ready_controller):
        ready_controller.execute_step(0)
        ready_controller.reset()
        assert ready_controller._sites_chain == [None] * 8

    def test_reset_clears_output_sites_and_messages(self, ready_controller):
        ready_controller.execute_step(0)
        step0 = ready_controller.steps[0]
        assert step0.output_sites is not None
        ready_controller.reset()
        assert step0.output_sites is None
        assert step0.result_info == ""
        assert step0.error_msg == ""
        assert step0.elapsed_s == 0.0

    def test_reset_preserves_method_choice(self, ready_controller):
        ready_controller.steps[3].active_method = "loess"
        ready_controller.reset()
        assert ready_controller.steps[3].active_method == "loess"

    def test_reset_preserves_params(self, ready_controller):
        ready_controller.steps[3].params["half_window"] = 7
        ready_controller.reset()
        assert ready_controller.steps[3].params["half_window"] == 7

    def test_reset_on_fresh_controller_does_not_raise(self, controller):
        controller.reset()  # nothing has run yet


# ── to_config() / from_config() ──────────────────────────────────────────────


class TestConfigRoundTrip:
    def test_to_config_has_8_step_entries(self, controller):
        cfg = controller.to_config()
        assert "steps" in cfg
        assert len(cfg["steps"]) == 8

    def test_to_config_entry_shape(self, controller):
        cfg = controller.to_config()
        for entry, step in zip(cfg["steps"], controller.steps):
            assert entry["id"] == step.id
            assert entry["name"] == step.name
            assert entry["method"] == step.active_method
            assert entry["params"] == step.params

    def test_to_config_params_is_a_copy(self, controller):
        cfg = controller.to_config()
        cfg["steps"][0]["params"]["injected"] = True
        assert "injected" not in controller.steps[0].params

    def test_round_trip_preserves_method_and_params_when_method_unchanged(
        self, controller
    ):
        """Round trip works cleanly when the saved method equals the
        fresh controller's already-seeded default_method, since the
        param keys line up."""
        controller.steps[3].params["half_window"] = 9
        controller.steps[3].params["weights"] = "gauss"
        cfg = controller.to_config()

        fresh = PipelineController()
        fresh.from_config(cfg)

        assert fresh.steps[3].active_method == "ama"
        assert fresh.steps[3].params["half_window"] == 9
        assert fresh.steps[3].params["weights"] == "gauss"

    def test_round_trip_drops_params_for_non_default_method(self, controller):
        """KNOWN BUG #7 — see module docstring: from_config() sets
        active_method but never calls reset_params() first, so it applies
        saved params against the *stale* (previous-method) params dict.
        Any saved param key that doesn't happen to also exist in the
        stale dict is silently dropped, and the step is left with the
        wrong method's param values."""
        controller.steps[3].active_method = "bilateral"
        controller.steps[3].reset_params()
        controller.steps[3].params["sig_dist"] = 4.5
        controller.steps[3].params["sig_val"] = 0.75
        cfg = controller.to_config()

        fresh = PipelineController()
        fresh.from_config(cfg)

        # the method choice itself IS restored correctly...
        assert fresh.steps[3].active_method == "bilateral"
        # ...but its params are NOT: "sig_dist"/"sig_val" don't exist in
        # the stale "ama"-seeded params dict, so they're dropped, and the
        # step is left holding "ama"'s own default params instead.
        assert "sig_dist" not in fresh.steps[3].params
        assert "sig_val" not in fresh.steps[3].params
        assert fresh.steps[3].params == {
            "half_window": 3,
            "weights": "tri",
            "max_skew": 6.0,
        }

    def test_round_trip_all_steps_active_method_always_restored(
        self, controller
    ):
        """active_method is always restored correctly regardless of the
        params bug above."""
        controller.steps[1].active_method = "flag_only"
        controller.steps[5].active_method = "sweep"
        controller.steps[5].reset_params()
        controller.steps[6].active_method = "none"
        cfg = controller.to_config()

        fresh = PipelineController()
        fresh.from_config(cfg)
        for a, b in zip(controller.steps, fresh.steps):
            assert a.active_method == b.active_method

    def test_from_config_ignores_unknown_param_keys(self, controller):
        cfg = {
            "steps": [
                {"id": 0, "method": "current", "params": {"bogus": 1}}
            ]
        }
        controller.from_config(cfg)
        assert "bogus" not in controller.steps[0].params

    def test_from_config_skips_out_of_range_id(self, controller):
        cfg = {"steps": [{"id": 99, "method": "x", "params": {}}]}
        controller.from_config(cfg)  # must not raise

    def test_from_config_skips_entry_with_missing_id(self, controller):
        cfg = {"steps": [{"method": "x", "params": {}}]}
        controller.from_config(cfg)  # id is None -> skipped, must not raise

    def test_from_config_defaults_missing_method_to_default_method(
        self, controller
    ):
        controller.steps[3].active_method = "loess"  # perturb first
        cfg = {"steps": [{"id": 3, "params": {}}]}
        controller.from_config(cfg)
        assert controller.steps[3].active_method == (
            controller.steps[3].default_method
        )

    def test_from_config_empty_steps_list_is_noop(self, controller):
        before = controller.to_config()
        controller.from_config({"steps": []})
        after = controller.to_config()
        assert before == after

    def test_from_config_missing_steps_key_does_not_raise(self, controller):
        controller.from_config({})


# ── get_output_sites() ────────────────────────────────────────────────────────


class TestGetOutputSites:
    def test_no_input_no_steps_returns_none(self, controller):
        assert controller.get_output_sites() is None

    def test_before_any_step_returns_raw_input(self, ready_controller, willy_subset):
        assert ready_controller.get_output_sites() is willy_subset

    def test_after_one_step_returns_its_output(self, ready_controller):
        result = ready_controller.execute_step(0)
        assert ready_controller.get_output_sites() is result

    def test_returns_latest_completed_step_output(self, ready_controller):
        ready_controller.execute_step(0)
        r1 = ready_controller.execute_step(1)
        assert ready_controller.get_output_sites() is r1


# ── _n() staticmethod ─────────────────────────────────────────────────────────


class TestNStaticmethod:
    def test_real_sites(self, willy_subset):
        assert PipelineController._n(willy_subset) == len(willy_subset)
        assert PipelineController._n(willy_subset) == 5

    def test_plain_list(self):
        assert PipelineController._n([1, 2, 3]) == 3

    def test_empty_list(self):
        assert PipelineController._n([]) == 0

    def test_none_returns_zero(self):
        assert PipelineController._n(None) == 0

    def test_non_iterable_returns_zero(self):
        assert PipelineController._n(42) == 0

    def test_string_counts_characters(self):
        # str is iterable, so list("abc") == ["a", "b", "c"]
        assert PipelineController._n("abc") == 3


# ── _coerce_sites() staticmethod ─────────────────────────────────────────────


class TestCoerceSitesStaticmethod:
    def test_none_result_returns_fallback(self):
        fallback = object()
        assert PipelineController._coerce_sites(None, fallback) is fallback

    def test_sites_like_result_returned_as_is(self, willy_subset):
        fallback = object()
        assert (
            PipelineController._coerce_sites(willy_subset, fallback)
            is willy_subset
        )

    def test_plain_list_result_falls_back(self):
        fallback = object()
        assert (
            PipelineController._coerce_sites([1, 2, 3], fallback)
            is fallback
        )

    def test_result_wrapped_in_sites_attr(self, willy_subset):
        class Wrapper:
            sites = willy_subset

        assert (
            PipelineController._coerce_sites(Wrapper(), "fb")
            is willy_subset
        )

    def test_result_wrapped_in_data_attr(self, willy_subset):
        class Wrapper:
            data = willy_subset

        assert (
            PipelineController._coerce_sites(Wrapper(), "fb")
            is willy_subset
        )

    def test_sites_attr_none_falls_through_to_data_attr(self, willy_subset):
        class Wrapper:
            sites = None
            data = willy_subset

        assert (
            PipelineController._coerce_sites(Wrapper(), "fb")
            is willy_subset
        )

    def test_wrapped_attr_not_sites_like_falls_back(self):
        class Wrapper:
            sites = "not-sites-like"

        fallback = object()
        assert (
            PipelineController._coerce_sites(Wrapper(), fallback)
            is fallback
        )

    def test_no_relevant_attrs_falls_back(self):
        class Wrapper:
            pass

        fallback = object()
        assert (
            PipelineController._coerce_sites(Wrapper(), fallback)
            is fallback
        )


# ── execute_step(): Load (step 0) ────────────────────────────────────────────


class TestExecuteStepLoad:
    def test_current_method_returns_input_unchanged(
        self, ready_controller, willy_subset
    ):
        result = ready_controller.execute_step(0)
        assert result is willy_subset

    def test_current_method_sets_result_info(self, ready_controller):
        ready_controller.execute_step(0)
        assert "stations" in ready_controller.steps[0].result_info

    def test_current_method_marks_step_done(self, ready_controller):
        ready_controller.execute_step(0)
        assert ready_controller.steps[0].status == StepStatus.DONE

    def test_folder_method_no_folder_raises_value_error(self, controller):
        controller.steps[0].active_method = "folder"
        controller.steps[0].reset_params()
        with pytest.raises(ValueError, match="No folder selected"):
            controller.execute_step(0)
        assert controller.steps[0].status == StepStatus.ERROR
        assert controller.steps[0].error_msg

    def test_folder_method_empty_folder_raises_no_edi_files(
        self, controller, tmp_path
    ):
        controller.steps[0].active_method = "folder"
        controller.steps[0].reset_params()
        controller.steps[0].params["folder"] = str(tmp_path)
        with pytest.raises(ValueError, match="No EDI files found"):
            controller.execute_step(0)

    @pytest.mark.skipif(not _HAS_TIPPER, reason="TIPPER data not available")
    def test_folder_method_loads_real_edis(self, controller):
        controller.steps[0].active_method = "folder"
        controller.steps[0].reset_params()
        controller.steps[0].params["folder"] = str(_TIPPER)
        result = controller.execute_step(0)
        assert result is not None
        assert len(list(result)) > 0
        assert controller.steps[0].status == StepStatus.DONE
        # side effect noted in source: also updates _sites_input
        assert controller._sites_input is result


# ── execute_step(): QC Screening (step 1) ────────────────────────────────────


class TestExecuteStepQC:
    def test_drop_low_conf_default(self, ready_controller):
        result = ready_controller.execute_step(1)
        assert result is not None
        assert ready_controller.steps[1].status == StepStatus.DONE
        assert "stations after QC" in ready_controller.steps[1].result_info

    def test_flag_only_returns_input_unchanged(
        self, ready_controller, willy_subset
    ):
        ready_controller.steps[1].active_method = "flag_only"
        ready_controller.steps[1].reset_params()
        result = ready_controller.execute_step(1)
        assert result is willy_subset
        assert "unchanged" in ready_controller.steps[1].result_info

    def test_qc_score_method_param_never_forwarded(
        self, ready_controller, monkeypatch
    ):
        """KNOWN BUG #4 — see module docstring: the 'method' ParamSpec is
        defined in the catalogue but silently dropped before calling
        drop_low_confidence_frequencies()."""
        import pycsamt.emtools as et

        captured = {}
        orig = et.drop_low_confidence_frequencies

        def spy(*args, **kwargs):
            captured.update(kwargs)
            return orig(*args, **kwargs)

        monkeypatch.setattr(et, "drop_low_confidence_frequencies", spy)
        ready_controller.steps[1].params["method"] = "snr"
        ready_controller.execute_step(1)
        assert "method" not in captured


# ── execute_step(): Frequency Edit (step 2) ──────────────────────────────────


class TestExecuteStepFrequencyEdit:
    def test_edit_by_conf_default(self, ready_controller):
        result = ready_controller.execute_step(2)
        assert result is not None
        assert ready_controller.steps[2].status == StepStatus.DONE
        assert "recover" in ready_controller.steps[2].result_info

    def test_close_gaps(self, ready_controller):
        ready_controller.steps[2].active_method = "close_gaps"
        ready_controller.steps[2].reset_params()
        result = ready_controller.execute_step(2)
        assert result is not None
        assert ready_controller.steps[2].status == StepStatus.DONE
        assert "Gaps closed" in ready_controller.steps[2].result_info

    def test_close_gaps_max_gap_param_never_forwarded(
        self, ready_controller, monkeypatch
    ):
        """KNOWN BUG #5 — see module docstring: the 'max_gap' ParamSpec is
        defined in the catalogue but never forwarded to
        close_skew_gaps()."""
        import pycsamt.emtools as et

        captured = {}
        orig = et.close_skew_gaps

        def spy(*args, **kwargs):
            captured.update(kwargs)
            return orig(*args, **kwargs)

        monkeypatch.setattr(et, "close_skew_gaps", spy)
        ready_controller.steps[2].active_method = "close_gaps"
        ready_controller.steps[2].reset_params()
        ready_controller.steps[2].params["max_gap"] = 15
        ready_controller.execute_step(2)
        assert "max_gap" not in captured

    def test_regrid_raises_type_error(self, ready_controller):
        """KNOWN BUG #1 — see module docstring: regrid_logspace() has no
        n_freq kwarg, so this method variant always raises TypeError."""
        ready_controller.steps[2].active_method = "regrid"
        ready_controller.steps[2].reset_params()
        with pytest.raises(TypeError, match="n_freq"):
            ready_controller.execute_step(2)
        assert ready_controller.steps[2].status == StepStatus.ERROR
        assert ready_controller.steps[2].error_msg

    def test_regrid_success_path_when_patched(
        self, ready_controller, monkeypatch, willy_subset
    ):
        """Isolates the controller's own success-path logic (result_info
        / return value) from bug #1 by monkeypatching regrid_logspace()
        to accept the n_freq kwarg the controller actually passes."""
        import pycsamt.emtools as et

        monkeypatch.setattr(
            et, "regrid_logspace", lambda sites, n_freq=30, **kw: sites
        )
        ready_controller.steps[2].active_method = "regrid"
        ready_controller.steps[2].reset_params()
        result = ready_controller.execute_step(2)
        assert result is willy_subset
        assert ready_controller.steps[2].result_info == (
            "Regridded to 30 frequencies."
        )
        assert ready_controller.steps[2].status == StepStatus.DONE


# ── execute_step(): Static Shift (step 3) — every method variant ────────────


class TestExecuteStepStaticShift:
    @pytest.mark.parametrize("method", ["ama", "loess", "bilateral", "refmedian"])
    def test_variant_does_not_raise_and_marks_done(
        self, ready_controller, method
    ):
        ready_controller.steps[3].active_method = method
        ready_controller.steps[3].reset_params()
        result = ready_controller.execute_step(3)
        assert result is not None
        assert ready_controller.steps[3].status == StepStatus.DONE
        assert "SS corrected" in ready_controller.steps[3].result_info


# ── execute_step(): Noise Removal (step 4) ───────────────────────────────────


class TestExecuteStepNoiseRemoval:
    def test_pipeline_default(self, ready_controller):
        result = ready_controller.execute_step(4)
        assert result is not None
        assert ready_controller.steps[4].status == StepStatus.DONE
        assert "Noise removed" in ready_controller.steps[4].result_info

    def test_rpca(self, ready_controller):
        ready_controller.steps[4].active_method = "rpca"
        ready_controller.steps[4].reset_params()
        result = ready_controller.execute_step(4)
        assert result is not None
        assert ready_controller.steps[4].status == StepStatus.DONE

    def test_rpca_failure_falls_back_to_noise_pipeline(
        self, ready_controller, monkeypatch, willy_subset
    ):
        """Covers the internal try/except fallback: if
        rpca_offdiag_denoise() raises, _run_step must not propagate the
        error — it logs and falls back to remove_noise_pipeline()."""
        import pycsamt.emtools as et

        def boom(*a, **kw):
            raise RuntimeError("synthetic RPCA failure")

        monkeypatch.setattr(et, "rpca_offdiag_denoise", boom)
        ready_controller.steps[4].active_method = "rpca"
        ready_controller.steps[4].reset_params()
        logs = []
        result = ready_controller.execute_step(4, logs.append)
        assert result is not None
        assert ready_controller.steps[4].status == StepStatus.DONE
        assert any("RPCA failed" in m for m in logs)

    @pytest.mark.parametrize("method", ["median", "mean", "gaussian"])
    def test_emap_options_all_raise_value_error(
        self, ready_controller, method
    ):
        """KNOWN BUG #2 — see module docstring: every UI-offered 'method'
        option for the EMAP filter is invalid; apply_emap_filter() only
        accepts 'ama'/'flma'/'tma'."""
        ready_controller.steps[4].active_method = "emap"
        ready_controller.steps[4].reset_params()
        ready_controller.steps[4].params["method"] = method
        with pytest.raises(ValueError, match="ama.*flma.*tma"):
            ready_controller.execute_step(4)
        assert ready_controller.steps[4].status == StepStatus.ERROR


# ── execute_step(): Strike Analysis (step 5) — every method variant ─────────


class TestExecuteStepStrike:
    @pytest.mark.parametrize(
        "method", ["consensus", "phase_tensor", "sweep"]
    )
    def test_variant_does_not_raise_and_returns_input_unchanged(
        self, ready_controller, willy_subset, method
    ):
        ready_controller.steps[5].active_method = method
        ready_controller.steps[5].reset_params()
        result = ready_controller.execute_step(5)
        # strike step never transforms sites
        assert result is willy_subset
        assert ready_controller.steps[5].status == StepStatus.DONE

    def test_result_info_never_shows_estimated_angle(self, ready_controller):
        """KNOWN BUG #6 — see module docstring: the angle-column
        autodetection looks for 'strike'/'angle'/'theta'/'azimuth', but
        the real strike tables use column name 'ang', so the angle is
        never rendered into result_info."""
        ready_controller.execute_step(5)
        assert (
            ready_controller.steps[5].result_info == "Strike table computed."
        )

    def test_result_info_shows_angle_when_column_name_matches(
        self, ready_controller, monkeypatch
    ):
        """Exercises the (currently unreachable in practice, see bug #6)
        success branch directly: if the strike table did carry a
        recognised angle column, result_info would render the median
        angle."""
        import pycsamt.emtools as et
        import pandas as pd

        fake_tbl = pd.DataFrame({"station": ["s1", "s2"], "strike": [10.0, 20.0]})
        monkeypatch.setattr(
            et, "estimate_strike_consensus", lambda *a, **kw: fake_tbl
        )
        ready_controller.execute_step(5)
        assert "Strike" in ready_controller.steps[5].result_info
        assert "15.0" in ready_controller.steps[5].result_info

    def test_result_info_generic_message_for_empty_table(
        self, ready_controller, monkeypatch
    ):
        """Covers the 'tbl is empty / not a DataFrame' branch."""
        import pycsamt.emtools as et

        monkeypatch.setattr(
            et, "estimate_strike_consensus", lambda *a, **kw: None
        )
        ready_controller.execute_step(5)
        assert (
            ready_controller.steps[5].result_info
            == "Strike analysis complete."
        )


# ── execute_step(): Rotation (step 6) — every method variant ────────────────


class TestExecuteStepRotation:
    def test_none_skips_rotation(self, ready_controller, willy_subset):
        ready_controller.steps[6].active_method = "none"
        ready_controller.steps[6].reset_params()
        result = ready_controller.execute_step(6)
        assert result is willy_subset
        assert "No rotation" in ready_controller.steps[6].result_info

    @pytest.mark.parametrize("method", ["to_strike", "z_rotation"])
    def test_rotation_variants_do_not_raise(self, ready_controller, method):
        ready_controller.steps[6].active_method = method
        ready_controller.steps[6].reset_params()
        result = ready_controller.execute_step(6)
        assert result is not None
        assert ready_controller.steps[6].status == StepStatus.DONE
        assert "Rotated" in ready_controller.steps[6].result_info


# ── execute_step(): Export (step 7) ──────────────────────────────────────────


class TestExecuteStepExport:
    def test_no_folder_raises_value_error(self, ready_controller):
        with pytest.raises(ValueError, match="No output folder selected"):
            ready_controller.execute_step(7)
        assert ready_controller.steps[7].status == StepStatus.ERROR

    def test_with_folder_raises_type_error_missing_new_z(
        self, ready_controller, tmp_path
    ):
        """KNOWN BUG #3 — see module docstring: export_edis() requires a
        positional `new_z` argument that pipeline_controller never
        supplies, so the Export step is currently always broken."""
        out_dir = tmp_path / "export_out"
        ready_controller.steps[7].params["folder"] = str(out_dir)
        with pytest.raises(TypeError, match="new_z"):
            ready_controller.execute_step(7)
        assert ready_controller.steps[7].status == StepStatus.ERROR
        # side effect: the output folder is still created before the crash
        assert out_dir.exists()

    def test_success_path_when_export_edis_is_patched(
        self, ready_controller, monkeypatch, tmp_path
    ):
        """Isolates the controller's own success-path logic (folder
        creation, edi-list extraction, result_info/log messages, return
        value) from bug #3 by monkeypatching et.export_edis() to accept
        the call the controller actually makes."""
        import pycsamt.emtools as et

        calls = {}

        def fake_export_edis(edi_objs, savepath=None, **kws):
            calls["edi_objs"] = edi_objs
            calls["savepath"] = savepath

        monkeypatch.setattr(et, "export_edis", fake_export_edis)
        out_dir = tmp_path / "export_ok"
        ready_controller.steps[7].params["folder"] = str(out_dir)
        logs = []
        result = ready_controller.execute_step(7, logs.append)
        assert ready_controller.steps[7].status == StepStatus.DONE
        assert out_dir.exists()
        assert calls["savepath"] == str(out_dir)
        assert len(calls["edi_objs"]) == len(list(ready_controller._sites_input))
        assert "Exported" in ready_controller.steps[7].result_info
        assert result is ready_controller._sites_input

    def test_success_path_falls_back_to_as_list_when_items_lack_edi(
        self, controller, monkeypatch, tmp_path
    ):
        """Covers the ``if not edis: edis = list(sites_in.as_list())``
        fallback: when iterating input sites yields objects without an
        ``.edi`` attribute, the export step falls back to
        ``sites_in.as_list()``."""
        import pycsamt.emtools as et

        class FakeSitesNoEdiAttr:
            """Iterates plain objects (no .edi) but still exposes
            as_list()."""

            def __init__(self, items):
                self._items = items

            def __iter__(self):
                return iter(self._items)

            def as_list(self):
                return self._items

        fake_edis = ["edi-object-1", "edi-object-2"]
        fake_sites = FakeSitesNoEdiAttr(fake_edis)

        calls = {}
        monkeypatch.setattr(
            et,
            "export_edis",
            lambda edi_objs, savepath=None, **kw: calls.update(
                edi_objs=edi_objs, savepath=savepath
            ),
        )
        controller.set_input_sites(fake_sites)
        out_dir = tmp_path / "export_fallback"
        controller.steps[7].params["folder"] = str(out_dir)
        controller.execute_step(7)

        assert calls["edi_objs"] == fake_edis
        assert controller.steps[7].status == StepStatus.DONE


# ── execute_step(): full catalogue sweep ─────────────────────────────────────
#
# Every method of every step is exercised at least once. Steps 2/"regrid",
# 4/"emap" and all of step 7 are known-broken (see module docstring) and are
# asserted to raise their documented exception instead of being lumped into
# the generic must-not-raise expectation.


def _catalogue_cases():
    cases = []
    for step in _build_steps():
        if step.id in (0, 7):
            continue  # covered by their own dedicated test classes above
        for m in step.methods:
            cases.append((step.id, m.name))
    return cases


_KNOWN_BROKEN = {
    (2, "regrid"): TypeError,
    (4, "emap"): ValueError,
}

_CASES = _catalogue_cases()
_IDS = [f"step{sid}-{name}" for sid, name in _CASES]


class TestExecuteStepCatalogueSweep:
    @pytest.mark.parametrize("step_id,method_name", _CASES, ids=_IDS)
    def test_sweep(self, ready_controller, step_id, method_name):
        step = ready_controller.steps[step_id]
        step.active_method = method_name
        step.reset_params()

        logs = []
        expected_exc = _KNOWN_BROKEN.get((step_id, method_name))
        if expected_exc is not None:
            with pytest.raises(expected_exc):
                ready_controller.execute_step(step_id, logs.append)
            assert step.status == StepStatus.ERROR
            assert step.error_msg
            return

        result = ready_controller.execute_step(step_id, logs.append)
        assert step.status == StepStatus.DONE
        assert result is not None
        assert ready_controller._sites_chain[step_id] is result
        assert logs, "log_cb should receive at least one message"


# ── Sites-chain walking logic ─────────────────────────────────────────────────


class TestSitesChainWalking:
    def test_each_step_records_its_own_output(self, ready_controller):
        r0 = ready_controller.execute_step(0)
        r1 = ready_controller.execute_step(1)
        assert ready_controller._sites_chain[0] is r0
        assert ready_controller._sites_chain[1] is r1

    def test_walks_back_past_skipped_steps(self, ready_controller):
        """Running step 0 then jumping straight to step 3 must not raise —
        _run_step's input-resolution loop walks back through the None
        chain entries for steps 1 and 2 down to step 0's output."""
        ready_controller.execute_step(0)
        assert ready_controller._sites_chain[1] is None
        assert ready_controller._sites_chain[2] is None
        ready_controller.execute_step(3)
        assert ready_controller.steps[3].status == StepStatus.DONE

    def test_falls_back_to_sites_input_when_chain_entirely_none(
        self, ready_controller
    ):
        """Running a mid-pipeline step directly (nothing executed yet)
        must fall back to the raw input sites."""
        ready_controller.execute_step(3)
        assert ready_controller.steps[3].status == StepStatus.DONE


# ── Error handling / misc execute_step contract ──────────────────────────────


class TestExecuteStepErrorHandling:
    def test_error_sets_status_message_and_reraises(self, controller):
        controller.steps[0].active_method = "folder"
        controller.steps[0].reset_params()  # folder == ""
        with pytest.raises(ValueError):
            controller.execute_step(0)
        step0 = controller.steps[0]
        assert step0.status == StepStatus.ERROR
        assert step0.error_msg
        assert step0.elapsed_s >= 0

    def test_unknown_step_id_raises_index_error(self, controller):
        with pytest.raises(IndexError):
            controller.execute_step(99)

    def test_log_callback_receives_messages(self, ready_controller):
        logs = []
        ready_controller.execute_step(0, logs.append)
        assert len(logs) >= 1
        assert all(isinstance(m, str) for m in logs)

    def test_no_log_callback_does_not_raise(self, ready_controller):
        ready_controller.execute_step(0)  # log_cb defaults to None

    def test_successful_step_marks_done_with_elapsed_time(
        self, ready_controller
    ):
        ready_controller.execute_step(0)
        step0 = ready_controller.steps[0]
        assert step0.status == StepStatus.DONE
        assert step0.elapsed_s >= 0
        assert step0.error_msg == ""
