"""Comprehensive tests for the pyCSAMT pipeline.

Public entry point: :mod:`pycsamt.pipeline`
Implementation:     :mod:`pycsamt.pipeline` (this package)
Config:             :mod:`pycsamt.api.pipe`

Coverage
--------
* :mod:`~pycsamt.pipeline._registry`  — StepSpec, STEP_REGISTRY, lookup helpers
* :mod:`~pycsamt.pipeline._steps`     — Step, StepResult
* :mod:`~pycsamt.pipeline._presets`   — Preset, PRESETS, helpers
* :mod:`~pycsamt.api.pipe.config`          — PipelineAPIConfig singleton
* :mod:`~pycsamt.pipeline._output`    — OutputDir
* :mod:`~pycsamt.pipeline._pipeline`  — Pipeline construction + Pipeline.run
* :mod:`~pycsamt.pipeline._config`    — YAML / JSON / .py config loaders
* :class:`~pycsamt.pipeline.PipelineResult`

Testing strategy
----------------
Pipeline mechanics (chaining, error handling, output writing) are tested with
a lightweight :class:`_IdentityStep` subclass whose ``transform`` returns sites
unchanged and whose ``generate_qc_plots`` returns ``[]``.  This isolates the
pipeline engine from the emtools processing functions, which are already
covered in their own test suites.

A small set of *registry-only* tests verify that every registered step can
resolve a callable function without actually calling it.
"""

from __future__ import annotations

import json
import warnings
from typing import Any

import matplotlib

matplotlib.use("Agg")
import pytest

# ─── pipeline imports — all via the public namespace ────────────────────────
from pycsamt.pipeline import (
    PRESETS,
    # config — also re-exported from pycsamt.pipeline
    PYCSAMT_PIPE,
    STEP_REGISTRY,
    Pipeline,
    PipelineAPIConfig,
    PipelineResult,
    Preset,
    Step,
    StepResult,
    categories,
    configure_pipe,
    get_preset,
    list_presets,
    list_steps,
    lookup_step,
    preset_catalogue,
    reset_pipe,
    step_codes,
    step_names,
)

# internal implementation modules — imported directly for unit-testing internals
from pycsamt.pipeline._output import OutputDir
from pycsamt.pipeline._report import (
    make_html_report,
    make_text_report,
)

# ─────────────────────────────────────────────────────────────────────────────
# Test doubles
# ─────────────────────────────────────────────────────────────────────────────


class _IdentityStep(Step):
    """Step whose transform returns sites unchanged, with no QC plots.

    Used throughout run-mechanics tests to avoid real emtools function calls.
    """

    def transform(self, sites: Any) -> Any:
        return sites

    def generate_qc_plots(self, sites: Any) -> list:
        return []


def _identity(code: str, **params: Any) -> _IdentityStep:
    step = _IdentityStep.__new__(_IdentityStep)
    step.spec = lookup_step(code)
    step.params = {**step.spec.defaults, **params}
    return step


class _CountingSites:
    """Minimal sites-like object that records how many steps touched it."""

    def __init__(self, n: int = 3, count: int = 0):
        self._n = n
        self.count = count

    def __len__(self) -> int:
        return self._n


class _MutatingStep(Step):
    """Step that increments ``sites.count`` — proves chaining works."""

    def transform(self, sites: _CountingSites) -> _CountingSites:
        return _CountingSites(n=sites._n, count=sites.count + 1)

    def generate_qc_plots(self, sites: Any) -> list:
        return []


def _mutating(code: str) -> _MutatingStep:
    step = _MutatingStep.__new__(_MutatingStep)
    step.spec = lookup_step(code)
    step.params = dict(step.spec.defaults)
    return step


class _RaisingStep(Step):
    """Step that always raises RuntimeError — for error-handling tests."""

    def transform(self, sites: Any) -> Any:
        raise RuntimeError("deliberate step failure")

    def generate_qc_plots(self, sites: Any) -> list:
        return []


def _raising(code: str) -> _RaisingStep:
    step = _RaisingStep.__new__(_RaisingStep)
    step.spec = lookup_step(code)
    step.params = dict(step.spec.defaults)
    return step


# ─────────────────────────────────────────────────────────────────────────────
# Fixtures
# ─────────────────────────────────────────────────────────────────────────────


@pytest.fixture(autouse=True)
def _reset_pipe_cfg():
    """Guarantee PYCSAMT_PIPE is at defaults for every test."""
    reset_pipe()
    yield
    reset_pipe()


@pytest.fixture()
def simple_pipe() -> Pipeline:
    return Pipeline(
        [
            ("notch", _identity("NR001")),
            ("band", _identity("FREQ001")),
            ("align", _identity("FREQ004")),
        ],
        name="test_pipe",
    )


@pytest.fixture()
def sites() -> _CountingSites:
    return _CountingSites(n=5, count=0)


# ─────────────────────────────────────────────────────────────────────────────
# 1 · Step registry
# ─────────────────────────────────────────────────────────────────────────────


class TestStepRegistry:
    def test_registry_has_55_entries(self):
        assert len(STEP_REGISTRY) == 55

    def test_lookup_by_code(self):
        spec = lookup_step("NR001")
        assert spec.code == "NR001"

    def test_lookup_by_name(self):
        spec = lookup_step("notch_powerline")
        assert spec.code == "NR001"

    def test_lookup_code_and_name_same_object(self):
        by_code = lookup_step("NR001")
        by_name = lookup_step("notch_powerline")
        assert by_code is by_name

    def test_lookup_unknown_raises_key_error(self):
        with pytest.raises(KeyError, match="not found|No pipeline step"):
            lookup_step("DOES_NOT_EXIST")

    def test_list_steps_returns_all(self):
        all_specs = list_steps()
        assert len(all_specs) == 55

    def test_list_steps_by_category_noise_removal(self):
        nr = list_steps("noise_removal")
        assert len(nr) == 14
        assert all(s.category == "noise_removal" for s in nr)

    def test_list_steps_by_category_frequency(self):
        freq = list_steps("frequency")
        assert len(freq) == 9

    def test_list_steps_by_category_export(self):
        export = list_steps("export")
        assert len(export) == 3
        assert {s.code for s in export} == {"EXPORT001", "EXPORT002", "EXPORT003"}
        assert all(s.origin == "builtin" for s in export)

    def test_categories_sorted(self):
        cats = categories()
        assert cats == sorted(cats)
        assert "noise_removal" in cats
        assert "static_shift" in cats
        assert "export" in cats

    def test_step_codes_sorted(self):
        codes = step_codes()
        assert codes == sorted(codes)
        assert "NR001" in codes
        assert "SS003" in codes

    def test_step_names_sorted(self):
        names = step_names()
        assert names == sorted(names)
        assert "notch_powerline" in names

    def test_all_specs_have_label(self):
        for spec in STEP_REGISTRY.values():
            assert spec.label, f"{spec.code} has empty label"

    def test_all_specs_have_name(self):
        for spec in STEP_REGISTRY.values():
            assert spec.name, f"{spec.code} has empty name"

    def test_special_case_override_fns_callable(self):
        for code in ("NR007", "SS002", "SS003"):
            spec = STEP_REGISTRY[code]
            assert callable(spec.override_fn), f"{code} override_fn not callable"

    def test_all_specs_resolvable(self):
        for code, spec in STEP_REGISTRY.items():
            fn = spec.get_fn()
            assert callable(fn), f"{code} did not resolve to a callable"

    def test_qc_defs_are_tuples_of_strings(self):
        for code, spec in STEP_REGISTRY.items():
            for entry in spec.qc_defs:
                assert len(entry) == 2, f"{code}: qc_def not a 2-tuple"
                assert all(
                    isinstance(s, str) for s in entry
                ), f"{code}: qc_def contains non-strings"


# ─────────────────────────────────────────────────────────────────────────────
# 2 · Step class
# ─────────────────────────────────────────────────────────────────────────────


class TestStepClass:
    def test_init_by_code(self):
        step = Step("NR001")
        assert step.spec.code == "NR001"

    def test_init_by_name(self):
        step = Step("notch_powerline")
        assert step.spec.code == "NR001"

    def test_defaults_are_merged(self):
        step = Step("NR001")
        assert "mains_hz" in step.params
        assert step.params["mains_hz"] == 50

    def test_user_params_override_defaults(self):
        step = Step("NR001", mains_hz=60)
        assert step.params["mains_hz"] == 60
        # defaults still present for keys the user did not override
        assert "n_harm" in step.params

    def test_unknown_code_raises_key_error(self):
        with pytest.raises(KeyError):
            Step("TOTALLY_UNKNOWN")

    def test_repr_contains_code_and_label(self):
        r = repr(Step("NR001"))
        assert "NR001" in r
        assert "Power-line" in r

    def test_to_dict_structure(self):
        step = Step("NR001", mains_hz=60)
        d = step.to_dict()
        assert d["code"] == "NR001"
        assert d["name"] == "notch_powerline"
        assert d["params"]["mains_hz"] == 60

    def test_diagnostic_step_returns_sites_unchanged(self):
        """DIM001 (returns_sites=False) must pass sites through."""
        step = _identity("DIM001")
        sentinel = object()
        assert step.transform(sentinel) is sentinel

    def test_qc_snapshot_passthrough(self):
        step = _identity("QC001")
        sentinel = [1, 2, 3]
        result = step.transform(sentinel)
        assert result is sentinel


# ─────────────────────────────────────────────────────────────────────────────
# 3 · Presets
# ─────────────────────────────────────────────────────────────────────────────


class TestPresets:
    _PRESET_NAMES = [
        "basic_qc",
        "noise_reduction",
        "full_processing",
        "tensor_analysis",
        "dimensionality_filter",
        "publication_ready",
        "stratagem_mt",
        "mt_qc",
        "amt_qc",
        "csamt_qc",
        "csumt_qc",
    ]

    def test_eleven_presets_registered(self):
        assert len(PRESETS) == 11

    @pytest.mark.parametrize("name", _PRESET_NAMES)
    def test_preset_loads_without_error(self, name):
        preset = get_preset(name)
        assert isinstance(preset, Preset)

    @pytest.mark.parametrize("name", _PRESET_NAMES)
    def test_preset_steps_not_empty(self, name):
        preset = get_preset(name)
        assert len(preset.steps) > 0

    @pytest.mark.parametrize("name", _PRESET_NAMES)
    def test_preset_steps_are_valid_tuples(self, name):
        preset = get_preset(name)
        for label, step in preset.steps:
            assert isinstance(label, str)
            assert isinstance(step, Step)

    def test_get_preset_invalid_raises_key_error(self):
        with pytest.raises(KeyError, match="Unknown preset"):
            get_preset("not_a_preset")

    def test_list_presets_returns_all(self):
        all_presets = list_presets()
        assert len(all_presets) == 11

    def test_preset_catalogue_is_non_empty_string(self):
        cat = preset_catalogue()
        assert isinstance(cat, str)
        assert "basic_qc" in cat
        assert "publication_ready" in cat

    def test_from_preset_pipeline(self):
        pipe = Pipeline.from_preset("basic_qc")
        assert isinstance(pipe, Pipeline)
        assert pipe.name == "basic_qc"
        assert len(pipe) == 5  # NR001, FREQ002, FREQ001, FREQ004, QC001


class TestMethodAwarePresets:
    """Structural checks for the mt_qc/amt_qc/csamt_qc/csumt_qc presets.

    Real end-to-end runs against MT/AMT/CSAMT data live in
    test_presets.py; these are fast, data-free structural assertions.
    """

    def _codes(self, name: str) -> list[str]:
        return [s.spec.code for _, s in get_preset(name).steps]

    def test_mt_qc_opens_and_closes_with_preview(self):
        codes = self._codes("mt_qc")
        assert codes[0] == "PRE001"
        assert codes[-1] == "PRE002"

    def test_all_method_presets_include_smart_qc_steps(self):
        for name in ("mt_qc", "amt_qc", "csamt_qc", "csumt_qc"):
            codes = self._codes(name)
            assert "QC005" in codes
            assert "QC006" in codes
            assert "QC007" in codes

    def test_mt_and_amt_qc_have_no_source_effects_steps(self):
        for name in ("mt_qc", "amt_qc"):
            codes = self._codes(name)
            assert "SRC001" not in codes
            assert "SRC002" not in codes
            assert "QC002" not in codes

    def test_csamt_and_csumt_qc_include_near_field_correction(self):
        for name in ("csamt_qc", "csumt_qc"):
            codes = self._codes(name)
            assert "SRC001" in codes
            assert "SRC002" in codes
            assert "QC002" in codes

    def test_csamt_near_field_step_defaults_to_auto_resolve(self):
        preset = get_preset("csamt_qc")
        src001 = next(s for _, s in preset.steps if s.spec.code == "SRC001")
        assert src001.params["source_offset"] is None

    def test_only_csumt_qc_includes_depth_section(self):
        assert "QC004" in self._codes("csumt_qc")
        for name in ("mt_qc", "amt_qc", "csamt_qc"):
            assert "QC004" not in self._codes(name)

    def test_amt_band_matches_stratagem_mt_band(self):
        amt = next(
            s for _, s in get_preset("amt_qc").steps if s.spec.code == "FREQ001"
        )
        strat = next(
            s
            for _, s in get_preset("stratagem_mt").steps
            if s.spec.code == "FREQ001"
        )
        assert amt.params["band_hz"] == strat.params["band_hz"]

    def test_csumt_band_uses_real_csumt_constants(self):
        from pycsamt.emtools.csumt import F_MAX_CSUMT, F_MIN_CSUMT

        csumt = next(
            s for _, s in get_preset("csumt_qc").steps if s.spec.code == "FREQ001"
        )
        assert csumt.params["band_hz"] == (F_MIN_CSUMT, F_MAX_CSUMT)

    def test_get_preset_for_method_maps_to_expected_presets(self):
        from pycsamt.pipeline import get_preset_for_method

        assert get_preset_for_method("MT").name == "mt_qc"
        assert get_preset_for_method("amt").name == "amt_qc"
        assert get_preset_for_method("CSAMT").name == "csamt_qc"
        assert get_preset_for_method("csumt").name == "csumt_qc"

    def test_get_preset_for_method_rejects_unknown_method(self):
        from pycsamt.pipeline import get_preset_for_method

        with pytest.raises(ValueError, match="not one of"):
            get_preset_for_method("bogus")

    def test_get_preset_for_method_rejects_method_without_preset(self):
        from pycsamt.pipeline import get_preset_for_method

        with pytest.raises(ValueError, match="no method-aware preset"):
            get_preset_for_method("CSEM")


# ─────────────────────────────────────────────────────────────────────────────
# 4 · PipelineAPIConfig
# ─────────────────────────────────────────────────────────────────────────────


class TestPipelineAPIConfig:
    def test_default_output_root(self):
        cfg = PipelineAPIConfig()
        assert cfg.output_root == "pipe_results"

    def test_default_processed_subdir(self):
        assert PYCSAMT_PIPE.processed_subdir == "processed"

    def test_default_plots_subdir(self):
        assert PYCSAMT_PIPE.plots_subdir == "plots"

    def test_configure_changes_value(self):
        configure_pipe(plot_dpi=300)
        assert PYCSAMT_PIPE.plot_dpi == 300

    def test_configure_invalid_key_raises(self):
        with pytest.raises(AttributeError, match="Unknown pipe config key"):
            configure_pipe(totally_unknown_key=True)

    def test_context_manager_restores_on_exit(self):
        original_dpi = PYCSAMT_PIPE.plot_dpi
        with PYCSAMT_PIPE.context(plot_dpi=999):
            assert PYCSAMT_PIPE.plot_dpi == 999
        assert PYCSAMT_PIPE.plot_dpi == original_dpi

    def test_context_manager_restores_on_exception(self):
        original_fmt = PYCSAMT_PIPE.plot_fmt
        try:
            with PYCSAMT_PIPE.context(plot_fmt="svg"):
                assert PYCSAMT_PIPE.plot_fmt == "svg"
                raise ValueError("boom")
        except ValueError:
            pass
        assert PYCSAMT_PIPE.plot_fmt == original_fmt

    def test_reset_restores_defaults(self):
        configure_pipe(plot_dpi=9999, plot_fmt="svg")
        reset_pipe()
        assert PYCSAMT_PIPE.plot_dpi == 150
        assert PYCSAMT_PIPE.plot_fmt == "png"

    def test_clone_is_independent(self):
        cloned = PYCSAMT_PIPE.clone()
        cloned.plot_dpi = 42
        assert PYCSAMT_PIPE.plot_dpi != 42

    def test_repr_includes_keys(self):
        r = repr(PYCSAMT_PIPE)
        assert "output_root" in r
        assert "plot_dpi" in r


# ─────────────────────────────────────────────────────────────────────────────
# 5 · OutputDir
# ─────────────────────────────────────────────────────────────────────────────


class TestOutputDir:
    def test_setup_creates_root(self, tmp_path):
        root = tmp_path / "pipe_out"
        od = OutputDir(root)
        od.setup()
        assert root.exists()

    def test_setup_creates_processed_and_plots(self, tmp_path):
        root = tmp_path / "pipe_out"
        od = OutputDir(root)
        od.setup()
        assert (root / "processed").exists()
        assert (root / "plots").exists()

    def test_step_plot_dir_creates_subdirectory(self, tmp_path):
        od = OutputDir(tmp_path / "out")
        od.setup()
        d = od.step_plot_dir(1, "notch_powerline")
        assert d.exists()
        assert d.name == "01_notch_powerline"

    def test_step_plot_dir_padding_two_digits(self, tmp_path):
        od = OutputDir(tmp_path / "out")
        od.setup()
        d = od.step_plot_dir(9, "step9")
        assert d.name.startswith("09_")

    def test_save_pipeline_config_writes_yaml(self, tmp_path):
        od = OutputDir(tmp_path / "out")
        od.setup()
        p = od.save_pipeline_config("name: test_pipe\nsteps: []\n")
        assert p is not None
        assert p.read_text() == "name: test_pipe\nsteps: []\n"

    def test_save_text_writes_file(self, tmp_path):
        od = OutputDir(tmp_path / "out")
        od.setup()
        p = od.save_text("hello\nworld\n", "summary.txt")
        assert p is not None
        assert "hello" in p.read_text()

    def test_repr_shows_status(self, tmp_path):
        od = OutputDir(tmp_path / "out")
        r = repr(od)
        assert "pending" in r
        od.setup()
        assert "ready" in repr(od)


# ─────────────────────────────────────────────────────────────────────────────
# 6 · Pipeline construction
# ─────────────────────────────────────────────────────────────────────────────


class TestPipelineConstruction:
    def test_from_step_list(self):
        pipe = Pipeline([_identity("NR001")])
        assert len(pipe) == 1
        assert pipe._steps[0][0] == "notch_powerline"  # auto-label from name

    def test_from_tuple_list(self):
        pipe = Pipeline([("my_notch", _identity("NR001"))])
        assert pipe._steps[0][0] == "my_notch"

    def test_empty_pipeline(self):
        pipe = Pipeline([])
        assert len(pipe) == 0

    def test_bad_input_raises_type_error(self):
        with pytest.raises(TypeError):
            Pipeline(["NR001"])  # strings are not Steps

    def test_append_adds_to_end(self):
        pipe = Pipeline([("notch", _identity("NR001"))])
        pipe.append("band", _identity("FREQ001"))
        assert len(pipe) == 2
        assert pipe._steps[-1][0] == "band"

    def test_append_returns_self_for_chaining(self):
        pipe = Pipeline([])
        result = pipe.append("a", _identity("NR001"))
        assert result is pipe

    def test_remove_found(self):
        pipe = Pipeline([("notch", _identity("NR001")), ("band", _identity("FREQ001"))])
        pipe.remove("notch")
        assert len(pipe) == 1
        assert pipe._steps[0][0] == "band"

    def test_remove_missing_raises_key_error(self):
        pipe = Pipeline([("notch", _identity("NR001"))])
        with pytest.raises(KeyError, match="notch_missing"):
            pipe.remove("notch_missing")

    def test_insert_at_position(self):
        pipe = Pipeline([("a", _identity("NR001")), ("c", _identity("FREQ001"))])
        pipe.insert(1, "b", _identity("FREQ004"))
        assert pipe._steps[1][0] == "b"
        assert len(pipe) == 3

    def test_replace_step(self):
        pipe = Pipeline([("notch", _identity("NR001"))])
        pipe.replace("notch", _identity("NR002"))
        assert pipe._steps[0][1].spec.code == "NR002"

    def test_replace_missing_raises_key_error(self):
        pipe = Pipeline([("notch", _identity("NR001"))])
        with pytest.raises(KeyError):
            pipe.replace("ghost", _identity("NR002"))

    def test_clone_is_independent(self):
        pipe = Pipeline([("notch", _identity("NR001"))])
        clone = pipe.clone()
        clone.append("extra", _identity("FREQ004"))
        assert len(pipe) == 1

    def test_len(self, simple_pipe):
        assert len(simple_pipe) == 3

    def test_getitem_by_label(self, simple_pipe):
        step = simple_pipe["notch"]
        assert step.spec.code == "NR001"

    def test_getitem_by_index(self, simple_pipe):
        label, step = simple_pipe[0]
        assert label == "notch"

    def test_getitem_missing_label_raises(self, simple_pipe):
        with pytest.raises(KeyError):
            simple_pipe["ghost_step"]

    def test_from_preset_returns_pipeline(self):
        pipe = Pipeline.from_preset("noise_reduction")
        assert isinstance(pipe, Pipeline)
        assert len(pipe) > 0

    def test_from_preset_codes_match(self):
        pipe = Pipeline.from_preset("tensor_analysis")
        codes = [s.spec.code for _, s in pipe._steps]
        assert "TZ001" in codes
        assert "TZ002" in codes

    def test_from_preset_invalid_raises(self):
        with pytest.raises(KeyError):
            Pipeline.from_preset("not_a_preset")

    def test_modify_during_run_raises(self, simple_pipe, sites):
        """Mutating a running pipeline must raise RuntimeError."""
        # We simulate `_running=True` directly
        simple_pipe._running = True
        with pytest.raises(RuntimeError, match="currently running"):
            simple_pipe.append("extra", _identity("FREQ002"))
        simple_pipe._running = False  # cleanup

    def test_iter_yields_tuples(self, simple_pipe):
        for item in simple_pipe:
            assert isinstance(item, tuple)
            label, step = item
            assert isinstance(label, str)
            assert isinstance(step, Step)


# ─────────────────────────────────────────────────────────────────────────────
# 7 · Pipeline repr and describe
# ─────────────────────────────────────────────────────────────────────────────


class TestPipelineRepr:
    def test_repr_contains_step_count(self, simple_pipe):
        r = repr(simple_pipe)
        assert "3 steps" in r

    def test_repr_contains_pipeline_name(self, simple_pipe):
        r = repr(simple_pipe)
        assert "test_pipe" in r

    def test_repr_contains_all_codes(self, simple_pipe):
        r = repr(simple_pipe)
        assert "NR001" in r
        assert "FREQ001" in r
        assert "FREQ004" in r

    def test_repr_html_is_string(self, simple_pipe):
        html = simple_pipe._repr_html_()
        assert isinstance(html, str)
        assert "<table" in html

    def test_repr_html_contains_codes(self, simple_pipe):
        html = simple_pipe._repr_html_()
        assert "NR001" in html

    def test_describe_returns_dataframe_or_str(self, simple_pipe):
        result = simple_pipe.describe()
        # pandas may or may not be available; accept both
        assert hasattr(result, "__len__") or isinstance(result, str)

    def test_repr_empty_pipeline(self):
        pipe = Pipeline([], name="empty")
        r = repr(pipe)
        assert "0 step" in r or "empty" in r


# ─────────────────────────────────────────────────────────────────────────────
# 8 · Config file loaders (YAML / JSON / .py)
# ─────────────────────────────────────────────────────────────────────────────

_YAML_CFG = """\
name: my_workflow
steps:
  - name: notch
    code: NR001
    params:
      mains_hz: 60
  - code: FREQ004
"""

_JSON_CFG = json.dumps(
    {
        "name": "json_wf",
        "steps": [
            {"name": "align", "code": "FREQ004"},
            {"name": "ss", "code": "SS001"},
        ],
    }
)

_PY_CFG = """\
pipeline_config = {
    "name": "py_wf",
    "steps": [
        {"name": "notch", "code": "NR001", "params": {"mains_hz": 50}},
        {"code": "FREQ002"},
    ],
}
"""

_YAML_WITH_PRESET = """\
name: extended
preset: basic_qc
steps:
  - name: extra_ss
    code: SS001
"""


class TestPipelineConfigFiles:
    def test_yaml_roundtrip_name(self, tmp_path):
        pipe = Pipeline([("notch", _identity("NR001"))], name="roundtrip_pipe")
        yaml_path = tmp_path / "pipe.yaml"
        pipe.to_yaml(yaml_path)
        loaded = Pipeline.from_yaml(yaml_path)
        assert loaded.name == "roundtrip_pipe"

    def test_yaml_roundtrip_step_count(self, tmp_path):
        pipe = Pipeline(
            [("a", _identity("NR001")), ("b", _identity("FREQ001"))],
            name="n2",
        )
        yaml_path = tmp_path / "pipe.yaml"
        pipe.to_yaml(yaml_path)
        loaded = Pipeline.from_yaml(yaml_path)
        assert len(loaded) == 2

    def test_yaml_roundtrip_codes(self, tmp_path):
        pipe = Pipeline(
            [("notch", _identity("NR001")), ("align", _identity("FREQ004"))],
            name="codes_test",
        )
        yaml_path = tmp_path / "pipe.yaml"
        pipe.to_yaml(yaml_path)
        loaded = Pipeline.from_yaml(yaml_path)
        codes = [s.spec.code for _, s in loaded._steps]
        assert codes == ["NR001", "FREQ004"]

    def test_yaml_param_preserved(self, tmp_path):
        yaml_path = tmp_path / "cfg.yaml"
        yaml_path.write_text(_YAML_CFG, encoding="utf-8")
        pipe = Pipeline.from_yaml(yaml_path)
        notch_step = pipe["notch"]
        assert notch_step.params["mains_hz"] == 60

    def test_yaml_with_preset_expands_steps(self, tmp_path):
        yaml_path = tmp_path / "preset.yaml"
        yaml_path.write_text(_YAML_WITH_PRESET, encoding="utf-8")
        pipe = Pipeline.from_yaml(yaml_path)
        codes = [s.spec.code for _, s in pipe._steps]
        # basic_qc has 5 steps + 1 extra SS001
        assert "SS001" in codes
        assert "NR001" in codes

    def test_json_roundtrip(self, tmp_path):
        json_path = tmp_path / "cfg.json"
        json_path.write_text(_JSON_CFG, encoding="utf-8")
        pipe = Pipeline.from_json(json_path)
        assert pipe.name == "json_wf"
        assert len(pipe) == 2

    def test_py_config_loads(self, tmp_path):
        py_path = tmp_path / "cfg.py"
        py_path.write_text(_PY_CFG, encoding="utf-8")
        pipe = Pipeline.from_py(py_path)
        assert pipe.name == "py_wf"
        assert len(pipe) == 2

    def test_py_config_missing_variable_raises(self, tmp_path):
        py_path = tmp_path / "bad.py"
        py_path.write_text("x = 1\n", encoding="utf-8")
        with pytest.raises(AttributeError, match="pipeline_config"):
            Pipeline.from_py(py_path)

    def test_py_config_wrong_type_raises(self, tmp_path):
        py_path = tmp_path / "bad2.py"
        py_path.write_text("pipeline_config = [1, 2, 3]\n", encoding="utf-8")
        with pytest.raises(TypeError):
            Pipeline.from_py(py_path)

    def test_to_yaml_string_is_loadable(self, tmp_path):
        pipe = Pipeline(
            [("a", _identity("FREQ001", fmin=(0.001, 10000)))],
            name="tuple_test",
        )
        yaml_str = pipe.to_yaml_string()
        # tuple converted to list — no !!python/tuple tags
        assert "!!python" not in yaml_str

    def test_to_yaml_string_roundtrip_with_tuples(self, tmp_path):
        """fmin tuple must survive yaml → reload without errors."""
        pipe = Pipeline(
            [("band", _identity("FREQ001", fmin=(0.001, 10000)))],
            name="tup",
        )
        yaml_path = tmp_path / "tup.yaml"
        pipe.to_yaml(yaml_path)
        pipe2 = Pipeline.from_yaml(yaml_path)
        step = pipe2["band"]
        assert step.params["fmin"] == [0.001, 10000]  # list after reload


# ─────────────────────────────────────────────────────────────────────────────
# 9 · Pipeline.run — mechanics
# ─────────────────────────────────────────────────────────────────────────────


class TestPipelineRun:
    def test_run_returns_pipeline_result(self, simple_pipe, sites):
        result = simple_pipe.run(sites, save_edis=False, save_report=False)
        assert isinstance(result, PipelineResult)

    def test_run_sites_in_is_input(self, simple_pipe, sites):
        result = simple_pipe.run(sites, save_edis=False, save_report=False)
        assert result.sites_in is sites

    def test_run_step_results_count_matches(self, simple_pipe, sites):
        result = simple_pipe.run(sites, save_edis=False, save_report=False)
        assert len(result.step_results) == len(simple_pipe)

    def test_run_elapsed_is_positive(self, simple_pipe, sites):
        result = simple_pipe.run(sites, save_edis=False, save_report=False)
        assert result.elapsed_sec > 0

    def test_run_ok_when_no_errors(self, simple_pipe, sites):
        result = simple_pipe.run(sites, save_edis=False, save_report=False)
        assert result.ok

    def test_run_n_errors_zero_when_clean(self, simple_pipe, sites):
        result = simple_pipe.run(sites, save_edis=False, save_report=False)
        assert result.n_errors == 0

    def test_run_step_results_are_ok(self, simple_pipe, sites):
        result = simple_pipe.run(sites, save_edis=False, save_report=False)
        for sr in result.step_results:
            assert sr.ok

    def test_run_step_results_have_elapsed(self, simple_pipe, sites):
        result = simple_pipe.run(sites, save_edis=False, save_report=False)
        for sr in result.step_results:
            assert sr.elapsed_sec >= 0

    def test_run_step_results_codes_match(self, simple_pipe, sites):
        result = simple_pipe.run(sites, save_edis=False, save_report=False)
        expected = ["NR001", "FREQ001", "FREQ004"]
        actual = [sr.step_code for sr in result.step_results]
        assert actual == expected

    def test_run_chaining_increments_count(self):
        """Each _MutatingStep increases sites.count by 1 — proves chaining."""
        pipe = Pipeline(
            [
                ("a", _mutating("NR001")),
                ("b", _mutating("FREQ001")),
                ("c", _mutating("FREQ004")),
            ]
        )
        start = _CountingSites(n=2, count=0)
        result = pipe.run(start, save_edis=False, save_report=False)
        assert result.sites_out.count == 3

    def test_run_pipeline_name_in_result(self, simple_pipe, sites):
        result = simple_pipe.run(sites, save_edis=False, save_report=False)
        assert result.pipeline_name == "test_pipe"

    def test_run_no_outdir_does_not_write_files(self, simple_pipe, sites):
        result = simple_pipe.run(sites, outdir=None, save_edis=False, save_report=False)
        assert result.outdir is None
        assert len(result.processed_paths) == 0

    # ── Error handling ──────────────────────────────────────────────────

    def test_run_error_warn_continues(self, sites):
        """on_step_error='warn': error recorded, pipeline continues."""
        configure_pipe(on_step_error="warn")
        pipe = Pipeline(
            [
                ("good_before", _identity("NR001")),
                ("bad_step", _raising("FREQ001")),
                ("good_after", _identity("FREQ004")),
            ]
        )
        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter("always")
            result = pipe.run(sites, save_edis=False, save_report=False)

        assert len(result.step_results) == 3
        assert result.step_results[1].error is not None
        assert result.step_results[0].ok
        assert result.step_results[2].ok
        assert result.n_errors == 1
        assert not result.ok
        assert any("deliberate step failure" in str(w.message) for w in ws)

    def test_run_error_skip_silent(self, sites):
        configure_pipe(on_step_error="skip")
        pipe = Pipeline(
            [
                ("bad", _raising("NR001")),
                ("ok", _identity("FREQ001")),
            ]
        )
        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter("always")
            result = pipe.run(sites, save_edis=False, save_report=False)

        assert result.step_results[0].error is not None
        # no user-visible warning with "skip"
        user_warns = [w for w in ws if "deliberate step failure" in str(w.message)]
        assert len(user_warns) == 0

    def test_run_error_raise_propagates(self, sites):
        configure_pipe(on_step_error="raise")
        pipe = Pipeline(
            [
                ("bad", _raising("NR001")),
                ("ok", _identity("FREQ001")),
            ]
        )
        with pytest.raises(RuntimeError, match="deliberate step failure"):
            pipe.run(sites, save_edis=False, save_report=False)

    def test_run_with_outdir_creates_tree(self, simple_pipe, sites, tmp_path):
        root = tmp_path / "out"
        simple_pipe.run(sites, outdir=root, save_edis=False, save_report=False)
        assert root.exists()
        assert (root / "processed").exists()
        assert (root / "plots").exists()

    def test_run_with_outdir_saves_yaml(self, simple_pipe, sites, tmp_path):
        root = tmp_path / "out"
        simple_pipe.run(sites, outdir=root, save_edis=False, save_report=False)
        assert (root / "pipeline.yaml").exists()

    def test_run_with_outdir_saves_reports(self, simple_pipe, sites, tmp_path):
        root = tmp_path / "out"
        simple_pipe.run(sites, outdir=root, save_edis=False, save_report=True)
        assert (root / "summary.txt").exists()
        assert (root / "report.html").exists()

    def test_run_plots_dir_has_step_subdirs(self, simple_pipe, sites, tmp_path):
        root = tmp_path / "out"
        simple_pipe.run(
            sites,
            outdir=root,
            save_plots=True,
            save_edis=False,
            save_report=False,
        )
        # _IdentityStep.generate_qc_plots returns [] so no figures, but
        # the run should still have created the plots/ tree structure.
        plots_dir = root / "plots"
        assert plots_dir.exists()

    def test_run_step_results_indices_are_1_based(self, simple_pipe, sites):
        result = simple_pipe.run(sites, save_edis=False, save_report=False)
        for i, sr in enumerate(result.step_results, start=1):
            assert sr.step_idx == i


# ─────────────────────────────────────────────────────────────────────────────
# 9b · Pipeline.run(cache=...)
# ─────────────────────────────────────────────────────────────────────────────


class _CountingMutatingStep(_MutatingStep):
    """_MutatingStep that also records how many times transform() actually ran."""

    def transform(self, sites):
        self._calls.append(1)
        return super().transform(sites)


def _counting_mutating(code: str, calls: list) -> _CountingMutatingStep:
    """Same bypass-__init__ construction _mutating() uses, plus a call counter."""
    step = _CountingMutatingStep.__new__(_CountingMutatingStep)
    step.spec = lookup_step(code)
    step.params = dict(step.spec.defaults)
    step._calls = calls
    return step


class TestPipelineCaching:
    def test_disabled_by_default_no_behavior_change(self, simple_pipe, sites):
        result = simple_pipe.run(sites, save_edis=False, save_report=False)
        assert all(not sr.cached for sr in result.step_results)

    def test_cache_false_is_identical_to_omitting_it(self, simple_pipe, sites):
        result = simple_pipe.run(
            sites, save_edis=False, save_report=False, cache=False
        )
        assert all(not sr.cached for sr in result.step_results)

    def test_second_run_hits_cache_and_skips_transform(self, tmp_path):
        calls: list[int] = []
        pipe = Pipeline(
            [
                ("a", _counting_mutating("NR001", calls)),
                ("b", _counting_mutating("FREQ001", calls)),
            ]
        )
        start = _CountingSites(n=3, count=0)

        r1 = pipe.run(start, outdir=None, save_report=False, cache=tmp_path)
        assert len(calls) == 2
        assert all(not sr.cached for sr in r1.step_results)

        r2 = pipe.run(start, outdir=None, save_report=False, cache=tmp_path)
        assert len(calls) == 2  # no new transform() calls on the second run
        assert all(sr.cached for sr in r2.step_results)
        assert r2.sites_out.count == r1.sites_out.count == 2

    def test_changed_upstream_param_invalidates_downstream_cache(self, tmp_path):
        calls: list[int] = []
        start = _CountingSites(n=3, count=0)

        pipe_a = Pipeline([("a", _counting_mutating("NR001", calls))])
        pipe_a.run(start, outdir=None, save_report=False, cache=tmp_path)
        assert len(calls) == 1

        # A different step at the same position -> different upstream state
        # for anything chained after it; a *second*, distinct pipeline command
        # naturally gets a fresh cache key chain, proving there's no false hit.
        pipe_b = Pipeline([("a", _counting_mutating("FREQ001", calls))])
        pipe_b.run(start, outdir=None, save_report=False, cache=tmp_path)
        assert len(calls) == 2  # ran again — different step code, different key

    def test_resume_after_simulated_crash(self, tmp_path):
        """A pipeline that raises on step 2 the first time, then succeeds on
        a rerun, must not recompute step 1 the second time — this is the
        "resumability" property this whole cache mechanism exists to give."""
        calls: list[int] = []
        start = _CountingSites(n=3, count=0)

        attempt_1 = Pipeline(
            [
                ("a", _counting_mutating("NR001", calls)),
                ("b", _raising("FREQ001")),
            ]
        )
        configure_pipe(on_step_error="raise")
        with pytest.raises(RuntimeError, match="deliberate step failure"):
            attempt_1.run(start, outdir=None, save_report=False, cache=tmp_path)
        assert len(calls) == 1  # step "a" completed and was cached before the crash

        # Rerun with step 2 fixed this time.
        attempt_2 = Pipeline(
            [
                ("a", _counting_mutating("NR001", calls)),
                ("b", _counting_mutating("FREQ001", calls)),
            ]
        )
        result = attempt_2.run(start, outdir=None, save_report=False, cache=tmp_path)
        assert len(calls) == 2  # step "a" replayed from cache, NOT recomputed
        assert result.step_results[0].cached is True
        assert result.step_results[1].cached is False
        assert result.ok

    def test_diagnostic_step_cache_hit_skips_recompute(self, tmp_path):
        calls: list[int] = []

        class _CountingDiagnosticStep(Step):
            def transform(self, sites):
                calls.append(1)
                return sites  # diagnostic: unchanged, but still "ran"

            def generate_qc_plots(self, sites):
                return []

        pipe = Pipeline([("qc", _CountingDiagnosticStep("QC001"))])
        start = _CountingSites(n=3, count=0)

        pipe.run(start, outdir=None, save_report=False, cache=tmp_path)
        assert len(calls) == 1

        pipe.run(start, outdir=None, save_report=False, cache=tmp_path)
        assert len(calls) == 1  # second call skipped entirely


# ─────────────────────────────────────────────────────────────────────────────
# 10 · PipelineResult
# ─────────────────────────────────────────────────────────────────────────────


class TestPipelineResult:
    @pytest.fixture()
    def clean_result(self, simple_pipe, sites):
        return simple_pipe.run(sites, outdir=None, save_edis=False, save_report=False)

    def test_plots_aggregates_all_step_plots(self, clean_result):
        # _IdentityStep returns no plots so this should be empty
        assert isinstance(clean_result.plots, list)

    def test_ok_true_when_no_errors(self, clean_result):
        assert clean_result.ok is True

    def test_n_errors_zero(self, clean_result):
        assert clean_result.n_errors == 0

    def test_summary_returns_string(self, clean_result):
        s = clean_result.summary()
        assert isinstance(s, str)
        assert "test_pipe" in s

    def test_repr_returns_string(self, clean_result):
        r = repr(clean_result)
        assert "PipelineResult" in r
        assert "elapsed" in r

    def test_processed_paths_empty_without_edi_write(self, clean_result):
        assert clean_result.processed_paths == []

    def test_outdir_none_without_outdir(self, clean_result):
        assert clean_result.outdir is None


# ─────────────────────────────────────────────────────────────────────────────
# 11 · Report generators
# ─────────────────────────────────────────────────────────────────────────────


def _make_step_result(idx: int, *, ok: bool = True, cached: bool = False) -> StepResult:
    return StepResult(
        step_idx=idx,
        step_name=f"step_{idx}",
        step_code=f"NR00{idx}",
        step_label=f"Step {idx} Label",
        params={"k": idx},
        elapsed_sec=0.1 * idx,
        plots=[],
        n_sites_in=5,
        n_sites_out=5 if ok else 0,
        error=None if ok else RuntimeError("test error"),
        cached=cached,
    )


class TestReportGenerators:
    def test_text_report_is_string(self):
        sr = [_make_step_result(i) for i in range(1, 4)]
        txt = make_text_report("my_pipe", sr, 1.23, None, 5, 5)
        assert isinstance(txt, str)

    def test_text_report_contains_pipeline_name(self):
        sr = [_make_step_result(1)]
        txt = make_text_report("willy_pipe", sr, 0.5, None, 5, 5)
        assert "willy_pipe" in txt

    def test_text_report_shows_step_codes(self):
        sr = [_make_step_result(1)]
        txt = make_text_report("pipe", sr, 0.5, None, 5, 5)
        assert "NR001" in txt

    def test_text_report_shows_error_section_when_errors(self):
        sr = [_make_step_result(1, ok=False)]
        txt = make_text_report("pipe", sr, 0.5, None, 5, 0)
        assert "Errors" in txt or "ERR" in txt

    def test_html_report_is_valid_html(self):
        sr = [_make_step_result(i) for i in range(1, 3)]
        html = make_html_report("html_pipe", sr, 1.0, None, 5, 5)
        assert "<html" in html
        assert "</html>" in html

    def test_html_report_contains_step_names(self):
        sr = [_make_step_result(1)]
        html = make_html_report("pipe", sr, 0.1, None, 5, 5)
        assert "step_1" in html

    def test_html_report_shows_error_span_on_bad_step(self):
        sr = [_make_step_result(1, ok=False)]
        html = make_html_report("pipe", sr, 0.1, None, 5, 0)
        assert 'class="err"' in html

    def test_text_report_shows_cached_marker(self):
        sr = [_make_step_result(1, cached=True), _make_step_result(2, cached=False)]
        txt = make_text_report("pipe", sr, 0.5, None, 5, 5)
        lines = txt.splitlines()
        cached_line = next(l for l in lines if "step_1" in l)
        uncached_line = next(l for l in lines if "step_2" in l)
        assert "yes" in cached_line
        assert "yes" not in uncached_line

    def test_html_report_shows_cached_badge(self):
        sr = [_make_step_result(1, cached=True)]
        html = make_html_report("pipe", sr, 0.1, None, 5, 5)
        assert "<b>Cached:</b> yes" in html

    def test_html_report_omits_cached_badge_when_not_cached(self):
        sr = [_make_step_result(1, cached=False)]
        html = make_html_report("pipe", sr, 0.1, None, 5, 5)
        assert "<b>Cached:</b>" not in html

    def test_html_report_uses_real_pycsamt_brand_hex(self):
        # _CSS was swapped from generic blue/green/red to the real
        # pyCSAMT brand tokens (docs/source/_static/css/custom.css).
        sr = [_make_step_result(1)]
        html = make_html_report("pipe", sr, 0.1, None, 5, 5)
        assert "#3e65b0" in html
        assert "#4a6fa5" not in html

    def test_html_report_still_renders_step_cards_after_palette_swap(self):
        sr = [_make_step_result(i) for i in range(1, 3)]
        html = make_html_report("pipe", sr, 0.1, None, 5, 5)
        assert 'class="step-card"' in html
        assert html.count('class="step-card"') == 2
