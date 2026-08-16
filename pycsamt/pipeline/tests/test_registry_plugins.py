"""Tests for the pipeline plugin/extensibility API.

Public entry point: :mod:`pycsamt.pipeline`
Implementation:     :mod:`pycsamt.pipeline._registry` (register_step /
                     unregister_step), :mod:`pycsamt.pipeline._plugins`
                     (discover_plugins)

These tests never touch the 47 built-in steps — see
:class:`~pycsamt.pipeline.tests.test_pipeline.TestStepRegistry` for the
built-in registry invariants, which must stay unaffected by anything here.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Any, Callable

import pytest

from pycsamt.pipeline import (
    STEP_REGISTRY,
    PluginLoadResult,
    StepSpec,
    discover_plugins,
    lookup_step,
    register_step,
    unregister_step,
)
from pycsamt.pipeline import _plugins as plugins_mod


def _demo_fn(sites: Any, **kw: Any) -> Any:
    return sites


def _other_demo_fn(sites: Any, **kw: Any) -> Any:
    return sites


@pytest.fixture(autouse=True)
def _cleanup_registry():
    """Remove any test-registered step, even if a test fails mid-way."""
    before = set(STEP_REGISTRY)
    yield
    for code in set(STEP_REGISTRY) - before:
        unregister_step(code, missing_ok=True)


def _demo_spec(code: str = "DEMO001", name: str = "demo_step") -> StepSpec:
    return StepSpec(
        code=code,
        name=name,
        label="Demo Plugin Step",
        category="plugin_demo",
        override_fn=_demo_fn,
    )


# ─────────────────────────────────────────────────────────────────────────────
# register_step / unregister_step
# ─────────────────────────────────────────────────────────────────────────────


class TestRegisterStep:
    def test_register_new_step_is_visible_everywhere(self):
        registered = register_step(_demo_spec())

        assert registered.origin == "plugin"
        assert lookup_step("DEMO001") is registered
        assert lookup_step("demo_step") is registered
        assert "DEMO001" in STEP_REGISTRY
        from pycsamt.pipeline import categories, list_steps  # noqa: PLC0415

        assert "plugin_demo" in categories()
        assert registered in list_steps("plugin_demo")

    def test_origin_is_always_stamped_plugin(self):
        spec = StepSpec(
            code="DEMO002",
            name="demo_step_2",
            label="Demo 2",
            category="plugin_demo",
            override_fn=_demo_fn,
            origin="builtin",  # caller tries to lie — must be overwritten
        )
        registered = register_step(spec)
        assert registered.origin == "plugin"

    def test_duplicate_code_without_replace_raises(self):
        register_step(_demo_spec())
        with pytest.raises(ValueError, match="already registered"):
            register_step(_demo_spec())

    def test_duplicate_name_without_replace_raises(self):
        register_step(_demo_spec(code="DEMO001", name="demo_step"))
        with pytest.raises(ValueError, match="already registered"):
            register_step(_demo_spec(code="DEMO003", name="demo_step"))

    def test_colliding_with_builtin_code_without_replace_raises(self):
        with pytest.raises(ValueError, match="already registered"):
            register_step(_demo_spec(code="NR001", name="not_notch_powerline"))

    def test_replace_existing_true_overwrites_same_entry(self):
        register_step(_demo_spec())
        updated = register_step(
            StepSpec(
                code="DEMO001",
                name="demo_step",
                label="Demo Plugin Step v2",
                category="plugin_demo",
                override_fn=_other_demo_fn,
            ),
            replace_existing=True,
        )
        assert lookup_step("DEMO001").label == "Demo Plugin Step v2"
        assert updated.override_fn is _other_demo_fn

    def test_replace_existing_true_can_patch_a_builtin(self):
        original = lookup_step("NR001")
        try:
            register_step(
                StepSpec(
                    code="NR001",
                    name="notch_powerline",
                    label="Patched notch",
                    category="noise_removal",
                    override_fn=_demo_fn,
                ),
                replace_existing=True,
            )
            patched = lookup_step("NR001")
            assert patched.label == "Patched notch"
            assert patched.origin == "plugin"
        finally:
            STEP_REGISTRY["NR001"] = original
            from pycsamt.pipeline._registry import (  # noqa: PLC0415
                _NAME_INDEX,
            )

            _NAME_INDEX["notch_powerline"] = original

    def test_replace_existing_true_but_code_and_name_point_elsewhere_raises(self):
        register_step(_demo_spec(code="DEMO001", name="demo_step"))
        # DEMO001 exists (name demo_step); NR001's name is notch_powerline.
        # Trying to register code=DEMO001 with name=notch_powerline is a
        # genuine conflict even with replace_existing=True.
        with pytest.raises(ValueError, match="already used by a different step"):
            register_step(
                StepSpec(
                    code="DEMO001",
                    name="notch_powerline",
                    label="Conflicting",
                    category="plugin_demo",
                    override_fn=_demo_fn,
                ),
                replace_existing=True,
            )

    def test_validate_true_rejects_unresolvable_module(self):
        bad = StepSpec(
            code="DEMO_BAD",
            name="demo_bad",
            label="Bad",
            category="plugin_demo",
            mod="pycsamt.pipeline._does_not_exist_xyz",
            fn_name="nope",
        )
        with pytest.raises((ModuleNotFoundError, AttributeError, RuntimeError)):
            register_step(bad)
        # Registry must be left untouched on failure.
        assert "DEMO_BAD" not in STEP_REGISTRY
        with pytest.raises(KeyError):
            lookup_step("demo_bad")

    def test_validate_false_skips_resolution_check(self):
        bad = StepSpec(
            code="DEMO_BAD",
            name="demo_bad",
            label="Bad",
            category="plugin_demo",
            mod="pycsamt.pipeline._does_not_exist_xyz",
            fn_name="nope",
        )
        registered = register_step(bad, validate=False)
        assert registered.code == "DEMO_BAD"


class TestUnregisterStep:
    def test_unregister_removes_code_and_name(self):
        register_step(_demo_spec())
        unregister_step("DEMO001")
        assert "DEMO001" not in STEP_REGISTRY
        with pytest.raises(KeyError):
            lookup_step("demo_step")

    def test_unregister_by_name_also_works(self):
        register_step(_demo_spec())
        unregister_step("demo_step")
        assert "DEMO001" not in STEP_REGISTRY

    def test_unregister_unknown_raises_key_error(self):
        with pytest.raises(KeyError):
            unregister_step("NOT_A_REAL_CODE")

    def test_unregister_unknown_missing_ok_noops(self):
        unregister_step("NOT_A_REAL_CODE", missing_ok=True)  # must not raise

    def test_reregister_after_unregister(self):
        register_step(_demo_spec())
        unregister_step("DEMO001")
        registered = register_step(_demo_spec())  # must not raise ValueError
        assert registered.code == "DEMO001"


# ─────────────────────────────────────────────────────────────────────────────
# discover_plugins
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class _FakeEntryPoint:
    name: str
    _target: Callable[[], None]

    def load(self) -> Callable[[], None]:
        return self._target


def _good_plugin_register() -> None:
    register_step(_demo_spec())


def _broken_plugin_register() -> None:
    raise RuntimeError("deliberately broken plugin")


class TestDiscoverPlugins:
    def test_discovers_and_runs_a_well_behaved_plugin(self, monkeypatch):
        monkeypatch.setattr(
            plugins_mod,
            "_iter_entry_points",
            lambda group: [_FakeEntryPoint("good", _good_plugin_register)],
        )
        results = discover_plugins()
        assert results == [PluginLoadResult(name="good", ok=True, error=None)]
        assert lookup_step("DEMO001").origin == "plugin"

    def test_broken_plugin_does_not_block_the_others_by_default(self, monkeypatch):
        monkeypatch.setattr(
            plugins_mod,
            "_iter_entry_points",
            lambda group: [
                _FakeEntryPoint("broken", _broken_plugin_register),
                _FakeEntryPoint("good", _good_plugin_register),
            ],
        )
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            results = discover_plugins(on_error="warn")

        assert [r.ok for r in results] == [False, True]
        assert results[0].name == "broken"
        assert "deliberately broken plugin" in results[0].error
        assert lookup_step("DEMO001").origin == "plugin"
        assert any("broken" in str(w.message) for w in caught)

    def test_strict_mode_raises_on_first_failure(self, monkeypatch):
        monkeypatch.setattr(
            plugins_mod,
            "_iter_entry_points",
            lambda group: [_FakeEntryPoint("broken", _broken_plugin_register)],
        )
        with pytest.raises(RuntimeError, match="deliberately broken plugin"):
            discover_plugins(on_error="raise")

    def test_invalid_on_error_value_rejected(self):
        with pytest.raises(ValueError, match="on_error"):
            discover_plugins(on_error="ignore")

    def test_no_entry_points_returns_empty_list(self, monkeypatch):
        monkeypatch.setattr(plugins_mod, "_iter_entry_points", lambda group: [])
        assert discover_plugins() == []

    def test_pipeline_init_never_calls_discover_plugins(self):
        """Guard the opt-in guarantee: pycsamt/pipeline/__init__.py must
        only *import* discover_plugins, never call it — plugin code must
        not run just because `import pycsamt.pipeline` ran."""
        import inspect

        import pycsamt.pipeline as pl

        source = inspect.getsource(pl)
        assert "discover_plugins(" not in source
