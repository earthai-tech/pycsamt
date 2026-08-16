"""Tests for the opt-in AI001 (audit_survey) pipeline step.

Public entry point: :mod:`pycsamt.pipeline` (``register_ai_steps``)
Implementation:     :mod:`pycsamt.pipeline.ai_steps`

Resolving this step imports ``pycsamt.ai`` (and, transitively, ``torch`` —
see the module docstring in ``pycsamt/pipeline/ai_steps.py``), so the
real-execution tests are skipped when no DL backend is importable, mirroring
the skip guard already used in ``pycsamt/ai/tests/test_inversion.py``.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from pycsamt.pipeline import (
    STEP_REGISTRY,
    Step,
    lookup_step,
    register_ai_steps,
    unregister_step,
)


def _has_backend() -> bool:
    from pycsamt.backends import get_backend

    return get_backend() != "none"


_PROJECT_ROOT = Path(__file__).resolve().parents[3]
_DATA_DIR = _PROJECT_ROOT / "data" / "AMT" / "WILLY_DATA" / "L22PLT"


@pytest.fixture(autouse=True)
def _cleanup_ai001():
    unregister_step("AI001", missing_ok=True)
    yield
    unregister_step("AI001", missing_ok=True)


class TestRegisterAiSteps:
    def test_not_registered_by_default(self):
        assert "AI001" not in STEP_REGISTRY

    def test_register_ai_steps_adds_ai001(self):
        registered = register_ai_steps()
        assert [s.code for s in registered] == ["AI001"]
        spec = lookup_step("AI001")
        assert spec.name == "audit_survey"
        assert spec.category == "ai"
        assert spec.origin == "plugin"
        assert spec.returns_sites is False

    def test_replace_existing_true_allows_recall(self):
        register_ai_steps()
        register_ai_steps(replace_existing=True)  # must not raise
        assert lookup_step("AI001").code == "AI001"

    def test_second_call_without_replace_raises(self):
        register_ai_steps()
        with pytest.raises(ValueError, match="already registered"):
            register_ai_steps()

    def test_pipeline_init_never_calls_register_ai_steps(self):
        """Guard the opt-in guarantee statically: importing pycsamt.pipeline
        alone must never force the torch import register_ai_steps() causes.
        A dynamic sys.modules check would be unreliable here since other
        tests/collection-time probes in the same process may have already
        imported torch for unrelated reasons (e.g. the backend skip-guard
        above) — the invariant that actually matters is that __init__.py's
        source only imports the name, never calls it."""
        import inspect

        import pycsamt.pipeline as pl

        source = inspect.getsource(pl)
        assert "register_ai_steps()" not in source


@pytest.mark.skipif(
    not _DATA_DIR.exists(), reason=f"WILLY EDI data not found at {_DATA_DIR}"
)
@pytest.mark.skipif(not _has_backend(), reason="no DL backend available")
class TestAuditSurveyStepExecution:
    def test_real_audit_run_against_willy_data(self, capsys):
        from pycsamt.emtools._core import ensure_sites

        register_ai_steps()
        paths = sorted(_DATA_DIR.glob("*.edi"))
        sites = ensure_sites([str(p) for p in paths])

        step = Step("AI001")
        result = step.transform(sites)

        assert result is sites  # diagnostic step: sites pass through unchanged
        printed = capsys.readouterr().out
        assert "stations" in printed.lower() or "frequency" in printed.lower()

    def test_report_path_writes_json(self, tmp_path):
        from pycsamt.emtools._core import ensure_sites

        register_ai_steps()
        paths = sorted(_DATA_DIR.glob("*.edi"))
        sites = ensure_sites([str(p) for p in paths])

        report_path = tmp_path / "audit.json"
        Step("AI001", report_path=str(report_path)).transform(sites)

        assert report_path.exists()
