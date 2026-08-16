"""Opt-in AI-backed pipeline step: the domain-gap survey audit.

Registering this step is deliberately **not automatic**. Resolving it forces
a real ``torch`` import (``pycsamt.ai``'s package ``__init__`` eagerly pulls
in ``pycsamt.ai.nets.drcnn``, which imports ``torch`` at module level so its
classes are picklable — see that module's own comment), so it must not be a
cost every pipeline user pays just for ``import pycsamt.pipeline``. Call
:func:`register_ai_steps` explicitly, or pass ``pycsamt pipe --with-ai-steps``
on the CLI, exactly the same opt-in shape as :func:`~pycsamt.pipeline.discover_plugins`.
"""

from __future__ import annotations

from typing import Any

__all__ = ["qc_audit_survey", "register_ai_steps"]


def qc_audit_survey(sites: Any, **kw: Any) -> None:
    """Run :func:`pycsamt.ai.domain_gap.audit.audit_survey` as a QC step.

    Diagnostic only: prints the report summary and, when ``report_path`` is
    supplied, writes the full report as JSON. ``sites`` is never modified.

    Parameters
    ----------
    sites:
        Anything accepted by ``ensure_sites`` (path, EDI collection, Sites).
    report_path:
        Optional destination for :meth:`SurveyAuditReport.write_json`.
    **kw:
        Forwarded to ``audit_survey`` (``recursive``, ``on_dup``, ``verbose``,
        ``freq_rtol``, ``band``, ``skew_th``, ``ellipt_th``,
        ``station_spacing_fallback``, ``metadata``).
    """
    from pycsamt.ai.domain_gap.audit import audit_survey

    report_path = kw.pop("report_path", None)
    report = audit_survey(sites, **kw)
    print(report.summary())
    if report_path:
        report.write_json(report_path)


def register_ai_steps(*, replace_existing: bool = False) -> list:
    """Register the opt-in AI step(s) into the pipeline :term:`step registry`.

    Currently just ``AI001`` / ``audit_survey``. Never called automatically —
    see the module docstring for why.

    Parameters
    ----------
    replace_existing:
        Forwarded to :func:`~pycsamt.pipeline.register_step`.

    Returns
    -------
    list[StepSpec]
        The registered spec(s), each stamped ``origin="plugin"`` by
        ``register_step``.
    """
    from ._registry import StepSpec, register_step

    spec = StepSpec(
        code="AI001",
        name="audit_survey",
        label="AI Domain-Gap Survey Audit",
        category="ai",
        returns_sites=False,
        mod="pycsamt.pipeline.ai_steps",
        fn_name="qc_audit_survey",
        defaults={},
    )
    return [register_step(spec, replace_existing=replace_existing)]
