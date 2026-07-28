# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""One-line front door to the pyCSAMT agent stack.

:class:`AgentMaster` is the programmatic master router: one object, one
``run()`` call from a plain-language request to the right specialist
agent.  ``run()`` first classifies the request with
:class:`~pycsamt.agents.router.IntentRouter` — the same classifier the
desktop app uses — and dispatches *question* requests to
:class:`~pycsamt.agents.package_qa.PackageQAAgent`, *code* requests to
:class:`~pycsamt.agents.code_gen.CodeGenerationAgent`, *metrics*
requests to :class:`~pycsamt.agents.metrics.MetricsAgent`, and
*workflow*/*plot* requests to
:class:`~pycsamt.agents.orchestrator.WorkflowOrchestratorAgent`.  So the
friendly entry point advertised across the docs works verbatim::

    from pycsamt.agents import AgentMaster

    master = AgentMaster(provider="anthropic")
    report = master.run(
        "Load data/edi/, flag stations with RMS > 2, build an Occam2D "
        "input for profile L22, launch inversion, and produce a PDF report."
    )

Everything heavier (the orchestrator, LLM clients) is imported lazily
on first use, so ``from pycsamt.agents import AgentMaster`` stays cheap
and free of circular imports.
"""

from __future__ import annotations

import time
from typing import TYPE_CHECKING, Any

from ..api.property import PyCSAMTObject

if TYPE_CHECKING:  # pragma: no cover - typing only
    from ._base import AgentResult
    from .orchestrator import WorkflowOrchestratorAgent

__all__ = ["AgentMaster"]

# Friendly names accepted in addition to the canonical BaseAgent
# providers ({"claude", "openai", "gemini", "deepseek", "minimax"}).
_PROVIDER_ALIASES = {
    "anthropic": "claude",
    "claude": "claude",
    "openai": "openai",
    "gpt": "openai",
    "gemini": "gemini",
    "google": "gemini",
    "deepseek": "deepseek",
    "minimax": "minimax",
}

_CAPABILITY_TEXT = (
    "I'm the pyCSAMT v2 agent stack. Ask a question about the package, "
    "request a standalone script, ask for a computed survey value, or "
    "describe a processing/inversion workflow to run — for example: "
    "'What does StaticShiftAgent do?', 'write me a script to load EDIs "
    "and run QC', 'what's the mean RMS for L22PLT?', or 'load data/edi/, "
    "run QC, and prepare an Occam2D input for profile L22'."
)


class AgentMaster(PyCSAMTObject):
    """Plain-language entry point to the agent workflows.

    Parameters
    ----------
    provider : str, default ``"claude"``
        LLM provider.  Friendly aliases are accepted:
        ``"anthropic"`` → ``"claude"`` and ``"google"`` → ``"gemini"``.
        Without an API key the agents fall back to the rule-based
        (regex/keyword) path, so ``AgentMaster()`` works offline.
    api_key : str, optional
        Provider API key.  When omitted, the provider's environment
        variable is used if set; otherwise the rule-based fallback
        runs at zero cost.
    model : str, optional
        Provider model override (defaults per provider).
    default_workflow : str, default ``"qc"``
        Workflow used when a request cannot be classified.

    Examples
    --------
    Plan first (no files touched), then execute:

    >>> from pycsamt.agents import AgentMaster
    >>> master = AgentMaster(provider="anthropic")
    >>> plan = master.plan(
    ...     "QC the EDI files and prepare a short report",
    ...     data_path="data/edi/",
    ... )
    >>> plan["workflow_type"]  # doctest: +SKIP
    'qc'
    >>> report = master.run(
    ...     "Load data/edi/, flag stations with RMS > 2, build an Occam2D "
    ...     "input for profile L22, launch inversion, and produce a PDF "
    ...     "report."
    ... )  # doctest: +SKIP

    See Also
    --------
    pycsamt.agents.WorkflowOrchestratorAgent :
        The dispatcher this façade drives; use it directly for
        structured workflow configurations.
    """

    __repr_fields__ = ("provider", "default_workflow")

    def __init__(
        self,
        provider: str = "claude",
        *,
        api_key: str | None = None,
        model: str | None = None,
        default_workflow: str = "qc",
    ) -> None:
        key = str(provider).strip().lower()
        if key not in _PROVIDER_ALIASES:
            raise ValueError(
                "provider must be one of "
                f"{sorted(set(_PROVIDER_ALIASES))}, got {provider!r}"
            )
        self.provider = _PROVIDER_ALIASES[key]
        self.default_workflow = default_workflow
        self._api_key = api_key
        self._model = model
        self._orchestrator: WorkflowOrchestratorAgent | None = None

    # ------------------------------------------------------------------
    @property
    def orchestrator(self) -> WorkflowOrchestratorAgent:
        """The lazily-built orchestrator behind this façade."""
        if self._orchestrator is None:
            # Local import: keeps module import light and cycle-free.
            from .orchestrator import (
                WorkflowOrchestratorAgent,
            )

            self._orchestrator = WorkflowOrchestratorAgent(
                api_key=self._api_key,
                model=self._model,
                llm_provider=self.provider,
                default_workflow=self.default_workflow,
            )
        return self._orchestrator

    # ------------------------------------------------------------------
    def run(
        self,
        request: str,
        *,
        data_path: str | None = None,
        output_dir: str | None = None,
        dry_run: bool = False,
        **extra: Any,
    ) -> AgentResult:
        """Route *request* to the right specialist agent and run it.

        ``request`` is first classified by
        :class:`~pycsamt.agents.router.IntentRouter` — a *question* is
        answered by :class:`~pycsamt.agents.package_qa.PackageQAAgent`,
        a *code* request generates a script via
        :class:`~pycsamt.agents.code_gen.CodeGenerationAgent`, a
        *metrics* request is computed by
        :class:`~pycsamt.agents.metrics.MetricsAgent`, and a
        *workflow* / *plot* request runs the full pipeline through
        :class:`~pycsamt.agents.orchestrator.WorkflowOrchestratorAgent`.
        A *meta* (capability) or *clarify* (ambiguous) request returns
        immediately with no data path required.

        Parameters
        ----------
        request : str
            Plain-language description of what to do.  Paths mentioned
            in the text are extracted when possible; pass ``data_path``
            / ``output_dir`` explicitly for scripts and CI.
        data_path : str, optional
            Survey input (EDI/AVG/J directory or file).
        output_dir : str, optional
            Where products (figures, inputs, reports) are written.
        dry_run : bool, default False
            For *workflow* / *plot* requests only: preview the selected
            chain without reading or writing. Ignored for questions,
            code, metrics, meta, and clarify requests, which never
            touch disk regardless.
        **extra
            Additional orchestrator payload fields, passed through for
            *workflow* / *plot* requests only.

        Returns
        -------
        AgentResult
            Status, per-step outputs, reasoning, and cost tracking.
        """
        from .router import (
            CLARIFY,
            CODE,
            META,
            METRICS,
            QUESTION,
            IntentRouter,
        )

        text = str(request)
        t0 = time.time()
        router = IntentRouter(
            api_key=self._api_key,
            model=self._model,
            llm_provider=self.provider,
        )
        decision = router.route(text)

        if decision.intent == META:
            return self._agent_result(
                status="success",
                summary="Capability summary.",
                data={"intent": META, "decision": decision},
                llm_interpretation=_CAPABILITY_TEXT,
                elapsed=time.time() - t0,
            )
        if decision.intent == CLARIFY:
            return self._agent_result(
                status="needs_review",
                summary=(
                    decision.clarification
                    or "Could you clarify what you'd like me to do — "
                    "answer a question, generate code, or run a "
                    "workflow?"
                ),
                data={"intent": CLARIFY, "decision": decision},
                elapsed=time.time() - t0,
            )
        if decision.intent == QUESTION:
            from .package_qa import PackageQAAgent

            qa = PackageQAAgent(
                api_key=self._api_key,
                model=self._model,
                llm_provider=self.provider,
            )
            result = qa.execute({"question": text})
            result.data.setdefault("intent", QUESTION)
            return result
        if decision.intent == CODE:
            return self._dispatch_code(text, data_path, output_dir)
        if decision.intent == METRICS:
            return self._dispatch_metrics(text, data_path)

        # WORKFLOW / PLOT (and any future intent this router version
        # doesn't classify above) → run the pipeline as before.
        payload: dict[str, Any] = {
            "request": text,
            "dry_run": bool(dry_run),
            **extra,
        }
        if data_path is not None:
            payload["data_path"] = str(data_path)
        if output_dir is not None:
            payload["output_dir"] = str(output_dir)
        result = self.orchestrator.execute(payload)
        result.data.setdefault("intent", decision.intent)
        return result

    def plan(self, request: str, **kwargs: Any) -> AgentResult:
        """Shortcut for :meth:`run` with ``dry_run=True``."""
        kwargs["dry_run"] = True
        return self.run(request, **kwargs)

    # ------------------------------------------------------------------
    # intent-specific dispatch helpers
    # ------------------------------------------------------------------

    def _dispatch_code(
        self,
        text: str,
        data_path: str | None,
        output_dir: str | None,
    ) -> AgentResult:
        """Generate a standalone script for *text* via CodeGenerationAgent."""
        from .code_gen import CodeGenerationAgent
        from .context import ContextInputAgent

        ctx = ContextInputAgent(
            api_key=self._api_key,
            model=self._model,
            llm_provider=self.provider,
        )
        ctx_result = ctx.execute({"request": text})
        cfg = ctx_result.data.get("config", {})
        if data_path is not None:
            cfg["data_path"] = str(data_path)

        cg = CodeGenerationAgent(
            api_key=self._api_key,
            model=self._model,
            llm_provider=self.provider,
        )
        result = cg.execute(
            {
                "workflow_config": cfg,
                "results": {},
                "output_dir": output_dir or "pycsamt_agent_output",
            }
        )
        result.data.setdefault("intent", "code")
        return result

    def _dispatch_metrics(
        self,
        text: str,
        data_path: str | None,
    ) -> AgentResult:
        """Compute the requested survey value(s) via MetricsAgent."""
        from .metrics import MetricsAgent, parse_metric_request

        kinds, _all_lines = parse_metric_request(text)
        ma = MetricsAgent()
        result = ma.execute(
            {"kinds": kinds or ["summary"], "data_path": data_path or ""}
        )
        result.data.setdefault("intent", "metrics")
        return result

    @staticmethod
    def _agent_result(
        *,
        status: str,
        summary: str,
        data: dict[str, Any],
        elapsed: float,
        llm_interpretation: str | None = None,
    ) -> AgentResult:
        from ._base import AgentResult

        return AgentResult(
            status=status,
            summary=summary,
            data=data,
            llm_interpretation=llm_interpretation,
            elapsed_seconds=elapsed,
        )
