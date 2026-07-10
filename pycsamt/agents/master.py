# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""One-line front door to the pyCSAMT agent stack.

:class:`AgentMaster` is a thin façade over
:class:`~pycsamt.agents.WorkflowOrchestratorAgent`: one object, one
``run()`` call from a plain-language request to an executed agent
chain.  It exists so the friendly entry point advertised across the
docs works verbatim::

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
    >>> plan = master.plan("QC the EDI files and prepare a short report",
    ...                    data_path="data/edi/")
    >>> plan["workflow_type"]          # doctest: +SKIP
    'qc'
    >>> report = master.run(
    ...     "Load data/edi/, flag stations with RMS > 2, build an Occam2D "
    ...     "input for profile L22, launch inversion, and produce a PDF "
    ...     "report."
    ... )                              # doctest: +SKIP

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
        """Classify *request* into a workflow and execute the chain.

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
            Preview the selected chain without reading or writing.
        **extra
            Additional orchestrator payload fields, passed through.

        Returns
        -------
        AgentResult
            Status, per-step outputs, reasoning, and cost tracking.
        """
        payload: dict[str, Any] = {
            "request": str(request),
            "dry_run": bool(dry_run),
            **extra,
        }
        if data_path is not None:
            payload["data_path"] = str(data_path)
        if output_dir is not None:
            payload["output_dir"] = str(output_dir)
        return self.orchestrator.execute(payload)

    def plan(self, request: str, **kwargs: Any) -> AgentResult:
        """Shortcut for :meth:`run` with ``dry_run=True``."""
        kwargs["dry_run"] = True
        return self.run(request, **kwargs)
