"""
pycsamt.agents._base
====================

Abstract foundation for all pycsamt agents.

Every agent:
  - Inherits :class:`BaseAgent` and implements :meth:`execute`.
  - Returns an :class:`AgentResult` dataclass.
  - Has access to pycsamt's full API (PYCSAMT_SECTION, PYCSAMT_STYLE,
    PYCSAMT_STATION_RENDERING, PYCSAMT_CONTROL, PLOT_CONFIG) so all
    figures produced by agents are indistinguishable from hand-crafted
    pycsamt plots.
  - Supports five LLM providers: Anthropic Claude (default), OpenAI,
    Google Gemini, DeepSeek, MiniMax.  When no API key is supplied every
    LLM-dependent method degrades gracefully.
"""

from __future__ import annotations

import json
import logging
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any

from ..api.agents import AGENT_CONFIG
from ..api.control import PYCSAMT_CONTROL
from ..api.plot import PLOT_CONFIG, add_colorbar

# ── pycsamt API singletons ────────────────────────────────────────────────────
from ..api.section import PYCSAMT_SECTION, SectionStyle
from ..api.station import PYCSAMT_STATION_RENDERING
from ..api.style import PYCSAMT_STYLE

logger = logging.getLogger(__name__)

# ── constants ─────────────────────────────────────────────────────────────────
_PROVIDERS = {"claude", "openai", "gemini", "deepseek", "minimax"}
_STATUS    = {"success", "failed", "needs_review"}

_DEFAULT_MODELS = {
    "claude":   "claude-sonnet-4-6",
    "openai":   "gpt-4o",
    "gemini":   "gemini-2.0-flash",
    "deepseek": "deepseek-chat",
    "minimax":  "MiniMax-M3",
}

_RETRY_DELAYS = (1.0, 2.0, 4.0)   # seconds between LLM retries


# ═══════════════════════════════════════════════════════════════════════════════
# AgentResult
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class AgentResult:
    """Standardised output returned by every pycsamt agent.

    Attributes
    ----------
    status : str
        ``"success"`` | ``"failed"`` | ``"needs_review"``
    summary : str
        One-sentence human-readable description of what happened.
    data : dict
        Agent-specific outputs (arrays, paths, dataframes, figures …).
    warnings : list of str
        Non-fatal issues encountered during execution.
    llm_interpretation : str or None
        Free-text interpretation written by the LLM, when available.
    elapsed_seconds : float
    cost_estimate_usd : float
        Estimated LLM API cost for this execution.
    error : str or None
        Exception message if ``status == "failed"``.
    error_fix_hint : str or None
        Suggested remediation for the error.

    Examples
    --------
    >>> result = some_agent.execute({"path": "/data/EDIs"})
    >>> result.status
    'success'
    >>> result["sites"]          # dict-like access to data
    <Sites 25 stations>
    >>> result.get("n_stations", 0)
    25
    """

    status: str
    summary: str
    data: dict[str, Any]          = field(default_factory=dict)
    warnings: list[str]           = field(default_factory=list)
    llm_interpretation: str | None = None
    elapsed_seconds: float        = 0.0
    cost_estimate_usd: float      = 0.0
    error: str | None             = None
    error_fix_hint: str | None    = None

    # dict-like helpers so callers can do result["key"] or result.get("key")
    def __getitem__(self, key: str) -> Any:
        return self.data[key]

    def __contains__(self, key: str) -> bool:
        return key in self.data

    def get(self, key: str, default: Any = None) -> Any:
        return self.data.get(key, default)

    def __bool__(self) -> bool:
        return self.status != "failed"

    def __repr__(self) -> str:
        warn_str = f", {len(self.warnings)} warnings" if self.warnings else ""
        cost_str = f", ${self.cost_estimate_usd:.4f}" if self.cost_estimate_usd else ""
        return (
            f"AgentResult(status={self.status!r}, "
            f"elapsed={self.elapsed_seconds:.1f}s{cost_str}{warn_str})"
        )

    @classmethod
    def failed(
        cls,
        error: str,
        *,
        hint: str | None = None,
        elapsed: float = 0.0,
    ) -> AgentResult:
        """Convenience constructor for failure results."""
        return cls(
            status="failed",
            summary=f"Agent failed: {error}",
            error=error,
            error_fix_hint=hint,
            elapsed_seconds=elapsed,
        )


# ═══════════════════════════════════════════════════════════════════════════════
# BaseAgent
# ═══════════════════════════════════════════════════════════════════════════════

class BaseAgent(ABC):
    """Abstract base class for all pycsamt agents.

    Parameters
    ----------
    name : str
        Human-readable agent name used in logs and reports.
    api_key : str or None
        LLM API key.  When ``None`` the agent runs without LLM support and
        every ``llm_interpretation`` field will be ``None``.
    model : str or None
        LLM model identifier.  Defaults to the provider's recommended model.
    llm_provider : {"claude", "openai", "gemini", "deepseek", "minimax"}
        Which provider to use.  Default ``"claude"``.
    section_preset : str
        Which :data:`~pycsamt.api.section.PYCSAMT_SECTION` preset governs
        figures produced by this agent.  Default ``"pseudosection"``.

    Examples
    --------
    Subclass and implement :meth:`execute`::

        class MyAgent(BaseAgent):
            SYSTEM_PROMPT = "You are an expert MT data analyst."

            def execute(self, input_data):
                t0 = time.time()
                # ... do work ...
                interp = self.query_llm("Interpret these results: ...")
                return AgentResult(
                    status="success",
                    summary="Analysis complete.",
                    data={"result": 42},
                    llm_interpretation=interp,
                    elapsed_seconds=time.time() - t0,
                    cost_estimate_usd=self._last_cost,
                )
    """

    #: Override in subclasses to give the LLM its domain expertise.
    SYSTEM_PROMPT: str = (
        "You are a geophysics expert specialising in "
        "magnetotelluric (MT/AMT/CSAMT) data processing and interpretation."
    )

    def __init__(
        self,
        name: str,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        section_preset: str = "pseudosection",
    ) -> None:
        llm_provider = llm_provider.lower()
        if llm_provider not in _PROVIDERS:
            raise ValueError(
                f"llm_provider must be one of {sorted(_PROVIDERS)}, "
                f"got {llm_provider!r}."
            )

        # Resolve effective (provider, key, model) — explicit args win;
        # AGENT_CONFIG fills in anything that was left as None/default.
        llm_provider, api_key, model = AGENT_CONFIG.resolve(
            llm_provider, api_key, model
        )

        self.name         = name
        self.api_key      = api_key
        self.llm_provider = llm_provider
        self.model        = model or _DEFAULT_MODELS[llm_provider]

        # ── pycsamt API injection ──────────────────────────────────────────
        self._section : SectionStyle = PYCSAMT_SECTION.style_for(section_preset)
        self._style                  = PYCSAMT_STYLE
        self._station                = PYCSAMT_STATION_RENDERING
        self._control                = PYCSAMT_CONTROL
        self._plot_cfg               = PLOT_CONFIG
        self._add_colorbar           = add_colorbar

        # cost accumulator reset each execute() call
        self._last_cost: float = 0.0

        self._log = logging.getLogger(f"pycsamt.agents.{name}")

    # ── abstract interface ────────────────────────────────────────────────────

    @abstractmethod
    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        """Run this agent on *input_data* and return an :class:`AgentResult`.

        Subclasses must implement this method.  The contract:

        * Reset ``self._last_cost = 0.0`` at the top.
        * Record wall-clock time with ``t0 = time.time()``.
        * Return ``AgentResult(elapsed_seconds=time.time()-t0,
          cost_estimate_usd=self._last_cost, ...)``.
        """

    # ── LLM interface ─────────────────────────────────────────────────────────

    def query_llm(
        self,
        prompt: str,
        system_message: str | None = None,
        *,
        temperature: float = 0.2,
        max_tokens: int = 1024,
    ) -> str | None:
        """Send *prompt* to the configured LLM and return the response text.

        Returns ``None`` when no API key is configured or all retries fail.
        Accumulates token cost into ``self._last_cost``.

        Parameters
        ----------
        prompt : str
        system_message : str or None
            Overrides :attr:`SYSTEM_PROMPT` for this call.
        temperature : float
        max_tokens : int

        Returns
        -------
        str or None
        """
        if not self.api_key:
            self._log.debug("No API key — LLM query skipped.")
            return None

        # Raise before the API call if the session budget is already exhausted.
        AGENT_CONFIG._check_budget()

        sys_msg = system_message or self.SYSTEM_PROMPT
        dispatch = {
            "claude":   self._query_claude,
            "openai":   self._query_openai,
            "gemini":   self._query_gemini,
            "deepseek": self._query_deepseek,
            "minimax":  self._query_minimax,
        }
        fn = dispatch[self.llm_provider]

        last_exc: Exception | None = None
        for delay in (*_RETRY_DELAYS, None):
            try:
                text, cost = fn(prompt, sys_msg, temperature, max_tokens)
                self._last_cost += cost
                AGENT_CONFIG._add_spend(cost)   # update session counter
                return text
            except Exception as exc:  # noqa: BLE001
                last_exc = exc
                is_rate = "rate" in str(exc).lower() or "429" in str(exc)
                if delay is None or not is_rate:
                    break
                self._log.warning("LLM rate limit, retrying in %ss…", delay)
                time.sleep(delay)

        self._log.error("LLM query failed: %s", last_exc)
        return None

    # ── provider backends ─────────────────────────────────────────────────────

    def _query_claude(
        self,
        prompt: str,
        system_message: str,
        temperature: float,
        max_tokens: int,
    ) -> tuple[str, float]:
        """Call Anthropic Claude and return (response_text, cost_usd)."""
        import anthropic  # lazy import

        client = anthropic.Anthropic(api_key=self.api_key)
        msg = client.messages.create(
            model=self.model,
            max_tokens=max_tokens,
            temperature=temperature,
            system=system_message,
            messages=[{"role": "user", "content": prompt}],
        )
        text = msg.content[0].text
        cost = AGENT_CONFIG.estimate_cost(
            self.llm_provider, self.model,
            msg.usage.input_tokens, msg.usage.output_tokens,
        )
        return text, cost

    def _query_openai(
        self,
        prompt: str,
        system_message: str,
        temperature: float,
        max_tokens: int,
    ) -> tuple[str, float]:
        """Call OpenAI and return (response_text, cost_usd)."""
        import openai  # lazy import

        client = openai.OpenAI(api_key=self.api_key)
        resp = client.chat.completions.create(
            model=self.model,
            max_tokens=max_tokens,
            temperature=temperature,
            messages=[
                {"role": "system", "content": system_message},
                {"role": "user",   "content": prompt},
            ],
        )
        text = resp.choices[0].message.content
        cost = AGENT_CONFIG.estimate_cost(
            self.llm_provider, self.model,
            resp.usage.prompt_tokens, resp.usage.completion_tokens,
        )
        return text, cost

    def _query_gemini(
        self,
        prompt: str,
        system_message: str,
        temperature: float,
        max_tokens: int,
    ) -> tuple[str, float]:
        """Call Google Gemini and return (response_text, cost_usd)."""
        import google.generativeai as genai  # lazy import

        genai.configure(api_key=self.api_key)
        gen_cfg = genai.GenerationConfig(
            temperature=temperature,
            max_output_tokens=max_tokens,
        )
        gemini_model = genai.GenerativeModel(
            self.model,
            system_instruction=system_message,
            generation_config=gen_cfg,
        )
        resp = gemini_model.generate_content(prompt)
        text = resp.text
        # Gemini doesn't expose exact token counts via the basic API; estimate
        n_in  = len(prompt.split()) * 4 // 3
        n_out = len(text.split())   * 4 // 3
        cost  = AGENT_CONFIG.estimate_cost(self.llm_provider, self.model, n_in, n_out)
        return text, cost

    def _query_deepseek(
        self,
        prompt: str,
        system_message: str,
        temperature: float,
        max_tokens: int,
    ) -> tuple[str, float]:
        r"""Call DeepSeek (OpenAI-compatible) and return (text, cost_usd).

        DeepSeek exposes an OpenAI-compatible REST API at
        ``https://api.deepseek.com``.  The ``openai`` Python package is
        reused with a custom ``base_url``.
        """
        import openai  # lazy import

        client = openai.OpenAI(
            api_key=self.api_key,
            base_url="https://api.deepseek.com",
        )
        resp = client.chat.completions.create(
            model=self.model,
            max_tokens=max_tokens,
            temperature=temperature,
            messages=[
                {
                    "role": "system",
                    "content": system_message,
                },
                {"role": "user", "content": prompt},
            ],
        )
        text = resp.choices[0].message.content
        cost = AGENT_CONFIG.estimate_cost(
            self.llm_provider,
            self.model,
            resp.usage.prompt_tokens,
            resp.usage.completion_tokens,
        )
        return text, cost

    def _query_minimax(
        self,
        prompt: str,
        system_message: str,
        temperature: float,
        max_tokens: int,
    ) -> tuple[str, float]:
        r"""Call MiniMax (OpenAI-compatible) and return (text, cost_usd).

        MiniMax exposes an OpenAI-compatible REST API at
        ``https://api.minimax.io/v1``.  The ``openai`` Python package is
        reused with a custom ``base_url``.
        """
        import openai  # lazy import

        client = openai.OpenAI(
            api_key=self.api_key,
            base_url="https://api.minimax.io/v1",
        )
        resp = client.chat.completions.create(
            model=self.model,
            max_tokens=max_tokens,
            temperature=temperature,
            messages=[
                {"role": "system", "content": system_message},
                {"role": "user",   "content": prompt},
            ],
        )
        text = resp.choices[0].message.content
        cost = AGENT_CONFIG.estimate_cost(
            self.llm_provider,
            self.model,
            resp.usage.prompt_tokens,
            resp.usage.completion_tokens,
        )
        return text, cost

    # ── JSON extraction helper ────────────────────────────────────────────────

    def extract_json(self, text: str) -> dict | list | None:
        """Extract the first JSON object / array from *text*.

        Returns ``None`` when no valid JSON is found.
        """
        if not text:
            return None
        # try whole string first
        try:
            return json.loads(text)
        except json.JSONDecodeError:
            pass
        # locate first { or [
        for start_char, end_char in (("{", "}"), ("[", "]")):
            start = text.find(start_char)
            if start == -1:
                continue
            # scan for balanced end
            depth = 0
            for i, ch in enumerate(text[start:], start):
                if ch == start_char:
                    depth += 1
                elif ch == end_char:
                    depth -= 1
                    if depth == 0:
                        try:
                            return json.loads(text[start:i + 1])
                        except json.JSONDecodeError:
                            break
        return None

    # ── input validation helpers ──────────────────────────────────────────────

    def require_keys(
        self,
        input_data: dict,
        *keys: str,
        agent_name: str | None = None,
    ) -> list[str]:
        """Return a list of missing required keys."""
        missing = [k for k in keys if k not in input_data]
        if missing:
            who = agent_name or self.name
            self._log.error("%s: missing required keys %s", who, missing)
        return missing

    # ── figure saving helper ──────────────────────────────────────────────────

    def _save_figure(
        self,
        fig: Any,
        output_dir: str | None,
        name: str,
        *,
        warnings_list: list[str] | None = None,
    ) -> str | None:
        """Save *fig* to ``<output_dir>/<name>`` using :attr:`_plot_cfg`.

        Returns the saved path string, or ``None`` when *output_dir* is falsy
        or saving fails.  Appends a warning string to *warnings_list* on error.
        """
        if not output_dir:
            return None
        import os
        os.makedirs(output_dir, exist_ok=True)
        try:
            paths = self._plot_cfg.save(fig, os.path.join(output_dir, name))
            return str(paths[0]) if paths else None
        except Exception as exc:  # noqa: BLE001
            msg = f"Could not save figure {name!r}: {exc}"
            self._log.warning(msg)
            if warnings_list is not None:
                warnings_list.append(msg)
            return None

    # ── repr ──────────────────────────────────────────────────────────────────

    def __repr__(self) -> str:
        llm = f"{self.llm_provider}/{self.model}" if self.api_key else "no-LLM"
        return f"{type(self).__name__}(name={self.name!r}, llm={llm!r})"


__all__ = ["AgentResult", "BaseAgent"]
