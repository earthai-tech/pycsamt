"""
pycsamt.api.agents
==================

Global LLM provider configuration for all pycsamt agents.

Configure once — every agent inherits the active provider, API key, model,
pricing overrides, and budget cap automatically.

Quick start
-----------
::

    from pycsamt.api.agents import AGENT_CONFIG

    AGENT_CONFIG.configure(provider="claude", api_key="sk-ant-…")

    from pycsamt.agents import DataQCAgent, PhaseAnalysisAgent

    qc = DataQCAgent()  # inherits provider + key automatically
    pt = PhaseAnalysisAgent()  # same

Or use the top-level convenience function::

    from pycsamt.agents import configure_agents

    configure_agents(provider="openai", api_key="sk-…")

Multi-provider workflow
-----------------------
Store keys for several providers upfront, then switch freely::

    AGENT_CONFIG.set_key("claude", "sk-ant-…")
    AGENT_CONFIG.set_key("openai", "sk-…")
    AGENT_CONFIG.set_key("gemini", "AIza…")

    AGENT_CONFIG.switch("claude")  # all agents use Claude
    # … some work …
    AGENT_CONFIG.switch("openai")  # now all agents use OpenAI, no re-supply

Temporary override (context manager)
-------------------------------------
::

    with AGENT_CONFIG.using(provider="gemini", api_key="AIza…"):
        r = DataQCAgent().execute({"path": "/data"})
    # original config is automatically restored here

Pricing — custom rates
-----------------------
Override any rate or add a model not in the built-in table::

    # update a rate that changed (USD per 1 M tokens)
    AGENT_CONFIG.set_rate("claude", "claude-opus-4-8", input=12.0, output=60.0)

    # register a model that isn't in the table yet
    AGENT_CONFIG.set_rate("openai", "gpt-5-turbo", input=5.0, output=20.0)

    # inspect the resolved rate for any provider/model
    print(AGENT_CONFIG.get_rate("claude", "claude-sonnet-4-6"))
    # → {"input": 3.0, "output": 15.0}

Pricing — session budget cap
------------------------------
::

    AGENT_CONFIG.set_budget(usd=2.0)  # raise BudgetExceededError above $2

    AGENT_CONFIG.spent_usd  # float — accumulated so far
    AGENT_CONFIG.remaining_usd  # float or None (None = no cap)

    AGENT_CONFIG.reset_budget()  # zero the spend counter; keeps the cap
    AGENT_CONFIG.reset_budget(cap=True)  # also remove the cap

Environment-variable fallback
------------------------------
Keys are resolved in this order:

1. Explicit key stored via :meth:`AgentConfig.configure` /
   :meth:`AgentConfig.set_key`.
2. Environment variable (first match wins per provider):

   =============== ===========================================================
   Provider        Environment variables checked
   =============== ===========================================================
   ``"claude"``    ``ANTHROPIC_API_KEY``, ``PYCSAMT_CLAUDE_API_KEY``
   ``"openai"``    ``OPENAI_API_KEY``, ``PYCSAMT_OPENAI_API_KEY``
   ``"gemini"``    ``GEMINI_API_KEY``, ``GOOGLE_API_KEY``,
                   ``GOOGLE_GENERATIVEAI_API_KEY``, ``PYCSAMT_GEMINI_API_KEY``
   ``"deepseek"``  ``DEEPSEEK_API_KEY``, ``PYCSAMT_DEEPSEEK_API_KEY``
   ``"minimax"``   ``MINIMAX_API_KEY``, ``MINIMAX_M3``, ``PYCSAMT_MINIMAX_API_KEY``
   =============== ===========================================================

Keys are also loaded automatically from ``.env.local`` at the repo root
(without overwriting variables already set in the OS environment).

3. ``None`` — agent runs without LLM support (regex/rule-based fallback).
"""

from __future__ import annotations

import os
import threading
from collections.abc import Generator
from contextlib import contextmanager
from pathlib import Path

# Thread-local flag: when True, _resolve_key()
# skips env-var lookup so offline mode is
# truly offline even when API keys exist in
# os.environ (e.g. loaded from .env.local).
_TLS = threading.local()


# ── Auto-load .env.local from repo root ───────────────────────────────────────
# Populates os.environ before any key resolution so that
# ANTHROPIC_API_KEY / OPENAI_API_KEY / DEEPSEEK_API_KEY etc.
# are available without requiring python-dotenv.
def _load_env_local() -> None:
    here = Path(__file__).resolve()
    # Walk up to find repo root (contains .env.local)
    for parent in here.parents:
        env_file = parent / ".env.local"
        if env_file.exists():
            for raw in env_file.read_text().splitlines():
                line = raw.strip()
                if not line or line.startswith("#") or "=" not in line:
                    continue
                k, _, v = line.partition("=")
                k = k.strip()
                v = v.strip()
                # never overwrite a key already set
                # in the real environment
                if k and v and k not in os.environ:
                    os.environ[k] = v
            break


_load_env_local()

__all__ = [
    "AgentConfig",
    "AGENT_CONFIG",
    "BudgetExceededError",
    "configure_agents",
    "reset_agents",
]

# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

_PROVIDERS: frozenset[str] = frozenset(
    {"claude", "openai", "gemini", "deepseek", "minimax"}
)

_DEFAULT_MODELS: dict[str, str] = {
    "claude": "claude-sonnet-4-6",
    "openai": "gpt-4o",
    "gemini": "gemini-2.0-flash",
    "deepseek": "deepseek-chat",
    "minimax": "MiniMax-M3",
}

_ENV_KEYS: dict[str, list[str]] = {
    "claude": [
        "ANTHROPIC_API_KEY",
        "PYCSAMT_CLAUDE_API_KEY",
    ],
    "openai": [
        "OPENAI_API_KEY",
        "PYCSAMT_OPENAI_API_KEY",
    ],
    "gemini": [
        "GEMINI_API_KEY",
        "GOOGLE_API_KEY",
        "GOOGLE_GENERATIVEAI_API_KEY",
        "PYCSAMT_GEMINI_API_KEY",
    ],
    "deepseek": [
        "DEEPSEEK_API_KEY",
        "PYCSAMT_DEEPSEEK_API_KEY",
    ],
    "minimax": [
        "MINIMAX_API_KEY",
        "MINIMAX_M3",
        "PYCSAMT_MINIMAX_API_KEY",
    ],
}

# Built-in rates (USD per 1 M tokens).  AgentConfig._custom_rates takes
# precedence; these are used only as fallback.
_BUILTIN_RATES: dict[str, dict[str, dict[str, float]]] = {
    "claude": {
        "claude-opus-4-8": {"input": 15.00, "output": 75.00},
        "claude-opus-4-7": {"input": 15.00, "output": 75.00},
        "claude-opus-4-6": {"input": 15.00, "output": 75.00},
        "claude-sonnet-4-6": {"input": 3.00, "output": 15.00},
        "claude-haiku-4-5": {"input": 0.80, "output": 4.00},
        "claude-haiku-4-5-20251001": {"input": 0.80, "output": 4.00},
        "claude-3-opus-20240229": {"input": 15.00, "output": 75.00},
        "claude-3-5-sonnet-20241022": {"input": 3.00, "output": 15.00},
        "claude-3-haiku-20240307": {"input": 0.25, "output": 1.25},
    },
    "openai": {
        "gpt-4o": {"input": 2.50, "output": 10.00},
        "gpt-4o-mini": {"input": 0.15, "output": 0.60},
        "gpt-4-turbo": {"input": 10.00, "output": 30.00},
        "gpt-4": {"input": 30.00, "output": 60.00},
        "gpt-3.5-turbo": {"input": 0.50, "output": 1.50},
    },
    "gemini": {
        "gemini-1.5-pro": {"input": 3.50, "output": 10.50},
        "gemini-1.5-flash": {"input": 0.075, "output": 0.30},
        "gemini-1.5-flash-8b": {"input": 0.0375, "output": 0.15},
        "gemini-2.0-flash": {"input": 0.10, "output": 0.40},
    },
    "deepseek": {
        # DeepSeek-V3 / deepseek-chat
        "deepseek-chat": {"input": 0.27, "output": 1.10},
        # DeepSeek-R1 (reasoning model)
        "deepseek-reasoner": {"input": 0.55, "output": 2.19},
    },
    "minimax": {
        # MiniMax-M3 (<=512K input tokens)
        "MiniMax-M3": {"input": 0.30, "output": 1.20},
    },
}

# Provider-level defaults when the exact model is not found anywhere
_PROVIDER_DEFAULTS: dict[str, dict[str, float]] = {
    "claude": {"input": 3.00, "output": 15.00},
    "openai": {"input": 2.50, "output": 10.00},
    "gemini": {"input": 3.50, "output": 10.50},
    "deepseek": {"input": 0.27, "output": 1.10},
    "minimax": {"input": 0.30, "output": 1.20},
}


# ---------------------------------------------------------------------------
# BudgetExceededError
# ---------------------------------------------------------------------------


class BudgetExceededError(RuntimeError):
    """Raised when an LLM call would be made after the session budget is used.

    Attributes
    ----------
    spent_usd : float
        Accumulated spend so far this session.
    budget_usd : float
        The cap that was set.
    """

    def __init__(self, spent: float, budget: float) -> None:
        self.spent_usd = spent
        self.budget_usd = budget
        super().__init__(
            f"Session budget of ${budget:.4f} exceeded "
            f"(spent so far: ${spent:.4f}). "
            f"Call AGENT_CONFIG.reset_budget() to clear the counter "
            f"or AGENT_CONFIG.set_budget(usd=<new_limit>) to raise the cap."
        )


# ---------------------------------------------------------------------------
# AgentConfig
# ---------------------------------------------------------------------------


class AgentConfig:
    """Global LLM configuration singleton for all pycsamt agents.

    :class:`~pycsamt.agents.BaseAgent` calls :meth:`resolve` on every
    instantiation and :meth:`get_rate` + :meth:`_add_spend` on every LLM
    call, so keys, rates, and the budget cap set here apply automatically to
    every agent in the session.

    Parameters
    ----------
    (none — use :meth:`configure` after construction)

    Attributes
    ----------
    provider : str or None
        Currently active LLM provider.
    model : str or None
        Resolved model (explicit override or provider default).
    api_key : str or None
        Resolved key for the active provider (stored key → env var → None).
    is_configured : bool
        ``True`` when a provider *and* a resolvable key are both present.
    spent_usd : float
        Accumulated LLM spend this session.
    remaining_usd : float or None
        Remaining budget, or ``None`` when no cap is set.

    Examples
    --------
    ::

        from pycsamt.api.agents import AGENT_CONFIG

        # one-call setup
        AGENT_CONFIG.configure(provider="claude", api_key="sk-ant-…")

        # pricing — override a rate (USD / 1 M tokens)
        AGENT_CONFIG.set_rate(
            "claude", "claude-opus-4-8", input=12.0, output=60.0
        )

        # pricing — add a model not in the built-in table
        AGENT_CONFIG.set_rate("openai", "gpt-5-turbo", input=5.0, output=20.0)

        # budget cap
        AGENT_CONFIG.set_budget(usd=5.0)
        print(AGENT_CONFIG.remaining_usd)  # 5.0

        # inspect
        print(AGENT_CONFIG.info())

        # reset
        AGENT_CONFIG.reset()
    """

    def __init__(self) -> None:
        self._provider: str | None = None
        self._model: str | None = None
        self._keys: dict[str, str] = {}  # provider → explicit key
        self._stack: list[dict] = []  # context-manager save-stack

        # pricing
        # structure: {provider: {model: {"input": float, "output": float}}}
        self._custom_rates: dict[str, dict[str, dict[str, float]]] = {}

        # budget
        self._budget_usd: float | None = None  # None = no cap
        self._spent_usd: float = 0.0

    # ------------------------------------------------------------------
    # Primary configuration API
    # ------------------------------------------------------------------

    def configure(
        self,
        *,
        provider: str,
        api_key: str,
        model: str | None = None,
    ) -> AgentConfig:
        """Set the active provider, key, and optional model in one call.

        Parameters
        ----------
        provider : {"claude", "openai", "gemini"}
            LLM provider to activate.
        api_key : str
            API key for *provider*.  Stored per-provider so switching away
            and back does not require re-supplying the key.
        model : str or None
            Model identifier override.  ``None`` uses the provider default.

        Returns
        -------
        AgentConfig
            *self* — allows chaining.

        Examples
        --------
        ::

            AGENT_CONFIG.configure(provider="claude", api_key="sk-ant-…")
            AGENT_CONFIG.configure(
                provider="openai",
                api_key="sk-…",
                model="gpt-4o-mini",
            )
        """
        self._validate_provider(provider)
        self._provider = provider.lower()
        self._keys[self._provider] = api_key
        self._model = model
        return self

    def set_key(self, provider: str, api_key: str) -> AgentConfig:
        """Store an API key for *provider* without changing the active provider.

        Useful for pre-loading keys for multiple providers so you can
        :meth:`switch` between them without re-supplying credentials.

        Parameters
        ----------
        provider : str
        api_key : str

        Returns
        -------
        AgentConfig
            *self*.

        Examples
        --------
        ::

            AGENT_CONFIG.set_key("claude", "sk-ant-…")
            AGENT_CONFIG.set_key("openai", "sk-…")
            AGENT_CONFIG.set_key("gemini", "AIza…")
            AGENT_CONFIG.switch("claude")
        """
        self._validate_provider(provider)
        self._keys[provider.lower()] = api_key
        return self

    def switch(
        self,
        provider: str,
        *,
        model: str | None = None,
    ) -> AgentConfig:
        """Switch the active provider.

        The key for *provider* must have been stored via :meth:`configure`
        or :meth:`set_key`, or be available as an environment variable.

        Parameters
        ----------
        provider : str
        model : str or None
            Override the model for this provider.  ``None`` keeps the
            current model override (or uses the provider default).

        Returns
        -------
        AgentConfig
            *self*.

        Examples
        --------
        ::

            AGENT_CONFIG.switch("openai")
            AGENT_CONFIG.switch("gemini", model="gemini-2.0-flash")
        """
        self._validate_provider(provider)
        self._provider = provider.lower()
        if model is not None:
            self._model = model
        return self

    def reset(self, *, keys: bool = True) -> AgentConfig:
        """Clear the active configuration, custom rates, and budget.

        Parameters
        ----------
        keys : bool
            When ``True`` (default) also wipe all stored per-provider keys.
            Pass ``False`` to clear only the active provider / model while
            keeping stored keys available for future :meth:`switch` calls.

        Returns
        -------
        AgentConfig
            *self*.

        Notes
        -----
        Custom rates and the budget cap/counter are always reset.
        Call :meth:`reset_budget` to zero only the spend counter.
        """
        self._provider = None
        self._model = None
        self._custom_rates.clear()
        self._budget_usd = None
        self._spent_usd = 0.0
        if keys:
            self._keys.clear()
        return self

    # ------------------------------------------------------------------
    # Pricing API
    # ------------------------------------------------------------------

    def set_rate(
        self,
        provider: str,
        model: str,
        *,
        input: float,  # noqa: A002  — mirrors industry terminology
        output: float,
    ) -> AgentConfig:
        """Override or add the cost rate for a specific provider + model.

        Custom rates take precedence over the built-in table for every cost
        estimate made by agents in this session.

        Parameters
        ----------
        provider : str
            ``"claude"`` | ``"openai"`` | ``"gemini"``
        model : str
            Model identifier exactly as passed to the LLM client, e.g.
            ``"claude-opus-4-8"`` or ``"gpt-5-turbo"``.
        input : float
            Cost per **1 000 000 input tokens** in USD.
        output : float
            Cost per **1 000 000 output tokens** in USD.

        Returns
        -------
        AgentConfig
            *self*.

        Examples
        --------
        ::

            # provider lowered their price
            AGENT_CONFIG.set_rate(
                "claude", "claude-opus-4-8", input=12.0, output=60.0
            )

            # a new model not yet in the built-in table
            AGENT_CONFIG.set_rate(
                "openai", "gpt-5-turbo", input=5.0, output=20.0
            )
        """
        self._validate_provider(provider)
        if input < 0 or output < 0:
            raise ValueError("Rates must be non-negative (USD / 1 M tokens).")
        p = provider.lower()
        self._custom_rates.setdefault(p, {})[model] = {
            "input": float(input),
            "output": float(output),
        }
        return self

    def get_rate(self, provider: str, model: str) -> dict[str, float]:
        """Return the resolved ``{"input": …, "output": …}`` rate.

        Resolution order: custom override → built-in table (exact match →
        prefix match) → provider default.

        Parameters
        ----------
        provider : str
        model : str

        Returns
        -------
        dict with keys ``"input"`` and ``"output"`` (USD / 1 M tokens).

        Examples
        --------
        ::

            rate = AGENT_CONFIG.get_rate("claude", "claude-sonnet-4-6")
            # {"input": 3.0, "output": 15.0}
        """
        p = provider.lower()

        # 1 — custom override (exact match)
        custom = self._custom_rates.get(p, {})
        if model in custom:
            return custom[model]

        # 2 — built-in table (exact match)
        builtin = _BUILTIN_RATES.get(p, {})
        if model in builtin:
            return builtin[model]

        # 3 — prefix match in custom, then built-in
        for table in (custom, builtin):
            for key, rate in table.items():
                if model.startswith(key) or key.startswith(model):
                    return rate

        # 4 — provider default
        return _PROVIDER_DEFAULTS.get(p, {"input": 3.00, "output": 15.00})

    def estimate_cost(
        self,
        provider: str,
        model: str,
        input_tokens: int,
        output_tokens: int,
    ) -> float:
        """Compute the USD cost for one LLM call using the resolved rate.

        Respects any rate overrides set via :meth:`set_rate`.

        Parameters
        ----------
        provider, model : str
        input_tokens, output_tokens : int

        Returns
        -------
        float — estimated cost in USD.

        Examples
        --------
        ::

            cost = AGENT_CONFIG.estimate_cost(
                "claude", "claude-sonnet-4-6", 500, 200
            )
        """
        rate = self.get_rate(provider, model)
        return (
            input_tokens * rate["input"] + output_tokens * rate["output"]
        ) / 1_000_000

    def list_rates(self, provider: str | None = None) -> dict:
        """Return the effective rate table (custom overrides merged with built-in).

        Parameters
        ----------
        provider : str or None
            When given, return only that provider's table.
            When ``None``, return all providers.

        Returns
        -------
        dict — same structure as the built-in table.

        Examples
        --------
        ::

            import pprint

            pprint.pprint(AGENT_CONFIG.list_rates("claude"))
        """
        providers = [provider.lower()] if provider else list(_BUILTIN_RATES)
        result: dict[str, dict[str, dict[str, float]]] = {}
        for p in providers:
            merged = dict(_BUILTIN_RATES.get(p, {}))
            merged.update(self._custom_rates.get(p, {}))  # custom wins
            result[p] = merged
        return result if len(providers) > 1 else result.get(providers[0], {})

    # ------------------------------------------------------------------
    # Budget API
    # ------------------------------------------------------------------

    def set_budget(self, *, usd: float) -> AgentConfig:
        """Set a session spend cap.

        Once :attr:`spent_usd` reaches *usd*, any subsequent LLM call
        will raise :class:`BudgetExceededError` before the API call is made.

        Parameters
        ----------
        usd : float
            Maximum spend in USD for this session.

        Returns
        -------
        AgentConfig
            *self*.

        Examples
        --------
        ::

            AGENT_CONFIG.set_budget(usd=2.0)
            # … after several agent calls …
            print(
                f"${AGENT_CONFIG.spent_usd:.4f} used, "
                f"${AGENT_CONFIG.remaining_usd:.4f} left"
            )
        """
        if usd <= 0:
            raise ValueError("Budget cap must be a positive USD amount.")
        self._budget_usd = float(usd)
        return self

    def reset_budget(self, *, cap: bool = False) -> AgentConfig:
        """Reset the session spend counter.

        Parameters
        ----------
        cap : bool
            When ``True`` also remove the budget cap (``set_budget`` must be
            called again to re-enable it).  Default ``False`` — zeroes the
            counter but keeps the existing cap.

        Returns
        -------
        AgentConfig
            *self*.

        Examples
        --------
        ::

            AGENT_CONFIG.reset_budget()  # zero counter, keep cap
            AGENT_CONFIG.reset_budget(cap=True)  # zero counter and remove cap
        """
        self._spent_usd = 0.0
        if cap:
            self._budget_usd = None
        return self

    @property
    def spent_usd(self) -> float:
        """Accumulated LLM spend this session (USD)."""
        return self._spent_usd

    @property
    def remaining_usd(self) -> float | None:
        """Remaining budget (USD), or ``None`` when no cap is set."""
        if self._budget_usd is None:
            return None
        return max(0.0, self._budget_usd - self._spent_usd)

    def _add_spend(self, cost: float) -> None:
        """Accumulate *cost* into the session counter.

        Called internally by :meth:`~pycsamt.agents.BaseAgent.query_llm`
        after every successful LLM response.
        """
        self._spent_usd += cost

    def _check_budget(self) -> None:
        """Raise :class:`BudgetExceededError` if the cap is already hit.

        Called internally by :meth:`~pycsamt.agents.BaseAgent.query_llm`
        *before* each API call.
        """
        if self._budget_usd is not None and self._spent_usd >= self._budget_usd:
            raise BudgetExceededError(self._spent_usd, self._budget_usd)

    # ------------------------------------------------------------------
    # Resolved-value properties
    # ------------------------------------------------------------------

    @property
    def provider(self) -> str | None:
        """Currently active provider, or ``None`` if unconfigured."""
        return self._provider

    @property
    def api_key(self) -> str | None:
        """Resolved key for the active provider.

        Resolution order: stored key → environment variable → ``None``.
        """
        return self._resolve_key(self._provider) if self._provider else None

    @property
    def model(self) -> str | None:
        """Active model (explicit override, or provider default)."""
        if self._model:
            return self._model
        return _DEFAULT_MODELS.get(self._provider) if self._provider else None

    @property
    def is_configured(self) -> bool:
        """``True`` when a provider and a resolvable key are both present."""
        return self._provider is not None and self.api_key is not None

    # ------------------------------------------------------------------
    # Resolution used by BaseAgent
    # ------------------------------------------------------------------

    def resolve(
        self,
        provider: str,
        api_key: str | None,
        model: str | None,
    ) -> tuple[str, str | None, str | None]:
        """Resolve the effective ``(provider, api_key, model)`` for an agent.

        Called automatically by :class:`~pycsamt.agents.BaseAgent.__init__`.
        Users should not call this directly.

        Resolution rules
        ----------------
        *If* ``api_key`` is explicitly given:
            Use ``provider``, ``api_key``, and ``model`` (or provider default)
            exactly as supplied.  The global config is ignored.

        *If* ``api_key`` is ``None``:
            - If ``provider`` equals the default ``"claude"`` **and** the
              global config has an active provider, inherit that provider.
            - Look up the key: stored key → env var.
            - Inherit the model from the global config when the effective
              provider matches the globally active provider.
        """
        provider = provider.lower()

        if api_key is not None:
            return provider, api_key, model or _DEFAULT_MODELS.get(provider)

        effective = provider
        if effective == "claude" and self._provider is not None:
            effective = self._provider

        resolved_key = self._resolve_key(effective)

        resolved_model = model
        if resolved_model is None and self._provider == effective:
            resolved_model = self._model
        resolved_model = resolved_model or _DEFAULT_MODELS.get(effective)

        return effective, resolved_key, resolved_model

    # ------------------------------------------------------------------
    # Inspection
    # ------------------------------------------------------------------

    def info(self) -> dict:
        """Return a summary dict suitable for display or logging.

        The API key is masked to its last four characters for safety.

        Returns
        -------
        dict
            Keys: ``provider``, ``model``, ``has_key``, ``key_source``,
            ``stored_providers``, ``custom_rate_models``,
            ``budget_usd``, ``spent_usd``, ``remaining_usd``.
        """
        key = self.api_key
        if key:
            source = (
                "explicit"
                if (self._provider and self._provider in self._keys)
                else "env"
            )
        else:
            source = "none"

        custom_models: dict[str, list[str]] = {
            p: sorted(m.keys()) for p, m in self._custom_rates.items() if m
        }

        return {
            "provider": self._provider,
            "model": self.model,
            "has_key": key is not None,
            "api_key_masked": f"…{key[-4:]}" if key else None,
            "key_source": source,
            "stored_providers": sorted(self._keys.keys()),
            "custom_rate_models": custom_models,
            "budget_usd": self._budget_usd,
            "spent_usd": round(self._spent_usd, 6),
            "remaining_usd": (
                round(self.remaining_usd, 6) if self.remaining_usd is not None else None
            ),
        }

    # ------------------------------------------------------------------
    # Context-manager: temporary override
    # ------------------------------------------------------------------

    @contextmanager
    def using(
        self,
        *,
        provider: str | None = None,
        api_key: str | None = None,
        model: str | None = None,
    ) -> Generator[AgentConfig, None, None]:
        """Temporarily override the global config inside a ``with`` block.

        All arguments are optional; omitted values are left unchanged.
        The original config (including custom rates and budget) is restored
        on exit even if an exception occurs.

        Parameters
        ----------
        provider : str or None
        api_key : str or None
        model : str or None

        Yields
        ------
        AgentConfig
            *self* with the override applied.

        Examples
        --------
        ::

            with AGENT_CONFIG.using(provider="gemini", api_key="AIza…"):
                r = DataQCAgent().execute(data)
            # original provider / key / model restored here
        """
        saved = {
            "_provider": self._provider,
            "_model": self._model,
            "_keys": dict(self._keys),
            "_custom_rates": {p: dict(m) for p, m in self._custom_rates.items()},
            "_budget_usd": self._budget_usd,
            "_spent_usd": self._spent_usd,
        }
        self._stack.append(saved)
        try:
            if provider is not None:
                self._validate_provider(provider)
                self._provider = provider.lower()
            if api_key is not None and self._provider is not None:
                self._keys[self._provider] = api_key
            if model is not None:
                self._model = model
            yield self
        finally:
            restored = self._stack.pop()
            self._provider = restored["_provider"]
            self._model = restored["_model"]
            self._keys = restored["_keys"]
            self._custom_rates = restored["_custom_rates"]
            self._budget_usd = restored["_budget_usd"]
            self._spent_usd = restored["_spent_usd"]

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _resolve_key(self, provider: str | None) -> str | None:
        """Return the API key for *provider*: stored key → env var → None.

        When :meth:`offline` context is active for the current thread,
        env-var lookup is skipped so that a key in the OS environment
        (e.g. from ``.env.local``) does not cause unexpected LLM calls.
        """
        if provider is None:
            return None
        explicit = self._keys.get(provider)
        if explicit:
            return explicit
        # Skip env resolution in offline mode.
        if getattr(_TLS, "force_offline", False):
            return None
        for env_var in _ENV_KEYS.get(provider, []):
            val = os.environ.get(env_var)
            if val:
                return val
        return None

    @contextmanager
    def offline(self) -> Generator[AgentConfig, None, None]:
        """Context manager: force truly offline mode in the current thread.

        While active, :meth:`_resolve_key` will not inspect environment
        variables, so agents created inside the block receive ``api_key=None``
        even when ``ANTHROPIC_API_KEY`` (or similar) is set in the OS
        environment.  Uses :mod:`threading.local` so concurrent requests
        in other threads are unaffected.

        Examples
        --------
        ::

            with AGENT_CONFIG.offline():
                agent = DataQCAgent()  # api_key will be None
        """
        old = getattr(_TLS, "force_offline", False)
        _TLS.force_offline = True
        try:
            yield self
        finally:
            _TLS.force_offline = old

    @staticmethod
    def _validate_provider(provider: str) -> None:
        if provider.lower() not in _PROVIDERS:
            raise ValueError(
                f"Unknown LLM provider {provider!r}. "
                f"Supported: {sorted(_PROVIDERS)}."
            )

    # ------------------------------------------------------------------
    # Dunder
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        if self._provider is None:
            return "AgentConfig(unconfigured)"
        key = self.api_key
        key_hint = f"…{key[-4:]}" if key else "no-key"
        budget_hint = (
            f", budget=${self._budget_usd:.2f} (spent=${self._spent_usd:.4f})"
            if self._budget_usd is not None
            else ""
        )
        return (
            f"AgentConfig(provider={self._provider!r}, "
            f"model={self.model!r}, key={key_hint!r}{budget_hint})"
        )

    def __bool__(self) -> bool:
        return self.is_configured


# ---------------------------------------------------------------------------
# Module-level singleton  (mirrors PYCSAMT_SECTION, PYCSAMT_STYLE, etc.)
# ---------------------------------------------------------------------------

#: Global singleton.  Import and configure this object once per session.
AGENT_CONFIG: AgentConfig = AgentConfig()


# ---------------------------------------------------------------------------
# Convenience wrappers  (mirrors configure_section(), configure_style(), …)
# ---------------------------------------------------------------------------


def configure_agents(
    *,
    provider: str,
    api_key: str,
    model: str | None = None,
) -> AgentConfig:
    """Configure :data:`AGENT_CONFIG` in one call and return it.

    Equivalent to ``AGENT_CONFIG.configure(...)``.

    Parameters
    ----------
    provider : {"claude", "openai", "gemini"}
    api_key : str
    model : str or None

    Returns
    -------
    AgentConfig

    Examples
    --------
    ::

        from pycsamt.agents import configure_agents

        configure_agents(provider="claude", api_key="sk-ant-…")
    """
    return AGENT_CONFIG.configure(provider=provider, api_key=api_key, model=model)


def reset_agents(*, keys: bool = True) -> AgentConfig:
    """Reset :data:`AGENT_CONFIG` to its unconfigured state.

    Equivalent to ``AGENT_CONFIG.reset(...)``.

    Parameters
    ----------
    keys : bool
        When ``True`` (default) also wipe stored per-provider keys.

    Returns
    -------
    AgentConfig
    """
    return AGENT_CONFIG.reset(keys=keys)
