"""
pycsamt.agents._pricing
=======================

Token-cost tables and estimation utilities for all supported LLM providers.

.. note::

   Rate lookups (:func:`get_rate`, :func:`estimate_cost`) delegate to the
   global :data:`~pycsamt.api.agents.AGENT_CONFIG` singleton, so any
   per-model overrides set via ``AGENT_CONFIG.set_rate(...)`` are
   automatically reflected here.

   :data:`PROVIDER_RATES` exposes the *built-in* baseline table (read-only
   reference); it does **not** include runtime overrides.

Last updated: 2025-06
"""

from __future__ import annotations

# Lazy import to avoid a hard dependency at module load time.
# AGENT_CONFIG is a singleton — importing it here is safe because
# api/agents.py has no circular dependency back into _pricing.py.
from ..api.agents import _BUILTIN_RATES, AGENT_CONFIG

# ---------------------------------------------------------------------------
# Public rate table (built-in baseline, no runtime overrides)
# ---------------------------------------------------------------------------

#: Built-in USD/1 M-token rates.  Use ``AGENT_CONFIG.list_rates()`` to
#: see the *effective* table including any :meth:`~AgentConfig.set_rate`
#: overrides.
PROVIDER_RATES: dict[str, dict[str, dict[str, float]]] = _BUILTIN_RATES


# ---------------------------------------------------------------------------
# Public helpers — delegate to AGENT_CONFIG so custom overrides apply
# ---------------------------------------------------------------------------

def get_rate(provider: str, model: str) -> dict[str, float]:
    """Return the resolved ``{"input": …, "output": …}`` rate.

    Delegates to :meth:`~pycsamt.api.agents.AgentConfig.get_rate`, so any
    per-model overrides set via ``AGENT_CONFIG.set_rate(...)`` take effect.

    Parameters
    ----------
    provider : str
        One of ``"claude"``, ``"openai"``,
        ``"gemini"``, ``"deepseek"``, ``"minimax"``.
    model : str
        Model identifier, e.g. ``"claude-sonnet-4-6"``.

    Returns
    -------
    dict with keys ``"input"`` and ``"output"`` (USD / 1 M tokens).
    """
    return AGENT_CONFIG.get_rate(provider, model)


def estimate_cost(
    provider: str,
    model: str,
    input_tokens: int,
    output_tokens: int,
) -> float:
    """Compute the USD cost for one LLM call.

    Delegates to :meth:`~pycsamt.api.agents.AgentConfig.estimate_cost`, so
    any per-model overrides set via ``AGENT_CONFIG.set_rate(...)`` are used.

    Parameters
    ----------
    provider, model : str
    input_tokens, output_tokens : int

    Returns
    -------
    float
        Estimated cost in USD.

    Examples
    --------
    >>> estimate_cost("claude", "claude-sonnet-4-6", 500, 200)
    0.00451...
    """
    return AGENT_CONFIG.estimate_cost(provider, model, input_tokens, output_tokens)


def format_cost(usd: float) -> str:
    """Return a human-readable cost string.

    Examples
    --------
    >>> format_cost(0.00012)
    '$0.000120'
    >>> format_cost(1.5)
    '$1.5000'
    """
    if usd < 0.001:
        return f"${usd:.6f}"
    if usd < 0.01:
        return f"${usd:.4f}"
    return f"${usd:.4f}"


def estimate_tokens(text: str) -> int:
    """Rough token count estimate (4 characters ≈ 1 token)."""
    return max(1, len(text) // 4)


__all__ = [
    "PROVIDER_RATES",
    "estimate_cost",
    "estimate_tokens",
    "format_cost",
    "get_rate",
]
