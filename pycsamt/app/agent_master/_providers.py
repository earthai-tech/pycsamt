# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Single source of truth for the Agent Master's LLM providers.

The settings panel renders one dropdown plus one contextual credential
panel, so a new provider is a single entry here — no extra ids, hidden
blocks or callback outputs.
"""

from __future__ import annotations

__all__ = [
    "OFFLINE",
    "PROVIDER_MODELS",
    "PROVIDER_META",
    "PROVIDER_OPTIONS",
    "is_llm",
    "label_for",
    "env_var",
    "models_for",
    "default_model",
]

#: The deterministic, key-free engine. Not an LLM provider.
OFFLINE = "offline"

PROVIDER_MODELS: dict[str, list[str]] = {
    "claude": [
        "claude-sonnet-4-6",
        "claude-opus-4-8",
        "claude-haiku-4-5-20251001",
    ],
    "openai": [
        "gpt-4o",
        "gpt-4o-mini",
        "gpt-4-turbo",
    ],
    "gemini": [
        "gemini-2.0-flash",
        "gemini-1.5-pro",
        "gemini-1.5-flash",
    ],
    "deepseek": [
        "deepseek-chat",
        "deepseek-reasoner",
    ],
    "minimax": [
        "MiniMax-M3",
    ],
}

#: ``provider -> (display name, environment variable fallback)``
PROVIDER_META: dict[str, tuple[str, str]] = {
    "claude":   ("Claude (Anthropic)", "ANTHROPIC_API_KEY"),
    "openai":   ("OpenAI",             "OPENAI_API_KEY"),
    "gemini":   ("Gemini (Google)",    "GOOGLE_API_KEY"),
    "deepseek": ("DeepSeek",           "DEEPSEEK_API_KEY"),
    "minimax":  ("MiniMax",            "MINIMAX_API_KEY"),
}

PROVIDER_OPTIONS: list[dict[str, str]] = [
    {"label": label, "value": key}
    for key, (label, _env) in PROVIDER_META.items()
] + [{"label": "Offline (rule-based, no key)", "value": OFFLINE}]


def is_llm(provider: str | None) -> bool:
    """True when *provider* needs an API key (i.e. is not offline)."""
    return bool(provider) and provider in PROVIDER_META


def label_for(provider: str | None) -> str:
    if provider == OFFLINE:
        return "Offline"
    meta = PROVIDER_META.get(provider or "")
    return meta[0] if meta else (provider or "")


def env_var(provider: str | None) -> str | None:
    meta = PROVIDER_META.get(provider or "")
    return meta[1] if meta else None


def models_for(provider: str | None) -> list[str]:
    return list(PROVIDER_MODELS.get(provider or "", []))


def default_model(provider: str | None) -> str | None:
    models = models_for(provider)
    return models[0] if models else None
