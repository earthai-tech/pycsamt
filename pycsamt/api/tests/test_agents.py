from __future__ import annotations

import pytest

from pycsamt.api import agents as agents_module
from pycsamt.api.agents import (
    AGENT_CONFIG,
    AgentConfig,
    BudgetExceededError,
    configure_agents,
    reset_agents,
)

_ALL_ENV_VARS = sorted(
    {var for names in agents_module._ENV_KEYS.values() for var in names}
)


@pytest.fixture(autouse=True)
def _clean_provider_env(monkeypatch):
    """Isolate every test from real provider keys loaded from .env.local."""
    for var in _ALL_ENV_VARS:
        monkeypatch.delenv(var, raising=False)
    yield


@pytest.fixture
def cfg() -> AgentConfig:
    return AgentConfig()


# ─────────────────────────────────────────────────────────────────────────
# _load_env_local
# ─────────────────────────────────────────────────────────────────────────


def test_load_env_local_reads_file_without_overwriting_existing(
    tmp_path, monkeypatch
):
    env_file = tmp_path / ".env.local"
    env_file.write_text(
        "\n".join(
            [
                "# a comment line",
                "",
                "NO_EQUALS_SIGN_HERE",
                "FAKE_TEST_VAR=hello",
                "ALREADY_SET=should_not_overwrite",
            ]
        ),
        encoding="utf-8",
    )
    fake_module_file = tmp_path / "sub" / "agents.py"
    monkeypatch.setattr(agents_module, "__file__", str(fake_module_file))
    monkeypatch.setenv("ALREADY_SET", "original")
    monkeypatch.delenv("FAKE_TEST_VAR", raising=False)

    agents_module._load_env_local()

    assert agents_module.os.environ["FAKE_TEST_VAR"] == "hello"
    assert agents_module.os.environ["ALREADY_SET"] == "original"


# ─────────────────────────────────────────────────────────────────────────
# BudgetExceededError
# ─────────────────────────────────────────────────────────────────────────


def test_budget_exceeded_error_attrs_and_message():
    err = BudgetExceededError(2.5, 2.0)
    assert err.spent_usd == 2.5
    assert err.budget_usd == 2.0
    assert "2.0000" in str(err)
    assert "2.5000" in str(err)


# ─────────────────────────────────────────────────────────────────────────
# configure / set_key / switch / reset
# ─────────────────────────────────────────────────────────────────────────


def test_configure_sets_provider_key_and_model(cfg):
    result = cfg.configure(provider="Claude", api_key="sk-1", model="m1")
    assert result is cfg
    assert cfg.provider == "claude"
    assert cfg.api_key == "sk-1"
    assert cfg.model == "m1"


def test_configure_rejects_unknown_provider(cfg):
    with pytest.raises(ValueError):
        cfg.configure(provider="bogus", api_key="x")


def test_set_key_stores_without_changing_active_provider(cfg):
    cfg.configure(provider="claude", api_key="sk-1")
    cfg.set_key("openai", "sk-2")
    assert cfg.provider == "claude"
    cfg.switch("openai")
    assert cfg.api_key == "sk-2"


def test_switch_changes_provider_and_optionally_model(cfg):
    cfg.set_key("claude", "sk-1")
    cfg.set_key("openai", "sk-2")
    cfg.switch("claude")
    assert cfg.provider == "claude"
    cfg.switch("openai", model="gpt-4o-mini")
    assert cfg.provider == "openai"
    assert cfg.model == "gpt-4o-mini"


def test_reset_clears_everything_by_default(cfg):
    cfg.configure(provider="claude", api_key="sk-1")
    cfg.set_rate("claude", "m1", input=1.0, output=2.0)
    cfg.set_budget(usd=5.0)
    cfg._add_spend(1.0)
    cfg.reset()
    assert cfg.provider is None
    assert cfg.model is None
    assert cfg.list_rates("claude").get("m1") is None
    assert cfg.spent_usd == 0.0
    assert cfg.remaining_usd is None
    assert cfg.api_key is None  # keys cleared too


def test_reset_keys_false_keeps_stored_keys(cfg):
    cfg.configure(provider="claude", api_key="sk-1")
    cfg.reset(keys=False)
    assert cfg.provider is None
    cfg.switch("claude")
    assert cfg.api_key == "sk-1"


# ─────────────────────────────────────────────────────────────────────────
# Pricing: set_rate / get_rate / estimate_cost / list_rates
# ─────────────────────────────────────────────────────────────────────────


def test_set_rate_rejects_negative_values(cfg):
    with pytest.raises(ValueError):
        cfg.set_rate("claude", "m1", input=-1.0, output=1.0)
    with pytest.raises(ValueError):
        cfg.set_rate("claude", "m1", input=1.0, output=-1.0)


def test_get_rate_custom_override_wins(cfg):
    cfg.set_rate("claude", "claude-sonnet-4-6", input=1.0, output=2.0)
    assert cfg.get_rate("claude", "claude-sonnet-4-6") == {
        "input": 1.0,
        "output": 2.0,
    }


def test_get_rate_builtin_exact_match(cfg):
    rate = cfg.get_rate("claude", "claude-sonnet-4-6")
    assert rate == {"input": 3.00, "output": 15.00}


def test_get_rate_prefix_match(cfg):
    rate = cfg.get_rate("claude", "claude-sonnet-4-6-20260101")
    assert rate == {"input": 3.00, "output": 15.00}


def test_get_rate_falls_back_to_provider_default(cfg):
    rate = cfg.get_rate("claude", "totally-unknown-model-xyz")
    assert rate == {"input": 3.00, "output": 15.00}


def test_get_rate_falls_back_to_hardcoded_default_for_unknown_provider(cfg):
    # provider not validated by get_rate (only configure/set_key/switch do)
    rate = cfg.get_rate("unknownprovider", "whatever")
    assert rate == {"input": 3.00, "output": 15.00}


def test_estimate_cost_uses_resolved_rate(cfg):
    cost = cfg.estimate_cost("claude", "claude-sonnet-4-6", 1_000_000, 0)
    assert cost == pytest.approx(3.00)


def test_list_rates_single_provider_merges_custom(cfg):
    cfg.set_rate("claude", "new-model", input=1.0, output=2.0)
    rates = cfg.list_rates("claude")
    assert "new-model" in rates
    assert "claude-sonnet-4-6" in rates


def test_list_rates_all_providers(cfg):
    rates = cfg.list_rates()
    assert set(rates.keys()) >= {"claude", "openai", "gemini"}


# ─────────────────────────────────────────────────────────────────────────
# Budget
# ─────────────────────────────────────────────────────────────────────────


def test_set_budget_rejects_non_positive(cfg):
    with pytest.raises(ValueError):
        cfg.set_budget(usd=0.0)
    with pytest.raises(ValueError):
        cfg.set_budget(usd=-1.0)


def test_budget_tracking_and_check(cfg):
    cfg.set_budget(usd=1.0)
    assert cfg.remaining_usd == 1.0
    cfg._add_spend(0.4)
    assert cfg.spent_usd == 0.4
    assert cfg.remaining_usd == pytest.approx(0.6)
    cfg._check_budget()  # should not raise yet
    cfg._add_spend(0.6)
    with pytest.raises(BudgetExceededError):
        cfg._check_budget()


def test_reset_budget_keeps_cap_by_default(cfg):
    cfg.set_budget(usd=1.0)
    cfg._add_spend(0.5)
    cfg.reset_budget()
    assert cfg.spent_usd == 0.0
    assert cfg.remaining_usd == 1.0


def test_reset_budget_cap_true_removes_cap(cfg):
    cfg.set_budget(usd=1.0)
    cfg.reset_budget(cap=True)
    assert cfg.remaining_usd is None


def test_remaining_usd_none_when_no_cap(cfg):
    assert cfg.remaining_usd is None


# ─────────────────────────────────────────────────────────────────────────
# Resolved-value properties
# ─────────────────────────────────────────────────────────────────────────


def test_provider_none_when_unconfigured(cfg):
    assert cfg.provider is None


def test_api_key_none_when_no_provider(cfg):
    assert cfg.api_key is None


def test_model_falls_back_to_default_for_provider(cfg):
    cfg.switch("openai")
    assert cfg.model == "gpt-4o"


def test_model_none_when_no_provider_and_no_override(cfg):
    assert cfg.model is None


def test_is_configured_true_and_false(cfg):
    assert cfg.is_configured is False
    cfg.configure(provider="claude", api_key="sk-1")
    assert cfg.is_configured is True


# ─────────────────────────────────────────────────────────────────────────
# resolve()
# ─────────────────────────────────────────────────────────────────────────


def test_resolve_explicit_api_key_ignores_global_config(cfg):
    cfg.configure(provider="openai", api_key="global-key")
    provider, key, model = cfg.resolve("claude", "explicit-key", None)
    assert provider == "claude"
    assert key == "explicit-key"
    assert model == "claude-sonnet-4-6"


def test_resolve_explicit_api_key_and_model(cfg):
    provider, key, model = cfg.resolve("openai", "k", "gpt-4o-mini")
    assert (provider, key, model) == ("openai", "k", "gpt-4o-mini")


def test_resolve_inherits_global_provider_when_default_claude(cfg):
    cfg.configure(provider="openai", api_key="sk-openai")
    provider, key, model = cfg.resolve("claude", None, None)
    assert provider == "openai"
    assert key == "sk-openai"
    assert model == "gpt-4o"


def test_resolve_inherits_model_only_when_provider_matches(cfg):
    cfg.configure(provider="openai", api_key="sk-openai", model="gpt-4o-mini")
    provider, key, model = cfg.resolve("claude", None, None)
    assert provider == "openai"
    assert model == "gpt-4o-mini"


def test_resolve_non_claude_default_does_not_inherit_global_provider(cfg):
    cfg.configure(provider="claude", api_key="sk-claude")
    cfg.set_key("gemini", "sk-gemini")
    provider, key, model = cfg.resolve("gemini", None, None)
    assert provider == "gemini"
    assert key == "sk-gemini"
    assert model == "gemini-2.0-flash"


def test_resolve_returns_none_key_when_nothing_available(cfg):
    provider, key, model = cfg.resolve("claude", None, None)
    assert provider == "claude"
    assert key is None
    assert model == "claude-sonnet-4-6"


# ─────────────────────────────────────────────────────────────────────────
# info()
# ─────────────────────────────────────────────────────────────────────────


def test_info_unconfigured(cfg):
    info = cfg.info()
    assert info["provider"] is None
    assert info["has_key"] is False
    assert info["api_key_masked"] is None
    assert info["key_source"] == "none"
    assert info["stored_providers"] == []
    assert info["custom_rate_models"] == {}
    assert info["remaining_usd"] is None


def test_info_explicit_key_source_and_masking(cfg):
    cfg.configure(provider="claude", api_key="sk-ant-abcd1234")
    info = cfg.info()
    assert info["key_source"] == "explicit"
    assert info["has_key"] is True
    assert info["api_key_masked"] == "…1234"


def test_info_env_key_source(cfg, monkeypatch):
    monkeypatch.setenv("ANTHROPIC_API_KEY", "env-key-value")
    cfg.switch("claude")
    info = cfg.info()
    assert info["key_source"] == "env"


def test_info_reports_custom_rate_models_and_budget(cfg):
    cfg.set_rate("claude", "m1", input=1.0, output=2.0)
    cfg.set_budget(usd=5.0)
    cfg._add_spend(1.0)
    info = cfg.info()
    assert info["custom_rate_models"] == {"claude": ["m1"]}
    assert info["budget_usd"] == 5.0
    assert info["spent_usd"] == 1.0
    assert info["remaining_usd"] == 4.0


# ─────────────────────────────────────────────────────────────────────────
# using() context manager
# ─────────────────────────────────────────────────────────────────────────


def test_using_applies_and_restores_override(cfg):
    cfg.configure(provider="claude", api_key="sk-claude")
    with cfg.using(provider="openai", api_key="sk-openai", model="gpt-4o-mini"):
        assert cfg.provider == "openai"
        assert cfg.api_key == "sk-openai"
        assert cfg.model == "gpt-4o-mini"
    assert cfg.provider == "claude"
    assert cfg.api_key == "sk-claude"


def test_using_restores_state_even_on_exception(cfg):
    cfg.configure(provider="claude", api_key="sk-claude")
    with pytest.raises(RuntimeError):
        with cfg.using(provider="openai", api_key="sk-openai"):
            assert cfg.provider == "openai"
            raise RuntimeError("boom")
    assert cfg.provider == "claude"
    assert cfg.api_key == "sk-claude"


def test_using_with_no_arguments_changes_nothing(cfg):
    cfg.configure(provider="claude", api_key="sk-claude")
    with cfg.using():
        assert cfg.provider == "claude"
    assert cfg.provider == "claude"


def test_using_api_key_ignored_when_no_active_provider(cfg):
    with cfg.using(api_key="sk-orphan"):
        assert cfg.provider is None
        assert cfg.api_key is None


# ─────────────────────────────────────────────────────────────────────────
# _resolve_key / offline()
# ─────────────────────────────────────────────────────────────────────────


def test_resolve_key_none_provider_returns_none(cfg):
    assert cfg._resolve_key(None) is None


def test_resolve_key_prefers_explicit_over_env(cfg, monkeypatch):
    monkeypatch.setenv("ANTHROPIC_API_KEY", "env-value")
    cfg.set_key("claude", "explicit-value")
    assert cfg._resolve_key("claude") == "explicit-value"


def test_resolve_key_falls_back_to_env(cfg, monkeypatch):
    monkeypatch.setenv("ANTHROPIC_API_KEY", "env-value")
    assert cfg._resolve_key("claude") == "env-value"


def test_resolve_key_none_when_nothing_set(cfg):
    assert cfg._resolve_key("claude") is None


def test_offline_skips_env_lookup(cfg, monkeypatch):
    monkeypatch.setenv("ANTHROPIC_API_KEY", "env-value")
    with cfg.offline():
        assert cfg._resolve_key("claude") is None
    assert cfg._resolve_key("claude") == "env-value"


def test_offline_restores_previous_flag_on_nested_use(cfg, monkeypatch):
    monkeypatch.setenv("ANTHROPIC_API_KEY", "env-value")
    with cfg.offline():
        with cfg.offline():
            assert cfg._resolve_key("claude") is None
        assert cfg._resolve_key("claude") is None
    assert cfg._resolve_key("claude") == "env-value"


# ─────────────────────────────────────────────────────────────────────────
# __repr__ / __bool__
# ─────────────────────────────────────────────────────────────────────────


def test_repr_unconfigured(cfg):
    assert repr(cfg) == "AgentConfig(unconfigured)"


def test_repr_configured_without_budget(cfg):
    cfg.configure(provider="claude", api_key="sk-ant-abcd1234")
    text = repr(cfg)
    assert "provider='claude'" in text
    assert "…1234" in text
    assert "budget=" not in text


def test_repr_configured_with_budget_and_no_key(cfg):
    cfg.switch("claude")
    cfg.set_budget(usd=2.0)
    text = repr(cfg)
    assert "no-key" in text
    assert "budget=$2.00" in text


def test_bool_reflects_is_configured(cfg):
    assert bool(cfg) is False
    cfg.configure(provider="claude", api_key="sk-1")
    assert bool(cfg) is True


# ─────────────────────────────────────────────────────────────────────────
# Module-level singleton + convenience wrappers
# ─────────────────────────────────────────────────────────────────────────


def test_configure_agents_and_reset_agents_use_global_singleton():
    try:
        result = configure_agents(provider="claude", api_key="sk-global")
        assert result is AGENT_CONFIG
        assert AGENT_CONFIG.provider == "claude"
        assert AGENT_CONFIG.api_key == "sk-global"

        reset_agents()
        assert AGENT_CONFIG.provider is None
        assert AGENT_CONFIG.api_key is None
    finally:
        AGENT_CONFIG.reset()
