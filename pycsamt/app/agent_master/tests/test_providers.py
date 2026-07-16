# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.agent_master._providers."""

from __future__ import annotations


class TestIsLlm:
    def test_offline_is_not_llm(self):
        from pycsamt.app.agent_master._providers import OFFLINE, is_llm

        assert is_llm(OFFLINE) is False

    def test_known_provider_is_llm(self):
        from pycsamt.app.agent_master._providers import is_llm

        assert is_llm("claude") is True

    def test_none_or_unknown_is_not_llm(self):
        from pycsamt.app.agent_master._providers import is_llm

        assert is_llm(None) is False
        assert is_llm("nope") is False


class TestLabelFor:
    def test_offline_label(self):
        from pycsamt.app.agent_master._providers import label_for

        assert label_for("offline") == "Offline"

    def test_known_provider_label(self):
        from pycsamt.app.agent_master._providers import label_for

        assert label_for("claude") == "Claude (Anthropic)"

    def test_unknown_provider_returns_itself(self):
        from pycsamt.app.agent_master._providers import label_for

        assert label_for("mystery") == "mystery"

    def test_none_returns_empty_string(self):
        from pycsamt.app.agent_master._providers import label_for

        assert label_for(None) == ""


class TestEnvVar:
    def test_known_provider_env_var(self):
        from pycsamt.app.agent_master._providers import env_var

        assert env_var("claude") == "ANTHROPIC_API_KEY"
        assert env_var("openai") == "OPENAI_API_KEY"

    def test_unknown_provider_returns_none(self):
        from pycsamt.app.agent_master._providers import env_var

        assert env_var("mystery") is None
        assert env_var(None) is None


class TestModelsFor:
    def test_known_provider_returns_models(self):
        from pycsamt.app.agent_master._providers import models_for

        models = models_for("claude")
        assert "claude-sonnet-4-6" in models

    def test_unknown_provider_returns_empty_list(self):
        from pycsamt.app.agent_master._providers import models_for

        assert models_for("mystery") == []
        assert models_for(None) == []

    def test_returned_list_is_a_copy(self):
        from pycsamt.app.agent_master._providers import (
            PROVIDER_MODELS,
            models_for,
        )

        models = models_for("claude")
        models.append("bogus-model")
        assert "bogus-model" not in PROVIDER_MODELS["claude"]


class TestDefaultModel:
    def test_returns_first_model(self):
        from pycsamt.app.agent_master._providers import (
            PROVIDER_MODELS,
            default_model,
        )

        assert default_model("claude") == PROVIDER_MODELS["claude"][0]

    def test_unknown_provider_returns_none(self):
        from pycsamt.app.agent_master._providers import default_model

        assert default_model("mystery") is None


class TestProviderOptions:
    def test_includes_offline_option(self):
        from pycsamt.app.agent_master._providers import (
            OFFLINE,
            PROVIDER_OPTIONS,
        )

        values = [opt["value"] for opt in PROVIDER_OPTIONS]
        assert OFFLINE in values

    def test_includes_all_meta_providers(self):
        from pycsamt.app.agent_master._providers import (
            PROVIDER_META,
            PROVIDER_OPTIONS,
        )

        values = {opt["value"] for opt in PROVIDER_OPTIONS}
        assert set(PROVIDER_META) <= values
