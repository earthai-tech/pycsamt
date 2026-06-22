# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for agent_registry (Phase 4) — no Qt required."""

from __future__ import annotations

import pytest
from pycsamt.app.desktop.agent_registry import (
    AGENT_REGISTRY,
    agent_names,
    default_params,
    get_entry,
    llm_agents,
    processing_agents,
)


# ── AGENT_REGISTRY structure ──────────────────────────────────────────────

def test_registry_is_non_empty():
    assert len(AGENT_REGISTRY) > 0


def test_all_entries_have_type():
    for name, entry in AGENT_REGISTRY.items():
        assert "type" in entry, f"'{name}' missing 'type'"
        assert entry["type"] in ("llm", "processing"), f"'{name}' bad type"


def test_llm_entries_have_class_name():
    for name, entry in AGENT_REGISTRY.items():
        if entry["type"] == "llm":
            assert "class_name" in entry, f"LLM agent '{name}' missing 'class_name'"


def test_processing_entries_have_fn_name():
    for name, entry in AGENT_REGISTRY.items():
        if entry["type"] == "processing":
            assert "fn_name" in entry, f"Processing agent '{name}' missing 'fn_name'"


def test_all_entries_have_description():
    for name, entry in AGENT_REGISTRY.items():
        assert "description" in entry, f"'{name}' missing 'description'"


def test_all_param_specs_have_type():
    for agent_name, entry in AGENT_REGISTRY.items():
        for p_name, spec in entry.get("params", {}).items():
            assert "type" in spec, (
                f"'{agent_name}'.params.'{p_name}' missing 'type'"
            )


def test_all_param_specs_have_default():
    for agent_name, entry in AGENT_REGISTRY.items():
        for p_name, spec in entry.get("params", {}).items():
            assert "default" in spec, (
                f"'{agent_name}'.params.'{p_name}' missing 'default'"
            )


def test_combo_params_have_options():
    for agent_name, entry in AGENT_REGISTRY.items():
        for p_name, spec in entry.get("params", {}).items():
            if spec["type"] == "combo":
                assert "options" in spec and len(spec["options"]) > 0, (
                    f"'{agent_name}'.params.'{p_name}' combo has no options"
                )


# ── Helper functions ──────────────────────────────────────────────────────

def test_agent_names_returns_list():
    names = agent_names()
    assert isinstance(names, list)
    assert len(names) == len(AGENT_REGISTRY)


def test_llm_agents_all_have_llm_type():
    for name in llm_agents():
        assert AGENT_REGISTRY[name]["type"] == "llm"


def test_processing_agents_all_have_processing_type():
    for name in processing_agents():
        assert AGENT_REGISTRY[name]["type"] == "processing"


def test_llm_and_processing_cover_all():
    all_names = set(agent_names())
    covered   = set(llm_agents()) | set(processing_agents())
    assert all_names == covered


def test_get_entry_known():
    assert get_entry("QC Quicklook") is not None


def test_get_entry_unknown():
    assert get_entry("DOES_NOT_EXIST") is None


def test_default_params_returns_dict():
    params = default_params("QC Quicklook")
    assert isinstance(params, dict)


def test_default_params_values_match_spec():
    for agent_name in agent_names():
        defaults = default_params(agent_name)
        entry    = get_entry(agent_name)
        for p_name, spec in entry.get("params", {}).items():
            assert p_name in defaults
            assert defaults[p_name] == spec["default"]


def test_default_params_unknown_agent():
    assert default_params("UNKNOWN") == {}
