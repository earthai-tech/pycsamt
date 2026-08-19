from __future__ import annotations

import logging
import uuid

import pytest

from pycsamt.core import config as cfg


def test_public_api_all():
    expected = {
        "StationNamePolicy",
        "CoreConfig",
        "get_config",
        "configure",
        "reset_config",
        "config_context",
        "to_dict",
        "register_adapter",
        "get_adapter",
        "list_adapters",
    }
    assert set(cfg.__all__) == expected


def test_station_name_policy_validate_synthesize():
    p = cfg.StationNamePolicy()
    assert p.validate("  A#1  ") == "A1"
    assert p.synthesize(7) == "S007"
    assert p.ensure(None, 12) == "S012"
    assert p.ensure("Site-1", None) == "Site-1"


def test_get_config_and_reset_defaults():
    cfg.reset_config()
    c = cfg.get_config()
    assert isinstance(c, cfg.CoreConfig)
    assert c.empty == pytest.approx(1.0e32)
    assert c.freq_order in {"asc", "desc"}
    assert c.target_format == "edi"


def test_configure_updates_and_logging_level():
    cfg.reset_config()
    cfg.configure(log_level="INFO", freq_order="asc")
    c = cfg.get_config()
    assert c.log_level == "INFO"
    logger = logging.getLogger("pycsamt")
    assert logger.level in (logging.INFO, 0)
    assert c.freq_order == "asc"


def test_config_context_restores_values():
    cfg.reset_config()
    before = cfg.get_config().strict
    with cfg.config_context(strict=True):
        assert cfg.get_config().strict is True
    assert cfg.get_config().strict == before


def test_configure_validation_errors():
    cfg.reset_config()
    with pytest.raises(AttributeError):
        cfg.configure(unknown_field=1)
    with pytest.raises(ValueError):
        cfg.configure(freq_order="up")
    with pytest.raises(ValueError):
        cfg.configure(on_duplicate_station="bad")


def test_adapter_registry_register_get_list():
    key = f"fmt_{uuid.uuid4().hex[:8]}"

    def dummy_factory(x):
        return x

    cfg.register_adapter(key, dummy_factory)
    got = cfg.get_adapter(key)
    assert got is dummy_factory

    listing = cfg.list_adapters()
    assert key in listing
    assert "dummy" in listing[key]


def test_to_dict_has_core_fields():
    cfg.reset_config()
    d = cfg.to_dict()
    assert isinstance(d, dict)
    for k in (
        "empty",
        "strict",
        "freq_order",
        "station_policy",
    ):
        assert k in d
