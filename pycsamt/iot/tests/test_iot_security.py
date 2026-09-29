"""Tests for :mod:`pycsamt.iot.security` (TLS/credential configuration).

Covers validation, header generation per auth scheme, secret redaction
in ``as_dict``/``repr`` vs. raw values in ``reveal``, ``SecurityConfig``
composition/coercion, and ``SecurityConfig.from_env``.
"""

from __future__ import annotations

import pytest

from pycsamt.iot.security import (
    AuthScheme,
    Credential,
    SecurityConfig,
    TLSConfig,
    redact_secret,
)


# ---------------------------------------------------------------------------
# redact_secret
# ---------------------------------------------------------------------------
def test_redact_secret_none_and_empty_pass_through():
    assert redact_secret(None) is None
    assert redact_secret("") is None


def test_redact_secret_masks_nonempty_value():
    assert redact_secret("hunter2") == "***REDACTED***"


# ---------------------------------------------------------------------------
# TLSConfig
# ---------------------------------------------------------------------------
def test_tlsconfig_defaults():
    tls = TLSConfig()
    assert tls.enabled is False
    assert tls.verify is True
    assert tls.ca_cert is None


def test_tlsconfig_certfile_keyfile_must_be_paired():
    with pytest.raises(ValueError, match="together"):
        TLSConfig(certfile="cert.pem")
    with pytest.raises(ValueError, match="together"):
        TLSConfig(keyfile="key.pem")
    # both given is fine
    tls = TLSConfig(certfile="cert.pem", keyfile="key.pem")
    assert tls.certfile == "cert.pem"


def test_tlsconfig_as_dict_contains_all_fields():
    tls = TLSConfig(enabled=True, ca_cert="ca.pem", min_version="TLSv1.2")
    d = tls.as_dict()
    assert d == dict(
        enabled=True,
        ca_cert="ca.pem",
        certfile=None,
        keyfile=None,
        verify=True,
        min_version="TLSv1.2",
    )


def test_tlsconfig_coerces_bool_like_strings():
    tls = TLSConfig(enabled="true", verify="no")
    assert tls.enabled is True
    assert tls.verify is False


# ---------------------------------------------------------------------------
# Credential
# ---------------------------------------------------------------------------
def test_credential_default_scheme_none_requires_nothing():
    cred = Credential()
    assert cred.scheme is AuthScheme.NONE
    assert cred.headers() == {}


def test_credential_bearer_requires_token():
    with pytest.raises(ValueError, match="bearer scheme requires a token"):
        Credential(scheme=AuthScheme.BEARER)
    cred = Credential(scheme=AuthScheme.BEARER, token="abc123")
    assert cred.headers() == {"Authorization": "Bearer abc123"}


def test_credential_basic_requires_username_and_password():
    with pytest.raises(ValueError, match="basic scheme requires"):
        Credential(scheme=AuthScheme.BASIC, username="bob")
    with pytest.raises(ValueError, match="basic scheme requires"):
        Credential(scheme=AuthScheme.BASIC, password="secret")
    cred = Credential(scheme=AuthScheme.BASIC, username="bob", password="secret")
    headers = cred.headers()
    assert headers["Authorization"].startswith("Basic ")


def test_credential_api_key_requires_api_key():
    with pytest.raises(ValueError, match="api_key scheme requires"):
        Credential(scheme=AuthScheme.API_KEY)
    cred = Credential(scheme=AuthScheme.API_KEY, api_key="xyz")
    assert cred.headers() == {"X-API-Key": "xyz"}


def test_credential_api_key_custom_header_name():
    cred = Credential(
        scheme=AuthScheme.API_KEY, api_key="xyz", api_key_header="X-Custom"
    )
    assert cred.headers() == {"X-Custom": "xyz"}


def test_credential_scheme_accepts_string():
    cred = Credential(scheme="bearer", token="abc")
    assert cred.scheme is AuthScheme.BEARER


def test_credential_none_scheme_with_stray_fields_returns_no_headers():
    # scheme=NONE never emits headers even if token/username set
    cred = Credential(scheme=AuthScheme.NONE, token="unused")
    assert cred.headers() == {}


def test_credential_as_dict_redacts_secrets_but_reveal_shows_raw():
    cred = Credential(
        scheme=AuthScheme.BASIC, username="bob", password="secret"
    )
    d = cred.as_dict()
    assert d["password"] == "***REDACTED***"
    assert d["username"] == "bob"
    assert d["scheme"] == "basic"

    raw = cred.reveal()
    assert raw["password"] == "secret"
    assert raw["username"] == "bob"
    assert raw["scheme"] == "basic"


def test_credential_repr_excludes_secret_fields():
    cred = Credential(scheme=AuthScheme.BEARER, token="topsecret")
    r = repr(cred)
    assert "topsecret" not in r


# ---------------------------------------------------------------------------
# SecurityConfig
# ---------------------------------------------------------------------------
def test_securityconfig_defaults():
    cfg = SecurityConfig()
    assert isinstance(cfg.tls, TLSConfig)
    assert isinstance(cfg.credential, Credential)
    assert cfg.require_tls is False
    assert cfg.allowed_protocols is None


def test_securityconfig_coerces_dict_tls_and_credential():
    cfg = SecurityConfig(
        tls={"enabled": True},
        credential={"scheme": "api_key", "api_key": "k1"},
    )
    assert isinstance(cfg.tls, TLSConfig)
    assert cfg.tls.enabled is True
    assert isinstance(cfg.credential, Credential)
    assert cfg.credential.api_key == "k1"


def test_securityconfig_rejects_wrong_types():
    with pytest.raises(TypeError, match="TLSConfig"):
        SecurityConfig(tls=123)
    with pytest.raises(TypeError, match="Credential"):
        SecurityConfig(credential=123)


def test_securityconfig_require_tls_without_enabled_raises():
    with pytest.raises(ValueError, match="require_tls"):
        SecurityConfig(require_tls=True)
    # fine when tls.enabled is True
    cfg = SecurityConfig(tls=TLSConfig(enabled=True), require_tls=True)
    assert cfg.require_tls is True


def test_securityconfig_allowed_protocols_normalised_lowercase():
    cfg = SecurityConfig(allowed_protocols=["MQTT", "Http"])
    assert cfg.allowed_protocols == ["mqtt", "http"]
    assert cfg.allows("MQTT") is True
    assert cfg.allows("serial") is False


def test_securityconfig_allows_everything_when_unset():
    cfg = SecurityConfig()
    assert cfg.allows("anything") is True


def test_securityconfig_client_options_tls_and_headers():
    cfg = SecurityConfig(
        tls=TLSConfig(
            enabled=True, certfile="c.pem", keyfile="k.pem", ca_cert="ca.pem"
        ),
        credential=Credential(scheme=AuthScheme.BEARER, token="tok"),
    )
    opts = cfg.client_options()
    assert opts["tls"] is True
    assert opts["certfile"] == "c.pem"
    assert opts["keyfile"] == "k.pem"
    assert opts["ca_cert"] == "ca.pem"
    assert opts["headers"] == {"Authorization": "Bearer tok"}
    assert opts["token"] == "tok"
    assert "username" not in opts
    assert "password" not in opts


def test_securityconfig_client_options_basic_username_password():
    cfg = SecurityConfig(
        credential=Credential(
            scheme=AuthScheme.BASIC, username="bob", password="secret"
        )
    )
    opts = cfg.client_options()
    assert opts["username"] == "bob"
    assert opts["password"] == "secret"
    assert "Authorization" in opts["headers"]
    assert "tls" not in opts


def test_securityconfig_client_options_no_tls_no_credential():
    cfg = SecurityConfig()
    assert cfg.client_options() == {}


def test_securityconfig_as_dict_redacts_and_serialises():
    cfg = SecurityConfig(
        tls=TLSConfig(enabled=True),
        credential=Credential(scheme=AuthScheme.API_KEY, api_key="secretkey"),
        allowed_protocols=["mqtt"],
    )
    d = cfg.as_dict()
    assert d["tls"]["enabled"] is True
    assert d["credential"]["api_key"] == "***REDACTED***"
    assert d["require_tls"] is False
    assert d["allowed_protocols"] == ["mqtt"]


def test_securityconfig_as_dict_allowed_protocols_none():
    cfg = SecurityConfig()
    assert cfg.as_dict()["allowed_protocols"] is None


# ---------------------------------------------------------------------------
# SecurityConfig.from_env
# ---------------------------------------------------------------------------
def test_from_env_token_wins_over_others(monkeypatch):
    monkeypatch.setenv("PYCSAMT_IOT_TOKEN", "tok123")
    monkeypatch.setenv("PYCSAMT_IOT_API_KEY", "shouldnotuse")
    cfg = SecurityConfig.from_env()
    assert cfg.credential.scheme is AuthScheme.BEARER
    assert cfg.credential.token == "tok123"


def test_from_env_api_key_used_when_no_token(monkeypatch):
    monkeypatch.delenv("PYCSAMT_IOT_TOKEN", raising=False)
    monkeypatch.setenv("PYCSAMT_IOT_API_KEY", "apikey123")
    cfg = SecurityConfig.from_env()
    assert cfg.credential.scheme is AuthScheme.API_KEY
    assert cfg.credential.api_key == "apikey123"


def test_from_env_basic_used_when_username_and_password(monkeypatch):
    monkeypatch.delenv("PYCSAMT_IOT_TOKEN", raising=False)
    monkeypatch.delenv("PYCSAMT_IOT_API_KEY", raising=False)
    monkeypatch.setenv("PYCSAMT_IOT_USERNAME", "bob")
    monkeypatch.setenv("PYCSAMT_IOT_PASSWORD", "secret")
    cfg = SecurityConfig.from_env()
    assert cfg.credential.scheme is AuthScheme.BASIC
    assert cfg.credential.username == "bob"


def test_from_env_no_credentials_gives_none_scheme(monkeypatch):
    for name in ("TOKEN", "API_KEY", "USERNAME", "PASSWORD", "TLS"):
        monkeypatch.delenv(f"PYCSAMT_IOT_{name}", raising=False)
    cfg = SecurityConfig.from_env()
    assert cfg.credential.scheme is AuthScheme.NONE


def test_from_env_tls_flag_and_paths(monkeypatch):
    for name in ("TOKEN", "API_KEY", "USERNAME", "PASSWORD"):
        monkeypatch.delenv(f"PYCSAMT_IOT_{name}", raising=False)
    monkeypatch.setenv("PYCSAMT_IOT_TLS", "true")
    monkeypatch.setenv("PYCSAMT_IOT_CA_CERT", "ca.pem")
    monkeypatch.setenv("PYCSAMT_IOT_CERTFILE", "c.pem")
    monkeypatch.setenv("PYCSAMT_IOT_KEYFILE", "k.pem")
    cfg = SecurityConfig.from_env()
    assert cfg.tls.enabled is True
    assert cfg.tls.ca_cert == "ca.pem"
    assert cfg.tls.certfile == "c.pem"
    assert cfg.tls.keyfile == "k.pem"


def test_from_env_invalid_tls_flag_falls_back_to_false(monkeypatch):
    for name in ("TOKEN", "API_KEY", "USERNAME", "PASSWORD"):
        monkeypatch.delenv(f"PYCSAMT_IOT_{name}", raising=False)
    monkeypatch.setenv("PYCSAMT_IOT_TLS", "not-a-bool")
    cfg = SecurityConfig.from_env()
    assert cfg.tls.enabled is False


def test_from_env_custom_prefix(monkeypatch):
    monkeypatch.setenv("MYAPP_TOKEN", "tok-custom")
    cfg = SecurityConfig.from_env(prefix="MYAPP_")
    assert cfg.credential.token == "tok-custom"


if __name__ == "__main__":  # pragma: no-cover
    pytest.main([__file__])
