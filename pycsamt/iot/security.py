"""Security configuration for IoT telemetry transports.

Field IoT reviewers expect a credible security story: TLS for transport
encryption and a defined credential scheme (bearer token, basic auth, or
API key). This module provides small, serialisable configuration objects
that:

* carry TLS material and credentials,
* **redact secrets** from ``repr``/``as_dict`` so they never leak into
  logs, manifests, or notebooks, and
* project into the ``**options`` accepted by
  :func:`~pycsamt.iot.protocols.build_telemetry_client`.

The objects hold configuration only; they do not implement cryptography
themselves — encryption is delegated to the underlying transport (TLS via
the OS/OpenSSL, tokens via HTTP headers).
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from enum import Enum
from typing import Any

from ..api.property import PyCSAMTObject
from . import _common as _c

__all__ = [
    "AuthScheme",
    "TLSConfig",
    "Credential",
    "SecurityConfig",
    "redact_secret",
]

_REDACTED = "***REDACTED***"


def redact_secret(value: Any) -> str | None:
    """Return a redaction placeholder when *value* is a non-empty secret."""
    if value is None or value == "":
        return None
    return _REDACTED


class AuthScheme(str, Enum):
    """Supported authentication schemes."""

    NONE = "none"
    BEARER = "bearer"
    BASIC = "basic"
    API_KEY = "api_key"


@dataclass(repr=False)
class TLSConfig(PyCSAMTObject):
    """Transport-layer security material for a telemetry client."""

    enabled: bool = False
    ca_cert: str | None = None
    certfile: str | None = None
    keyfile: str | None = None
    verify: bool = True
    min_version: str | None = None

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        self.enabled = _c.as_bool(self.enabled)
        self.ca_cert = _c.as_optional_str(self.ca_cert, "ca_cert")
        self.certfile = _c.as_optional_str(self.certfile, "certfile")
        self.keyfile = _c.as_optional_str(self.keyfile, "keyfile")
        self.verify = _c.as_bool(self.verify)
        self.min_version = _c.as_optional_str(self.min_version, "min_version")
        if (self.certfile is None) != (self.keyfile is None):
            raise ValueError("certfile and keyfile must be provided together.")

    def as_dict(self) -> dict[str, Any]:
        """Return TLS configuration (paths are not secrets)."""
        return dict(
            enabled=self.enabled,
            ca_cert=self.ca_cert,
            certfile=self.certfile,
            keyfile=self.keyfile,
            verify=self.verify,
            min_version=self.min_version,
        )


@dataclass(repr=False)
class Credential(PyCSAMTObject):
    """Authentication credential with secret redaction.

    Secrets (``token``, ``password``, ``api_key``) are never shown by
    ``repr`` or :meth:`as_dict`; use :meth:`reveal` explicitly when a real
    transport needs the raw values.
    """

    __repr_exclude__ = {"token", "password", "api_key", "logger", "verbose"}

    scheme: AuthScheme | str = AuthScheme.NONE
    token: str | None = None
    username: str | None = None
    password: str | None = None
    api_key: str | None = None
    api_key_header: str = "X-API-Key"

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        self.scheme = _c.normalise_enum(self.scheme, AuthScheme, "scheme")
        self.token = _c.as_optional_str(self.token, "token")
        self.username = _c.as_optional_str(self.username, "username")
        self.password = _c.as_optional_str(self.password, "password")
        self.api_key = _c.as_optional_str(self.api_key, "api_key")
        self.api_key_header = _c.as_nonempty_str(self.api_key_header, "api_key_header")
        if self.scheme is AuthScheme.BEARER and not self.token:
            raise ValueError("bearer scheme requires a token.")
        if self.scheme is AuthScheme.BASIC and not (self.username and self.password):
            raise ValueError("basic scheme requires username and password.")
        if self.scheme is AuthScheme.API_KEY and not self.api_key:
            raise ValueError("api_key scheme requires an api_key.")

    def headers(self) -> dict[str, str]:
        """Return HTTP headers implementing the configured scheme."""
        if self.scheme is AuthScheme.BEARER and self.token:
            return {"Authorization": f"Bearer {self.token}"}
        if self.scheme is AuthScheme.BASIC and self.username:
            import base64

            raw = f"{self.username}:{self.password or ''}".encode()
            token = base64.b64encode(raw).decode("ascii")
            return {"Authorization": f"Basic {token}"}
        if self.scheme is AuthScheme.API_KEY and self.api_key:
            return {self.api_key_header: self.api_key}
        return {}

    def reveal(self) -> dict[str, Any]:
        """Return the raw credential values (handle with care)."""
        return dict(
            scheme=self.scheme.value
            if isinstance(self.scheme, AuthScheme)
            else str(self.scheme),
            token=self.token,
            username=self.username,
            password=self.password,
            api_key=self.api_key,
            api_key_header=self.api_key_header,
        )

    def as_dict(self) -> dict[str, Any]:
        """Return a redacted, serialisable credential summary."""
        return dict(
            scheme=self.scheme.value
            if isinstance(self.scheme, AuthScheme)
            else str(self.scheme),
            token=redact_secret(self.token),
            username=self.username,
            password=redact_secret(self.password),
            api_key=redact_secret(self.api_key),
            api_key_header=self.api_key_header,
        )


@dataclass(repr=False)
class SecurityConfig(PyCSAMTObject):
    """Combined TLS + credential policy for telemetry transports."""

    tls: TLSConfig = field(default_factory=TLSConfig)
    credential: Credential = field(default_factory=Credential)
    require_tls: bool = False
    allowed_protocols: list[str] | None = None

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        if isinstance(self.tls, dict):
            self.tls = TLSConfig(**self.tls)
        if not isinstance(self.tls, TLSConfig):
            raise TypeError("tls must be a TLSConfig.")
        self.tls.validate()
        if isinstance(self.credential, dict):
            self.credential = Credential(**self.credential)
        if not isinstance(self.credential, Credential):
            raise TypeError("credential must be a Credential.")
        self.credential.validate()
        self.require_tls = _c.as_bool(self.require_tls)
        if self.allowed_protocols is not None:
            self.allowed_protocols = [
                _c.as_nonempty_str(p, "protocol").lower()
                for p in self.allowed_protocols
            ]
        if self.require_tls and not self.tls.enabled:
            raise ValueError("require_tls is set but TLS is not enabled.")

    def allows(self, protocol: str) -> bool:
        """Return whether *protocol* is permitted by policy."""
        if self.allowed_protocols is None:
            return True
        return str(protocol).lower() in self.allowed_protocols

    def client_options(self) -> dict[str, Any]:
        """Return options for :func:`build_telemetry_client`.

        Includes TLS flags, credential headers, and (for MQTT) username /
        password so a transport can authenticate. Raw secrets are included
        because this feeds the live client, not a log.
        """
        options: dict[str, Any] = {}
        if self.tls.enabled:
            options["tls"] = True
            if self.tls.certfile:
                options["certfile"] = self.tls.certfile
            if self.tls.keyfile:
                options["keyfile"] = self.tls.keyfile
            if self.tls.ca_cert:
                options["ca_cert"] = self.tls.ca_cert
        headers = self.credential.headers()
        if headers:
            options["headers"] = headers
        if self.credential.token:
            options["token"] = self.credential.token
        if self.credential.username:
            options["username"] = self.credential.username
        if self.credential.password:
            options["password"] = self.credential.password
        return options

    def as_dict(self) -> dict[str, Any]:
        """Return a redacted, serialisable security summary."""
        return dict(
            tls=self.tls.as_dict(),
            credential=self.credential.as_dict(),
            require_tls=self.require_tls,
            allowed_protocols=(
                list(self.allowed_protocols)
                if self.allowed_protocols is not None
                else None
            ),
        )

    @classmethod
    def from_env(cls, prefix: str = "PYCSAMT_IOT_") -> SecurityConfig:
        """Build a security config from environment variables.

        Recognised variables (with *prefix*): ``TOKEN``, ``API_KEY``,
        ``USERNAME``, ``PASSWORD``, ``TLS`` (truthy), ``CA_CERT``,
        ``CERTFILE``, ``KEYFILE``.
        """

        def _env(name: str) -> str | None:
            return os.environ.get(prefix + name)

        token = _env("TOKEN")
        api_key = _env("API_KEY")
        username = _env("USERNAME")
        password = _env("PASSWORD")
        if token:
            credential = Credential(scheme=AuthScheme.BEARER, token=token)
        elif api_key:
            credential = Credential(scheme=AuthScheme.API_KEY, api_key=api_key)
        elif username and password:
            credential = Credential(
                scheme=AuthScheme.BASIC, username=username, password=password
            )
        else:
            credential = Credential()
        tls_enabled = False
        tls_raw = _env("TLS")
        if tls_raw is not None:
            try:
                tls_enabled = _c.as_bool(tls_raw)
            except ValueError:
                tls_enabled = False
        tls = TLSConfig(
            enabled=tls_enabled,
            ca_cert=_env("CA_CERT"),
            certfile=_env("CERTFILE"),
            keyfile=_env("KEYFILE"),
        )
        return cls(tls=tls, credential=credential, require_tls=False)
