"""HTTP(S) telemetry transport using the standard library.

:class:`HTTPTelemetryClient` POSTs each packet as a JSON body to
``endpoint``. It relies only on :mod:`urllib`, so no third-party HTTP
dependency is required. Bearer tokens and extra headers may be supplied
through options or a :class:`~pycsamt.iot.security.SecurityConfig`.
"""

from __future__ import annotations

from typing import Any

from ..core import TelemetryPacket
from .base import (
    BaseTelemetryClient,
    IoTProtocol,
    TelemetryError,
)

__all__ = ["HTTPTelemetryClient"]


class HTTPTelemetryClient(BaseTelemetryClient):
    """POST telemetry packets to an HTTP(S) endpoint."""

    protocol = IoTProtocol.HTTP

    def __init__(
        self,
        endpoint: str | None = None,
        *,
        dry_run: bool = False,
        timeout: float = 10.0,
        token: str | None = None,
        headers: dict[str, str] | None = None,
        method: str = "POST",
        **options: Any,
    ) -> None:
        super().__init__(
            endpoint,
            dry_run=dry_run,
            timeout=timeout,
            token=token,
            headers=headers,
            method=method,
            **options,
        )

    def _connect(self) -> None:
        # HTTP is connectionless; validate configuration eagerly instead.
        url = self._require_endpoint()
        if not str(url).lower().startswith(("http://", "https://")):
            raise TelemetryError(
                f"HTTP endpoint must start with http:// or https:// (got {url!r})."
            )

    def _headers(self) -> dict[str, str]:
        headers = {"Content-Type": "application/json"}
        extra = self.options.get("headers") or {}
        headers.update({str(k): str(v) for k, v in extra.items()})
        token = self.options.get("token")
        if token:
            headers.setdefault("Authorization", f"Bearer {token}")
        return headers

    def _transport_send(self, packet: TelemetryPacket) -> str:
        import urllib.error
        import urllib.request

        url = self._require_endpoint()
        request = urllib.request.Request(
            url,
            data=self._payload_bytes(packet),
            headers=self._headers(),
            method=str(self.options.get("method", "POST")),
        )
        timeout = float(self.options.get("timeout", 10.0))
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                status = getattr(response, "status", response.getcode())
        except urllib.error.HTTPError as exc:
            raise TelemetryError(
                f"HTTP {exc.code} from {url}: {exc}"
            ) from exc
        except urllib.error.URLError as exc:
            raise TelemetryError(
                f"HTTP request to {url} failed: {exc}"
            ) from exc
        if not 200 <= int(status) < 300:
            raise TelemetryError(f"HTTP endpoint returned status {status}.")
        return f"http {status}"

    def _transport_healthcheck(self) -> bool:
        # Reachability is only known on a real request; assume config-valid.
        return bool(self.endpoint)
