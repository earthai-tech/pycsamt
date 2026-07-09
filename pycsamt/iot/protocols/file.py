"""File-backed telemetry transport (newline-delimited JSON).

A :class:`FileTelemetryClient` appends each packet as one JSON object per
line (JSONL) to ``endpoint``. It doubles as a source: :meth:`receive`
walks previously written lines. This transport is dependency-free and is
the natural default for offline logging, replay, and testing.
"""

from __future__ import annotations

import json
import os
from typing import Any, Dict, List, Optional

from ..core import TelemetryPacket
from .base import BaseTelemetryClient, IoTProtocol, TelemetryError

__all__ = ["FileTelemetryClient"]


class FileTelemetryClient(BaseTelemetryClient):
    """Append/replay telemetry packets to a JSON-lines file."""

    protocol = IoTProtocol.FILE

    def __init__(
        self,
        endpoint: Optional[str] = None,
        *,
        dry_run: bool = False,
        append: bool = True,
        **options: Any,
    ) -> None:
        super().__init__(endpoint, dry_run=dry_run, append=append, **options)
        self._read_cursor = 0
        self._read_buffer: List[Dict[str, Any]] = []

    def _connect(self) -> None:
        path = self._require_endpoint()
        parent = os.path.dirname(os.path.abspath(path))
        if parent:
            os.makedirs(parent, exist_ok=True)
        mode = "a" if self.options.get("append", True) else "w"
        if mode == "w" and os.path.exists(path):
            open(path, "w", encoding="utf-8").close()

    def _disconnect(self) -> None:
        self._read_buffer = []
        self._read_cursor = 0

    def _transport_send(self, packet: TelemetryPacket) -> str:
        path = self._require_endpoint()
        line = json.dumps(packet.as_dict(), default=str)
        with open(path, "a", encoding="utf-8") as handle:
            handle.write(line + "\n")
        return f"appended to {path}"

    def _load_buffer(self) -> None:
        path = self._require_endpoint()
        if not os.path.exists(path):
            self._read_buffer = []
            return
        with open(path, "r", encoding="utf-8") as handle:
            rows: List[Dict[str, Any]] = []
            for raw in handle:
                raw = raw.strip()
                if not raw:
                    continue
                try:
                    rows.append(json.loads(raw))
                except json.JSONDecodeError as exc:
                    raise TelemetryError(
                        f"malformed telemetry line in {path}: {exc}"
                    ) from exc
        self._read_buffer = rows

    def _transport_receive(
        self, *, timeout: Optional[float] = None
    ) -> Optional[Dict[str, Any]]:
        if self._read_cursor == 0 or not self._read_buffer:
            self._load_buffer()
        if self._read_cursor >= len(self._read_buffer):
            return None
        row = self._read_buffer[self._read_cursor]
        self._read_cursor += 1
        return row

    def read_all(self) -> List[Dict[str, Any]]:
        """Return every packet dictionary written to the file."""
        self._load_buffer()
        return list(self._read_buffer)

    def _transport_healthcheck(self) -> bool:
        path = self._require_endpoint()
        parent = os.path.dirname(os.path.abspath(path)) or "."
        return os.access(parent, os.W_OK)
