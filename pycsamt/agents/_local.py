# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Local Ollama transport and request-scoped limits, without cloud credentials.

Uses the native chat API. No proxy, redirects, remote endpoints, automatic
model downloads, or cloud fallback are permitted by this provider.
"""

from __future__ import annotations

import http.client
import ipaddress
import json
import math
import queue
import socket
import threading
import time
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field
from urllib.parse import urlsplit


class LocalModelError(RuntimeError):
    """Actionable local-provider failure; never triggers cloud fallback."""


class LocalCancelled(LocalModelError):
    """The user cancelled the active request."""


@dataclass(frozen=True)
class LocalSettings:
    endpoint: str = "http://127.0.0.1:11434"
    model: str = "qwen2.5-coder:1.5b"
    timeout: float = 60.0
    context_tokens: int = 8192
    output_tokens: int = 1024
    temperature: float = 0.2
    max_calls: int = 4

    def __post_init__(self):
        url = urlsplit(self.endpoint)
        host = url.hostname or ""
        try:
            loopback = (
                host == "localhost" or ipaddress.ip_address(host).is_loopback
            )
        except ValueError:
            loopback = False
        if (
            url.scheme != "http"
            or not loopback
            or url.username
            or url.password
            or url.path not in ("", "/")
            or url.query
            or url.fragment
        ):
            raise ValueError(
                "Local Ollama endpoint must be a loopback HTTP URL, e.g. http://127.0.0.1:11434."
            )
        # Access port here to reject invalid port strings before any I/O.
        if url.port is not None and not 1 <= url.port <= 65535:
            raise ValueError("Invalid local endpoint port.")
        if not self.model.strip() or "cloud" in self.model.lower():
            raise ValueError(
                "Choose an installed local model, not an Ollama cloud model."
            )
        for name, low, high in (
            ("timeout", 1, 600),
            ("context_tokens", 512, 131072),
            ("output_tokens", 1, 16384),
            ("temperature", 0, 2),
            ("max_calls", 1, 20),
        ):
            value = getattr(self, name)
            if not math.isfinite(value) or not low <= value <= high:
                raise ValueError(
                    f"Local {name} must be between {low} and {high}."
                )
        if self.output_tokens >= self.context_tokens:
            raise ValueError(
                "Local output limit must be smaller than the context window."
            )

    @classmethod
    def from_mapping(cls, cfg: dict) -> LocalSettings:
        return cls(
            endpoint=str(cfg.get("ollama_endpoint") or cls.endpoint).strip(),
            model=str(cfg.get("model_ollama") or cls.model).strip(),
            timeout=float(cfg.get("ollama_timeout", cls.timeout)),
            context_tokens=int(cfg.get("ollama_context", cls.context_tokens)),
            output_tokens=int(cfg.get("ollama_output", cls.output_tokens)),
            temperature=float(cfg.get("ollama_temperature", cls.temperature)),
            max_calls=int(cfg.get("ollama_max_calls", cls.max_calls)),
        )


@dataclass
class LocalRequest:
    settings: LocalSettings
    cancelled: Callable[[], bool] = lambda: False
    started: float = field(default_factory=time.monotonic)
    calls: int = 0
    usage: list[dict] = field(default_factory=list)
    verified_models: set[str] = field(default_factory=set)

    def check(self):
        if self.cancelled():
            raise LocalCancelled("Local generation stopped by user.")
        if time.monotonic() - self.started >= self.settings.timeout:
            raise LocalModelError(
                "Local request timed out. Try a smaller model, less context, or a longer time limit."
            )


_REQUEST: ContextVar[LocalRequest | None] = ContextVar(
    "pycsamt_local_request", default=None
)


def current_request() -> LocalRequest | None:
    return _REQUEST.get()


def local_only() -> bool:
    return current_request() is not None


@contextmanager
def local_session(
    settings: LocalSettings, cancelled: Callable[[], bool] = lambda: False
) -> Iterator[LocalRequest]:
    request = LocalRequest(settings, cancelled)
    token = _REQUEST.set(request)
    try:
        yield request
    finally:
        _REQUEST.reset(token)


def _exchange(
    request: LocalRequest,
    path: str,
    payload: dict | None = None,
    *,
    stream: bool = False,
):
    """Run blocking HTTP in a bounded worker so Stop also interrupts loading."""
    request.check()
    url = urlsplit(request.settings.endpoint)
    host = "127.0.0.1" if url.hostname == "localhost" else url.hostname
    connection = http.client.HTTPConnection(
        host, url.port or 80, timeout=min(2, request.settings.timeout)
    )
    result: queue.Queue = queue.Queue(maxsize=1)

    def work():
        try:
            connection.connect()
            if connection.sock:
                connection.sock.settimeout(request.settings.timeout)
            connection.request(
                "POST" if payload is not None else "GET",
                path,
                body=json.dumps(payload) if payload is not None else None,
                headers={"Content-Type": "application/json"},
            )
            response = connection.getresponse()
            if response.status != 200:
                detail = response.read(16384).decode("utf-8", errors="replace")
                if response.status == 404:
                    raise LocalModelError(
                        "Local model not found. Install the selected model in Ollama and check its exact name."
                    )
                raise LocalModelError(
                    f"Ollama returned HTTP {response.status}: {detail[:500]}. Check server logs and available memory."
                )
            if stream:
                text, final, size = [], None, 0
                for line in response:
                    request.check()
                    size += len(line)
                    if size > 8 * 1024 * 1024:
                        raise LocalModelError(
                            "Local response exceeded the size limit."
                        )
                    if not line.strip():
                        continue
                    part = json.loads(line)
                    if part.get("error"):
                        raise LocalModelError(f"Ollama: {part['error']}")
                    text.append(part.get("message", {}).get("content", ""))
                    if part.get("done"):
                        final = part
                        break
                if final is None:
                    raise LocalModelError(
                        "Ollama stream ended before completion; no complete answer received."
                    )
                if final.get("done_reason") == "length":
                    raise LocalModelError(
                        "Local output limit reached. Increase the output limit or request a shorter answer."
                    )
                value = ("".join(text), final)
            else:
                raw = response.read(8 * 1024 * 1024 + 1)
                if len(raw) > 8 * 1024 * 1024:
                    raise LocalModelError(
                        "Local metadata response exceeded the size limit."
                    )
                value = json.loads(raw)
            result.put((value, None))
        except Exception as exc:
            result.put((None, exc))
        finally:
            connection.close()

    worker = threading.Thread(target=work, daemon=True, name="pycsamt-ollama")
    worker.start()
    try:
        while True:
            request.check()
            try:
                value, error = result.get(timeout=0.05)
                break
            except queue.Empty:
                continue
        if isinstance(error, LocalModelError):
            raise error
        if isinstance(error, (TimeoutError, socket.timeout)):
            raise LocalModelError(
                "Ollama timed out. Try a smaller model or increase the request time limit."
            ) from error
        if isinstance(error, OSError):
            raise LocalModelError(
                "Cannot reach local Ollama. Start 'ollama serve' and check the endpoint in Settings."
            ) from error
        if error:
            raise LocalModelError(
                f"Invalid Ollama response: {error}"
            ) from error
        return value
    finally:
        # Shutdown unblocks a read held by the worker, including model load.
        active_socket = connection.sock
        if active_socket:
            try:
                active_socket.shutdown(socket.SHUT_RDWR)
            except OSError:
                pass
        connection.close()
        worker.join(timeout=0.1)


def verify_model(request: LocalRequest, model: str) -> dict:
    info = _exchange(request, "/api/show", {"model": model})
    if info.get("remote_host") or info.get("remote_model"):
        raise LocalModelError(
            "This Ollama model is hosted remotely. Choose a downloaded local model."
        )
    return info


def model_status(settings: LocalSettings) -> dict:
    request = LocalRequest(settings)
    models = _exchange(request, "/api/tags").get("models", [])
    info = verify_model(request, settings.model)
    return {
        "models": models,
        "selected": settings.model,
        "details": info.get("details", {}),
        "status": "ready",
        "note": "Installed locally; generation speed has not been tested.",
    }


def generate(
    prompt: str,
    system: str,
    *,
    model: str,
    max_tokens: int,
    temperature: float,
) -> tuple[str, dict]:
    request = current_request() or LocalRequest(
        LocalSettings(model=model, temperature=temperature)
    )
    request.check()
    if request.calls >= request.settings.max_calls:
        raise LocalModelError(
            "Local request call budget exhausted; no further generation attempted."
        )
    request.calls += 1
    if model not in request.verified_models:
        verify_model(request, model)
        request.verified_models.add(model)
    settings = request.settings
    output = min(max_tokens, settings.output_tokens)
    # Conservative UTF-8 byte bound, not a tokenizer claim. Never silently
    # truncate grounded context or ask Ollama to drop earlier instructions.
    if (
        len((system + prompt).encode("utf-8")) + output + 128
        > settings.context_tokens
    ):
        raise LocalModelError(
            "Retrieved prompt exceeds the conservative local context budget. Increase context or shorten the request/evidence."
        )
    text, final = _exchange(
        request,
        "/api/chat",
        {
            "model": model,
            "messages": [
                {"role": "system", "content": system},
                {"role": "user", "content": prompt},
            ],
            "stream": True,
            "keep_alive": "5m",
            "options": {
                "num_ctx": settings.context_tokens,
                "num_predict": output,
                "temperature": 0.0
                if temperature == 0
                else settings.temperature,
            },
        },
        stream=True,
    )
    if not text.strip():
        raise LocalModelError("Ollama returned no answer text.")
    usage = {
        k: final.get(k)
        for k in (
            "prompt_eval_count",
            "eval_count",
            "total_duration",
            "load_duration",
        )
    }
    usage.update(
        model=model,
        api_cost_usd=0.0,
        compute_cost="local hardware; not estimated",
    )
    request.usage.append(usage)
    return text, usage
