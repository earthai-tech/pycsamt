# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Request-scoped observations of language calls, gates and tool dispatch.

Records are observations, not certificates of scientific validity. Unknown
metadata is retained as unknown; opaque scientific objects are not stringified.
"""
from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar
from datetime import datetime, timezone
from functools import wraps
import json
import hashlib
import math
from pathlib import Path
import time
from uuid import uuid4

from ._checkpoint_identity import file_hash, input_identities

_CURRENT = ContextVar("pycsamt_provenance", default=None)
_LANGUAGE_CALL = ContextVar("pycsamt_language_call", default=None)
_PRIVATE = {"api_key", "authorization", "password", "access_token", "secret"}
BOOKKEEPING = {"workflow_plan.json", "agent_trace.json", "environment.json",
               "output_manifest.json", "run_record.json"}


def safe_value(value):
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else {"nonfinite": str(value)}
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(k): safe_value(v) for k, v in value.items() if str(k).lower() not in _PRIVATE}
    if isinstance(value, (list, tuple)):
        return [safe_value(v) for v in value]
    # Numerical arrays live in scientific artifacts; do not substitute repr()
    # for data, or silently claim it can reconstruct arbitrary Python objects.
    return {"python_type": f"{type(value).__module__}.{type(value).__qualname__}",
            "serialized": False}


class RunRecord:
    def __init__(self, payload):
        self.run_id = uuid4().hex
        self.events = []
        self.secrets = set()
        self.paths = set()
        self.errors = []
        self.status = "running"
        self.emit("request_received", request=payload.get("request", payload.get("question")),
                  history=payload.get("history"), history_provided="history" in payload)
        try:
            self.emit("input_manifest", inputs=input_identities(payload))
        except (TypeError, ValueError, OSError) as exc:
            self.emit("input_manifest", inputs=None, error=str(exc))

    def emit(self, event, **fields):
        item = {"sequence": len(self.events), "event": event,
                "utc": datetime.now(timezone.utc).isoformat(), **safe_value(fields)}
        self.events.append(item)
        return item

    def snapshot(self):
        value = {"schema_version": 1, "run_id": self.run_id, "status": self.status,
                 "events": self.events, "recording_errors": self.errors,
                 "credential_redaction": bool(self.secrets)}
        # Known credentials are also removed from error text or echoed prompts.
        def redact(obj):
            if isinstance(obj, str):
                for secret in self.secrets:
                    obj = obj.replace(secret, "<credential-redacted>")
                return obj
            if isinstance(obj, dict):
                return {k: redact(v) for k, v in obj.items()}
            if isinstance(obj, list):
                return [redact(v) for v in obj]
            return obj
        return redact(value)

    def persist(self, path):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(path.suffix + ".tmp")
        temporary.write_text(json.dumps(self.snapshot(), indent=2, allow_nan=False) + "\n", encoding="utf-8")
        temporary.replace(path)
        self.paths.add(path)


def current_record():
    return _CURRENT.get()


def emit(event, **fields):
    record = current_record()
    return record.emit(event, **fields) if record is not None else None


@contextmanager
def provenance_scope(payload):
    existing = current_record()
    if existing is not None:
        yield existing
        return
    record = RunRecord(payload)
    token = _CURRENT.set(record)
    try:
        yield record
    except BaseException as exc:
        record.status = "failed"
        record.emit("request_exception", error=str(exc), exception_type=type(exc).__name__)
        raise
    finally:
        record.emit("request_finished", status=record.status)
        for path in list(record.paths):
            try:
                record.persist(path)
            except OSError as exc:
                record.errors.append(f"Final record persistence failed: {exc}")
        _CURRENT.reset(token)


def recorded_request(function):
    @wraps(function)
    def wrapped(self, *args, **kwargs):
        if current_record() is not None:
            return function(self, *args, **kwargs)
        request = args[0] if args else kwargs.get("input_data", kwargs.get("request", {}))
        payload = request if isinstance(request, dict) else {"request": request, **kwargs}
        with provenance_scope(payload) as record:
            result = function(self, *args, **kwargs)
            record.status = result.status
        if record.errors:
            if result.status == "success":
                result.status = "needs_review"
            result.warnings.extend(record.errors)
            record.status = result.status
            record.emit("recording_incomplete", status=result.status)
        result.data["provenance"] = record.snapshot()
        result.data["recording_complete"] = not record.errors
        return result
    return wrapped


def recorded_language_call(function):
    @wraps(function)
    def wrapped(self, prompt, system_message=None, *, temperature=0.2, max_tokens=1024):
        record = current_record()
        if record is None:
            return function(self, prompt, system_message, temperature=temperature, max_tokens=max_tokens)
        if getattr(self, "api_key", None):
            record.secrets.add(self.api_key)
        call_id = uuid4().hex
        record.emit("language_call_started", call_id=call_id, agent=self.name,
                    provider=self.llm_provider, requested_model=self.model,
                    provider_model_version=None, version_status="not_reported",
                    prompt=prompt, system_message=system_message or self.SYSTEM_PROMPT,
                    temperature=temperature, max_tokens=max_tokens,
                    available=self.llm_available)
        started = time.perf_counter()
        self.last_usage = {}
        token = _LANGUAGE_CALL.set(call_id)
        try:
            text = function(self, prompt, system_message, temperature=temperature, max_tokens=max_tokens)
        except BaseException as exc:
            record.emit("language_call_finished", call_id=call_id, status="failed", error=str(exc),
                        elapsed_seconds=time.perf_counter() - started)
            raise
        finally:
            _LANGUAGE_CALL.reset(token)
        record.emit("language_call_finished", call_id=call_id,
                    status="returned" if text is not None else "no_response",
                    response=text, usage=self.last_usage, elapsed_seconds=time.perf_counter() - started)
        return text
    return wrapped


def transport_attempt(agent):
    emit("language_transport_attempt", call_id=_LANGUAGE_CALL.get(), agent=agent.name,
         provider=agent.llm_provider, model=agent.model)


def provider_response(agent, response):
    """Keep actual response identifiers when the installed SDK exposes them."""
    agent.last_usage = {"returned_model_id": getattr(response, "model", None),
                        "response_id": getattr(response, "id", None)}
    emit("language_provider_response", call_id=_LANGUAGE_CALL.get(), **agent.last_usage)


def artifact_snapshot(output_dir):
    """Observe bytes before/after a step; do not infer producer from filename."""
    if not output_dir:
        return {}
    root = Path(output_dir)
    if not root.exists():
        return {}
    entries = {}
    for path in sorted(root.rglob("*")):
        if path.is_symlink() or not path.is_file() or ".checkpoints" in path.parts or path.name in BOOKKEEPING:
            continue
        entries[path.relative_to(root).as_posix()] = {"sha256": file_hash(path), "size_bytes": path.stat().st_size}
    return entries


def observe_artifacts(output_dir):
    if current_record() is None:
        return {}
    try:
        return artifact_snapshot(output_dir)
    except OSError as exc:
        current_record().errors.append(f"Artifact observation failed: {exc}")
        emit("artifact_observation_failed", error=str(exc))
        return {}


def record_retrieval(agent, query, context, retriever=None):
    if current_record() is None:
        return
    from ..assistant.rag.index_store import read_manifest

    chunks = getattr(context, "chunks", []) if context is not None else []
    corpus = getattr(retriever, "chunks", None)
    corpus_hash = None
    if corpus is not None:
        content = hashlib.sha256()
        for chunk in corpus:
            content.update(json.dumps(chunk.to_dict(), sort_keys=True).encode("utf-8"))
            content.update(b"\n")
        corpus_hash = content.hexdigest()
    emit("retrieval", agent=agent, query=query,
         context=getattr(context, "context_text", None),
         citations=getattr(context, "citations", []),
         ranked_chunks=[{"id": c.id, "source_path": c.source_path, "symbol": c.symbol,
                         "text": c.text} for c in chunks],
         active_corpus_sha256=corpus_hash,
         persisted_index_manifest=read_manifest(),
         persisted_manifest_is_active_corpus_proof=False)
