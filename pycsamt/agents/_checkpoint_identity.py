# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Conservative content identities for workflow approvals and local resume.

Opaque inputs disable reuse. Pickle bytes are never used to infer input identity.
"""
from __future__ import annotations

import hashlib
import inspect
import json
import math
from pathlib import Path
from typing import Any


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(64 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def data_identity(path: str | Path) -> dict:
    """Identify an input file/directory by names and bytes, not timestamps."""
    path = Path(path).resolve(strict=True)
    if path.is_file():
        return {"path": str(path), "sha256": file_hash(path)}
    if not path.is_dir():
        raise ValueError(f"Not a regular data file or directory: {path}")
    files = []
    for child in sorted(path.rglob("*")):
        if child.is_symlink():
            raise ValueError("Symbolic links require an explicit input manifest.")
        if child.is_file():
            files.append([child.relative_to(path).as_posix(), file_hash(child)])
    return {"path": str(path), "files": files}


def canonical(value: Any) -> Any:
    """Accept only explicit deterministic values; never fall back to repr()."""
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float) and math.isfinite(value):
        return value
    if isinstance(value, Path):
        return str(value.resolve())
    if isinstance(value, (list, tuple)):
        return [canonical(item) for item in value]
    if isinstance(value, dict) and all(isinstance(k, str) for k in value):
        return {k: canonical(v) for k, v in sorted(value.items())}
    raise ValueError(f"No stable checkpoint identity for {type(value).__name__}")


def digest(value: Any) -> str:
    payload = json.dumps(canonical(value), sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def input_identities(value: Any, prefix: str = "") -> dict:
    """Find declared input paths, including nested step/solver settings."""
    found = {}
    if isinstance(value, dict):
        for key, item in value.items():
            label = f"{prefix}.{key}" if prefix else key
            if key in {"path", "data_path", "edi_dir", "checkpoint", "project_dir", "binary_path"} and item:
                found[label] = data_identity(item)
            elif isinstance(item, (dict, list, tuple)):
                found.update(input_identities(item, label))
    elif isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            found.update(input_identities(item, f"{prefix}[{index}]"))
    return found


def _callable_identity(fn) -> dict | None:
    if fn is None:
        return None
    if not inspect.isfunction(fn):
        raise ValueError("Opaque step adapters cannot be resumed.")
    source = inspect.getsource(fn)
    closure = inspect.getclosurevars(fn)
    # Global scientific functions are covered by their module source below;
    # a mutable custom global value cannot be treated as a stable adapter.
    globals_used = {}
    for key, value in closure.globals.items():
        if inspect.ismodule(value):
            module_file = getattr(value, "__file__", None)
            if not module_file:
                raise ValueError("Adapter module has no source identity.")
            globals_used[key] = file_hash(Path(module_file))
        elif inspect.isfunction(value) or inspect.isclass(value):
            globals_used[key] = file_hash(Path(inspect.getfile(value)))
        else:
            globals_used[key] = canonical(value)
    return {"source": source, "defaults": canonical(fn.__defaults__),
            "kwdefaults": canonical(fn.__kwdefaults__),
            "closure": canonical(closure.nonlocals), "globals": globals_used}


def workflow_identity(config: dict, steps: list, workflow: str) -> str | None:
    """Return a reuse identity, or None when anything is unverifiable.

    Includes package runtime Python sources, constructor settings we can recover,
    step adapters, run configuration and content of declared input paths.
    Missing constructor attributes and opaque values intentionally disable reuse.
    """
    try:
        import sys
        from importlib.metadata import version, PackageNotFoundError

        canonical(config)  # reject opaque configurations before scanning sources
        declared_inputs = input_identities(config)
        chain = []
        for step in steps:
            agent = step.agent
            settings = {}
            for name, parameter in inspect.signature(type(agent).__init__).parameters.items():
                if name in {"self", "api_key"}:
                    continue
                if parameter.kind in {parameter.VAR_POSITIONAL, parameter.VAR_KEYWORD}:
                    raise ValueError("Variadic agent configuration needs an explicit contract.")
                if not hasattr(agent, name):
                    raise ValueError(f"Agent constructor setting {name} is not inspectable.")
                settings[name] = canonical(getattr(agent, name))
            agent_file = Path(inspect.getfile(type(agent)))
            chain.append({"step": step.name, "required": step.required,
                          "class": f"{type(agent).__module__}.{type(agent).__qualname__}",
                          "source": file_hash(agent_file), "settings": settings,
                          "adapter": _callable_identity(step.input_fn)})
        package = Path(__file__).resolve().parents[1]
        # Test code is not an execution dependency; custom agent/adapter source
        # files are included separately above even when defined in tests.
        sources = [[p.relative_to(package).as_posix(), file_hash(p)]
                   for p in sorted(package.rglob("*.py")) if "tests" not in p.parts]
        versions = {}
        for name in ("pycsamt", "numpy", "scipy", "torch", "tensorflow"):
            try:
                versions[name] = version(name)
            except PackageNotFoundError:
                versions[name] = None
        return digest({"workflow": workflow, "config": config, "inputs": declared_inputs,
                       "chain": chain, "package_sources": sources,
                       "python": sys.version, "platform": sys.platform, "versions": versions})
    except (ValueError, TypeError, OSError):
        return None
