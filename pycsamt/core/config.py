# -*- coding: utf-8 -*-
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any, Callable, Dict, Iterator, Mapping
from typing import Optional
from contextlib import contextmanager
from pathlib import Path
import os
import re
import copy
import json
import logging
import warnings

try:  # 3.11+
    import tomllib  # type: ignore[attr-defined]
except Exception:  # pragma: no cover
    tomllib = None  # type: ignore

__all__ = [
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
]


@dataclass
class StationNamePolicy:
    allow_pattern: str = r"A-Za-z0-9_\-"
    maxlen: int = 32
    prefix: str = "S"
    pad: int = 3
    strip: bool = True
    custom_normalize: Callable[[str], str] = staticmethod(
        lambda s: s
    )

    def validate(self, name: Optional[str]) -> Optional[str]:
        if not name:
            return None
        s = name.strip() if self.strip else name
        s = self.custom_normalize(s)
        s = re.sub(fr"[^ {self.allow_pattern}]", "", s)
        s = s[: self.maxlen]
        return s or None

    def synthesize(self, station_id: Optional[Any]) -> str:
        if station_id is None:
            return f"{self.prefix}UNK"
        sid = str(station_id).strip()
        if sid.isdigit():
            return f"{self.prefix}{int(sid):0{self.pad}d}"
        token = re.sub(r"\W+", "", sid)[: self.maxlen]
        return token or f"{self.prefix}UNK"

    def ensure(
        self,
        name: Optional[str],
        station_id: Optional[Any],
    ) -> str:
        return self.validate(name) or self.synthesize(station_id)


@dataclass
class CoreConfig:
    empty: float = 1.0e32
    strict: bool = False
    case_sensitive_sections: bool = False
    on_duplicate_station: str = "replace"
    target_format: str = "edi"
    log_level: str = "WARNING"
    freq_order: str = "desc"
    freq_tol: float = 1e-9
    compute_res_from_z: bool = True
    compute_z_from_res: bool = True
    load_spectra: bool = True
    load_time_series: bool = False
    station_policy: StationNamePolicy = field(
        default_factory=StationNamePolicy
    )
    error_fill_value: float = float("nan")
    infer_errors: bool = True
    encoding: str = "utf-8"
    newline: str = "\n"
    backend: Dict[str, Any] = field(default_factory=dict)

    def copy(self) -> "CoreConfig":
        return copy.deepcopy(self)


_CFG: CoreConfig = CoreConfig()
_ADAPTERS: Dict[str, Callable[..., Any]] = {}


def get_config() -> CoreConfig:
    return _CFG


def to_dict() -> Dict[str, Any]:
    return asdict(_CFG)


def configure(**kwargs: Any) -> CoreConfig:
    global _CFG
    for key, value in kwargs.items():
        if not hasattr(_CFG, key):
            raise AttributeError(
                f"Unknown config field: {key!r}"
            )
        if (
            key == "on_duplicate_station"
            and value not in {"replace", "keep", "error"}
        ):
            raise ValueError(
                "on_duplicate_station must be one of "
                "'replace', 'keep', 'error'"
            )
        if key == "freq_order" and value not in {"asc", "desc"}:
            raise ValueError("freq_order must be 'asc' or 'desc'")
        if key == "target_format" and value != "edi":
            warnings.warn(
                "Only 'edi' is supported as target_format"
            )
        setattr(_CFG, key, value)

    try:
        lvl = getattr(
            logging,
            _CFG.log_level.upper(),
            logging.WARNING,
        )
        logging.getLogger("pycsamt").setLevel(lvl)
    except Exception:
        pass
    return _CFG


def reset_config() -> None:
    global _CFG
    _CFG = CoreConfig()


@contextmanager
def config_context(**overrides: Any) -> Iterator[CoreConfig]:
    global _CFG
    old = _CFG.copy()
    try:
        configure(**overrides)
        yield _CFG
    finally:
        _CFG = old


def register_adapter(
    key: str,
    factory: Callable[..., Any],
) -> None:
    if not key or not isinstance(key, str):
        raise ValueError("Adapter key must be a non‑empty string")
    _ADAPTERS[key.lower()] = factory


def get_adapter(key: str) -> Optional[Callable[..., Any]]:
    return _ADAPTERS.get(key.lower())


def list_adapters() -> Dict[str, str]:
    out: Dict[str, str] = {}
    for k, f in _ADAPTERS.items():
        try:
            out[k] = getattr(
                f,
                "__qualname__",
                getattr(f, "__name__", repr(f)),
            )
        except Exception:
            out[k] = repr(f)
    return out


_DEF_PATHS = (
    Path(os.environ.get("PYCSAMT_CONFIG", "")),
    Path.home() / ".config" / "pycsamt.toml",
    Path.home() / ".pycsamt.toml",
    Path.home() / ".pycsamt.json",
)


def _load_user_config() -> None:
    for p in _DEF_PATHS:
        if not p or not str(p):
            continue
        if p.exists() and p.is_file():
            try:
                if p.suffix.lower() == ".toml" and tomllib is not None:
                    data = tomllib.loads(
                        p.read_text(encoding="utf-8")
                    )
                else:
                    data = json.loads(
                        p.read_text(encoding="utf-8")
                    )
                if isinstance(data, Mapping):
                    payload = data.get("core", data)
                    if isinstance(payload, Mapping):
                        configure(**dict(payload))
                return
            except Exception as exc:
                warnings.warn(
                    f"Failed to load config from {p}: {exc}"
                )
                return


_load_user_config()
