"""Runtime configuration, adapter registration, and configuration contexts."""

from __future__ import annotations

import copy
import json
import logging
import math
import os
import re
import warnings
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Callable

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
    r"""
    Rules to validate and synthesize station names.

    This policy is applied during inter-format conversion (AVG or
    Jones -> EDI). When a provided name is missing or invalid, a
    synthetic name is derived from the station id.

    Parameters
    ----------
    allow_pattern : str, optional
        Regex character class (without brackets) that defines the
        whitelist of allowed characters. The default accepts ASCII
        letters, digits, underscore and dash.
    maxlen : int, optional
        Maximum length of a station name after normalization.
    prefix : str, optional
        Prefix used when creating synthetic names.
    pad : int, optional
        Zero-padding width for numeric station ids.
    strip : bool, optional
        If True, strip leading and trailing whitespace first.
    custom_normalize : Callable[[str], str], optional
        Hook called before validation. Can perform additional
        transliteration or case normalization.

    Notes
    -----
    Normalization runs as: strip → custom_normalize → filter by
    ``allow_pattern`` → truncate to ``maxlen``. Empty results are
    treated as invalid.

    Examples
    --------
    >>> from pycsamt.core.config import StationNamePolicy
    >>> pol = StationNamePolicy(prefix=\"S\", pad=4)
    >>> pol.ensure(\" Site-1  \", station_id=None)
    'Site-1'
    >>> pol.ensure(None, station_id=7)
    'S0007'
    """

    allow_pattern: str = r"A-Za-z0-9_\-"
    maxlen: int = 32
    prefix: str = "S"
    pad: int = 3
    strip: bool = True
    custom_normalize: Callable[[str], str] = staticmethod(lambda s: s)

    def _try_float(self, v) -> float | None:
        try:
            x = float(v)
            if math.isfinite(x):
                return x
        except Exception:
            pass
        return None

    def validate(self, name: str | None) -> str | None:
        r"""
        Validate and normalize a name.

        Parameters
        ----------
        name : str or None
            Candidate station name.

        Returns
        -------
        str or None
            Normalized name if valid; otherwise ``None``.

        Notes
        -----
        The method does not synthesize names. Use ``ensure`` for a
        name-or-fallback behavior.
        """

        if not name:
            return None
        s = name.strip() if self.strip else name
        s = self.custom_normalize(s)
        s = re.sub(rf"[^ {self.allow_pattern}]", "", s)
        s = s[: self.maxlen]
        return s or None

    def synthesize(self, station_id: Any | None) -> str:
        r"""
        Create a deterministic synthetic name from a station id.

        Parameters
        ----------
        station_id : Any or None
            Station identifier. Numeric ids are zero-padded. Non-
            numeric ids are compacted by removing non-word chars.

        Returns
        -------
        str
            Synthetic station name.

        Examples
        --------
        >>> StationNamePolicy().synthesize(12)
        'S012'
        >>> StationNamePolicy(prefix='X').synthesize('AB-01')
        'AB01'
        """
        if station_id is None:
            return f"{self.prefix}UNK"

        # numeric-friendly path
        x = self._try_float(station_id)
        if x is not None:
            iv = int(round(x))
            return f"{self.prefix}{iv:0{self.pad}d}"

        # non-numeric fallback (compact token)
        sid = str(station_id).strip()
        token = re.sub(r"\W+", "", sid)[: self.maxlen]
        return token or f"{self.prefix}UNK"

    def ensure(
        self,
        name: str | None,
        station_id: Any | None,
    ) -> str:
        r"""
        Return a valid station name, validating or synthesizing.

        Parameters
        ----------
        name : str or None
            Provided station name.
        station_id : Any or None
            Station identifier used if ``name`` is invalid.

        Returns
        -------
        str
            Validated name or a synthetic fallback.

        See Also
        --------
        validate : Validate only, without fallback.
        synthesize : Create a name from the station id only.
        """
        s = self.validate(name)

        # if name is clearly non-numeric, trust it
        if s and not s.isdigit():
            return s

        # prefer numeric id when available (avoids 150.0->"1500")
        x = self._try_float(station_id)
        if x is not None:
            iv = int(round(x))
            # if user passed a numeric name, keep when identical
            if s and s.isdigit() and int(s) == iv:
                return s
            return f"{self.prefix}{iv:0{self.pad}d}"

        # last resort
        return s or self.synthesize(station_id)


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
    backend: dict[str, Any] = field(default_factory=dict)

    def copy(self) -> CoreConfig:
        r"""
        Return a deep copy of the configuration.

        Returns
        -------
        CoreConfig
            A detached copy that can be safely mutated.

        Notes
        -----
        Used internally by :func:`config_context` to restore the
        previous state atomically.
        """
        return copy.deepcopy(self)


_CFG: CoreConfig = CoreConfig()
_ADAPTERS: dict[str, Callable[..., Any]] = {}


def get_config() -> CoreConfig:
    r"""
    Return the live :class:`CoreConfig` singleton.

    Returns
    -------
    CoreConfig
        The active configuration object.

    Notes
    -----
    The returned instance is mutable. Prefer :func:`configure`
    or :func:`config_context` for controlled updates.
    """
    return _CFG


def to_dict() -> dict[str, Any]:
    r"""
    Serialize the current configuration to a plain dict.

    Returns
    -------
    dict
        A JSON-serializable mapping of all fields.

    Examples
    --------
    >>> from pycsamt.core.config import to_dict
    >>> d = to_dict()
    >>> 'freq_order' in d
    True
    """
    return asdict(_CFG)


def configure(**kwargs: Any) -> CoreConfig:
    r"""
    Update configuration fields with validation.

    Parameters
    ----------
    **kwargs : Any
        Field names and values to set on the global config.

    Returns
    -------
    CoreConfig
        The updated configuration object.

    Raises
    ------
    AttributeError
        If an unknown field name is provided.
    ValueError
        If a value is invalid for a known field.

    Notes
    -----
    Also aligns the package logger level to ``log_level``,
    if that logger has been created.

    Examples
    --------
    >>> from pycsamt.core.config import configure
    >>> _ = configure(freq_order='asc', infer_errors=False)
    """
    global _CFG
    for key, value in kwargs.items():
        if not hasattr(_CFG, key):
            raise AttributeError(f"Unknown config field: {key!r}")
        if key == "on_duplicate_station" and value not in {
            "replace",
            "keep",
            "error",
        }:
            raise ValueError(
                "on_duplicate_station must be one of "
                "'replace', 'keep', 'error'"
            )
        if key == "freq_order" and value not in {"asc", "desc"}:
            raise ValueError("freq_order must be 'asc' or 'desc'")
        if key == "target_format" and value != "edi":
            warnings.warn(
                "Only 'edi' is supported as target_format", stacklevel=2
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
    r"""
    Reset the global configuration to factory defaults.

    Notes
    -----
    A new :class:`CoreConfig` instance replaces the existing
    singleton. Any references to the previous instance will not
    reflect subsequent changes.
    """
    global _CFG
    _CFG = CoreConfig()


@contextmanager
def config_context(**overrides: Any) -> Iterator[CoreConfig]:
    r"""
    Context manager for temporary configuration overrides.

    Parameters
    ----------
    **overrides : Any
        Field names and temporary values.

    Yields
    ------
    CoreConfig
        The active configuration during the context.

    Notes
    -----
    On exit, the previous configuration is restored atomically,
    even if an exception is raised.

    Examples
    --------
    >>> from pycsamt.core.config import config_context, get_config
    >>> with config_context(strict=True):
    ...     assert get_config().strict is True
    >>> assert get_config().strict is False
    """
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
    r"""
    Register a format adapter that yields an EDI object or
    collection.

    Parameters
    ----------
    key : str
        Format key, e.g. ``'avg'`` or ``'j'``.
    factory : Callable[..., Any]
        Callable that accepts a source object and returns an EDI
        object or an EDICollection.

    Raises
    ------
    ValueError
        If ``key`` is not a non-empty string.

    Notes
    -----
    Adapters are resolved lazily at call time, which avoids heavy
    imports and circular dependencies.

    See Also
    --------
    get_adapter : Retrieve a registered adapter.
    list_adapters : Inspect the registry.
    """

    if not key or not isinstance(key, str):
        raise ValueError("Adapter key must be a non‑empty string")
    _ADAPTERS[key.lower()] = factory


def get_adapter(key: str) -> Callable[..., Any] | None:
    r"""
    Return the adapter factory for a key.

    Parameters
    ----------
    key : str
        Format key used during registration.

    Returns
    -------
    Callable or None
        The factory callable if found, else ``None``.
    """
    return _ADAPTERS.get(key.lower())


def list_adapters() -> dict[str, str]:
    r"""
    List registered adapters with readable names.

    Returns
    -------
    dict
        Mapping ``{key: 'qualname'}`` for registered factories.

    Examples
    --------
    >>> from pycsamt.core.config import list_adapters
    >>> isinstance(list_adapters(), dict)
    True
    """
    out: dict[str, str] = {}
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
    r"""
    Best-effort user config loader (TOML or JSON).

    Search order is:

    1. ``$PYCSAMT_CONFIG`` if set,
    2. ``~/.config/pycsamt.toml``,
    3. ``~/.pycsamt.toml``,
    4. ``~/.pycsamt.json``.

    Notes
    -----
    Only keys matching :class:`CoreConfig` are applied via
    :func:`configure`. If ``{\"core\": {...}}`` exists, that
    sub-mapping is used.

    Errors are downgraded to warnings and do not stop import.
    """

    for p in _DEF_PATHS:
        if not p or not str(p):
            continue
        if p.exists() and p.is_file():
            try:
                if p.suffix.lower() == ".toml" and tomllib is not None:
                    data = tomllib.loads(p.read_text(encoding="utf-8"))
                else:
                    data = json.loads(p.read_text(encoding="utf-8"))
                if isinstance(data, Mapping):
                    payload = data.get("core", data)
                    if isinstance(payload, Mapping):
                        configure(**dict(payload))
                return
            except Exception as exc:
                warnings.warn(
                    f"Failed to load config from {p}: {exc}", stacklevel=2
                )
                return


_load_user_config()

CoreConfig.__doc__ = r"""
Global configuration container for :mod:`pycsamt`.

The object stores defaults for parsing, conversion, and
population of EDI data structures. Use :func:`configure` for
updates or :func:`config_context` for temporary overrides.

Parameters
----------
empty : float
    Sentinel used in legacy files to mark missing values.
strict : bool
    If True, raise on recoverable parse issues; otherwise,
    degrade gracefully when possible.
case_sensitive_sections : bool
    Treat section names as case-sensitive when reading.
on_duplicate_station : {'replace', 'keep', 'error'}
    Policy when collection loaders encounter duplicate station
    ids.
target_format : {'edi'}
    Preferred internal representation for processing. Keep
    as 'edi' to normalize all inputs to EDI objects.
log_level : str
    Default logging level for package loggers.
freq_order : {'asc', 'desc'}
    Preferred order for frequency vectors after normalization.
freq_tol : float
    Relative tolerance for de-duplicating near-equal
    frequencies.
compute_res_from_z : bool
    Compute (rho, phase) from Z when absent.
compute_z_from_res : bool
    Reconstruct Z from (rho, phase) when Z is missing.
load_spectra : bool
    If True, preserve spectra sections when available.
load_time_series : bool
    If True, preserve time-series sections when available.
station_policy : StationNamePolicy
    Validation and synthesis rules for station names.
error_fill_value : float
    Fill value used when error arrays are required by a
    consumer but not available.
infer_errors : bool
    If True, allow simple heuristics to infer missing errors
    from weights or coherency measures when present.
encoding : str
    Default text encoding for reading legacy files.
newline : str
    Line terminator to use when writing files.
backend : dict
    Reserved space for per-backend knobs.

Notes
-----
The instance stored in this module is mutable by design.
Use :func:`configure` or :func:`config_context` instead of
assigning attributes directly, so validation is applied.

Examples
--------
>>> from pycsamt.core.config import configure, get_config
>>> _ = configure(freq_order='desc', strict=False)
>>> get_config().freq_order
'desc'
"""
