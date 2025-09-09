# -*- coding: utf-8 -*-
"""
Lightweight parsing utilities for A.G. Jones J‑format files.

This module exposes small, focused helpers used by higher‑level
components (e.g., ``j.py``, ``components.py``) to lex and parse
J files. It does **not** assemble site objects; it only handles
normalization, tokenization, and row conversion with robust error
messages.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence
import io

from ..exceptions import JParseError 

from .config import (
    ENCODING_DEFAULT,
    INFO_KEY_ALIASES,
    KIND_COMPLEX_TF,
    KIND_RHO_PHI,
    RE_BLANK,
    RE_COMMENT,
    RE_DATATYPE_UNITS,
    RE_INFO,
    RE_NPOINTS,
    RE_ROW_R,
    RE_ROW_TF,
    RE_STATION,
    TE_TM_TO_TENSOR,
    UNITS_CANONICAL,
)



@dataclass(frozen=True)
class DataType:
    """Decoded data‑type descriptor.

    Attributes
    ----------
    kind : str
        One of ``R, S, Z, Q, C, T`` (upper‑case).
    comp : str
        Two‑letter component code, e.g., ``XY`` or mode ``TE``.
    units : str | None
        Canonical units string (``"ohms"`` or ``"mV/km/nT"``) or
        ``None`` if not specified.
    tensor_hint : str | None
        For ``TE``/``TM``, the canonical tensor entry (e.g., ``RXY``),
        else ``None``.
    """

    kind: str
    comp: str
    units: str | None
    tensor_hint: str | None = None


@dataclass(frozen=True)
class ParsedRow:
    """A normalized row from an R/S or TF block.

    Attributes
    ----------
    period : float
        Period in seconds (always positive).
    freq : float
        Frequency in Hz (``1/period``). If input stored frequency,
        ``freq`` equals that value and ``period`` is inverted.
    values : Mapping[str, float]
        Other numeric columns by name.
    flags : Mapping[str, bool]
        Flags such as ``stored_as_freq`` or ``rejected``.
    """

    period: float
    freq: float
    values: Mapping[str, float]
    flags: Mapping[str, bool]


def iter_lines(obj: io.TextIOBase | str | Path | Sequence[str], *,
               encoding: str = ENCODING_DEFAULT,
               keepends: bool = False) -> Iterator[str]:
    """Yield lines from a path, file‑like, or a sequence of strings.

    Parameters
    ----------
    obj : file‑like | str | Path | sequence of str
        Source to read from. Paths are opened in text mode.
    encoding : str, default="utf-8"
        Text encoding used for paths.
    keepends : bool, default=False
        Whether to keep line endings in yielded strings.

    Yields
    ------
    str
        One line at a time.
    """
    if isinstance(obj, (str, Path)):
        with open(obj, "r", encoding=encoding) as f:
            for ln in f:
                yield ln if keepends else ln.rstrip("\r\n")
        return

    if isinstance(obj, io.TextIOBase):
        for ln in obj:
            yield ln if keepends else ln.rstrip("\r\n")
        return

    # Assume an iterable of strings
    for ln in obj:  # type: ignore[assignment]
        yield ln if keepends else str(ln).rstrip("\r\n")



def is_comment(s: str) -> bool:
    """Return True if ``s`` is a comment line."""
    return bool(RE_COMMENT.match(s))


def is_blank(s: str) -> bool:
    """Return True if ``s`` is blank or whitespace only."""
    return bool(RE_BLANK.match(s))


def parse_info(line: str, *, lineno: int | None = None) -> tuple[str, str]:
    """Parse an information line ``>KEY = value``.

    Returns the canonical attribute name and the raw value string.
    """
    m = RE_INFO.match(line)
    if not m:
        raise JParseError(_fmt_err("Malformed info line", lineno, line))
    key = m.group("key").upper()
    val = m.group("val").strip()
    return INFO_KEY_ALIASES.get(key, key.lower()), val


def parse_station(line: str, *, lineno: int | None = None) -> str:
    """Parse a station line (alnum up to 6 chars)."""
    m = RE_STATION.match(line)
    if not m:
        raise JParseError(_fmt_err("Malformed station id", lineno, line))
    return m.group("station").upper()


def parse_datatype_units(line: str, *, lineno: int | None = None) -> DataType:
    """Parse a data‑type line like ``ZXY SI`` or ``RTE``.

    The result contains the kind, component, optional units, and a
    ``tensor_hint`` for TE/TM modes.
    """
    line = " ".join(str(line).strip().split())
    
    m = RE_DATATYPE_UNITS.match(line)
    if not m:
        raise JParseError(_fmt_err("Malformed data‑type", lineno, line))

    kind = m.group("kind").upper()
    comp = m.group("comp").upper()

    unit_token = m.groupdict().get("unit")
    units = None
    if unit_token:
        units = _normalize_units(unit_token)

    tensor_hint = None
    if comp in TE_TM_TO_TENSOR:
        tensor_hint = TE_TM_TO_TENSOR[comp]

    return DataType(kind=kind, comp=comp, units=units,
                    tensor_hint=tensor_hint)


def parse_npoints(line: str, *, lineno: int | None = None) -> int:
    """Parse the number of rows to follow."""
    m = RE_NPOINTS.match(line)
    if not m:
        raise JParseError(_fmt_err("Expected integer count", lineno, line))
    return int(m.group("n"))


def parse_row(kind: str, line: str, *, lineno: int | None = None) -> ParsedRow:
    """Parse a single data row for ``kind``.

    Handles R/S (rho‑phi) and complex TF kinds (Z/Q/C/T). Returns a
    :class:`ParsedRow` with period/frequency normalized and rejection
    flags applied according to the spec.
    """
    k = kind.upper()
    if k in KIND_RHO_PHI:
        m = RE_ROW_R.match(line)
        if not m:
            raise JParseError(_fmt_err("Malformed R/S row", lineno, line))
        gd = {name: _to_float(m.group(name), name, lineno)
              for name in ("p", "rho", "pha", "rhomax", "rhomin",
                           "phamax", "phamin", "wrho", "wpha")}
        period, freq, stored_as_freq = _normalize_period(gd.pop("p"))
        flags = {
            "stored_as_freq": stored_as_freq,
            "rejected": (gd["rho"] < 0.0) or (gd["wrho"] < 0.0),
        }
        return ParsedRow(period=period, freq=freq,
                         values=gd, flags=flags)

    if k in KIND_COMPLEX_TF:
        m = RE_ROW_TF.match(line)
        if not m:
            raise JParseError(_fmt_err("Malformed TF row", lineno, line))
        gd = {name: _to_float(m.group(name), name, lineno)
              for name in ("p", "real", "imag", "error", "weight")}
        period, freq, stored_as_freq = _normalize_period(gd.pop("p"))
        flags = {"stored_as_freq": stored_as_freq}
        return ParsedRow(period=period, freq=freq,
                         values=gd, flags=flags)

    raise JParseError(
        _fmt_err(f"Unsupported kind '{kind}'", lineno, line)
    )


def _normalize_units(token: str) -> str:
    """Normalize a unit token to a canonical string.

    Accepts liberal variants (e.g., ``S.I.``, ``si``, ``FIELD``) and
    returns ``"ohms"`` or ``"mV/km/nT"``.
    """
    t = token.strip().upper().replace(".", "")
    key = "SI" if t.startswith("SI") else "FIELD"
    return UNITS_CANONICAL[key]


def _normalize_period(raw_p: float) -> tuple[float, float, bool]:
    """Return ``(period, freq, stored_as_freq)``.

    If ``raw_p < 0``, interpret input as frequency in Hz per spec.
    Otherwise, input is a period in seconds.
    """
    if raw_p < 0.0:
        f = abs(raw_p)
        if f == 0.0:
            raise JParseError("Frequency cannot be zero")
        return 1.0 / f, f, True
    # period
    p = raw_p
    if p <= 0.0:
        raise JParseError("Period must be positive when provided")
    return p, 1.0 / p, False


def _to_float(s: str, name: str, lineno: int | None) -> float:
    try:
        return float(s)
    except Exception as exc:  # noqa: BLE001
        raise JParseError(_fmt_err(f"Bad float for '{name}'",
                                   lineno, s)) from exc


def _fmt_err(msg: str, lineno: int | None, frag: str | None = None) -> str:
    loc = f" at line {lineno}" if lineno is not None else ""
    if frag is None:
        return f"{msg}{loc}."
    snip = str(frag).strip()
    if len(snip) > 60:
        snip = snip[:57] + "..."
    return f"{msg}{loc}: '{snip}'."


def strip_nondata(lines: Iterable[str]) -> Iterator[str]:
    """Yield only non‑comment lines (keep blanks for block logic).

    This keeps structure intact while removing comment noise.
    """
    for i, ln in enumerate(lines, 1):
        if is_comment(ln):
            continue
        yield ln


def iter_info(lines: Iterable[str]) -> Iterator[tuple[str, str]]:
    """Yield ``(key, value)`` pairs from information lines.

    Stops when the first non‑info, non‑comment, non‑blank line is
    encountered. Suitable for reading the header block.
    """
    for i, ln in enumerate(lines, 1):
        if is_blank(ln) or is_comment(ln):
            continue
        if RE_INFO.match(ln):
            yield parse_info(ln, lineno=i)
            continue
        # First non‑info line ends the info block
        break


def iter_rows(kind: str, lines: Iterable[str]) -> Iterator[ParsedRow]:
    """Yield :class:`ParsedRow` objects for rows of ``kind``.

    Parameters
    ----------
    kind : str
        One of ``R, S, Z, Q, C, T``.
    lines : iterable of str
        The raw lines containing only the row body.
    """
    for i, ln in enumerate(lines, 1):
        if is_blank(ln) or is_comment(ln):
            continue
        yield parse_row(kind, ln, lineno=i)


__all__ = [
    # errors & data holders
    "JParseError", "DataType", "ParsedRow",
    # IO
    "iter_lines",
    # classifiers / parsers
    "is_comment", "is_blank", "parse_info", "parse_station",
    "parse_datatype_units", "parse_npoints", "parse_row",
    # generators
    "strip_nondata", "iter_info", "iter_rows",
]
