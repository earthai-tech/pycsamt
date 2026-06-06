# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.tdem.reader
===================

Unified TEM file reader.

:class:`TEMReader` wraps every format-specific function in
:mod:`pycsamt.tdem.io` behind a single, consistent interface:

* Auto-detect the file format from its extension and magic bytes.
* Carry shared acquisition defaults (current, loop geometry, …) so
  you only set them once.
* Propagate ``verbose`` / ``logger`` to every call.
* Optionally cache results for later inspection.

The individual :func:`~pycsamt.tdem.io.read_xyz`,
:func:`~pycsamt.tdem.io.read_geosoft_dat`, … functions remain
available for direct use; :class:`TEMReader` just calls them —
no logic is duplicated.

Quick start
-----------
::

    from pycsamt.tdem import TEMReader

    # configure once
    reader = TEMReader(
        current=8.0,
        loop_side=200.0,
        data_unit="nV/Am2",
        verbose=1,
    )

    # auto-detect format
    soundings = reader.read("survey.dat")           # Geosoft DAT
    soundings = reader.read("profile.tem")          # AMIRA or WalkTEM
    soundings = reader.read("sounding.avg")         # Zonge GDP

    # explicit format
    soundings = reader.read("myfile.txt", fmt="xyz",
                             time_unit="ms", data_unit="nT/s")

    # named method (same as read with fmt="geosoft")
    soundings = reader.read_geosoft_dat("survey.dat")

    # inspect cached results
    print(reader.results)
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Union

from ..api.property import MetadataMixin, PyCSAMTObject
from ..log.logger import get_logger
from . import io
from ._base import TEMSounding

__all__ = ["TEMReader"]

Pathish = Union[str, Path]


# ─────────────────────────────────────────────────────────────────────────────
# Format registry
# ─────────────────────────────────────────────────────────────────────────────

# Each entry: format_name → (reader_callable, [file_extensions])
# Extensions are lower-cased without the leading dot.
_REGISTRY: dict[str, tuple] = {
    "xyz":      (io.read_xyz,          ["xyz", "csv", "txt"]),
    "temavg":   (io.read_temavg,       ["avg"]),
    "tem_z":    (io.read_tem_z,        ["z"]),
    "tem_log":  (io.read_tem_log,      ["log"]),
    "geosoft":  (io.read_geosoft_dat,  ["dat"]),
    "amira":    (io.read_amira,        ["tem"]),
    "zonge":    (io.read_zonge,        ["avg", "tem"]),
    "walkttem": (io.read_walkttem,     ["tem"]),
}

# ── Extension → candidate format names ───────────────────────────────────────
# When more than one format shares an extension (e.g. .avg, .tem), the
# auto-detector uses magic-byte inspection to disambiguate.
_EXT_MAP: dict[str, list[str]] = {}
for _fmt, (_fn, _exts) in _REGISTRY.items():
    for _ext in _exts:
        _EXT_MAP.setdefault(_ext, []).append(_fmt)


# ─────────────────────────────────────────────────────────────────────────────
# Format auto-detection
# ─────────────────────────────────────────────────────────────────────────────

def _first_lines(path: Path, n: int = 10) -> list[str]:
    """Return the first *n* non-blank, non-comment lines of *path*."""
    lines: list[str] = []
    try:
        with open(path, errors="replace") as fh:
            for raw in fh:
                stripped = raw.strip()
                if stripped:
                    lines.append(stripped)
                if len(lines) >= n:
                    break
    except OSError:
        pass
    return lines


def _detect_format(path: Path) -> str:
    """Infer the TEM file format from extension and magic bytes.

    Returns a key from :data:`_REGISTRY`.  Falls back to ``"xyz"`` if
    the extension is not recognised.
    """
    ext = path.suffix.lstrip(".").lower()
    candidates = _EXT_MAP.get(ext, [])

    # Unambiguous by extension
    if len(candidates) == 1:
        return candidates[0]

    # Nothing recognised → default
    if not candidates:
        return "xyz"

    # Disambiguation by magic bytes
    lines = _first_lines(path, n=15)
    joined = " ".join(lines).upper()

    if ext == "avg":
        # Zonge TEMAVG processed format starts with `\ TEMAVG`
        if any(ln.upper().startswith(("\\", "/ TEMAVG", "\\ TEMAVG")) or
               "TEMAVG" in ln.upper()
               for ln in lines[:3]):
            # Could still be Zonge GDP (starts with *)
            if any(ln.startswith("*") for ln in lines[:3]):
                return "zonge"
            return "temavg"
        return "zonge"

    if ext == "tem":
        # WalkTEM / Aarhus Workbench uses TXOPERATION or GATESET keywords
        if any(kw in joined
               for kw in ("TXOPERATION", "GATESET", "GENERALHEADER",
                          "DATASTATUS", "END_TXOPERATION")):
            return "walkttem"
        # AMIRA / EMIT uses TRANSMITTER_ keywords or `;` comments
        if any(kw in joined
               for kw in ("TRANSMITTER_AREA", "TRANSMITTER_CURRENT",
                          "GATE_TIMES", "DATA_UNIT")) or any(
               ln.startswith(";") for ln in lines[:5]):
            return "amira"
        # Zonge GDP uses `*` comment lines
        if any(ln.startswith("*") for ln in lines[:5]):
            return "zonge"
        return "amira"   # safest default for unknown .tem

    if ext in ("dat",):
        # Geosoft DAT: comment lines start with `/`
        if any(ln.startswith("/") for ln in lines[:5]):
            return "geosoft"
        return "geosoft"   # .dat with no `/` lines → still try geosoft

    if ext in ("z",):
        return "tem_z"

    if ext in ("log",):
        return "tem_log"

    return "xyz"


# ─────────────────────────────────────────────────────────────────────────────
# TEMReader
# ─────────────────────────────────────────────────────────────────────────────

class TEMReader(PyCSAMTObject, MetadataMixin):
    """Unified TEM file reader with format auto-detection.

    :class:`TEMReader` is a thin orchestration layer over the
    individual format functions in :mod:`pycsamt.tdem.io`.  It adds:

    * **Auto-detection** — ``read(path)`` infers the file format from
      the extension and magic bytes.
    * **Shared defaults** — acquisition parameters set on the instance
      (``current``, ``tx_area``, …) are injected into every call so
      you only specify them once.
    * **Per-call override** — keyword arguments passed to ``read()`` or
      the named ``read_*`` methods always take priority over instance
      defaults.
    * **Verbose / logger** — set ``verbose=1`` for INFO messages or
      ``verbose=2`` for DEBUG.  Pass a custom :class:`logging.Logger`
      via the ``logger`` parameter.
    * **Result cache** — set ``store=True`` to accumulate all results
      in :attr:`results` keyed by file path.

    Parameters
    ----------
    current : float or None
        Transmitter current in Amperes.
    tx_area : float or None
        Transmitter loop area in m².
    loop_side : float or None
        Side of a square transmitter loop in m.
    loop_radius : float or None
        Radius of a circular transmitter loop in m.
    rx_area : float
        Receiver coil area in m².  Default 1.0.
    rx_turns : int
        Receiver coil turns.  Default 1.
    data_unit : str
        Default data column units.  Default ``"nV/Am2"``.
    data_type : str
        Default data type.  Default ``"dBdt"``.
    gate_times_unit : str
        Default gate-time unit.  Default ``"ms"``.
    store : bool
        Accumulate results in :attr:`results` when ``True``.
        Default ``False``.
    verbose : int
        Verbosity level: 0 = silent, 1 = INFO, 2 = DEBUG.
    logger : logging.Logger or None
        Custom logger.  If ``None`` a module-level logger is used.

    Examples
    --------
    >>> from pycsamt.tdem import TEMReader
    >>> reader = TEMReader(current=8.0, loop_side=200.0, verbose=1)
    >>> soundings = reader.read("survey.dat")     # Geosoft DAT
    >>> soundings = reader.read("profile.tem")    # AMIRA / WalkTEM
    >>> soundings = reader.read_zonge("run.avg")  # explicit format
    """

    # Fields shown by PyCSAMTObject.__repr__
    __repr_fields__ = (
        "current", "tx_area", "loop_side", "loop_radius",
        "rx_area", "rx_turns", "data_unit", "data_type",
        "gate_times_unit", "store", "verbose",
    )
    __repr_exclude__ = {"logger", "_log", "_results"}

    def __init__(
        self,
        *,
        current: float | None = None,
        tx_area: float | None = None,
        loop_side: float | None = None,
        loop_radius: float | None = None,
        rx_area: float = 1.0,
        rx_turns: int = 1,
        data_unit: str = "nV/Am2",
        data_type: str = "dBdt",
        gate_times_unit: str = "ms",
        store: bool = False,
        verbose: int = 0,
        logger: logging.Logger | None = None,
    ) -> None:
        self.current         = current
        self.tx_area         = tx_area
        self.loop_side       = loop_side
        self.loop_radius     = loop_radius
        self.rx_area         = rx_area
        self.rx_turns        = rx_turns
        self.data_unit       = data_unit
        self.data_type       = data_type
        self.gate_times_unit = gate_times_unit
        self.store           = store
        self.verbose         = verbose
        self._log: logging.Logger = logger or get_logger(__name__)
        self._results: dict[str, Any] = {}

    # ── Logging helpers ───────────────────────────────────────────────────────

    def _info(self, msg: str) -> None:
        if self.verbose >= 1:
            self._log.info(msg)

    def _debug(self, msg: str) -> None:
        if self.verbose >= 2:
            self._log.debug(msg)

    def _warn(self, msg: str) -> None:
        self._log.warning(msg)

    # ── Defaults / kwargs merging ─────────────────────────────────────────────

    def _defaults(self) -> dict[str, Any]:
        """Return non-None instance-level acquisition defaults."""
        keys = (
            "current", "tx_area", "loop_side", "loop_radius",
            "rx_area", "rx_turns", "data_unit", "data_type",
            "gate_times_unit",
        )
        return {k: getattr(self, k)
                for k in keys
                if getattr(self, k, None) is not None}

    def _merge(self, **call_kwargs: Any) -> dict[str, Any]:
        """Merge instance defaults with per-call kwargs.

        Per-call values always win; ``None`` call values are ignored so
        they do not shadow a real instance default.
        """
        merged = self._defaults()
        merged.update(
            {k: v for k, v in call_kwargs.items() if v is not None}
        )
        return merged

    # ── Registry access ───────────────────────────────────────────────────────

    @property
    def formats(self) -> list[str]:
        """Return the list of recognised format names."""
        return list(_REGISTRY)

    # ── Result cache ──────────────────────────────────────────────────────────

    @property
    def results(self) -> dict[str, Any]:
        """Mapping of path → read result (populated when ``store=True``)."""
        return dict(self._results)

    def clear(self) -> None:
        """Clear the result cache."""
        self._results.clear()

    # ── Main dispatch ─────────────────────────────────────────────────────────

    def read(
        self,
        path: Pathish,
        *,
        fmt: str | None = None,
        **kwargs: Any,
    ) -> Any:
        """Read a TEM file, auto-detecting its format when *fmt* is omitted.

        Parameters
        ----------
        path : str or Path
            Path to the TEM data file.
        fmt : str or None
            Explicit format key.  One of ``"xyz"``, ``"temavg"``,
            ``"tem_z"``, ``"tem_log"``, ``"geosoft"``, ``"amira"``,
            ``"zonge"``, ``"walkttem"``.  When ``None`` the format is
            inferred automatically.
        **kwargs
            Passed to the underlying reader function after merging with
            instance defaults (per-call values win).

        Returns
        -------
        Any
            Whatever the underlying reader returns:
            a single :class:`~pycsamt.tdem.TEMSounding`,
            a :class:`list` of :class:`~pycsamt.tdem.TEMSounding`,
            a :class:`~pycsamt.tdem.avg.TEMAVG`, etc.

        Raises
        ------
        FileNotFoundError
        ValueError
            If *fmt* is given but not in :attr:`formats`.
        """
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"File not found: {p}")

        if fmt is not None and fmt not in _REGISTRY:
            raise ValueError(
                f"Unknown format {fmt!r}.  "
                f"Valid formats: {self.formats}"
            )

        resolved_fmt = fmt or _detect_format(p)
        reader_fn, _ = _REGISTRY[resolved_fmt]

        merged = self._merge(**kwargs)
        self._info(f"Reading {p.name!r} as {resolved_fmt!r}")
        self._debug(f"  full path: {p}  kwargs: {merged}")

        try:
            result = reader_fn(p, **merged)
        except TypeError as exc:
            # Some readers (temavg, tem_z, tem_log) do not accept
            # acquisition kwargs; retry without them.
            self._debug(
                f"  retrying {resolved_fmt!r} without acquisition kwargs "
                f"({exc})"
            )
            safe_keys = {"verbose", "logger"}
            safe_kw   = {k: v for k, v in merged.items() if k in safe_keys}
            result = reader_fn(p, **safe_kw)

        n = len(result) if hasattr(result, "__len__") else 1
        self._info(f"  {n} record(s) loaded from {p.name!r}")

        if self.store:
            self._results[str(p)] = result

        return result

    # ── Named format methods ──────────────────────────────────────────────────
    # Each delegates to read() with an explicit fmt= so the caller benefits
    # from defaults merging and caching without any duplicated logic.

    def read_xyz(
        self,
        path: Pathish,
        **kwargs: Any,
    ) -> TEMSounding:
        """Read a generic XYZ / CSV TEM sounding.

        See :func:`pycsamt.tdem.io.read_xyz` for parameter details.
        """
        return self.read(path, fmt="xyz", **kwargs)

    def read_temavg(
        self,
        path: Pathish,
        **kwargs: Any,
    ):
        """Read a Zonge TEMAVG processed ``.AVG`` file.

        See :func:`pycsamt.tdem.io.read_temavg` for parameter details.
        """
        return self.read(path, fmt="temavg", **kwargs)

    def read_tem_z(
        self,
        path: Pathish,
        **kwargs: Any,
    ):
        """Read a Zonge TEMAVG contour ``.Z`` file.

        See :func:`pycsamt.tdem.io.read_tem_z` for parameter details.
        """
        return self.read(path, fmt="tem_z", **kwargs)

    def read_tem_log(
        self,
        path: Pathish,
        **kwargs: Any,
    ):
        """Read a Zonge TEMAVG processing ``.LOG`` file.

        See :func:`pycsamt.tdem.io.read_tem_log` for parameter details.
        """
        return self.read(path, fmt="tem_log", **kwargs)

    def read_geosoft_dat(
        self,
        path: Pathish,
        **kwargs: Any,
    ) -> list[TEMSounding]:
        """Read a Geosoft Oasis Montaj ``.dat`` file.

        See :func:`pycsamt.tdem.io.read_geosoft_dat` for parameter details.
        """
        return self.read(path, fmt="geosoft", **kwargs)

    def read_amira(
        self,
        path: Pathish,
        **kwargs: Any,
    ) -> list[TEMSounding]:
        """Read an AMIRA / EMIT ``.tem`` file.

        See :func:`pycsamt.tdem.io.read_amira` for parameter details.
        """
        return self.read(path, fmt="amira", **kwargs)

    def read_zonge(
        self,
        path: Pathish,
        **kwargs: Any,
    ) -> list[TEMSounding]:
        """Read a Zonge GDP ``.avg`` / ``.tem`` sounding file.

        See :func:`pycsamt.tdem.io.read_zonge` for parameter details.
        """
        return self.read(path, fmt="zonge", **kwargs)

    def read_walkttem(
        self,
        path: Pathish,
        **kwargs: Any,
    ) -> list[TEMSounding]:
        """Read a WalkTEM / Aarhus Workbench ``.tem`` file.

        See :func:`pycsamt.tdem.io.read_walkttem` for parameter details.
        """
        return self.read(path, fmt="walkttem", **kwargs)

    # ── Batch read ────────────────────────────────────────────────────────────

    def read_batch(
        self,
        paths: list[Pathish],
        *,
        fmt: str | None = None,
        **kwargs: Any,
    ) -> dict[str, Any]:
        """Read multiple files and return a mapping of path → result.

        Parameters
        ----------
        paths : list of path-like
            File paths to read.
        fmt : str or None
            Force a single format for all files.  ``None`` → per-file
            auto-detection.
        **kwargs
            Forwarded to every :meth:`read` call.

        Returns
        -------
        dict[str, Any]
            Mapping from ``str(path)`` to the read result.
        """
        out: dict[str, Any] = {}
        for p in paths:
            key = str(p)
            self._info(f"Batch reading {key!r} …")
            out[key] = self.read(p, fmt=fmt, **kwargs)
        return out

    # ── Configure ─────────────────────────────────────────────────────────────

    def configure(self, **kwargs: Any) -> "TEMReader":
        """Update acquisition defaults in-place.

        Returns *self* for chaining::

            reader.configure(current=6.6, loop_side=40.0).read("file.tem")

        Parameters
        ----------
        **kwargs
            Any writable :class:`TEMReader` attribute.
        """
        _valid = {
            "current", "tx_area", "loop_side", "loop_radius",
            "rx_area", "rx_turns", "data_unit", "data_type",
            "gate_times_unit", "store", "verbose",
        }
        for k, v in kwargs.items():
            if k not in _valid:
                raise AttributeError(
                    f"TEMReader has no configurable attribute {k!r}.  "
                    f"Valid: {sorted(_valid)}"
                )
            setattr(self, k, v)
        return self

    # ── Summary ───────────────────────────────────────────────────────────────

    def summary(self, **kw) -> str:
        lines = [
            "TEMReader",
            f"  current         = {self.current!r}",
            f"  tx_area         = {self.tx_area!r}",
            f"  loop_side       = {self.loop_side!r}",
            f"  loop_radius     = {self.loop_radius!r}",
            f"  rx_area         = {self.rx_area}",
            f"  rx_turns        = {self.rx_turns}",
            f"  data_unit       = {self.data_unit!r}",
            f"  data_type       = {self.data_type!r}",
            f"  gate_times_unit = {self.gate_times_unit!r}",
            f"  store           = {self.store}",
            f"  verbose         = {self.verbose}",
            f"  formats         = {self.formats}",
            f"  cached results  = {len(self._results)}",
        ]
        return "\n".join(lines)
