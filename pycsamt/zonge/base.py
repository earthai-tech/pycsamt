#       Author: LKouadio <etanoyau@gmail.com>
#       License: LGPL-3.0-or-later
"""
pycsamt.zonge.base
==================

Foundational building blocks for the *Zonge* sub-package.

Goals
-----
- Provide a **single, flexible base** every component can
  inherit from (e.g., Station, Frequency, Resistivity…).
- Offer a robust **AVGFrame** holder for data + metadata.
- Define a clear **component contract**:
    • ``from_avg(...)``  → build from an AVG file/table
    • ``read(...)``       → parse & populate attributes
    • ``write()``         → reconstruct a textual AVG block
- Be agnostic to legacy/modern files.  Upstream loaders
  (e.g., ``utils.load_avg``) normalize columns; components
  consume the tidy frame.

Design notes
------------
- All public objects keep string widths short for readability.
- Columns are lower-case canonical names (e.g., ``'freq'``,
  ``'rho'``, ``'phase'``, ``'comp'``…).
- Components target *subsets* of the frame and may expose
  additional, derived attributes internally.
"""

from __future__ import annotations

import json
import warnings
from abc import ABC, abstractmethod
from collections.abc import Mapping, MutableMapping, Sequence
from dataclasses import asdict, field
from datetime import datetime, timezone
from pathlib import Path
from typing import (
    Any,
    Literal,
)

import numpy as np
import pandas as pd

from ..compat.python import dc
from ..exceptions import AvgDataError
from ..log.logger import get_logger
from ._transfer import LegacyAVGBase
from .schema import ALL_ALIASES
from .utils import (
    _standardise_columns,
    classify_avg_format,
    load_avg,
)

logger = get_logger(__name__)

__all__ = [
    "FieldAliases",
    "AvgRow",
    "AVGFrame",
    "AVGComponentBase",
    "guess_kind_from_df",
]


@dc(slots=True)
class AVGFrame:
    r"""A container for a tidy AVG table and its metadata.

    This dataclass acts as a standardized wrapper for the data
    and metadata parsed from a Zonge AVG file. It ensures that
    data passed between different parts of the package is
    consistent and self-contained.

    The `__post_init__` hook automatically standardizes the
    column names of the incoming DataFrame, ensuring that any
    instance of `AVGFrame` always contains a clean, canonical
    dataset.

    Parameters
    ----------
    data : pandas.DataFrame
        A DataFrame containing the tabular data from an AVG file.
        Column names are automatically normalized to the package's
        internal canonical schema upon initialization.
    meta : dict, optional
        A dictionary holding the header metadata (e.g., survey
        parameters, units) extracted from the ``$keyword=value``
        lines in the AVG file.
    source : pathlib.Path, optional
        The filesystem path to the original ``.avg`` file, used for
        provenance and logging.

    Examples
    --------
    >>> import pandas as pd
    >>> from pycsamt.zonge.base import AVGFrame
    >>> df = pd.DataFrame({"Resistivity": [100], "Freq": [1024]})
    >>> frame = AVGFrame(data=df, meta={"Survey.Type": "CSAMT"})
    >>> print(frame.data.columns)
    Index(['rho', 'freq'], dtype='object')
    >>> print(frame.meta["Survey.Type"])
    CSAMT

    See Also
    --------
    pycsamt.zonge.avg.AVG : The main user-facing data object.
    pycsamt.zonge.info.DataInfo : Aggregator that uses AVGFrame.
    """

    data: pd.DataFrame
    meta: dict[str, Any] = field(default_factory=dict)
    source: Path | None = None

    def __post_init__(self):
        self.data = _standardise_columns(self.data)

    @property
    def nrows(self) -> int:
        """Row count."""
        return int(len(self.data))

    @property
    def columns(self) -> tuple[str, ...]:
        """Column labels (tuple for immutability)."""
        return tuple(map(str, self.data.columns))

    def copy(self) -> AVGFrame:
        """Deep copy the frame (safe for mutation)."""
        return AVGFrame(
            data=self.data.copy(deep=True),
            meta=dict(self.meta),
            source=self.source,
        )

    def to_json(self, *, orient: str = "records", indent: int = 0) -> str:
        """Serialise *data* to JSON (metadata excluded)."""
        return self.data.to_json(orient=orient, indent=indent)

    def meta_as_json(self, *, indent: int = 0) -> str:
        """Serialise metadata to JSON."""
        # pandas 2.x removed pd.io.json.dumps/loads; use stdlib json
        # (default=str covers numpy scalars and paths in the metadata)
        return json.dumps(dict(self.meta), indent=indent or None, default=str)

    def asdict(self) -> dict[str, Any]:
        """Plain dict for diagnostics / logging."""
        return {
            "data": json.loads(self.to_json()),
            "meta": dict(self.meta),
            "source": str(self.source) if self.source else None,
        }

    def __str__(self) -> str:
        src = f", file='{self.source.name}'" if self.source else ""
        cols = ", ".join(self.columns[:6])
        tail = "…" if len(self.columns) > 6 else ""
        return (
            f"AVGFrame[{self.nrows}×{len(self.columns)}]" f"{src} cols=[{cols}{tail}]"
        )

    def __repr__(self) -> str:
        keys = list(self.meta)[:4]
        more = "…" if len(self.meta) > 4 else ""
        return (
            f"AVGFrame(nrows={self.nrows}, "
            f"ncols={len(self.columns)}, "
            f"meta_keys={keys}{more})"
        )


class FieldAliases:
    r"""Dynamically provides all known aliases for canonical names.

    This class serves as a convenient, centralized namespace for
    accessing all known variations of a given canonical column
    name. It is designed to be a single source of truth for
    developers who need to work with different possible column
    names.

    Notes
    -----
    The primary purpose of this class is to improve code
    readability and maintainability by avoiding the use of "magic
    strings" for column names. Instead of writing
    ``df[("rho", "Resistivity")]``, a developer can use the more
    explicit ``df[FieldAliases.rho]``.

    The class attributes are not hardcoded. They are generated
    dynamically by introspecting the ``ALL_ALIASES`` dictionary
    from the :mod:`~pycsamt.zonge.schema` module. This ensures
    that `FieldAliases` is always in sync with the master data
    schema and adheres to the DRY (Don't Repeat Yourself)
    principle.

    Examples
    --------
    You can access the tuple of aliases for any canonical name
    directly as a class attribute:

    >>> from pycsamt.zonge.base import FieldAliases
    >>> # Get all known names for apparent resistivity
    >>> FieldAliases.rho
    ('ARes.mag', 'Resistivity')

    >>> # Get all known names for H-field phase error
    >>> FieldAliases.s_hphz
    ('B.perr', 'H.perr', 'sHphz')

    See Also
    --------
    pycsamt.zonge.schema.ALL_ALIASES : The source dictionary used
        to build this class.
    """

    # Dynamically populate the class attributes from the schema
    for canon, aliases in ALL_ALIASES.items():
        # Ensure the attribute name is a valid Python identifier
        attr_name = canon.replace(".", "_").replace("%", "pct")
        locals()[attr_name] = aliases

    # Statically define a few key ones for type hinting/clarity
    station: tuple[str, ...] = ALL_ALIASES.get("station", ())
    freq: tuple[str, ...] = ALL_ALIASES.get("freq", ())
    comp: tuple[str, ...] = ALL_ALIASES.get("comp", ())
    rho: tuple[str, ...] = ALL_ALIASES.get("rho", ())
    phase: tuple[str, ...] = ALL_ALIASES.get("phase", ())


@dc(slots=True)
class AvgRow:
    r"""A structured, format-agnostic representation of a single row.

    This dataclass provides a type-safe container for the data
    from a single measurement record in an AVG file. It serves as
    a convenient Data Transfer Object (DTO) for operations that
    require iterating over individual data points, such as testing,
    serialization to JSON, or specific calculations.

    Attributes
    ----------
    station : int or float
        The station identifier for the measurement.
    freq : float
        The frequency at which the measurement was taken (Hz).
    comp : str
        The component label (e.g., 'ExHy'). Defaults to 'ExHy' if
        not provided.
    amps : float, optional
        The transmitter current amplitude (A).
    emag, ephz : float, optional
        The magnitude and phase of the electric field.
    hmag, hphz : float, optional
        The magnitude and phase of the magnetic field.
    rho, phase : float, optional
        The calculated apparent resistivity (ohm-m) and impedance
        phase (mrad).
    e_err, h_err, rho_err : float, optional
        The relative percent error for E-field, H-field, and
        resistivity magnitudes.
    e_perr, h_perr, z_perr : float, optional
        The standard deviation (error) for E-field phase, H-field
        phase, and impedance phase.
    """

    station: int | float
    freq: float
    comp: str
    amps: float | None = None
    emag: float | None = None
    ephz: float | None = None
    hmag: float | None = None
    hphz: float | None = None
    rho: float | None = None
    phase: float | None = None

    # Optional quality fields (kept simple)
    e_err: float | None = None  # relative %
    e_perr: float | None = None  # mrad
    h_err: float | None = None
    h_perr: float | None = None
    rho_err: float | None = None
    z_perr: float | None = None

    def __post_init__(self) -> None:
        self.comp = str(self.comp).strip() or "ExHy"

    def asdict(self) -> dict[str, Any]:
        d = asdict(self)
        return d

    def __str__(self) -> str:
        r_fmt = (
            f"{self.rho:.1f}"
            if isinstance(self.rho, (int, float)) and np.isfinite(self.rho)
            else "nan"
        )
        return (
            f"AvgRow(stn={self.station}, f={self.freq:g} Hz, "
            f"comp={self.comp}, rho={r_fmt})"
        )

    __repr__ = __str__


class AVGComponentBase(ABC):
    r"""Abstract base class for a single AVG data component.

    This class defines the fundamental "contract" that all data
    components (e.g., Station, Resistivity, Phase) must follow.
    It ensures that every component can be initialized from a
    standardized data source and can serialize itself back into a
    textual format.

    Subclasses are expected to be specialized containers that manage
    a specific subset of columns from the main AVG DataFrame.

    Attributes
    ----------
    required : set[str]
        A class attribute specifying the set of canonical column
        names that must be present in the source DataFrame for the
        component to be initialized.
    provides : set[str]
        A class attribute specifying the set of canonical column
        names that this component is responsible for managing or
        creating.
    _frame : pandas.DataFrame
        A protected attribute holding the component's data as a
        tidy DataFrame with canonical column names.
    _meta : dict
        A protected attribute holding relevant metadata.

    Methods
    -------
    read(source, meta=None)
        Abstract method to parse a source DataFrame and populate
        the component's internal state.
    write()
        Abstract method to serialize the component's data into a
        sequence of text lines suitable for an AVG file.
    from_avg(avg, meta=None, **kws)
        A classmethod factory that provides a convenient way to
        build a component from various sources, including file
        paths or in-memory DataFrames.

    Notes
    -----
    The design of this base class ensures that all components
    operate on a consistent, tidy data structure. The `from_avg`
    factory method is particularly important, as it handles the
    automatic transformation of legacy AVG files into the modern,
    canonical format before passing the data to the component's
    `read` method. This makes all components inherently agnostic
    to the original file format.
    """

    # Sets listing which columns must be present / will be produced.
    required: set[str] = set()
    provides: set[str] = set()

    def __init__(
        self,
        data: pd.DataFrame | None = None,
        meta: MutableMapping[str, Any] | None = None,
        *,
        name: str | None = None,
        verbose=0,
    ) -> None:
        self._frame: pd.DataFrame = (
            data.copy(deep=True) if data is not None else pd.DataFrame()
        )
        self._meta: dict[str, Any] = dict(meta or {})
        self._name: str = name or self.__class__.__name__
        self.verbose = verbose
        self._logger = get_logger(self._name)

    @classmethod
    def from_avg(
        cls,
        avg: str
        | Path
        | AVGFrame
        | pd.DataFrame
        | tuple[pd.DataFrame, Mapping[str, Any]],
        *,
        meta: Mapping[str, Any] | None = None,
        **kws,
    ) -> AVGComponentBase:
        """
        Build a component from a path / AVGFrame / dataframe.

        The method accepts:
          • path → uses ``utils.load_avg`` (modern or legacy)
          • ``AVGFrame`` → uses its data/meta directly
          • ``(df, meta)`` tuple
          • bare ``df`` + explicit ``meta`` kwarg
        """
        df: pd.DataFrame
        m: Mapping[str, Any]

        if isinstance(avg, (str, Path)):
            # When loading from a file, classify it first
            path = Path(avg)
            lines = path.read_text(errors="replace").splitlines()
            kind = classify_avg_format(lines)
            df, m = load_avg(path)

            # If it's legacy, transform it to modern structure
            if kind == 1:
                transformer = LegacyAVGBase()
                ds = transformer.from_dataframe(df, meta=m)
                df = ds.to_dataframe().reset_index()
                m = ds.attrs

            frame = AVGFrame(df, dict(m), path)

        elif isinstance(avg, AVGFrame):
            frame = avg
        elif isinstance(avg, tuple) and len(avg) == 2:
            df, m = avg
            frame = AVGFrame(df, dict(m))
        elif isinstance(avg, pd.DataFrame):
            frame = AVGFrame(avg, dict(meta or {}))
        else:
            raise TypeError(
                "from_avg expects Path|AVGFrame|DataFrame|" "(DataFrame, meta) tuple."
            )

        obj = cls()
        try:
            obj.read(frame.data, frame.meta)
        except TypeError:
            # Fallback for older components
            obj.read(frame.data)

        return obj

    @abstractmethod
    def read(
        self,
        source: pd.DataFrame,
        meta: Mapping[str, Any] | None = None,
    ) -> None:
        """
        Parse **source** and mutate component state.

        Implementations should:
          1) validate required columns using ``_require(...)``
          2) copy/select needed columns into ``self._frame``
          3) set/merge metadata into ``self._meta``
        """
        ...

    @abstractmethod
    def write(self) -> Sequence[str]:
        """
        Serialise the component to **text** lines.  Implementations
        may delegate to ``_write_csv_block`` for a consistent block
        format (header + CSV).  The base does not write to disk.
        """
        ...

    @property
    def frame(self) -> pd.DataFrame:
        """Read-only view of the component table."""
        return self._frame

    @property
    def meta(self) -> dict[str, Any]:
        """Free-form metadata."""
        return self._meta

    @property
    def name(self) -> str:
        """Short human-readable component name."""
        return self._name

    @property
    def shape(self) -> tuple[int, int]:
        """Table shape (rows, cols)."""
        return (int(len(self._frame)), int(self._frame.shape[1] or 0))

    def asdict(self, *, include_meta: bool = True) -> dict[str, Any]:
        """Plain dict with data (+ meta optionally)."""
        d: dict[str, Any] = {"data": self._frame.to_dict("list")}
        if include_meta:
            d["meta"] = dict(self._meta)
        return d

    def to_json(self, *, indent: int = 0) -> str:
        """JSON serialiser for diagnostics."""
        return pd.io.json.dumps(self.asdict(), indent=indent)

    def _require(self, *cols: str) -> None:
        """
        Ensure columns are present; raise a clear error otherwise.
        """
        missing = [c for c in cols if c not in self._frame.columns]
        if missing:
            raise AvgDataError(f"{self._name}: missing columns {missing}")

    def _select(self, *cols: str) -> pd.DataFrame:
        """
        Return a *copy* with only selected columns (silently drops
        absent ones).  Useful when composing write-blocks.
        """
        keep = [c for c in cols if c in self._frame.columns]
        return self._frame.loc[:, keep].copy()

    def _write_csv_block(
        self,
        *,
        cols: Sequence[str],
        title: str | None = None,
        float_fmt: str = "%.6g",
        na_rep: str = "",
        include_meta: bool = True,
        stamp: bool = True,
    ) -> Sequence[str]:
        """
        Build a *text* block with an optional title/meta preamble
        and a CSV table of the selected ``cols``.

        This does **not** write to disk; callers own the I/O.
        """
        lines: list[str] = []

        # Title / banner (optional but nice for humans)
        if title:
            lines.append(f"\\ {title}")

        # $key=value meta lines (optional)
        if include_meta and self._meta:
            for k, v in self._meta.items():
                # Avoid multi-line explosions in headers
                v_str = str(v).replace("\n", " ").strip()
                lines.append(f"${k}={v_str}")

        # UTC stamp for provenance if desired
        if stamp:
            ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
            lines.append(f"$Written={ts}")

        # Guard against empty data
        if self._frame.empty:
            lines.append("")  # keep a blank for clean joins
            return lines

        # Emit CSV header + rows
        table = self._select(*cols)
        # Ensure stable column order as requested
        table = table.loc[:, list(cols)]
        csv = table.to_csv(index=False, float_format=float_fmt, na_rep=na_rep)
        # Separate header from CSV with a blank line for clarity
        lines.append("")
        lines.extend(csv.splitlines())
        return lines

    def __str__(self) -> str:
        r, c = self.shape
        cols = ", ".join(self._frame.columns[:6])
        tail = "…" if self._frame.shape[1] > 6 else ""
        return f"{self._name}[{r}×{c}] cols=[{cols}{tail}]"

    __repr__ = __str__


def guess_kind_from_df(
    df_or_frame: pd.DataFrame | AVGFrame,
    meta: Mapping[str, Any] | None = None,
    *,
    transform: bool = False,
    mode: Literal["strict", "soft"] = "soft",
    error: Literal["raise", "warn", "ignore"] = "raise",
    verbose: bool = False,
) -> int | tuple[pd.DataFrame, Mapping[str, Any], int]:
    r"""Infers the AVG format kind from a DataFrame's columns.

    This helper inspects the column names of a DataFrame to make
    an educated guess about whether it represents a legacy (kind-1)
    or modern (kind-2) AVG file. It can optionally transform
    legacy data into the modern structure.

    Parameters
    ----------
    df_or_frame : pd.DataFrame or AVGFrame
        The DataFrame or AVGFrame object to inspect.
    meta : mapping, optional
        An optional metadata dictionary. Required if `transform` is
        ``True`` and the source is a legacy DataFrame.
    transform : bool, default False
        If ``True`` and a legacy (kind-1) format is detected, this
        function will use the :class:`~.LegacyAVGBase` transformer
        to convert the data into the modern structure.
    error : {'raise', 'warn', 'ignore'}, default 'raise'
        Determines the behavior when a format cannot be reliably
        determined.
        - 'raise': Raises an `AvgFileError`.
        - 'warn': Issues a warning and defaults to kind-2.
        - 'ignore': Silently defaults to kind-2.
    verbose : bool, default False
        Controls logging output during transformation.

    Returns
    -------
    int or tuple
        - If `transform` is ``False``, returns an integer: ``1`` for
          legacy, ``2`` for modern.
        - If `transform` is ``True``, returns a tuple:
          ``(df, meta, kind)``, where `df` and `meta` are the
          (potentially transformed) data and metadata.

    Raises
    ------
    TypeError
        If the input is not a pandas DataFrame or AVGFrame.
    AvgFileError
        If `error` is 'raise' and the format cannot be determined.

    Notes
    -----
    The detection logic uses a multi-pass heuristic:
    1.  It first checks for columns with dot-notation (e.g.,
        'ARes.mag'), which are exclusive to modern files.
    2.  If none are found, it checks for legacy-specific QC
        column names (e.g., '%Emag', 'sPhz').
    3.  If neither is found, it checks for the presence of
        internal canonical QC names (e.g., 'pc_emag', 's_phz'),
        treating them as modern.
    4.  Finally, it falls back to a default controlled by the
        `error` parameter.
    """

    if isinstance(df_or_frame, AVGFrame):
        df = df_or_frame.data
        meta = meta or df_or_frame.meta
    elif isinstance(df_or_frame, pd.DataFrame):
        df = df_or_frame
    else:
        raise AvgDataError("Input must be a pandas DataFrame or an AVGFrame.")

    cols = set(df.columns)
    kind = 0

    # 1. Check for modern indicators (most reliable)
    if any("." in str(c) for c in cols):
        kind = 2
    else:
        # 2. Check for legacy indicators
        legacy_indicators = {
            "%Emag",
            "sEphz",
            "%Hmag",
            "sHphz",
            "%Rho",
            "sPhz",
        }
        if any(indicator in cols for indicator in legacy_indicators):
            kind = 1
        else:
            # 3. Check for already-standardized canonical names
            canonical_indicators = {
                "pc_emag",
                "s_ephz",
                "pc_hmag",
                "s_hphz",
                "pc_rho",
                "s_phz",
            }
            if any(c in canonical_indicators for c in cols):
                kind = 2  # Treat as structurally modern

    # 4. Handle cases where no clear indicators are found
    if kind == 0:
        if mode == "strict":
            msg = "Could not determine AVG kind from DataFrame columns."
            if error == "raise":
                raise AvgDataError(msg)
            elif error == "warn":
                warnings.warn(msg + " Defaulting to modern (kind-2).", stacklevel=2)
                kind = 2
            # 'ignore'

        kind = 2

    if transform and kind == 1:
        if verbose:
            logger.info("Transforming legacy DataFrame to modern.")
        transformer = LegacyAVGBase()
        ds = transformer.from_dataframe(df, meta=meta)
        return ds.to_dataframe().reset_index(), ds.attrs, 2

    elif transform:
        return df, meta or {}, kind

    return kind


__all__.extend(["LegacyAVGBase"])
