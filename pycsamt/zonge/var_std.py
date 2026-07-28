# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Phase-variation quality metrics:

* :class:`SPhz`   – stdev of Z phase (mrad or deg)
* :class:`SEphz`  – stdev of E phase (mrad or deg)
* :class:`SHphz`  – stdev of H phase (mrad or deg)

These are **scalar** QC variables defined per (station, freq, comp)
row. They are **not** 2×2 tensors, but provide a
:meth:`to_tensor_like` adapter that reuses :class:`TensorBase`
gridding when you want a (station × freq × 2 × 2) view with values
placed at the matching component slot and the rest as NaN.

Notes
-----
* Accepted column names are both **legacy** (e.g. ``sPhz``) and
  **modern** (e.g. ``Z.perr``). Case-insensitive matching is used.
* Units default to ``mrad`` via ``Unit.Phase`` in dataset attrs.
  You can call :meth:`convert_unit` to flip between ``mrad`` and
  ``deg`` safely.

"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

# from dataclasses import dataclass
from typing import (
    Any,
    ClassVar,
)

import numpy as np
import pandas as pd

from ..exceptions import AvgDataError
from ..log.logger import get_logger
from .base import AVGComponentBase
from .utils import find_and_rename_column, to_numeric_if_possible

logger = get_logger(__name__)

__all__ = ["SPhz", "SEphz", "SHphz", "PhaseSigma"]


class PhaseStdBase(AVGComponentBase):
    """
    Common implementation for phase stdev QC variables.

    Parameters
    ----------
    data : pandas.DataFrame | sequence, optional
        If a tidy DataFrame is given, we look for
        ``station, freq, comp`` and a suitable phase-stdev column
        (legacy or modern). If a vector-like is given, we create a
        minimal tidy frame with ``station=NaN`` and ``comp='ExHy'``.
    meta : mapping, optional
        Metadata to stash; ``Unit.Phase`` defaults to ``mrad``.

    Attributes
    ----------
    VAR_NAME : str
        Canonical internal column name used in this object.
    KEY_CANDIDATES : list[str]
        Accepted source column names (legacy + modern).
    LABEL : str
        Short label used for banners and debug prints.
    """

    # Class-level constants; NOT instance fields
    VAR_NAME: ClassVar[str] = ""
    # KEY_CANDIDATES: ClassVar[Tuple[str, ...]] = ()
    LABEL: ClassVar[str] = ""

    def __init__(
        self,
        data: pd.DataFrame | None = None,
        meta: Mapping[str, Any] | None = None,
        *,
        name: str | None = None,
        verbose: bool = False,
    ) -> None:
        # Explicitly call the parent's initializer
        super().__init__(data=data, meta=meta, name=name, verbose=verbose)

    def read(  # noqa: D401 (docstring above)
        self,
        source: pd.DataFrame | Sequence[float] | np.ndarray | pd.Series,
        meta: Mapping[str, Any] | None = None,
        **kws: Any,
    ) -> None:
        """
        Parse *source* and build an internal tidy frame with the
        canonical variable name :pyattr:`VAR_NAME`.

        If *source* is vector-like, construct a minimal tidy frame
        with columns ``station, freq, comp`` when possible (missing
        coords become NaN / 'ExHy').
        """
        if not self.VAR_NAME:
            raise RuntimeError(
                f"{self.__class__.__name__}: subclass constants not set."
            )

        self._meta = dict(meta or {})

        # ensure Unit.Phase presence for roundtrip/export
        self._meta.setdefault("Unit.Phase", "mrad")

        # vector-like path -
        if isinstance(source, (list, tuple, np.ndarray, pd.Series)):
            vec = pd.to_numeric(pd.Series(source), errors="coerce")
            df = pd.DataFrame({self.VAR_NAME: vec})
            df["station"] = kws.get("station", np.nan)
            df["freq"] = kws.get("freq", np.nan)
            df["comp"] = _norm_comp(kws.get("comp", "ExHy"))
            self._frame = df[["station", "freq", "comp", self.VAR_NAME]]
            return

        # dataframe path
        if not isinstance(source, pd.DataFrame):
            raise TypeError(
                f"{self.__class__.__name__}.read expects " "DataFrame or vector-like"
            )

        df = source.copy()
        df = find_and_rename_column(df, self.VAR_NAME)

        # locate a compatible source column
        # col = _find_col(df, self.KEY_CANDIDATES)
        # if col is None:
        #     raise AvgDataError(
        #         f"{self.LABEL}: no compatible phase-stdev column "
        #         f"found among {self.KEY_CANDIDATES!r}"
        #     )

        if self.VAR_NAME not in df.columns:
            df[self.VAR_NAME] = np.nan
            if self.verbose:
                logger.debug(
                    f"'{self.VAR_NAME}' not found in source. "
                    "Creating as empty column."
                )

        # df = df.rename(columns={col: self.VAR_NAME})

        # coords – inject conservative defaults if missing
        if "comp" not in df.columns:
            df["comp"] = "ExHy"
        if "station" not in df.columns:
            df["station"] = np.nan
        if "freq" not in df.columns:
            df["freq"] = np.nan

        # coerce numerics
        df["station"] = to_numeric_if_possible(df["station"])
        df["freq"] = pd.to_numeric(df["freq"], errors="coerce")
        df[self.VAR_NAME] = df[self.VAR_NAME].map(_to_num)

        self._frame = df[["station", "freq", "comp", self.VAR_NAME]]
        self._meta = dict(meta or {})
        self._meta.setdefault("Unit.Phase", "mrad")

    def write(self) -> list[str]:
        """
        Serialise to a compact CSV block with a small meta preamble.

        Returns
        -------
        list[str]
            Lines ready to prepend to a ``.avg`` block.
        """
        title = f"$Phase-Stdev · {self.LABEL}"
        # export a tiny meta banner with Unit.Phase for clarity
        meta = {"Unit.Phase": self._meta.get("Unit.Phase", "mrad")}
        tmp = self.__class__()  # ephemeral holder for helpers
        tmp._frame = self._frame.copy()
        tmp._meta = meta

        return tmp._write_csv_block(
            cols=["station", "freq", "comp", self.VAR_NAME],
            title=title,
            include_meta=True,
            stamp=True,
        )

    def to_xarray(
        self,
        *,
        coords: Sequence[str] = ("station", "freq", "comp"),
        attrs: dict[str, Any] | None = None,
    ):
        r"""
        Convert the component table into an :class:`xarray.Dataset`
        with a single data variable named ``VAR_NAME`` arranged on
        the requested coordinate grid.

        The output uses the intersection of the requested *coords*
        that are present in the component frame, preserving the
        order of dimensions (default:
        :math:`\text{station} \rightarrow \text{freq} \rightarrow
        \text{comp}`).

        Duplicate rows with identical coordinates are averaged so
        each grid cell contains a unique value.

        Parameters
        ----------
        coords : sequence of str, optional
            Preferred coordinate columns in priority order. Only
            those present in the frame are used. By default:
            ``("station", "freq", "comp")``.
        attrs : mapping, optional
            Extra attributes to merge into dataset-level metadata.

        Returns
        -------
        xarray.Dataset
            Dataset with dimensions given by the present subset of
            *coords* and one data variable named ``self.VAR_NAME``.

        Notes
        -----
        If the variable represents phase standard deviation
        (i.e., ``VAR_NAME`` ends with ``"phz"``), the attribute
        ``"Unit.Phase"`` is ensured and defaults to ``"mrad"`` if
        missing. This mirrors common CSAMT/CSAVGW conventions.
        """
        if self._frame.empty:
            raise AvgDataError("empty frame; nothing to export.")

        # Work on a copy; we'll normalize a few columns below.
        df = self._frame.copy()

        # Ensure a 'comp' column exists so we always get a comp dim.
        if "comp" not in df.columns:
            df["comp"] = "ExHy"

        # Determine which coord columns we actually have.
        idx_cols = [c for c in coords if c in df.columns]
        if not idx_cols:
            raise AvgDataError(
                f"no coordinate columns found; expected any of {coords!r}"
            )

        # Light type normalization:
        # - station may be float or label → avoid coercion that would
        #   mangle labels; keep as-is when not numeric.
        # - freq should be numeric for sorting/gridding.
        if "station" in idx_cols:
            df["station"] = to_numeric_if_possible(df["station"])
        if "freq" in idx_cols:
            df["freq"] = pd.to_numeric(df["freq"], errors="coerce")

        # Provide a stable, interpretable order for 'comp'. We put a
        # canonical list first and append any unexpected labels.
        if "comp" in idx_cols:
            canon = [
                "ExHy",
                "ExHx",
                "EyHx",
                "EyHy",
                "Zxx",
                "Zxy",
                "Zyx",
                "Zyy",
                "Zvec",
                "Zdet",
            ]
            present = pd.Series(df["comp"].astype(str).unique()).tolist()
            extras = [c for c in present if c not in canon]
            cats = canon + extras
            df["comp"] = pd.Categorical(
                df["comp"].astype(str), categories=cats, ordered=True
            )

        # We will export only the single data variable of interest.
        var = self.VAR_NAME
        if var not in df.columns:
            raise AvgDataError(
                f"{var!r} not found in frame columns {list(df.columns)!r}"
            )

        # Deduplicate: average numeric values across identical coords.
        dup = df.duplicated(subset=idx_cols, keep=False)
        if bool(dup.any()):
            gb = df.groupby(idx_cols, sort=True, dropna=False)
            dfv = gb[[var]].mean()
            tidy = dfv.reset_index()
        else:
            tidy = df.sort_values(idx_cols, kind="mergesort")

        # Build the Dataset: MultiIndex → dense grid.
        ds = tidy.set_index(idx_cols)[[var]].to_xarray()

        # Order dimensions as requested, dropping any that are absent.
        dim_order = [d for d in coords if d in ds.dims]
        ds = ds.transpose(*dim_order)

        # Compose/merge attributes. Start from the component meta,
        # ensure unit hints for phase variables, then layer user attrs.
        merged = dict(self._meta)
        if var.lower().endswith("phz"):
            merged.setdefault("Unit.Phase", "mrad")
        if attrs:
            merged.update(attrs)
        ds.attrs.update(merged)

        return ds

    def to_tensor_like(
        self,
        *,
        align: str = "union",
        station: float | int | None = None,
        fill_value: float = np.nan,
        sort_freq: bool = True,
        agg: str | None = "mean",
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Place scalar values into a (·, ·, 2, 2) grid at the slot
        indicated by each row's component label.

        This is a **layout convenience**, not tensor algebra.

        Returns
        -------
        tensor : ndarray
            ``(n_freq, 2, 2)`` for a single *station* or
            ``(n_station, n_freq, 2, 2)`` otherwise.
        freqs : ndarray
            The frequency grid used.
        stations : ndarray
            Station axis (empty for single-station requests).
        """
        # local import to avoid a hard dependency / import cycles
        from .tensor import (
            TensorBase,  # pylint: disable=import-outside-toplevel
        )

        keep = ["station", "freq", "comp", self.VAR_NAME]
        df = self._frame.loc[:, [c for c in keep if c in self._frame]]

        tb = TensorBase.from_avg((df, {}))
        return tb.to_tensor(
            var=self.VAR_NAME,
            station=station,
            agg=agg,
            fill_value=fill_value,
            sort_freq=sort_freq,
            align=align,
        )

    def convert_unit(self, target: str = "mrad") -> None:
        r"""
        Convert the phase-stdev column **in place** between
        :math:`\mathrm{mrad}` and :math:`\mathrm{deg}` and update
        ``Unit.Phase`` in :pyattr:`meta`.

        The conversion factors come from:

        .. math::
           1~\mathrm{rad} = 1000~\mathrm{mrad},\qquad
           1~\mathrm{deg} = \frac{\pi}{180}~\mathrm{rad}

        Hence:

        .. math::
           1~\mathrm{mrad} =
           \frac{180}{\pi \cdot 1000}~\mathrm{deg}
           \approx 0.05729578~\mathrm{deg}

        Parameters
        ----------
        target : {'mrad', 'deg'}, default 'mrad'
            Desired output unit for the stored values.

        Notes
        -----
        * If the current unit already matches *target*, this is a
          no-op.
        * Non-numeric cells are coerced to ``NaN`` and remain so.
        * If the component's data column is absent, the method
          returns quietly (nothing to convert).
        """
        # Resolve current ↔ target units (case-insensitive).
        cur = str(self._meta.get("Unit.Phase", "mrad")).lower()
        tgt = str(target).lower()
        if cur == tgt:
            return
        if tgt not in {"mrad", "deg"}:
            raise ValueError("target must be 'mrad' or 'deg'")

        # Identify the canonical data column.
        col = self.VAR_NAME
        if col not in self._frame.columns:
            # Nothing to convert; keep meta unchanged.
            return

        # Coerce to float; preserve NaNs for missing / bad cells.
        x = pd.to_numeric(self._frame[col], errors="coerce")

        # Compute the scale once, using exact constants.
        if cur == "mrad" and tgt == "deg":
            # 1 mrad = 180 / (π * 1000) deg
            factor = 180.0 / (np.pi * 1000.0)
        elif cur == "deg" and tgt == "mrad":
            # 1 deg = (π / 180) * 1000 mrad
            factor = (np.pi / 180.0) * 1000.0
        else:
            # Future-proof against unexpected unit labels.
            raise ValueError(f"unsupported conversion {cur} → {tgt}")

        # Apply conversion and update unit metadata.
        self._frame[col] = x * factor
        self._meta["Unit.Phase"] = tgt

    def __str__(self) -> str:  # noqa: D401
        """Human-readable summary string."""
        r, c = self.shape
        return f"{self.__class__.__name__}[{r}×{c}] var={self.VAR_NAME!s}"


class SPhz(PhaseStdBase):
    """
    Standard deviation of *impedance phase* (``Z``).

    Recognised source columns
    -------------------------
    * ``sPhz``
    * ``Z.perr``

    Internal canonical column: ``'sphz'``.
    """

    VAR_NAME = "s_phz"
    # KEY_CANDIDATES = ("sPhz", "SPhz", "Z.perr", "z.perr", "sphz")
    LABEL = "Z phase σ (sPhz)"


class SEphz(PhaseStdBase):
    """
    Standard deviation of *E-phase*.

    Recognised source columns
    -------------------------
    * ``sEphz``
    * ``E.perr``

    Internal canonical column: ``'sephz'``.
    """

    VAR_NAME = "s_ephz"
    # KEY_CANDIDATES = ("sEphz", "SEphz", "E.perr", "e.perr", "sephz")
    LABEL = "E phase σ (sEphz)"


class SHphz(PhaseStdBase):
    """
    Standard deviation of *H-phase*.

    Recognised source columns
    -------------------------
    * ``sHphz``
    * ``H.perr``

    Internal canonical column: ``'shphz'``.
    """

    VAR_NAME = "s_hphz"
    # KEY_CANDIDATES = ("sHphz", "SHphz", "H.perr", "h.perr", "shphz")
    LABEL = "H phase σ (sHphz)"


def _find_col(
    df: pd.DataFrame,
    candidates: Sequence[str],
) -> str | None:
    """
    Return the first column name present in *df* among *candidates*.

    Matching is case-insensitive and ignores surrounding spaces.
    """
    low = {str(c).strip().lower(): c for c in df.columns}
    for want in candidates:
        key = str(want).strip().lower()
        if key in low:
            return low[key]
    return None


def _to_num(x: Any) -> float | np.floating | np.nan:
    """
    Robust numeric coercion:

    * empty/asterisk/'nan' → NaN
    * integral floats become ints when safe (not required here)
    """
    if x is None:
        return np.nan
    s = str(x).strip()
    if s in {"", "*", "nan", "NaN", "None", "null"}:
        return np.nan
    try:
        return float(s)
    except Exception:
        return np.nan


def _norm_comp(x: Any) -> str:
    """
    Normalise component label to a canonical upper-case token.
    """
    if x is None:
        return "ExHy"
    return str(x).strip() or "ExHy"


PhaseSigma = SPhz  # dedicated aliases per field
