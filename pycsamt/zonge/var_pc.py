# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
Percent-variation (QC) components for Zonge AVG tables.

This module provides small, format-agnostic containers for the
three "percent-error" columns commonly found in AMTAVG/CSAVGW-
style data:

* PcEmag – relative |E| error in percent (``%Emag`` / ``E.%err``)
* PcHmag – relative |H| error in percent (``%Hmag`` / ``B.%err``)
* PcRho  – relative apparent-resistivity error in percent
           (``%Rho`` / ``ARes.%err`` / ``rho.%err``)

All classes inherit from :class:`AVGComponentBase`, read tidy
tables (modern or legacy), normalize column names, and can export
to :class:`xarray.Dataset` for multi-dimensional workflows.
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
from .utils import _to_numeric_percent, find_and_rename_column
from .utils import to_xarray as _to_xr

logger = get_logger(__name__)

__all__ = [
    "PcEmag",
    "PcRho",
    "PcHmag",
    "PcHmag",
    "EmagPctErr",
    "HmagPctErr",
    "RhoPctErr",
]


class PercentVarBase(AVGComponentBase):
    """
    Abstract base for percent-variation QC columns.

    Subclasses specify:

    * ``VAR_NAME`` – canonical output name (e.g., ``"pc_emag"``)
    * ``ALIASES``  – ordered list of legacy/modern column aliases

    They inherit a consistent ``read()/write()/to_xarray()``.
    """

    # Class-level constants (not dataclass fields)

    # # canonical variable name written in the frame
    VAR_NAME: ClassVar[str] = ""  # e.g., "pc_emag"

    # # # ordered list of candidate column labels in source tables
    # ALIASES: ClassVar[Tuple[str, ...]] = ()   # e.g., ("%Emag", "E.%err", ...)

    # # default banner/title used by ``write()``
    TITLE: ClassVar[str] = "Percent Variation"

    # # default dataset attribute for units
    UNIT_ATTR: ClassVar[str] = "Unit.Percent"

    def __init__(
        self,
        data: pd.DataFrame | None = None,
        meta: Mapping[str, Any] | None = None,
        *,
        name: str | None = None,
        verbose: bool = False,
    ) -> None:
        super().__init__(data=data, meta=meta, name=name, verbose=verbose)

    def read(
        self,
        source: pd.DataFrame,
        meta: Mapping[str, Any] | None = None,
        **kws: Any,
    ) -> None:
        """
        Load the percent-variation column from a tidy table.

        The method:

        1. locates a column among ``ALIASES``,
        2. ensures ``station`` / ``freq`` / ``comp`` exist,
        3. normalizes the percent column into ``VAR_NAME``,
        4. stores a compact frame with those four columns, and
        5. merges header attributes into ``self._meta``.

        Parameters
        ----------
        source
            Tidy :class:`pandas.DataFrame`. Must contain at least
            one alias for the percent column and, preferably, the
            ``station`` and ``freq`` coordinates. If ``comp`` is
            absent, it is injected as ``"ExHy"``.
        meta
            Free-form header/keyword mapping to stash as attrs.
        """
        # Guard to catch misconfigured subclasses early
        if not self.VAR_NAME:
            raise RuntimeError(
                f"{self.__class__.__name__}: VAR_NAME/ALIASES"
                " must be set on the subclass."
            )

        if not isinstance(source, pd.DataFrame):
            raise TypeError("PercentVarBase.read expects a DataFrame.")

        df = source.copy()
        self._meta = dict(meta or {})
        self._meta.setdefault(self.UNIT_ATTR, "%")

        # var_col = _first_present(df, self.ALIASES)
        # if var_col is None:
        #     raise AvgDataError(
        #         f"{self.__class__.__name__}: none of aliases "
        #         f"{self.ALIASES!r} present in table."
        #     )
        # After standardization, we expect the canonical VAR_NAME.
        # If not present, create it with NaNs for consistency.

        # Use the new helper to standardize the column
        df = find_and_rename_column(df, self.VAR_NAME)

        if self.VAR_NAME not in df.columns:
            df[self.VAR_NAME] = np.nan
            if self.verbose:
                logger.debug(f"'{self.VAR_NAME}' not in source. Creating empty.")

        # ensure coords exist (inject conservative defaults)
        if "comp" not in df.columns:
            df["comp"] = "ExHy"
        if "station" not in df.columns:
            df["station"] = np.nan
        if "freq" not in df.columns:
            raise AvgDataError("frequency column 'freq' is required.")

        # normalize percent column → float
        df[self.VAR_NAME] = _to_numeric_percent(df[self.VAR_NAME].copy())

        keep = [
            c for c in ("station", "freq", "comp", self.VAR_NAME) if c in df.columns
        ]
        # store a compact, predictable layout
        # self._frame = df.loc[:, ["station", "freq", "comp", self.VAR_NAME]]
        self._frame = df.loc[:, keep].copy()
        self._meta = dict(meta or {})

        # ensure a stable units hint at dataset level
        self._meta.setdefault(self.UNIT_ATTR, "%")

        return self

    def write(
        self,
        *,
        float_fmt: str = "%.6g",
        na_rep: str = "",
    ) -> Sequence[str]:
        """
        Serialise the component to a human-friendly CSV block.

        A small banner is included, followed by a UTC timestamp
        and the minimal 4-column table:

        ``station, freq, comp, <VAR_NAME>``

        Parameters
        ----------
        float_fmt
            Format for floating-point export.
        na_rep
            Representation for missing values.

        Returns
        -------
        list[str]
            Lines ready to prepend/append when writing AVG files.
        """
        if self._frame.empty:
            return [f"\\ ${self.TITLE}", "$Written="]  # minimal stub

        return self._write_csv_block(
            cols=["station", "freq", "comp", self.VAR_NAME],
            title=f"$ {self.TITLE}",
            float_fmt=float_fmt,
            na_rep=na_rep,
            include_meta=True,
            stamp=True,
        )

    def to_xarray(
        self,
        *,
        coords: Sequence[str] = ("station", "freq", "comp"),
        attrs: dict[str, Any] | None = None,
    ):
        """
        Convert to an :class:`xarray.Dataset`.

        The resulting dataset has dimensions given by *coords*
        (subset of the columns present), and a single data
        variable with the canonical name ``VAR_NAME``.

        Parameters
        ----------
        coords
            Coordinate columns to grid against (default order
            ``station → freq → comp``).
        attrs
            Extra attributes to merge into the dataset-level
            metadata.  ``Unit.Percent='%'`` is ensured.

        Returns
        -------
        xarray.Dataset
            Dataset with one numeric variable named ``VAR_NAME``.
        """
        if self._frame.empty:
            raise AvgDataError("empty percent-variation frame.")

        df = self._frame.copy()
        merged = dict(self._meta)
        merged.setdefault(self.UNIT_ATTR, "%")
        if attrs:
            merged.update(attrs)

        return _to_xr(
            df,
            coords=coords,
            data_vars=[self.VAR_NAME],
            attrs=merged,
        )

    def to_tensor_like(
        self,
        *,
        align: str = "union",
        station: float | int | None = None,
        fill_value: float = np.nan,
        sort_freq: bool = True,
        agg: str | None = "mean",
    ):
        """
        Use TensorBase's gridding logic to place this scalar variable
        in a 2×2 container. Entries will be filled only where a row's
        `comp` exists; other slots are NaN. This is for alignment/
        reshaping convenience, not tensor algebra.
        """
        from .tensor import (
            TensorBase,  # local import to avoid hard coupling
        )

        df = self._frame.copy()
        # hand only the 4 required columns to TensorBase
        keep = ["station", "freq", "comp", self.VAR_NAME]
        df = df.loc[:, [c for c in keep if c in df.columns]]

        tb = TensorBase.from_avg((df, {}))
        return tb.to_tensor(
            var=self.VAR_NAME,
            station=station,
            agg=agg,
            fill_value=fill_value,
            sort_freq=sort_freq,
            align=align,
        )

    # friendly diagnostics
    def __str__(self) -> str:
        r, c = self.shape
        return f"{self.__class__.__name__}[{r}×{c}] var={self.VAR_NAME}"


class PcEmag(PercentVarBase):
    r"""
    Percent error on electric-field magnitude, :math:`|E|`.

    This component normalizes legacy aliases into a canonical
    variable named ``pc_emag``:

    * ``'%Emag'`` (legacy AMTAVG/MTEdit),
    * ``'E.%err'`` (CSAVGW / modern),
    * case-sensitive exact matches are used.

    Notes
    -----
    Values are stored in *percent*; the dataset attribute
    ``Unit.Percent`` is set to ``'%'`` when exporting to
    :class:`xarray.Dataset`.
    """

    VAR_NAME = "pc_emag"
    # ALIASES = get_aliases(VAR_NAME, kind ='qc') # ("%Emag", "E.%err")
    TITLE = "Percent |E| Variation"
    UNIT_ATTR = "Unit.Percent"


class PcHmag(PercentVarBase):
    r"""
    Percent error on magnetic-field magnitude, :math:`|H|`.

    This component normalizes legacy aliases into a canonical
    variable named ``pc_hmag``:

    * ``'%Hmag'`` (legacy AMTAVG/MTEdit),
    * ``'B.%err'`` (CSAVGW, where ``B`` denotes the H-field),
    * ``'H.%err'`` (rare but seen).

    Notes
    -----
    Values are stored in *percent*; the dataset attribute
    ``Unit.Percent`` is set to ``'%'`` when exporting to
    :class:`xarray.Dataset`.
    """

    VAR_NAME = "pc_hmag"
    # 'H' vs 'B' modern label differences covered; legacy %Hmag too
    # ALIASES = get_aliases(VAR_NAME, kind ='qc')#("%Hmag", "B.%err", "H.%err")
    TITLE = "Percent |H| Variation"
    UNIT_ATTR = "Unit.Percent"


class PcRho(PercentVarBase):
    r"""
    Percent error on apparent resistivity, :math:`\rho_a`.

    This component normalizes legacy aliases into a canonical
    variable named ``pc_rho``:

    * ``'%Rho'`` (legacy AMTAVG/MTEdit),
    * ``'ARes.%err'`` (CSAVGW),
    * ``'rho.%err'`` (modern lower-case variant).

    Notes
    -----
    Values are stored in *percent*; the dataset attribute
    ``Unit.Percent`` is set to ``'%'`` when exporting to
    :class:`xarray.Dataset`.
    """

    VAR_NAME = "pc_rho"
    # ALIASES = get_aliases(VAR_NAME, kind ='qc')# ("%Rho", "ARes.%err", "rho.%err")
    TITLE = "Percent ρa Variation"
    UNIT_ATTR = "Unit.Percent"


EmagPctErr = PcEmag
HmagPctErr = PcHmag
RhoPctErr = PcRho
