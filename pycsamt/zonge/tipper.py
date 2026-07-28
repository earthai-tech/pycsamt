# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
Tipper Transfer Function Component.

This module provides the `Tipper` class for managing the
magnetotelluric Tipper transfer function, which relates the
vertical magnetic field (Hz) to the horizontal components
(Hx, Hy).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import (
    Any,
)

import numpy as np
import pandas as pd

from ..exceptions import AvgDataError
from .base import AVGComponentBase
from .utils import _standardise_columns, _to_complex
from .utils import to_xarray as _to_xr

__all__ = ["Tipper"]


class Tipper(AVGComponentBase):
    r"""Geomagnetic transfer function (Tipper) component.

    The Tipper, also known as the induction vector, is a transfer
    function that relates the vertical component of the magnetic
    field (:math:`H_z`) to the horizontal components
    (:math:`H_x`, :math:`H_y`).

    .. math::
        H_z = T_x H_x + T_y H_y

    This class is designed to hold the complex Tipper components
    :math:`T_x` and :math:`T_y` after they have been calculated.

    Attributes
    ----------
    tx, ty : pd.Series
        The complex Tipper components for the x and y directions.

    See Also
    --------
    pycsamt.zonge.avg.AMTAVG.calculate_tipper : The method used to
        compute the Tipper values.
    """

    def __init__(
        self,
        data: pd.DataFrame | None = None,
        meta: Mapping[str, Any] | None = None,
        *,
        name: str | None = None,
        verbose: bool = False,
    ) -> None:
        """Initializes the Tipper component."""
        super().__init__(data=data, meta=meta, name=name or "Tipper", verbose=verbose)

    def read(
        self,
        source: pd.DataFrame,
        meta: Mapping[str, Any] | None = None,
        **kws: Any,
    ) -> None:
        """
        Populate the Tipper component from a DataFrame.

        This method expects a DataFrame that contains the
        calculated Tipper components, named 'tx' and 'ty'. If these
        columns are missing, they will be created and filled with
        NaNs for structural consistency.
        """
        if not isinstance(source, pd.DataFrame):
            raise TypeError("Tipper.read expects a pandas.DataFrame.")

        df = _standardise_columns(source.copy())
        self._meta = dict(meta or {})

        # Ensure required columns exist, creating them if necessary
        for col in ["tx", "ty"]:
            if col not in df.columns:
                df[col] = np.nan
                if self.verbose:
                    self._logger.debug(f"'{col}' not in source. Creating empty column.")

        # Ensure coordinates exist for consistency
        if "station" not in df.columns:
            df["station"] = np.nan
        if "freq" not in df.columns:
            raise AvgDataError("Tipper requires a 'freq' column.")

        # Normalize types
        from .utils import to_numeric_if_possible

        df["station"] = to_numeric_if_possible(df["station"])
        df["freq"] = pd.to_numeric(df["freq"], errors="coerce")
        df["tx"] = df["tx"].map(_to_complex)
        df["ty"] = df["ty"].map(_to_complex)

        keep_cols = ["station", "freq", "tx", "ty"]
        self._frame = df.loc[:, [c for c in keep_cols if c in df.columns]]

        return self

    def write(self) -> list[str]:
        """
        Serialise to a compact CSV block with a meta preamble.
        """
        if self._frame.empty:
            return []

        return self._write_csv_block(
            cols=["station", "freq", "tx", "ty"],
            title="$Tipper Block",
            include_meta=True,
            stamp=True,
        )

    @property
    def tx(self) -> pd.Series:
        """The complex Tipper component Tx."""
        return self._frame.get("tx", pd.Series(dtype="complex128"))

    @property
    def ty(self) -> pd.Series:
        """The complex Tipper component Ty."""
        return self._frame.get("ty", pd.Series(dtype="complex128"))

    def to_xarray(
        self,
        *,
        coords: Sequence[str] = ("station", "freq"),
        attrs: dict[str, Any] | None = None,
    ):
        """
        Convert the Tipper data into an xarray.Dataset.
        """
        if self._frame.empty:
            raise AvgDataError("Empty frame; nothing to export.")

        merged = dict(self._meta)
        if attrs:
            merged.update(attrs)

        return _to_xr(
            self._frame.copy(),
            coords=coords,
            data_vars=["tx", "ty"],
            attrs=merged,
        )

    def __str__(self) -> str:
        """Provide a concise, human-readable representation."""
        if self._frame.empty:
            return "Tipper(status=empty)"

        n_st = (
            self._frame["station"].nunique() if "station" in self._frame.columns else 0
        )
        n_frq = self._frame["freq"].nunique() if "freq" in self._frame.columns else 0
        return f"Tipper(rows={len(self._frame)}, stations={n_st}, freqs={n_frq})"

    def __repr__(self) -> str:
        """Provide an unambiguous developer representation."""
        return self.__str__()
