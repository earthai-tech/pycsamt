# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
Impedance tensor (Z) component for Zonge AVG data.

This module provides the `Z` class, which computes the complex
impedance tensor from apparent resistivity and phase data. It
inherits from TensorBase to provide powerful reshaping and
analysis capabilities.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import (
    Any,
)

import numpy as np
import pandas as pd

from ..constants import MU_0, PI
from ..exceptions import AvgDataError
from ..log.logger import get_logger
from ..utils._dependency import import_optional_dependency
from .tensor import TensorBase, _norm_comp
from .utils import _standardise_columns, _to_num

logger = get_logger(__name__)

__all__ = ["Z"]


class Z(TensorBase):
    r"""Complex impedance tensor (Z) component.

    This class reads a tidy AVG table and computes the complex
    impedance Z from apparent resistivity (:math:`\rho_a`) and
    impedance phase (:math:`\phi`). It provides properties to
    access the complex tensor, its real and imaginary parts, and
    the propagated error.

    The impedance is calculated using the standard formula [1]_:

    .. math::
        Z = \sqrt{\rho_a \cdot \omega \cdot \mu_0} \cdot
        e^{i \cdot \phi}

    where :math:`\omega = 2\pi f`.

    Attributes
    ----------
    z : pd.Series
        The complex impedance :math:`Z` in Ohms [Ω].
    z_real, z_imag : pd.Series
        The real and imaginary parts of the impedance tensor.
    z_err : pd.Series
        The propagated error in the magnitude of :math:`Z`,
        :math:`|dZ|`, in Ohms [Ω].
    z_xx, z_xy, z_yx, z_yy : pd.Series
        Convenience properties to access the individual components
        of the complex impedance tensor.

    Examples
    --------
    >>> from pycsamt.zonge import Z
    >>> from pycsamt.zonge.avg import AVG
    >>> avg = AVG.from_file("data/avg/K2.avg")
    >>> z_component = avg.z
    >>> # Get the complex impedance for all measurements
    >>> complex_z_values = z_component.z
    >>> # Get only the Z_xy component
    >>> z_xy = z_component.z_xy

    References
    ----------
    .. [1] Vozoff, K. (1972). The magnetotelluric method in the
           exploration of sedimentary basins. Geophysics, 37(1),
           98-141.

    See Also
    --------
    TensorBase : The parent class providing tensor-shaping logic.
    """

    def read(
        self,
        source: pd.DataFrame,
        meta: Mapping[str, Any] | None = None,
        **kws: Any,
    ) -> None:
        """
        Read and prepare data for impedance calculation.
        """
        if not isinstance(source, pd.DataFrame):
            raise TypeError("Z.read expects a pandas.DataFrame.")

        # df = source.copy()
        df = _standardise_columns(source.copy())
        self._meta = dict(meta or {})

        # # --- Required columns ---

        # missing = [k for k, v in required.items() if v is None]
        required = ["rho", "phase", "freq"]
        missing = [c for c in required if c not in df.columns]
        if missing:
            raise AvgDataError(f"Z: missing required columns: {missing}")
        # --- Optional error columns ---
        optional_qc = ["pc_rho", "s_phz"]
        for col in optional_qc:
            if col not in df.columns:
                df[col] = np.nan
                if self.verbose:
                    self._logger.debug(f"'{col}' not in source. Creating empty.")

        # Ensure coords exist
        if "station" not in df.columns:
            df["station"] = np.nan
        if "comp" not in df.columns:
            df["comp"] = "ExHy"

        # Normalize component labels to uppercase canonical form
        df["comp"] = df["comp"].map(_norm_comp)
        df.dropna(subset=["comp"], inplace=True)

        # Normalize types
        for col in ["rho", "phase", "freq", "pc_rho", "s_phz"]:
            if col in df.columns:
                df[col] = df[col].map(_to_num)

        keep_cols = [
            "station",
            "freq",
            "comp",
            "rho",
            "phase",
            "pc_rho",
            "s_phz",
        ]
        self._frame = df.loc[:, [c for c in keep_cols if c in df.columns]]

        return self

    def _get_component_series(
        self, comp_names: tuple[str, ...], series: pd.Series
    ) -> pd.Series:
        """Helper to filter a property series by component."""
        if self._frame.empty or "comp" not in self._frame.columns:
            return pd.Series(dtype=series.dtype)

        mask = self._frame["comp"].isin(comp_names)
        return series[mask]

    @property
    def z(self) -> pd.Series:
        """
        Complex impedance Z [Ω].
        """
        if self._frame.empty:
            return pd.Series(dtype="complex128")

        rho = self._frame["rho"]
        phase_mrad = self._frame["phase"]
        freq = self._frame["freq"]

        # Convert phase from milliradians to radians
        phase_rad = phase_mrad * 1e-3
        omega = 2 * PI * freq

        # Calculate magnitude of Z
        z_mag = np.sqrt(rho * omega * MU_0)

        # Calculate complex impedance
        return z_mag * np.exp(1j * phase_rad)

    @property
    def z_real(self) -> pd.Series:
        """
        Real part of the impedance tensor, Z' [Ω].
        """
        return self.z.apply(np.real)

    @property
    def z_imag(self) -> pd.Series:
        """
        Imaginary part of the impedance tensor, Z'' [Ω].
        """
        return self.z.apply(np.imag)

    @property
    def z_err(self) -> pd.Series:
        r"""
        Propagated error in the magnitude of Z, |dZ| [Ω].

        Calculated via standard error propagation from the
        relative error in resistivity (dρ/ρ) and the absolute
        error in phase (dφ).

        .. math::
            |dZ| \approx \sqrt{
            (\frac{\partial |Z|}{\partial \rho} d\rho)^2 +
            (|Z| d\phi)^2
            }

        Since phase errors are often dominant and uncorrelated,
        a simpler estimate is often used:

        .. math::
            |dZ| \approx \frac{1}{2} |Z| \frac{d\rho}{\rho}
        """
        if self._frame.empty or "rho" not in self._frame.columns:
            return pd.Series(dtype="float64")

        has_rho_err = "pc_rho" in self._frame.columns
        has_phi_err = "s_phz" in self._frame.columns

        if not has_rho_err and not has_phi_err:
            return pd.Series(dtype="float64", index=self._frame.index)

        z_mag = np.sqrt(self._frame["rho"] * (2 * PI * self._frame["freq"]) * MU_0)

        term_rho_sq = 0.0
        if has_rho_err:
            # Relative error drho/rho
            rel_err_rho = self._frame["pc_rho"] / 100.0
            term_rho_sq = (0.5 * rel_err_rho) ** 2

        term_phi_sq = 0.0
        if has_phi_err:
            dphi_rad = self._frame["s_phz"] * 1e-3
            term_phi_sq = dphi_rad**2

        # Propagated error in magnitude |Z|
        # d|Z| = 0.5 * |Z| * (drho/rho)

        return z_mag * np.sqrt(term_rho_sq + term_phi_sq)

    @property
    def z_xx(self) -> pd.Series:
        """Complex impedance for the Zxx component."""
        return self._get_component_series(("ZXX", "ExHX"), self.z)

    @property
    def z_xy(self) -> pd.Series:
        """Complex impedance for the Zxy component."""
        return self._get_component_series(("ZXY", "EXHY"), self.z)

    @property
    def z_yx(self) -> pd.Series:
        """Complex impedance for the Zyx component."""
        return self._get_component_series(("ZYX", "EYHX"), self.z)

    @property
    def z_yy(self) -> pd.Series:
        """Complex impedance for the Zyy component."""
        return self._get_component_series(("ZYY", "EYHY"), self.z)

    @property
    def z_xx_err(self) -> pd.Series:
        """Propagated error for the Zxx component."""
        return self._get_component_series(("ZXX", "EXHX"), self.z_err)

    @property
    def z_xy_err(self) -> pd.Series:
        """Propagated error for the Zxy component."""
        return self._get_component_series(("ZXY", "EXHY"), self.z_err)

    @property
    def z_yx_err(self) -> pd.Series:
        """Propagated error for the Zyx component."""
        return self._get_component_series(("ZYX", "EYHX"), self.z_err)

    @property
    def z_yy_err(self) -> pd.Series:
        """Propagated error for the Zyy component."""
        return self._get_component_series(("ZYY", "EYHY"), self.z_err)

    def to_tensor(
        self,
        *,
        var: str = "z",
        station: int | float | None = None,
        agg: str | None = "mean",
        fill_value: float = np.nan,
        sort_freq: bool = True,
        align: str = "union",
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Convert impedance data into a 2x2 tensor.
        """
        temp_frame = self._frame.copy()

        if var == "z":
            temp_frame["__real"] = self.z_real
            temp_frame["__imag"] = self.z_imag

            # Create a temporary instance for base method call
            tb = TensorBase()
            tb._frame = temp_frame

            T_real, freqs, stations = tb.to_tensor(
                var="__real",
                station=station,
                agg=agg,
                fill_value=fill_value,
                sort_freq=sort_freq,
                align=align,
            )
            T_imag, _, _ = tb.to_tensor(
                var="__imag",
                station=station,
                agg=agg,
                fill_value=fill_value,
                sort_freq=sort_freq,
                align=align,
            )
            return T_real + 1j * T_imag, freqs, stations

        elif var in ("z_real", "z_imag", "z_err"):
            temp_frame[var] = getattr(self, var)
            tb = TensorBase()
            tb._frame = temp_frame
            return tb.to_tensor(
                var=var,
                station=station,
                agg=agg,
                fill_value=fill_value,
                sort_freq=sort_freq,
                align=align,
            )
        else:
            return super().to_tensor(
                var=var,
                station=station,
                agg=agg,
                fill_value=fill_value,
                sort_freq=sort_freq,
                align=align,
            )

    def to_xarray(
        self,
        *,
        var: str = "z",
        station: int | float | None = None,
        agg: str | None = "mean",
        fill_value: float = np.nan,
        attrs: Mapping[str, Any] | None = None,
    ):
        """
        Return a 3D or 4D xarray.DataArray.
        """
        import_optional_dependency(
            "xarray",
            extra="xarray is required for to_xarray()",
            errors="raise",
        )

        import xarray as xr

        # Use the new to_tensor method to get the data
        T, freqs, stations = self.to_tensor(
            var=var, station=station, agg=agg, fill_value=fill_value
        )

        e_axis = np.array(["Ex", "Ey"])
        h_axis = np.array(["Hx", "Hy"])

        merged_attrs = dict(self._meta)
        if attrs:
            merged_attrs.update(attrs)

        if stations.size == 0:
            da = xr.DataArray(
                T,
                dims=("freq", "e", "h"),
                coords={"freq": freqs, "e": e_axis, "h": h_axis},
                attrs=merged_attrs,
                name=var,
            )
        else:
            da = xr.DataArray(
                T,
                dims=("station", "freq", "e", "h"),
                coords={
                    "station": stations,
                    "freq": freqs,
                    "e": e_axis,
                    "h": h_axis,
                },
                attrs=merged_attrs,
                name=var,
            )
        return da

    def write(self) -> Sequence[str]:
        """
        Serializes the core Z data to a CSV block.
        """
        if self._frame.empty:
            return ["\\ $Z (Impedance) Block", ""]

        # Create a temporary frame for writing
        df_write = self._frame[["station", "freq", "comp"]].copy()

        df_write["z_real"] = self.z_real
        df_write["z_imag"] = self.z_imag
        df_write["z_err"] = self.z_err

        return self._write_csv_block(
            cols=list(df_write.columns),
            title="$Z (Impedance) Block",
            include_meta=True,
            stamp=True,
        )

    def __str__(self) -> str:
        """Provide a concise, human-readable representation."""

        if self._frame.empty:
            return "Z(status=empty)"

        # Safely get unique counts for the summary
        n_st = (
            self._frame["station"].nunique() if "station" in self._frame.columns else 0
        )
        n_frq = self._frame["freq"].nunique() if "freq" in self._frame.columns else 0
        n_comp = self._frame["comp"].nunique() if "comp" in self._frame.columns else 0

        return (
            f"Z(rows={len(self._frame)}, "
            f"stations={n_st}, "
            f"freqs={n_frq}, "
            f"components={n_comp})"
        )

    def __repr__(self) -> str:
        """Provide an unambiguous developer representation."""
        # For this class, a detailed __str__ is also a good __repr__
        return self.__str__()
