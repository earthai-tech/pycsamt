# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
pycsamt.zonge.processing
------------------------

This module provides the ASTATIC class, which contains methods
for advanced data conditioning and analysis, such as static shift
correction, filtering, and data merging, mirroring the
capabilities of Zonge's ASTATIC program.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
)

import numpy as np
import pandas as pd

from ..constants import MU_0, PI
from ..exceptions import NotReadError, ProcessingError
from ..utils.validation import has_read
from .config import Zonge
from .proc_utils import ama, flma, tma

if TYPE_CHECKING:
    from .avg import AMTAVG, BaseAVG
    from .base import AVGFrame

__all__ = ["ASTATIC"]


class ASTATIC(Zonge):
    r"""A class for advanced processing of Zonge AVG data.

    This class provides a suite of methods for data conditioning
    and analysis that mirror the functionality of Zonge's ASTATIC
    software [1]_. It operates on a loaded
    :class:`~.avg.AMTAVG` object, allowing for complex
    operations like static shift correction, data filtering, and
    interpolation.

    Parameters
    ----------
    avg_data : :class:`~.avg.BaseAVG`, optional
        An initialized and loaded AVG object. If provided, the
        processor is immediately ready for use.
    verbose : bool, default False
        If ``True``, log messages will be printed to the console
        during processing operations, providing insight into the
        process.

    Attributes
    ----------
    avg : :class:`~.avg.AMTAVG` or None
        The AVG data object that the processor is operating on.
        It is populated by either passing it to the constructor or
        by using the :meth:`read` method.

    Methods
    -------
    read(source, meta=None)
        Loads a data source into the processor.
    correct_static_shift(...)
        Corrects for static shift using various spatial filters.
    correct_capacitive_coupling(...)
        Remediates high-frequency distortions from capacitive
        coupling.

    Notes
    -----
    The `ASTATIC` class is designed using a "composition over
    inheritance" approach. It is not a type of `AVG` file;
    rather, it is a tool that *contains* and *operates on* an
    `AVG` object.

    Most processing methods, such as `correct_static_shift`,
    modify the underlying `avg` object's DataFrame *in place* when
    the `update_components` parameter is set to ``True``. This
    ensures that all data components are consistently updated
    after a processing step.

    Examples
    --------
    The typical workflow involves loading data with an `AMTAVG`
    object and then passing that object to the `ASTATIC`
    processor.

    >>> from pycsamt.zonge import AMTAVG, ASTATIC
    >>> # 1. Load the data
    >>> avg = AMTAVG.from_file("data/avg/K2.avg")
    >>>
    >>> # 2. Create a processor and apply a correction
    >>> processor = ASTATIC().read(avg)
    >>> processor.correct_static_shift(
    ...     reference_freq=4096, filter_method="tma"
    ... )
    >>>
    >>> # 3. The original avg object is now updated and can be saved
    >>> avg.to_modern("K2_corrected.avg")

    References
    ----------
    .. [1] Zonge International, Inc. (2014). *ASTATIC v3.70
           User Manual*.

    See Also
    --------
    pycsamt.zonge.avg.AMTAVG : The primary data container class.
    pycsamt.zonge.proc_utils : The module containing the core
        filtering algorithms.
    """

    def __init__(
        self,
        avg_data: BaseAVG | AMTAVG | None = None,
        verbose: bool = False,
        **kws,
    ):
        super().__init__(verbose=verbose)
        self.avg: BaseAVG | AMTAVG | None = avg_data

        if self.avg is not None:
            self.read(avg_data)

    def read(
        self,
        source: str | Path | AVGFrame | pd.DataFrame | BaseAVG,
        meta: Mapping[str, Any] | None = None,
    ) -> ASTATIC:
        r"""Load a data source into the processor.

        This is the primary method for associating a dataset with
        the ASTATIC processor. It is designed to be flexible,
        accepting various input types and ensuring that the
        processor is initialized with a valid, fully-loaded AVG
        data object.

        Parameters
        ----------
        source : str, Path, AVGFrame, pd.DataFrame, or BaseAVG
            The data source to load. This can be:

            - A string or `pathlib.Path` pointing to a Zonge AVG
              file. A new :class:`~.avg.AMTAVG` instance will be
              created internally to read the file.
            - A `pandas.DataFrame`. A new `AMTAVG` instance will be
              created to read the DataFrame.
            - A pre-existing, loaded object that inherits from
              :class:`~.avg.BaseAVG` (e.g., an `AVG` or `AMTAVG`
              instance).

        meta : mapping, optional
            An optional dictionary of metadata. This is only used
            when `source` is a `pd.DataFrame`.

        Returns
        -------
        self : ASTATIC
            The method returns the instance of the class, allowing
            for convenient method chaining.

        Raises
        ------
        TypeError
            If the `source` is of an unsupported type.
        NotReadError
            If an `AVG`-like object is passed as the `source` but
            it has not been populated with data yet.

        Notes
        -----
        This method acts as a smart constructor for the processor.
        If the provided `source` is not already a loaded `AVG`
        object, the method takes on the responsibility of creating
        one. This ensures that the `self.avg` attribute is always
        a valid, ready-to-use data container.

        Examples
        --------
        >>> from pycsamt.zonge import AMTAVG, ASTATIC
        >>> # --- Reading from a pre-loaded AVG object ---
        >>> avg = AMTAVG.from_file("data/avg/K2.avg")
        >>> processor = ASTATIC().read(avg)
        >>>
        >>> # --- Reading directly from a file path ---
        >>> processor_from_file = ASTATIC().read("data/avg/K1.avg")

        See Also
        --------
        pycsamt.zonge.avg.AVG.from_file : The primary factory for
            creating an AVG data object.
        pycsamt.utils.validation.has_read : The underlying utility
            used to validate loaded data objects.
        """
        from .avg import AMTAVG

        if isinstance(source, (str, Path, pd.DataFrame)):
            # If source is a file or raw DataFrame, create a new
            # AMTAVG instance to read and process it.
            self.avg = AMTAVG(verbose=self.verbose).read(source, meta=meta)
        elif hasattr(source, "__has_read__"):
            # If source is an AVG-like object, check if it's loaded
            has_read(source)
            self.avg = source
        else:
            raise TypeError(f"Unsupported source type: {type(source).__name__}")

        if self.verbose:
            src_name = (
                self.avg._source_path.name if self.avg._source_path else "DataFrame"
            )
            self._logger.info(
                f"ASTATIC processor initialized with data from '{src_name}'."
            )
        return self

    def correct_capacitive_coupling(
        self,
        contact_resistance: float | pd.Series | str,
        setup_length: float | pd.Series | str,
        wire_capacitance: float | pd.Series | str = 15.0,
        update_components: bool = True,
    ) -> pd.DataFrame:
        r"""Correct for capacitive coupling effects in E-field data.

        This method remediates distortions in high-frequency data
        caused by capacitive coupling between receiver wires and
        the ground.

        Parameters
        ----------
        contact_resistance : float, pd.Series, or str
            The contact resistance at the electrodes, in Ohms. Can
            be a single value for all measurements, a Series, or
            the name of a column in the dataset.
        setup_length : float, pd.Series, or str
            The length of the setup wire, in meters. This is often
            the distance from the GDP to the dipole center.
        wire_capacitance : float, pd.Series, or str, default 15.0
            The wire-to-ground capacitance in picofarads per meter
            (pF/m). This is often an empirical tuning parameter.
        update_components : bool, default True
            If ``True``, the main DataFrame (`self.avg.info.df`) is
            updated in place with the corrected E-field values, and
            all dependent components (Z, rho, phase) are
            recalculated.

        Returns
        -------
        pandas.DataFrame
            A new DataFrame containing the corrected `emag` and
            `ephz` columns.

        Notes
        -----
        The correction is based on a simple circuit model where the
        capacitive admittance (:math:`Y_c`) shunts the true Earth
        impedance. The measured voltage is corrected to estimate
        the true voltage that would be measured without this effect.

        The correction factor is given by:

        .. math::
            E_{true} = \frac{E_{measured}}{1 + Z_c Y_c}

        where :math:`Z_c` is the contact impedance (resistance) and
        :math:`Y_c = i \omega C` is the capacitive admittance.
        """
        has_read(self, attributes="avg")
        df = self.avg.df.copy()

        required = ["emag", "ephz", "freq"]
        missing = [c for c in required if c not in df.columns]
        if missing:
            raise ProcessingError(
                "DataFrame is missing required columns for capacitive "
                f"coupling correction: {missing}"
            )

        # --- Standardize inputs into Series ---
        def _resolve_param(param: float | pd.Series | str, name: str) -> pd.Series:
            if isinstance(param, str):
                if param not in df.columns:
                    raise ProcessingError(f"Column '{param}' not found.")
                return df[param]
            elif isinstance(param, pd.Series):
                return param
            return pd.Series(param, index=df.index)

        Zc = _resolve_param(contact_resistance, "contact_resistance")
        L = _resolve_param(setup_length, "setup_length")
        C_per_m = (
            _resolve_param(wire_capacitance, "wire_capacitance") * 1e-12
        )  # pF/m to F/m

        # --- Perform Correction ---
        omega = 2 * PI * df["freq"]
        # Total capacitance C = C_per_meter * Length
        C_total = C_per_m * L
        # Capacitive admittance Yc = i * omega * C
        Yc = 1j * omega * C_total

        # Reconstruct the complex measured E-field
        e_phase_rad = df["ephz"] * 1e-3
        E_measured = df["emag"] * np.exp(1j * e_phase_rad)

        # Apply the correction
        E_corrected = E_measured / (1 + Zc * Yc)

        # Deconstruct back into magnitude and phase
        df["emag"] = np.abs(E_corrected)
        df["ephz"] = np.angle(E_corrected, deg=False) * 1000.0  # to mrad

        if update_components:
            # Update the main AVG object and trigger recalculation
            # of all dependent components (Z, rho, phase, etc.)
            self.avg.info.read(df, self.avg.info.meta)
            if self.avg.verbose:
                self._logger.info(
                    "Capacitive coupling correction applied. " "AVG object updated."
                )

        return df[["emag", "ephz"]]

    def correct_static_shift(
        self,
        reference_freq: float,
        *,
        filter_method: Literal["tma", "flma", "ama"] = "tma",
        update_components: bool = True,
        **kwargs,
    ) -> pd.DataFrame:
        r"""Correct for static shift using a spatial filter.

        This method applies a spatial filter to the data at a
        single reference frequency to estimate and correct for
        static shift effects.

        Parameters
        ----------
        reference_freq : float
            The frequency (in Hz) at which to perform the spatial
            filtering. This should typically be the highest
            frequency with clean data.
        filter_method : {'tma', 'flma', 'ama'}, default 'tma'
            The filtering algorithm to use:
            - 'tma': Trimmed Moving Average (works on resistivity).
            - 'flma': Fixed-Length Moving Average (works on impedance).
            - 'ama': Adaptive Moving Average (works on impedance).
        update_components : bool, default True
            If ``True``, the main DataFrame (`self.avg.info.df`) is
            updated in place with the static-shifted resistivity
            values, and all components are re-read.
        **kwargs
            Additional keyword arguments to be passed to the chosen
            filtering function (e.g., `window_size` for `tma`,
            `filter_width_dipoles` for `flma`).

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the station, the original data
            profile, the smoothed profile, and the shift factor.
        """
        has_read(self, attributes="avg")
        if self.avg.df is None:
            raise NotReadError("AVG object's DataFrame is not loaded.")

        df = self.avg.df
        required = ["station", "freq", "rho", "phase", "comp"]
        if not all(c in df.columns for c in required):
            raise ProcessingError(
                "DataFrame must contain 'station', 'freq', 'rho', "
                "'phase', and 'comp' columns."
            )

        stations = sorted(df["station"].unique())

        # Interpolate data to the reference frequency for all stations
        interp_data = {}
        for station in stations:
            st_data = df[df["station"] == station].sort_values("freq")
            if len(st_data) < 2:
                continue

            log_freq = np.log(st_data["freq"])
            interp_data[station] = {
                "rho": np.exp(
                    np.interp(
                        np.log(reference_freq),
                        log_freq,
                        np.log(st_data["rho"]),
                    )
                ),
                "phase": np.interp(np.log(reference_freq), log_freq, st_data["phase"]),
            }

        rho_profile = pd.Series({s: d["rho"] for s, d in interp_data.items()})

        # Dispatch to the appropriate filter
        if filter_method == "tma":
            smoothed_rho = tma(rho_profile, **kwargs)

        elif filter_method in ("flma", "ama"):
            omega = 2 * PI * reference_freq
            phase_profile = pd.Series({s: d["phase"] for s, d in interp_data.items()})
            # Calculate complex impedance at the reference frequency
            z_profile = np.sqrt(rho_profile * omega * MU_0) * np.exp(
                1j * phase_profile * 1e-3
            )

            if filter_method == "flma":
                smoothed_z = flma(z_profile, rho_profile.index, **kwargs)
            else:  # ama
                smoothed_z = ama(
                    z_profile,
                    rho_profile.index,
                    frequency=reference_freq,
                    **kwargs,
                )

            # Convert smoothed impedance back to resistivity
            smoothed_rho = (np.abs(smoothed_z) ** 2) / (omega * MU_0)
        else:
            raise ValueError(f"Unknown filter method: '{filter_method}'")

        shift_factors = smoothed_rho / rho_profile

        if update_components:
            df_copy = df.copy()
            df_copy["rho"] *= df_copy["station"].map(shift_factors)
            df_copy["rho_sc"] = df_copy["rho"]

            self.avg.info.read(df_copy, self.avg.info.meta)
            if self.avg.verbose:
                self._logger.info(
                    f"{filter_method.upper()} static shift applied. "
                    "AVG object updated."
                )

        return pd.DataFrame(
            {
                "station": stations,
                "rho_original": rho_profile,
                "rho_smoothed": smoothed_rho,
                "shift_factor": shift_factors,
            }
        )
