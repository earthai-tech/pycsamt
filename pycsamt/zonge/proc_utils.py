# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
pycsamt.zonge.proc_utils
------------------------

This module provides the core filtering algorithms used for
static shift correction, mirroring the methods described in the
Zonge ASTATIC program manual.
"""

from __future__ import annotations

import warnings
from collections.abc import Iterable
from typing import (
    Literal,
    overload,
)

import numpy as np
import pandas as pd
from scipy.interpolate import UnivariateSpline, interp1d
from scipy.signal import hilbert
from scipy.signal.windows import hann
from scipy.stats import trim_mean

from ..constants import MU_0, PI
from ..decorators import isdf
from ..exceptions import ProcessingError

__all__ = [
    "tma",
    "flma",
    "ama",
    "interpolate_to_log_space",
    "smooth_rho_from_phase",
    "get_reference_frequency",
    "get_skew",
    "get_strike",
    "prepare_strike_frame",
]


@overload
def tma(
    rho_profile: str,
    *,
    data: pd.DataFrame,
    window_size: int = 5,
    trim_proportion: float = 0.2,
) -> pd.DataFrame:
    ...


@overload
def tma(
    rho_profile: np.ndarray,
    *,
    data: None = None,
    window_size: int = 5,
    trim_proportion: float = 0.2,
) -> np.ndarray:
    ...


@overload
def tma(
    rho_profile: pd.Series,
    *,
    data: None = None,
    window_size: int = 5,
    trim_proportion: float = 0.2,
) -> pd.Series:
    ...


def tma(
    rho_profile: pd.Series | np.ndarray | str,
    *,
    data: pd.DataFrame | None = None,
    window_size: int = 5,
    trim_proportion: float = 0.2,
) -> pd.Series | np.ndarray | pd.DataFrame:
    r"""Apply a Trimmed Moving Average (TMA) filter.

    This filter is designed to remove single-station static
    offsets while preserving broader geological trends. It can
    operate on a pandas Series, a NumPy array, or a column
    within a pandas DataFrame.

    Parameters
    ----------
    rho_profile : pd.Series, np.ndarray, or str
        The apparent resistivity data to be filtered. This can be:
        - A pandas Series of resistivity values.
        - A 1D NumPy array of resistivity values.
        - The string name of the column to use if `data` is
          provided.
    data : pd.DataFrame, optional
        A DataFrame containing the resistivity data. Required if
        `rho_profile` is provided as a string.
    window_size : int, default 5
        The size of the moving window for the filter. Must be
        an odd integer.
    trim_proportion : float, default 0.2
        The proportion of observations to trim from each end of
        the window before computing the mean. For a window of 5,
        0.2 trims 1 value from each end (the min and max).

    Returns
    -------
    pd.Series, np.ndarray, or pd.DataFrame
        The smoothed resistivity profile, returned in the same
        format as the input (`rho_profile`). If a DataFrame was
        passed to `data`, the function returns a new DataFrame
        with an added column for the smoothed data.

    Raises
    ------
    ValueError
        If `rho_profile` is a string but `data` is not provided,
        or if `window_size` is not an odd integer.
    """
    if window_size % 2 == 0:
        raise ValueError("`window_size` must be an odd integer.")

    input_type = "series"  # Default
    if isinstance(rho_profile, str):
        if data is None:
            raise ValueError(
                "A DataFrame must be provided via the 'data' "
                "parameter when 'rho_profile' is a string."
            )
        series = data[rho_profile]
        input_type = "dataframe"
    elif isinstance(rho_profile, np.ndarray):
        series = pd.Series(rho_profile)
        input_type = "array"
    else:
        series = rho_profile

    # Perform the rolling trimmed mean calculation
    smoothed_series = series.rolling(
        window=window_size, center=True, min_periods=window_size // 2 + 1
    ).apply(lambda x: trim_mean(x.dropna(), trim_proportion))

    # Return the data in the same format as the input
    if input_type == "dataframe":
        df_out = data.copy()
        df_out[f"{rho_profile}_tma"] = smoothed_series
        return df_out
    elif input_type == "array":
        return smoothed_series.to_numpy()

    return smoothed_series


@overload
def flma(
    z_profile: str,
    stations: str,
    dipole_length: float,
    *,
    data: pd.DataFrame,
    filter_width_dipoles: float = 5.0,
) -> pd.DataFrame:
    ...


@overload
def flma(
    z_profile: np.ndarray,
    stations: np.ndarray,
    dipole_length: float,
    *,
    data: None = None,
    filter_width_dipoles: float = 5.0,
) -> np.ndarray:
    ...


@overload
def flma(
    z_profile: pd.Series,
    stations: pd.Series,
    dipole_length: float,
    *,
    data: None = None,
    filter_width_dipoles: float = 5.0,
) -> pd.Series:
    ...


def flma(
    z_profile: pd.Series | np.ndarray | str,
    stations: pd.Series | np.ndarray | str,
    dipole_length: float,
    *,
    data: pd.DataFrame | None = None,
    filter_width_dipoles: float = 5.0,
) -> pd.Series | np.ndarray | pd.DataFrame:
    r"""Apply a Fixed-Length Moving Average (FLMA) filter.

    This filter smooths complex impedance data using a spatial
    Hanning window whose width is a fixed multiple of the
    receiver dipole length.

    Parameters
    ----------
    z_profile : pd.Series, np.ndarray, or str
        The complex impedance data (Z) to be filtered. Can be:
        - A pandas Series of complex impedance values.
        - A 1D NumPy array of complex impedance values.
        - The string name of the column if `data` is provided.
    stations : pd.Series, np.ndarray, or str
        The station locations (coordinates). Must correspond to
        the `z_profile` data. Can be a Series, array, or column
        name if `data` is provided.
    dipole_length : float
        The length of the E-field receiver dipole in the same
        units as the station locations.
    data : pd.DataFrame, optional
        A DataFrame containing the impedance and station data.
        Required if `z_profile` or `stations` are strings.
    filter_width_dipoles : float, default 5.0
        The total width of the Hanning window, expressed in
        multiples of the `dipole_length`.

    Returns
    -------
    pd.Series, np.ndarray, or pd.DataFrame
        The smoothed complex impedance profile, returned in the
        same format as the input. If a DataFrame was passed to
        `data`, a new DataFrame with an added column for the
        smoothed data is returned.
    """
    # --- Input Validation and Standardization ---
    input_type = "series"
    z_series: pd.Series
    stn_series: pd.Series

    if isinstance(z_profile, str):
        if data is None or not isinstance(data, pd.DataFrame):
            raise ValueError(
                "A DataFrame must be provided via the 'data' "
                "parameter when 'z_profile' is a string."
            )
        z_series = data[z_profile]
        stn_series = data[stations]
        input_type = "dataframe"
    elif isinstance(z_profile, np.ndarray):
        z_series = pd.Series(z_profile)
        stn_series = pd.Series(stations)
        input_type = "array"
    else:
        z_series = z_profile
        stn_series = stations

    # --- Core Filtering Logic ---
    smoothed_z = pd.Series(index=z_series.index, dtype="complex128")
    filter_radius = (filter_width_dipoles / 2.0) * dipole_length

    for i, stn_center in enumerate(stn_series):
        window_min = stn_center - filter_radius
        window_max = stn_center + filter_radius

        window_mask = (stn_series >= window_min) & (stn_series <= window_max)
        z_in_window = z_series[window_mask]

        if z_in_window.empty or z_in_window.isnull().all():
            smoothed_z.iloc[i] = np.nan
            continue

        weights = hann(len(z_in_window))
        smoothed_z.iloc[i] = np.average(z_in_window.dropna(), weights=weights)

    # --- Return data in the original format ---
    if input_type == "dataframe":
        df_out = data.copy()
        df_out[f"{z_profile}_flma"] = smoothed_z
        return df_out
    elif input_type == "array":
        return smoothed_z.to_numpy()

    return smoothed_z


@overload
def ama(
    z_profile: str,
    stations: str,
    dipole_length: float,
    frequency: float,
    *,
    data: pd.DataFrame,
    skin_depth_factor: float = 2.0,
    iterations: int = 3,
) -> pd.DataFrame:
    ...


@overload
def ama(
    z_profile: np.ndarray,
    stations: np.ndarray,
    dipole_length: float,
    frequency: float,
    *,
    data: None = None,
    skin_depth_factor: float = 2.0,
    iterations: int = 3,
) -> np.ndarray:
    ...


@overload
def ama(
    z_profile: pd.Series,
    stations: pd.Series,
    dipole_length: float,
    frequency: float,
    *,
    data: None = None,
    skin_depth_factor: float = 2.0,
    iterations: int = 3,
) -> pd.Series:
    ...


def ama(
    z_profile: pd.Series | np.ndarray | str,
    stations: pd.Series | np.ndarray | str,
    dipole_length: float,
    frequency: float,
    *,
    data: pd.DataFrame | None = None,
    skin_depth_factor: float = 2.0,
    iterations: int = 3,
) -> pd.Series | np.ndarray | pd.DataFrame:
    r"""Apply an Adaptive Moving Average (AMA) filter.

    This filter smooths complex impedance data using a spatial
    Hanning window whose width is adapted iteratively based on the
    local skin depth.

    Parameters
    ----------
    z_profile : pd.Series, np.ndarray, or str
        The complex impedance data (Z) to be filtered.
    stations : pd.Series, np.ndarray, or str
        The station locations (coordinates).
    dipole_length : float
        The length of the E-field receiver dipole.
    frequency : float
        The frequency of the data in Hz.
    data : pd.DataFrame, optional
        A DataFrame containing the impedance and station data.
        Required if `z_profile` or `stations` are strings.
    skin_depth_factor : float, default 2.0
        The filter width, expressed as a multiple of the local
        skin depth.
    iterations : int, default 3
        The number of iterations to perform to allow the filter
        width to adapt to the smoothed data.

    Returns
    -------
    pd.Series, np.ndarray, or pd.DataFrame
        The smoothed complex impedance profile, returned in the
        same format as the input.
    """
    # --- Input Validation and Standardization ---
    input_type = "series"
    z_series: pd.Series
    stn_series: pd.Series

    if isinstance(z_profile, str):
        if data is None or not isinstance(data, pd.DataFrame):
            raise ValueError(
                "A DataFrame must be provided via 'data' when "
                "'z_profile' is a string."
            )
        z_series = data[z_profile]
        stn_series = data[stations]
        input_type = "dataframe"
    elif isinstance(z_profile, np.ndarray):
        z_series = pd.Series(z_profile)
        stn_series = pd.Series(stations)
        input_type = "array"
    else:
        z_series = z_profile
        stn_series = stations

    # --- Core Filtering Logic ---
    omega = 2 * PI * frequency
    smoothed_z = z_series.copy()

    for _ in range(iterations):
        rho_local = (np.abs(smoothed_z) ** 2) / (omega * MU_0)
        skin_depth = 503 * np.sqrt(rho_local / frequency)
        filter_radius = (skin_depth_factor / 2.0) * skin_depth

        current_pass_z = pd.Series(index=z_series.index, dtype="complex128")
        for i, stn_center in enumerate(stn_series):
            radius = filter_radius.iloc[i]
            window_min = stn_center - radius
            window_max = stn_center + radius

            window_mask = (stn_series >= window_min) & (stn_series <= window_max)
            z_in_window = smoothed_z[window_mask]

            if z_in_window.empty or z_in_window.isnull().all():
                current_pass_z.iloc[i] = np.nan
                continue

            weights = hann(len(z_in_window))
            current_pass_z.iloc[i] = np.average(z_in_window.dropna(), weights=weights)
        smoothed_z = current_pass_z

    # --- Return data in the original format ---
    if input_type == "dataframe":
        df_out = data.copy()
        df_out[f"{z_profile}_ama"] = smoothed_z
        return df_out
    elif input_type == "array":
        return smoothed_z.to_numpy()

    return smoothed_z


@isdf
def interpolate_to_log_space(
    df: pd.DataFrame,
    *,
    freq_min: float | None = None,
    freq_max: float | None = None,
    num_points: int = 50,
    interp_kind: str = "cubic",
) -> pd.DataFrame:
    r"""Interpolate AVG data onto a regular log-spaced grid.

    This function resamples the apparent resistivity and phase
    data for each station and component onto a new, logarithmically
    spaced frequency axis. This is a common preprocessing step for
    inversion and modeling.

    Parameters
    ----------
    df : pandas.DataFrame
        The input DataFrame containing the AVG data. Must include
        'station', 'freq', 'rho', and 'phase' columns.
    freq_min : float, optional
        The minimum frequency for the new grid. If ``None``, it is
        inferred from the minimum frequency in the data.
    freq_max : float, optional
        The maximum frequency for the new grid. If ``None``, it is
        inferred from the maximum frequency in the data.
    num_points : int, default 50
        The number of logarithmically spaced points to create
        between `freq_min` and `freq_max`.
    interp_kind : str, default 'cubic'
        The kind of interpolation to perform, passed to
        `scipy.interpolate.interp1d`. Common options are
        'slinear', 'quadratic', and 'cubic'.

    Returns
    -------
    pandas.DataFrame
        A new DataFrame with the data interpolated onto the
        specified logarithmic frequency grid.
    """
    required_cols = ["station", "freq", "rho", "phase", "comp"]
    if not all(col in df.columns for col in required_cols):
        raise ProcessingError(
            "Input DataFrame must contain 'station', 'freq', 'rho', "
            "'phase', and 'comp' columns."
        )

    # Determine frequency range for the new grid
    f_min = freq_min or df["freq"].min()
    f_max = freq_max or df["freq"].max()
    new_freq_log = np.logspace(np.log10(f_min), np.log10(f_max), num_points)

    interpolated_data = []
    # Group by each unique sounding curve
    for (station, comp), group in df.groupby(["station", "comp"]):
        # Sort by frequency to ensure correct interpolation order
        sounding = group.sort_values("freq").dropna(subset=["rho", "phase"])
        if len(sounding) < 2:
            continue  # Cannot interpolate with fewer than 2 points

        log_freq = np.log(sounding["freq"])
        log_rho = np.log(sounding["rho"])
        phase_rad = sounding["phase"] * 1e-3  # to radians

        # Create interpolation functions
        interp_rho = interp1d(
            log_freq,
            log_rho,
            kind=interp_kind,
            bounds_error=False,
            fill_value=np.nan,
        )
        interp_phase = interp1d(
            log_freq,
            phase_rad,
            kind=interp_kind,
            bounds_error=False,
            fill_value=np.nan,
        )

        # Evaluate at the new frequency points
        new_rho_log = interp_rho(np.log(new_freq_log))
        new_phase_rad = interp_phase(np.log(new_freq_log))

        # Assemble the new DataFrame for this sounding
        interp_sounding = pd.DataFrame(
            {
                "station": station,
                "comp": comp,
                "freq": new_freq_log,
                "rho": np.exp(new_rho_log),
                "phase": new_phase_rad * 1000.0,  # back to mrad
            }
        )
        interpolated_data.append(interp_sounding)

    if not interpolated_data:
        return pd.DataFrame()  # Return empty if no data was processed

    return pd.concat(interpolated_data, ignore_index=True)


@isdf
def smooth_rho_from_phase(
    df: pd.DataFrame,
    *,
    smoothing_factor: float = 0.1,
) -> pd.DataFrame:
    r"""Smooth apparent resistivity using the Hilbert transform.

    This function reconstructs a smoother apparent resistivity
    curve from the impedance phase data, leveraging the causal
    nature of magnetotelluric data.

    Parameters
    ----------
    df : pandas.DataFrame
        The input DataFrame. Must include 'station', 'freq',
        'rho', and 'phase' columns.
    smoothing_factor : float, default 0.1
        The smoothing factor for the `UnivariateSpline` applied
        to the phase data before the Hilbert transform. A smaller
        value results in less smoothing.

    Returns
    -------
    pandas.DataFrame
        A new DataFrame with a 'rho_smoothed' column containing
        the reconstructed, smoother apparent resistivity values.
    """
    required_cols = ["station", "freq", "rho", "phase", "comp"]
    if not all(col in df.columns for col in required_cols):
        raise ProcessingError(
            "Input DataFrame must contain 'station', 'freq', 'rho', "
            "'phase', and 'comp' columns."
        )

    df_out = df.copy()
    df_out["rho_smoothed"] = np.nan

    # Group by each unique sounding curve
    for (_station, _comp), group in df.groupby(["station", "comp"]):
        sounding = group.sort_values("freq").dropna(subset=["rho", "phase"])
        if len(sounding) < 5:  # Need enough points for spline
            continue

        log_freq = np.log(sounding["freq"])
        phase_rad = sounding["phase"] * 1e-3  # to radians

        # 1. Smooth the phase data with a spline
        spline = UnivariateSpline(log_freq, phase_rad, s=smoothing_factor)
        phase_smooth = spline(log_freq)

        # 2. Apply Hilbert transform to get the imaginary part of ln(Z)
        # The imaginary part of log(Z) is the phase
        imag_log_z = phase_smooth

        # 3. The real part of log(Z) is the Hilbert transform of the phase
        real_log_z = -np.imag(hilbert(imag_log_z))

        # 4. Reconstruct the impedance magnitude
        log_z_mag = real_log_z
        z_mag_smooth = np.exp(log_z_mag)

        # 5. Calculate the new, smoother apparent resistivity
        omega = 2 * PI * sounding["freq"]
        rho_smooth = (z_mag_smooth**2) / (omega * MU_0)

        # Place the smoothed values back into the output DataFrame
        df_out.loc[sounding.index, "rho_smoothed"] = rho_smooth

    return df_out


@isdf
def get_reference_frequency(
    df: pd.DataFrame,
    mode: str | float = "auto",
    *,
    qc_column: str = "pc_rho",
    qc_threshold: float = 20.0,
) -> float:
    r"""Determine a suitable reference frequency for static shift.

    This utility automatically selects a reference frequency based
    on the principle of using the "highest frequency with clean
    data," as recommended in the ASTATIC manual.

    Parameters
    ----------
    df : pandas.DataFrame
        The input DataFrame containing the AVG data. Must include
        'station', 'freq', and a quality control column.
    mode : {'auto'} or float, default 'auto'
        The method for determining the reference frequency.
        - 'auto': Automatically select the frequency based on QC
          metrics.
        - float: A specific frequency value to use, which will be
          returned directly.
    qc_column : str, default 'pc_rho'
        The canonical name of the quality control column to use for
        filtering "clean" data points.
    qc_threshold : float, default 20.0
        The threshold for the `qc_column`. Data points where the
        value in this column is *below* the threshold are
        considered clean.

    Returns
    -------
    float
        The determined reference frequency in Hz.

    Notes
    -----
    The 'auto' mode follows a robust, data-driven approach:
    1. It filters the dataset to include only "clean" data points
       (where `qc_column` < `qc_threshold`).
    2. For each station, it finds the highest frequency that has
       clean data.
    3. It then returns the **median** of these highest frequencies.
       The median is used to provide a stable estimate that is
       robust against outlier stations.

    If no data points meet the clean criteria, the function will
    issue a warning and fall back to using the absolute maximum
    frequency present in the entire dataset.
    """
    if isinstance(mode, (int, float)):
        return float(mode)

    required_cols = ["station", "freq", qc_column]
    if not all(col in df.columns for col in required_cols):
        raise ProcessingError(
            "Input DataFrame must contain 'station', 'freq', and "
            f"the specified qc_column ('{qc_column}')."
        )

    # 1. Filter for "clean" data points
    clean_df = df[df[qc_column] < qc_threshold]

    if clean_df.empty:
        warnings.warn(
            f"No data found with '{qc_column}' < {qc_threshold}. "
            "Falling back to the absolute maximum frequency.",
            stacklevel=2,
        )
        return df["freq"].max()

    # 2. Find the highest clean frequency for each station
    max_freq_per_station = clean_df.groupby("station")["freq"].max()

    # 3. Return the median of these frequencies for robustness
    return max_freq_per_station.median()


@isdf
def get_strike(df: pd.DataFrame) -> pd.DataFrame:
    r"""Calculate the geoelectric strike angle.

    This function determines the principal axis direction of the
    impedance tensor for each station and frequency. The strike is
    the angle to which the coordinate system must be rotated to
    minimize the diagonal components of the impedance tensor
    (:math:`Z_{xx}`, :math:`Z_{yy}`), which is a common assumption
    for 2D geological structures.

    Parameters
    ----------
    df : pandas.DataFrame
        A tidy DataFrame containing the impedance tensor components.
        Must include 'station', 'freq', 'comp', and a complex 'z'
        column from which the tensor elements can be pivoted.

    Returns
    -------
    pandas.DataFrame
        A DataFrame with columns for 'station', 'freq', and the
        calculated 'strike_angle' in degrees clockwise from North.

    Notes
    -----
    The calculation is performed using a common tensor
    decomposition method based on the real parts of the impedance
    tensor components [1]_. The formula used is:

    .. math::
        \theta_s = \frac{1}{2} \arctan\left(
        \frac{\text{Re}(Z_{xy} + Z_{yx})}
        {\text{Re}(Z_{xx} - Z_{yy})}
        \right)

    This angle represents the orientation of the primary 2D
    geological structure.

    Examples
    --------
    >>> from pycsamt.zonge import AMTAVG
    >>> from pycsamt.zonge.proc_utils import get_strike
    >>> avg = AMTAVG.from_file("data/avg/K2.avg")
    >>> # The Z component's frame already has the required columns
    >>> strike_df = get_strike(avg.z.frame)
    >>> print(strike_df.head())

    References
    ----------
    .. [1] Simpson, F., & Bahr, K. (2005). *Practical
           Magnetotellurics*. Cambridge University Press.

    See Also
    --------
    get_skew : Calculate the dimensionality indicator (skew).
    pycsamt.zonge.avg.AMTAVG.rotate : Rotate the tensor by a given
        angle.
    """
    required = ["station", "freq", "comp", "z"]
    if not all(c in df.columns for c in required):
        raise ProcessingError(
            f"Input DataFrame is missing one of required columns: {required}"
        )

    # Pivot to get tensor components for each measurement point
    # tensor_df = df.pivot_table(
    #     index=['station', 'freq'], columns='comp', values='z'
    # ).reset_index()

    # Split complex data before pivoting ---
    df_copy = df.copy()
    df_copy["z_real"] = np.real(df_copy["z"])
    df_copy["z_imag"] = np.imag(df_copy["z"])

    tensor_df_real = df_copy.pivot_table(
        index=["station", "freq"], columns="comp", values="z_real"
    )
    tensor_df_imag = df_copy.pivot_table(
        index=["station", "freq"], columns="comp", values="z_imag"
    )
    # Recombine into a complex tensor after pivoting
    tensor_df = (tensor_df_real + 1j * tensor_df_imag).reset_index()

    # Ensure all four components are present
    for comp in ["ExHx", "ExHy", "EyHx", "EyHy"]:
        if comp not in tensor_df.columns:
            tensor_df[comp] = 0 + 0j

    zxx, zxy = tensor_df["ExHx"], tensor_df["ExHy"]
    zyx, zyy = tensor_df["EyHx"], tensor_df["EyHy"]

    # Calculate strike angle using tensor decomposition
    real_a = np.real(zxy + zyx)
    real_b = np.real(zxx - zyy)
    strike_rad = 0.5 * np.arctan2(real_a, real_b)
    strike_deg = np.rad2deg(strike_rad)

    return pd.DataFrame(
        {
            "station": tensor_df["station"],
            "freq": tensor_df["freq"],
            "strike_angle": strike_deg,
        }
    )


@isdf
def get_skew(df: pd.DataFrame) -> pd.DataFrame:
    r"""Calculate Swift's skew for the impedance tensor.

    Skew is a dimensionless parameter that quantifies the
    three-dimensionality of the subsurface conductivity structure.
    It is a crucial diagnostic tool for determining whether a 1D
    or 2D interpretation is appropriate for the data.

    Parameters
    ----------
    df : pandas.DataFrame
        A tidy DataFrame containing the impedance tensor components.
        Must include 'station', 'freq', 'comp', and a complex 'z'
        column.

    Returns
    -------
    pandas.DataFrame
        A DataFrame with columns for 'station', 'freq', and the
        calculated 'skew'.

    Notes
    -----
    The function calculates the conventional skew definition
    proposed by Swift [1]_. The formula is:

    .. math::
        \text{skew} = \frac{|Z_{xx} + Z_{yy}|}{|Z_{xy} - Z_{yx}|}

    Generally, skew values are interpreted as follows:
    - **skew < 0.1**: Data can be considered 1D.
    - **0.1 <= skew <= 0.3**: Data is likely 2D.
    - **skew > 0.3**: Data is likely 3D or affected by significant
      local distortion.

    Examples
    --------
    >>> from pycsamt.zonge import AMTAVG
    >>> from pycsamt.zonge.proc_utils import get_skew
    >>> avg = AMTAVG.from_file("data/avg/K2.avg")
    >>> skew_df = get_skew(avg.z.frame)
    >>> print(skew_df.head())

    References
    ----------
    .. [1] Swift, C. M. (1967). A magnetotelluric investigation
           of an electrical conductivity anomaly in the southwestern
           United States. Ph.D. Thesis, MIT.

    See Also
    --------
    get_strike : Calculate the geoelectric strike angle.
    """
    required = ["station", "freq", "comp", "z"]
    if not all(c in df.columns for c in required):
        raise ProcessingError(
            f"Input DataFrame is missing one of required columns: {required}"
        )

    # Pivot to get tensor components for each measurement point
    df_copy = df.copy()
    df_copy["z_real"] = np.real(df_copy["z"])
    df_copy["z_imag"] = np.imag(df_copy["z"])

    tensor_df_real = df_copy.pivot_table(
        index=["station", "freq"], columns="comp", values="z_real"
    )
    tensor_df_imag = df_copy.pivot_table(
        index=["station", "freq"], columns="comp", values="z_imag"
    )
    # Recombine into a complex tensor after pivoting
    tensor_df = (tensor_df_real + 1j * tensor_df_imag).reset_index()

    # Ensure all four components are present
    for comp in ["ExHx", "ExHy", "EyHx", "EyHy"]:
        if comp not in tensor_df.columns:
            tensor_df[comp] = 0 + 0j

    zxx, zxy = tensor_df["ExHx"], tensor_df["ExHy"]
    zyx, zyy = tensor_df["EyHx"], tensor_df["EyHy"]

    # Calculate Swift's skew
    numerator = np.abs(zxx + zyy)
    denominator = np.abs(zxy - zyx)
    # Avoid division by zero
    skew = np.divide(
        numerator,
        denominator,
        out=np.full_like(numerator, np.nan, dtype=float),
        where=(denominator != 0),
    )

    return pd.DataFrame(
        {
            "station": tensor_df["station"],
            "freq": tensor_df["freq"],
            "skew": skew,
        }
    )


def prepare_strike_frame(
    *,
    z_frame: pd.DataFrame | None = None,
    df: pd.DataFrame | None = None,
    prefer: Literal["z", "df"] = "z",
    station_col: str = "station",
    freq_col: str = "freq",
    comp_col: str = "comp",
    z_col: str = "z",
    rho_col: str = "rho",
    phase_col: str = "phase",
    phase_unit: Literal["auto", "deg", "mrad", "rad"] = "auto",
    drop_na: bool = True,
    na_policy: Literal["any", "all"] = "any",
    components: Iterable[str] | None = None,
    ensure_sorted: bool = True,
    copy: bool = True,
    mu0: float = 4.0 * np.pi * 1e-7,
) -> pd.DataFrame:
    r"""
    Build a strike-ready frame with complex ``z``.

    Prefers an input frame that already contains a complex
    impedance column. If not available, derives ``z`` from
    ``rho``, ``phase`` and ``freq``:

    .. math::
        |Z| = \sqrt{\mu_0 \, \omega \, \rho_a},\quad
        Z = |Z|\,[\cos\phi + i\sin\phi]

    Parameters
    ----------
    z_frame : DataFrame, optional
        Frame that already holds a complex column ``z``.
    df : DataFrame, optional
        Long frame with columns for apparent resistivity,
        phase and frequency (``rho``, ``phase``, ``freq``).
    prefer : {"z","df"}, default "z"
        Source priority when both frames are given.
    station_col, freq_col, comp_col, z_col, rho_col, phase_col
        Column names in provided frames.
    phase_unit : {"auto","deg","mrad","rad"}, default "auto"
        Unit for ``phase``. ``"auto"`` guesses:
        - ``|phase|max <= π*1.5`` → rad
        - else if ``<= 180`` → deg
        - else → mrad
    drop_na : bool, default True
        Drop NA rows before computing.
    na_policy : {"any","all"}, default "any"
        Row-drop policy when ``drop_na=True``.
    components : iterable of str, optional
        If given, keep only these components.
    ensure_sorted : bool, default True
        Sort output by (station, freq, comp).
    copy : bool, default True
        Work on a copy of the source frame(s).
    mu0 : float, default 4π·1e−7
        Magnetic permeability of free space.

    Returns
    -------
    DataFrame
        Columns: ``station``, ``freq``, ``comp``, ``z``.

    Raises
    ------
    ProcessingError
        If required columns are missing or inputs are
        not provided.
    """
    need_base = {station_col, freq_col, comp_col}

    def _subset(fr: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
        miss = [c for c in cols if c not in fr.columns]
        if miss:
            raise ProcessingError(f"Missing columns: {miss}")
        out = fr.loc[:, cols]
        return out.copy() if copy else out

    def _mask_components(fr: pd.DataFrame) -> pd.DataFrame:
        if components is None:
            return fr
        return fr[fr[comp_col].isin(list(components))]

    # -------------------- try z_frame path -------------------- #
    if z_frame is not None and prefer == "z":
        cols = [station_col, freq_col, comp_col, z_col]
        try:
            out = _subset(z_frame, cols)
            out = _mask_components(out)
            if drop_na:
                out = out.dropna(how=na_policy)
            # try to ensure complex dtype
            if not np.issubdtype(out[z_col].dtype, np.complexfloating):
                out[z_col] = out[z_col].astype(complex)
            out = out.rename(
                columns={
                    station_col: "station",
                    freq_col: "freq",
                    comp_col: "comp",
                    z_col: "z",
                }
            )
            if ensure_sorted:
                out = out.sort_values(["station", "freq", "comp"])
            return out
        except ProcessingError:
            pass  # fall back to df path

    # --------------------- try df (derive) -------------------- #
    if df is not None:
        cols = list(need_base | {rho_col, phase_col, freq_col})
        out = _subset(df, cols)
        out = _mask_components(out)
        if drop_na:
            out = out.dropna(how=na_policy)

        # coerce numerics
        for c in (rho_col, phase_col, freq_col):
            out[c] = pd.to_numeric(out[c], errors="coerce")
        if drop_na:
            out = out.dropna(how=na_policy)

        if out.empty:
            raise ProcessingError("No data after NA drop.")

        # magnitude from apparent resistivity (SI)
        f = out[freq_col].to_numpy(float)
        r = out[rho_col].to_numpy(float)
        # avoid negative / zero
        f = np.clip(f, 1e-12, np.inf)
        r = np.clip(r, 0.0, np.inf)
        omega = 2.0 * np.pi * f
        zmag = np.sqrt(mu0 * omega * r)

        # phase → radians
        ph = out[phase_col].to_numpy(float)
        if phase_unit == "deg":
            ph_rad = np.deg2rad(ph)
        elif phase_unit == "mrad":
            ph_rad = ph / 1000.0
        elif phase_unit == "rad":
            ph_rad = ph
        else:
            # auto
            max_abs = np.nanmax(np.abs(ph)) if ph.size else 0.0
            if max_abs <= np.pi * 1.5:
                ph_rad = ph  # already radians
            elif max_abs <= 180.0:
                ph_rad = np.deg2rad(ph)
            else:
                ph_rad = ph / 1000.0

        zc = zmag * (np.cos(ph_rad) + 1j * np.sin(ph_rad))

        ret = out.loc[:, [station_col, freq_col, comp_col]]
        ret = ret.copy() if copy else ret
        ret["z"] = zc
        ret = ret.rename(
            columns={
                station_col: "station",
                freq_col: "freq",
                comp_col: "comp",
            }
        )
        if ensure_sorted:
            ret = ret.sort_values(["station", "freq", "comp"])
        return ret

    # ------------------------ fallback z_frame ----------------- #
    if z_frame is not None:
        cols = [station_col, freq_col, comp_col, z_col]
        out = _subset(z_frame, cols)
        out = _mask_components(out)
        if drop_na:
            out = out.dropna(how=na_policy)
        if not np.issubdtype(out[z_col].dtype, np.complexfloating):
            out[z_col] = out[z_col].astype(complex)
        out = out.rename(
            columns={
                station_col: "station",
                freq_col: "freq",
                comp_col: "comp",
                z_col: "z",
            }
        )
        if ensure_sorted:
            out = out.sort_values(["station", "freq", "comp"])
        return out

    raise ProcessingError(
        "Provide 'z_frame' with a complex 'z' column or "
        "'df' with 'rho','phase','freq' to derive z."
    )
