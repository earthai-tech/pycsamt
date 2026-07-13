# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
Statistical utilities.
"""

from collections.abc import Iterable, Sequence
from typing import Any, Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from ..exceptions import StatsError

__all__ = [
    "get_confidence_ratio",
    "drawn_boundaries",
    "remove_outliers",
    "scale_position",
]


def get_confidence_ratio(
    arr: Any,
    axis: Union[int, str] = 0,
    invalid: Optional[Union[Iterable[Any], Any]] = None,
    consider_nan: bool = True,
    as_percent: bool = False,
    average: bool = False,
) -> Union[np.ndarray, float]:
    """
    Compute confidence ratio along an axis by counting valid items.

    A value is "valid" if it is not in `invalid` and,
    optionally, not NaN.

    Parameters
    ----------
    arr : array-like
        Input data (numeric or object) of any dimensionality.
    axis : int or str, default 0
        Axis along which to compute ratio. Strings 'rows', 'cols'
        map to axes 1, 0 respectively for 2D inputs.
    invalid : scalar or iterable, optional
        Value(s) to treat as invalid in the data.
        Default: None (only NaNs if consider_nan=True).
    consider_nan : bool, default True
        If True, treat NaNs as invalid.
    as_percent : bool, default False
        If True, return ratio * 100.
    average : bool, default False
        If True, return mean of the ratio array.

    Returns
    -------
    ratio : ndarray or float
        Confidence ratio per slice or its mean.

    Raises
    ------
    StatsError
        If `axis` is invalid or data cannot be processed.
    """
    # convert to array
    try:
        a = np.asarray(arr)
    except Exception as e:
        raise StatsError(f"Cannot convert input to array: {e}")
    # map axis names
    if isinstance(axis, str):
        name = axis.lower()
        if name in ("cols", "columns", "x"):
            ax = 1
        elif name in ("rows", "y"):
            ax = 0
        else:
            raise StatsError(f"Unknown axis name: {axis}")
    else:
        ax = int(axis)
    nd = a.ndim
    if nd == 0:
        # scalar
        total = 1
        valid = 0
        v = a.item()
        if consider_nan and np.isnan(v):
            valid = 0
        elif invalid is not None and (
            v in invalid if hasattr(invalid, "__iter__") else v == invalid
        ):
            valid = 0
        else:
            valid = 1
        ratio = valid / total
        return ratio * 100 if as_percent else ratio
    # ensure axis in range
    if ax < 0 or ax >= nd:
        raise StatsError(f"Axis {axis} out of bounds for ndim={nd}")
    # prepare invalid mask
    mask = np.zeros_like(a, dtype=bool)
    if consider_nan and np.issubdtype(a.dtype, float):
        mask |= np.isnan(a)
    if invalid is not None:
        inv = np.array(invalid)
        # broadcast equality
        try:
            # elementwise
            mask |= np.isin(a, inv)
        except Exception:
            # fallback scalar compare
            mask |= a == invalid
    # count total and valid per slice
    total = a.shape[ax]
    # sum of invalids along axis
    n_invalid = np.sum(mask, axis=ax)
    n_valid = total - n_invalid
    ratio = n_valid / total
    if as_percent:
        ratio = ratio * 100
    if average:
        return float(np.mean(ratio))
    return ratio


def remove_outliers(
    data: Any,
    method: str = "IQR",
    threshold: float = 3.0,
    fill_value: Optional[float] = None,
    axis: int = 1,
    interpolate: bool = False,
    kind: str = "linear",
) -> Any:
    """Efficient strategy to remove outliers in the data.

    Indeed, an outlier is the data point of the given sample,
    observation, or distribution that shall lie outside the overall pattern.
    A commonly used rule says that one will consider a data point an
    outlier if it has more than 1.5 IQR below the first quartile or above
    the third.

    Two approaches are used to remove the outliers.

    - Inter Quartile Range (``IQR``)
      IQR is the most commonly used and most trusted approach used in
      the research field. Said differently, low outliers shall
      lie below Q1-1.5 IQR, and high outliers shall lie Q3+1.5IQR.
      One needs to calculate median, quartiles, including IQR, Q1,
      and Q3.

      .. math::

        Q1 = 1/4(n + 1)

        Q3 = 1/4 (n + 1)

        Q2 = Q3 – Q1

      To define the outlier base value is defined above and below
      datasets normal range namely Upper and Lower bounds, define the
      upper and the lower bound (1.5*IQR value is considered) :

      .. math::

         upper = Q3 +1.5*IQR

         lower = Q1 – 1.5*IQR

      In the above formula as according to statistics, the 0.5
      scale-up of :math:`IQR (new_IQR = IQR + 0.5*IQR)` is taken, to consider
      all the data between 2.7 standard deviations in the Gaussian
      Distribution

    - Z-score
      Is also called a standard score. This value/score helps to understand
      that how far is the data point from the mean. And after setting up
      a threshold value one can utilize z score values of data points
      to define the outliers.

      .. math::

          Zscore = (\text{data_point} -\text{mean}) / \text{std. deviation}

    Now to define an outlier threshold value is chosen which is
    generally 3.0. As 99.7% of the data points lie between +/- 3 standard
    deviation (using Gaussian Distribution approach).

    .. versionadded: 0.1.5

    Parameters
    -----------
    ar: Arraylike, pd.dataframe
       Arraylike  containing outliers to remove.

       .. versionadded:: 0.2.7
          Accepts dataframe and can remove outliers using the `z_score`.

    method: str, default='IQR'
      The selected approach to remove the outliers. It can be
      ['IQR'|'Z-score']. See Above for outlier explanations.  Note that
      when selecting ``"z-score"`` the threshold value greatly influence
      the quality of data considering as ooutliers.

    threshold: float, default=3
      Thershold values is useful for ``"z-score"`` as the value for considering
      data above as outliers.

    fill_value: float, optional
      Value to replace the outliers. If not given, outliers are suppressed
      in the array.

    axis: int, default=1
      axis from which to remove values. This is useful when two dimensional
      array is supplied. Default, delete outlier from the rows.

    interpolate: bool, default=False,
       If ``fill_value='NaN'``, interpolation can be triggered to get the
       closest value in array to replace missing values. Note that
       `fill_value` should be NaN for interpolation to be concise.

    kind: str, default='linear'
      kind of interpolation. It could be ['nearest'|'linear'|'cubic'].

    .. versionadded:: 0.2.8
       Interpolate NaN value after outliers removal.


    Returns
    --------
    arr: Array_like
        New array whith removed outliers.

    Examples
    ---------
    >>> import numpy as np
    >>> np.random.seed (42 )
    >>> from watex.utils.funcutils import remove_outliers
    >>> data = np.random.randn (7, 3 )
    >>> data_r = remove_outliers ( data )
    >>> data.shape , data_r.shape
    (7, 3) (5, 3)
    >>> remove_outliers ( data, fill_value =np.nan )
    array([[ 0.49671415, -0.1382643 ,  0.64768854],
           [ 1.52302986, -0.23415337, -0.23413696],
           [ 1.57921282,  0.76743473, -0.46947439],
           [ 0.54256004, -0.46341769, -0.46572975],
           [ 0.24196227,         nan,         nan],
           [-0.56228753, -1.01283112,  0.31424733],
           [-0.90802408,         nan,  1.46564877]])
    >>> # for one dimensional
    >>> remove_outliers ( data[:, 0] , fill_value =np.nan )
    array([ 0.49671415,  1.52302986,  1.57921282,  0.54256004,  0.24196227,
           -0.56228753,         nan])
    >>> remove_outliers ( data[:, 0] , fill_value =np.nan, interpolate=True  )
    >>> import matplotlib.pyplot as plt
    >>> plt.plot (np.arange (len(data ), data, ))
    """

    # DataFrame path
    if isinstance(data, pd.DataFrame):
        df = data.copy()
        if method.lower() == "iqr":
            Q1 = df.quantile(0.25)
            Q3 = df.quantile(0.75)
            IQR = Q3 - Q1
            lower = Q1 - 1.5 * IQR
            upper = Q3 + 1.5 * IQR
            mask = (df < lower) | (df > upper)
        elif method.lower() in ("zscore", "z-score"):
            z = (df - df.mean()) / df.std()
            mask = z.abs() > threshold
        else:
            raise StatsError(f"Unknown method: {method}")
        if fill_value is None:
            return df[~mask.any(axis=axis)]
        df = df.mask(mask, other=fill_value)
        if interpolate and pd.isna(fill_value):
            return df.interpolate(method=kind, axis=axis)
        return df
    # ndarray path
    try:
        arr = np.asarray(data, dtype=float)
    except Exception as e:
        raise StatsError(f"Cannot convert input: {e}")
    m = method.lower()
    if m == "iqr":
        flat = arr[~np.isnan(arr)]
        Q1 = np.percentile(flat, 25)
        Q3 = np.percentile(flat, 75)
        IQR = Q3 - Q1
        lower = Q1 - 1.5 * IQR
        upper = Q3 + 1.5 * IQR
        mask = (arr < lower) | (arr > upper)
    elif m in ("zscore", "z-score"):
        mean = np.nanmean(arr)
        std = np.nanstd(arr)
        z = (arr - mean) / std
        mask = np.abs(z) > threshold
    else:
        raise StatsError(f"Unknown method: {method}")
    # apply mask
    if fill_value is None:
        if arr.ndim == 1:
            clean = arr[~mask]
        else:
            clean = arr[~mask.any(axis=axis)]
    else:
        clean = arr.copy()
        clean[mask] = fill_value
        if interpolate and np.isnan(fill_value):
            # use pandas for interpolation
            df = pd.DataFrame(clean)
            df = df.interpolate(method=kind, axis=axis)
            clean = df.values
    return clean


def scale_position(
    y: Union[np.ndarray, pd.Series],
    x: Optional[Union[np.ndarray, pd.Series]] = None,
    func: Optional[callable] = None,
    column: Optional[Union[int, str]] = None,
    initial_params: Optional[Sequence[float]] = None,
    bounds: Optional[tuple] = None,
    dropna: bool = False,
    xy_numeric: bool = False,
    asarray: bool = True,
    show: bool = False,
    plot_kwargs: dict = None,
) -> Union[tuple, np.ndarray]:
    """
    Fit and rescale y vs. x using a model, return scaled y.

    Parameters
    ----------
    y : array-like or Series
        Dependent data to scale.
    x : array-like or Series, optional
        Independent variable. Defaults to range(len(y)).
    func : callable, optional
        Model f(x, *params). Defaults to linear.
    column : int or str, optional
        If y is DataFrame, select column.
    initial_params : sequence, optional
        Initial guess for parameters.
    bounds : tuple, optional
        Bounds for curve_fit.
    dropna : bool, default=False
        Drop NaNs before fitting.
    xy_numeric : bool, default=False
        Coerce x and y to numeric dtype.
    asarray : bool, default=True
        Return numpy arrays if True.
    show : bool, default=False
        Plot data and fit if True.
    plot_kwargs : dict, default {}
        Additional kwargs for plotting.

    Returns
    -------
    y_scaled, popt, pcov
    """

    if plot_kwargs is None:
        plot_kwargs = {}

    def _linear(x, a, b):
        return a * x + b

    f = func if callable(func) else _linear

    # handle DataFrame input
    if isinstance(y, pd.DataFrame):
        if column is None:
            raise ValueError("column required for DataFrame y")
        y = y.iloc[:, column] if isinstance(column, int) else y[column]

    # convert y to Series
    if isinstance(y, (np.ndarray, list)):
        y = pd.Series(y)
    if dropna:
        y = y.dropna()

    # prepare x
    if x is None:
        x = np.arange(len(y))
    if isinstance(x, (np.ndarray, list)):
        x = pd.Series(x)
    y, x = y.align(x, join="inner")

    if xy_numeric:
        y = pd.to_numeric(y, errors="raise")
        x = pd.to_numeric(x, errors="raise")

    popt, pcov = curve_fit(
        f,
        x.values,
        y.values,
        p0=initial_params,
        bounds=bounds or (-np.inf, np.inf),
    )
    y_fit = f(x.values, *popt)

    if show:
        plt.plot(x, y, ".", label="data", **plot_kwargs)
        plt.plot(x, y_fit, "-", label=f"fit: {popt}")
        plt.legend()
        plt.show()

    if asarray:
        return y_fit, popt, pcov
    return pd.Series(y_fit, index=y.index), popt, pcov


def drawn_boundaries(
    profile: Any, peak_value: float, peak_index: int
) -> tuple:
    """
    Identify anomaly boundaries around a peak in a 1D profile.

    Parameters
    ----------
    profile : array-like
        1D data profile.
    peak_value : float
        Value at the anomaly peak.
    peak_index : int
        Index of the anomaly peak in the profile.

    Returns
    -------
    peak_value : float
        The anomaly peak value.
    peak_index : int
        Index of the peak.
    boundaries : ndarray
        Sequence of values defining the anomaly boundary.
    """
    arr = np.asarray(profile, dtype=float)
    n = arr.size
    if not (0 <= peak_index < n):
        raise StatsError(
            f"Peak index {peak_index} out of bounds for length {n}"
        )
    peak = float(peak_value)
    # search left boundary
    left_vals = []
    max_diff = 0.0
    for i in range(peak_index - 1, -1, -1):
        diff = arr[i] - peak
        if diff > max_diff:
            max_diff = diff
            left_vals.append(arr[i])
        else:
            break
    left_vals.reverse()
    # search right boundary
    right_vals = []
    max_diff = 0.0
    for i in range(peak_index + 1, n):
        diff = arr[i] - peak
        if diff > max_diff:
            max_diff = diff
            right_vals.append(arr[i])
        else:
            break
    # assemble boundaries
    if left_vals:
        boundaries = np.concatenate([left_vals, [peak], right_vals])
    else:
        boundaries = np.concatenate([[peak], right_vals])
    return peak, peak_index, boundaries
