# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
Extended math utilities.
"""
import cmath
import warnings
from math import factorial, radians

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.spatial.distance import pdist, squareform

from ..api.typing import (
    ArrayLike,
    DataFrame,
    DType,
    F,
    List,
    NDArray,
    T,
)
from ..exceptions import ResistivityError
from ..gis.utils import utm_to_ll
from .arrayops import (
    is_iterable,
)
from .validation import (
    _assert_all_types,
    _is_arraylike_1d,
    _is_numeric_dtype,
    check_consistency_size,
)

mu0 = 4 * np.pi * 1e-7

__all__ = [
    'round_dipole_length',
    'compute_azimuth',
]



def round_dipole_length(value: float) -> float:
    """
    Round dipole length to the nearest multiple of 5 meters.
    Values within ±2.5 of a multiple of 5 round down, otherwise up.
    """
    if not isinstance(value, (int, float)):
        raise TypeError(f"Value must be numeric, got {type(value)}")
    # Compute nearest multiple of 5
    return float(5 * np.round(value / 5))

def compute_azimuth(
    easting: np.ndarray,
    northing: np.ndarray,
    utm_zone: str = '49N',
    extrapolate: bool = False
) -> np.ndarray:
    """
    Compute azimuth between successive points in UTM coords.

    Parameters
    ----------
    easting : np.ndarray
        Array of easting (m).
    northing : np.ndarray
        Array of northing (m).
    utm_zone : str, optional
        UTM zone for inverse projection (default '49N').
    extrapolate : bool, optional
        If True, prepends an extrapolated azimuth at index 0.
        Otherwise returns length-1 array of azimuths.

    Returns
    -------
    azimuths : np.ndarray
        Azimuth angles (degrees) between points.
    """
    # Validate inputs
    easting = np.asarray(easting, dtype=float)
    northing = np.asarray(northing, dtype=float)
    if easting.shape != northing.shape:
        raise ValueError("Easting and northing must have same shape.")

    n = easting.size
    if n < 2:
        return np.array([])

    # Convert UTM to lat/lon (radians)
    lat_deg, lon_deg = utm_to_ll(
        reference_ellipsoid=23,
        northing=northing,
        easting=easting,
        zone=utm_zone
    )
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)

    # Calculate deltas
    dlon = np.diff(lon)
    sin_dlon = np.sin(dlon)
    cos_dlon = np.cos(dlon)

    # Use haversine-based bearing formula
    numerator = sin_dlon * np.cos(lat[1:])
    denominator = (
        np.cos(lat[:-1]) * np.sin(lat[1:]) -
        np.sin(lat[:-1]) * np.cos(lat[1:]) * cos_dlon
    )
    az = np.arctan2(numerator, denominator)
    az_deg = np.rad2deg(az) % 360

    if extrapolate:
        # extrapolate first value
        f = interp1d(np.arange(1, n), az_deg, fill_value='extrapolate')
        first = f(0)
        az_out = np.concatenate(([first], az_deg))
    else:
        az_out = az_deg

    return np.around(az_out, 3)


def linkage_matrix(
    df: DataFrame ,
    columns:List[str] =None,
    kind:str ='design',
    metric:str ='euclidean',
    method:str ='complete',
    as_frame =False,
    optimal_ordering=False,
 )->NDArray:
    r""" Compute the distance matrix from the hierachical clustering algorithm

    Parameters
    ------------
    df: dataframe or NDArray of (n_samples, n_features)
        dataframe of Ndarray. If array is given , must specify the column names
        to much the array shape 1
    columns: list
        list of labels to name each columns of arrays of (n_samples, n_features)
        If dataframe is given, don't need to specify the columns.

    kind: str, ['squareform'|'condense'|'design'], default is {'design'}
        kind of approach to summing up the linkage matrix.
        Indeed, a condensed distance matrix is a flat array containing the
        upper triangular of the distance matrix. This is the form that ``pdist``
        returns. Alternatively, a collection of :math:`m` observation vectors
        in :math:`n` dimensions may be passed as  an :math:`m` by :math:`n`
        array. All elements of the condensed distance matrix must be finite,
        i.e., no NaNs or infs.
        Alternatively, we could used the ``squareform`` distance matrix to yield
        different distance values than expected.
        the ``design`` approach uses the complete inpout example matrix  also
        called 'design matrix' to lead correct linkage matrix similar to
        `squareform` and `condense``.

    metric : str or callable, default is {'euclidean'}
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string, it must be one of the options
        allowed by :func:`sklearn.metrics.pairwise.pairwise_distances`.
        If ``X`` is the distance array itself, use "precomputed" as the metric.
        Precomputed distance matrices must have 0 along the diagonal.

    method : str, optional, default is {'complete'}
        The linkage algorithm to use. See the ``Linkage Methods`` section below
        for full descriptions.

    optimal_ordering : bool, optional
        If True, the linkage matrix will be reordered so that the distance
        between successive leaves is minimal. This results in a more intuitive
        tree structure when the data are visualized. defaults to False, because
        this algorithm can be slow, particularly on large datasets. See
        also :func:`scipy.cluster.hierarchy.linkage`.


    Returns
    --------
    row_clusters: linkage matrix
        consist of several rows where each rw represents one merge. The first
        and second columns denotes the most dissimilar members of each cluster
        and the third columns reports the distance between those members


    Linkage Methods
    -----------------
    The following are methods for calculating the distance between the
    newly formed cluster :math:`u` and each :math:`v`.

    * method='single' assigns

      .. math::
         d(u,v) = \min(dist(u[i],v[j]))

      for all points :math:`i` in cluster :math:`u` and
      :math:`j` in cluster :math:`v`. This is also known as the
      Nearest Point Algorithm.

    * method='complete' assigns

      .. math::
         d(u, v) = \max(dist(u[i],v[j]))

      for all points :math:`i` in cluster u and :math:`j` in
      cluster :math:`v`. This is also known by the Farthest Point
      Algorithm or Voor Hees Algorithm.

    * method='average' assigns

      .. math::
         d(u,v) = \sum_{ij} \\frac{d(u[i], v[j])}{(|u|*|v|)}

      for all points :math:`i` and :math:`j` where :math:`|u|`
      and :math:`|v|` are the cardinalities of clusters :math:`u`
      and :math:`v`, respectively. This is also called the UPGMA
      algorithm.

    * method='weighted' assigns

      .. math::
         d(u,v) = (dist(s,v) + dist(t,v))/2

      where cluster u was formed with cluster s and t and v
      is a remaining cluster in the forest (also called WPGMA).

    * method='centroid' assigns

      .. math::
         dist(s,t) = ||c_s-c_t||_2

      where :math:`c_s` and :math:`c_t` are the centroids of
      clusters :math:`s` and :math:`t`, respectively. When two
      clusters :math:`s` and :math:`t` are combined into a new
      cluster :math:`u`, the new centroid is computed over all the
      original objects in clusters :math:`s` and :math:`t`. The
      distance then becomes the Euclidean distance between the
      centroid of :math:`u` and the centroid of a remaining cluster
      :math:`v` in the forest. This is also known as the UPGMC
      algorithm.

    * method='median' assigns :math:`d(s,t)` like the ``centroid``
      method. When two clusters :math:`s` and :math:`t` are combined
      into a new cluster :math:`u`, the average of centroids s and t
      give the new centroid :math:`u`. This is also known as the
      WPGMC algorithm.

    * method='ward' uses the Ward variance minimization algorithm.
      The new entry :math:`d(u,v)` is computed as follows,

      .. math::

         d(u,v) = \sqrt{\frac{|v|+|s|}{T}d(v,s)^2 \\
                      + \frac{|v|+|t|}{T}d(v,t)^2 \\
                      - \frac{|v|}{T}d(s,t)^2}

      where :math:`u` is the newly joined cluster consisting of
      clusters :math:`s` and :math:`t`, :math:`v` is an unused
      cluster in the forest, :math:`T=|v|+|s|+|t|`, and
      :math:`|*|` is the cardinality of its argument. This is also
      known as the incremental algorithm.

    Warning: When the minimum distance pair in the forest is chosen, there
    may be two or more pairs with the same minimum distance. This
    implementation may choose a different minimum than the MATLAB
    version.

    See Also
    --------
    scipy.spatial.distance.pdist : pairwise distance metrics

    References
    ----------
    .. [1] Daniel Mullner, "Modern hierarchical, agglomerative clustering
           algorithms", :arXiv:`1109.2378v1`.
    .. [2] Ziv Bar-Joseph, David K. Gifford, Tommi S. Jaakkola, "Fast optimal
           leaf ordering for hierarchical clustering", 2001. Bioinformatics
           :doi:`10.1093/bioinformatics/17.suppl_1.S22`

    """
    df = _assert_all_types(df, pd.DataFrame, np.ndarray)

    if columns is not None:
        if isinstance (columns , str):
            columns = [columns]
        if len(columns)!= df.shape [1]:
            raise TypeError("Number of columns must fit the shape of X."
                            f" got {len(columns)} instead of {df.shape [1]}"
                            )
        df = pd.DataFrame(data = df.values if hasattr(df, 'columns') else df ,
                          columns = columns )

    kind= str(kind).lower().strip()
    if kind not in ('squareform', 'condense', 'design'):
        raise ValueError(f"Unknown method {method!r}. Expect 'squareform',"
                         " 'condense' or 'design'.")

    labels = [f'ID_{i}' for i in range(len(df))]
    if kind =='squareform':
        row_dist = pd.DataFrame (squareform (
        pdist(df, metric= metric )), columns = labels  ,
        index = labels)
        row_clusters = linkage (row_dist, method =method, metric =metric
                                )
    if kind =='condense':
        row_clusters = linkage (pdist(df, metric =metric), method =method
                                )
    if kind =='design':
        row_clusters = linkage(df.values if hasattr (df, 'columns') else df,
                               method = method,
                               optimal_ordering=optimal_ordering )

    if as_frame:
        row_clusters = pd.DataFrame ( row_clusters,
                                     columns = [ 'row label 1',
                                                'row lable 2',
                                                'distance',
                                                'no. of items in clust.'
                                                ],
                                     index = ['cluster %d' % (i +1) for i in
                                              range(row_clusters.shape[0])
                                              ]
                                     )
    return row_clusters

def d_hanning_window(
        x: ArrayLike[DType[float]],
        xk: float ,
        W: int
        )-> F:
    """ Discrete hanning function.

    For futher details, please refer to https://doi.org/10.1190/1.2400625

    :param x: variable point along the window width
    :param xk: Center of the window `W`. It presumes to host the most weigth.
    :param W: int, window-size; preferably set to odd number. It must be less than
          the dipole length.
    :return: Anonymous function (x,xk, W) value
    """
    # x =check_y (x, input_name ='x')
    return  1/W * (1 + np.cos (
        2 * np.pi * (x-xk) /W)) if np.abs(x-xk) <= W/2 else  0.

def betaj (
        xj: int ,
        L: int ,
        W: int ,
        **kws
 )-> float :
    """ Weight factor function for convoluting at station/site j.

    The function deals with the discrete hanning window based on ideas presented
    in Torres-Verdin and Bostick, 1992, https://doi.org/10.1190/1.2400625.

    :param xj: int, position of the point to compute its weight.
    :param W: int, window size, presumes to be the number of dipole.
    :param L: int : length of dipole in meters
    :param kws: dict , additional :func:`scipy.intergate.quad` functions.

    :return: Weight value at the position `xj`, prefix-`x`is used to specify
        the direction. Commonly the survey direction is considered as `x`.

    :example:
        >>> from watex.exmath import betaj
        >>> # compute the weight point for window-size = 5 at position j =2
        >>> L= 1 ; W=5
        >>> betaj (xj = 2 , L=L, W=W )
        ... 0.35136534572813144
    """
    if W < L :
        raise ValueError("Window-size must be greater than the dipole length.")

    xk = W/2
    # vec_betaj = np.vectorize( betaj ) ; vec_betaj(0, 1, 5)
    return  quad (d_hanning_window, xj - L/2 , xj +L/2, args=( xk, W),
                  **kws)[0]

def rhoa2z (
        rhoa: NDArray[DType[T]],
        phs:ArrayLike,
        freq: ArrayLike
)-> NDArray[DType[T]]:
    r""" Convert apparent resistivity to impendance tensor z

    :param rhoa: Apparent resistivity in :math:`\Omega.m`
    :type rhoa: ndarray, shape (N, M)

    :param phs: Phase in degrees
    :type phs: ndarray, shape (N, M)
    :param freq: Frequency in Hertz
    :type freq: array-like , shape (N, )
    :
    :return: Impendance tensor; Tensor is a complex number in :math:`\Omega`.
    :rtype: ndarray, shape (N, M), dtype = 'complex'

    :example:
    >>> import numpy as np
    >>> rhoa = np.array([1623.73691735])
    >>> phz = np.array([45.])
    >>> f = np.array ([1014])
    >>> rhoa2z(rhoa, phz, f)
    ... array([[2.54950976+2.54950976j]])

    """

    rhoa = np.array(rhoa); freq = np.array(freq) ; phs = np.array(phs)

    if len(phs) != len(rhoa):
        raise ValueError ("Phase and rhoa must have the same length."
                          f" {len(phs)} & {len(rhoa)} are given.")

    if len(freq) != len(rhoa):
        raise ValueError("frequency and rhoa must have the same length."
                         "{len(freq} & {len(rhoa)} are given.")

    omega0 = 2 * np.pi * freq[:, None]
    z= np.sqrt(rhoa * omega0 * mu0 ) * (np.cos (
        np.deg2rad(phs)) + 1j * np.sin(np.deg2rad(phs)))

    return z

def rhophi2z(rho, phi, freq):
    r"""
    Convert impedance-style information given in Rho/Phi format
    into complex valued Z.

    Parameters
    -----------
    rho: ArrayLike 1D/2D
       Resistivity array in :math:`\Omega.m`. If array is two-dimensional,
       it should be 2x2 array (real).

    phi: ArrayLike 1D/2D
       Phase array in degree (:math:`\degree`). If array is two-dimensional,
       it should be 2x2 array (real).

    freq: float, arraylike 1d
      Frequency in Hz

    Returns
    ---------
    Z: Arraylike 1d or 2d , complex

      Z dimension depends to the inputs array `rho` and `phi`.

    Examples
    ---------
    >>> import numpy as np
    >>> from watex.utils.exmath import rhophi2z
    >>> rhophi2z (823 , 25 , 500 )
    array([1300.00682824+606.20313966j])
    >>> rho = np.array ([[823, 700], [723, 526]] )
    >>> phi = np.array ([[45, 50], [90, 180]])
    >>> rhophi2z (rho, phi , freq= 500  )
    array([[ 1.01427314e+03+1.01427314e+03j,  8.50328081e+02+1.01338154e+03j],
           [ 8.23227764e-14+1.34443297e+03j, -1.14673449e+03+1.40434473e-13j]])
    >>> rhophi2z (np.array ( [ 823, 700])  , np.array ([45, 50 ])  , [500, 700] )
    array([1014.27313876+1014.27313876j, 1006.12175325+1199.04921402j])
    >>> rho  = np.abs (np.random.randn (7, 3 ) * 100 )
    >>> phi = np.abs ( np.random.randn (7, 3 ) *180 % 90 )
    >>> freq = np.abs ( np.random.randn (7) * 100 )
    >>> rhophi2z (rho   , phi  , freq )

    """
    def _rhophi2z (r, p, f ):
        """ An isolated part of `rhophi2z """
        with warnings.catch_warnings():
            warnings.filterwarnings(action='ignore', category=RuntimeWarning)
            abs_z  = np.sqrt(5 * f * r)
        # `f` may arrive as a 1-element array: extract the scalar
        # explicitly (implicit conversion is removed in NumPy >= 2.3)
        return cmath.rect(np.asarray(abs_z, dtype=float).ravel()[0],
                          radians(p))

    is_array2x2 =False

    rho = np.array (
        is_iterable(rho, exclude_string= True ,
                    transform =True ))
    phi = np.array (
        is_iterable(phi, exclude_string= True ,
                    transform =True ))
    freq = np.array (
        is_iterable(freq, exclude_string= True ,
                    transform =True ))

    if ( rho.shape == (2,2) or  phi.shape == (2,2)):
        n=None
        if rho.shape != (2,2):
            n, t  ='Resistivity', rho
        elif phi.shape != (2,2):
            n , t ='Phase', phi
        if n is not None:
            raise ResistivityError ("Resistivity and Phase must be consistent."
                           f" Expect 2 x2 array for {n}. Got {t.shape}")

        is_array2x2 = True
    if not ( _is_numeric_dtype(rho) and _is_numeric_dtype(phi)):
        raise ResistivityError ('Resistivity and Phase arguments must be one (1D) or'
                       ' two dimensional (2x2) arrays (real)')

    if is_array2x2 :
        z = np.zeros((2,2),'complex')
        for i in range(2):
            for j in range(2):
                z[i, j ] = _rhophi2z(r = rho[i,j], p = phi[i,j], f = freq )
                # abs_z  = np.sqrt(5 * freq * rho[i,j])
                # z[i,j] = cmath.rect(abs_z , radians(phi[i,j]))
        return z

    check_consistency_size(rho, phi, freq )

    if _is_arraylike_1d (phi ):

        z = np.zeros_like ( phi , dtype ='complex')
        # when scalar is passed or 1d array is
        # given
        for ii in range ( len(phi)): #
            z [ii] = _rhophi2z ( rho[ii], phi[ii], freq[:, None ][ii] )
    else:
        # when non square matrix is given
        # range like freq and n_stations

        z = np.zeros(( len( freq), phi.shape [1]), dtype = 'complex')
        for i in range (len(freq)):
            for j in range(phi.shape[1]) :
                z[i, j ] =  _rhophi2z(rho[i, j], phi[i,j], freq[i] )

    return z

def z2rhoa (
        z:NDArray [DType[complex]],
        freq: ArrayLike[DType[float]]
)-> NDArray[DType[float]]:
    r""" Convert impendance tensor z  to apparent resistivity

    :param z: Impedance tensor  in :math:`\Omega`
    :type z: ndarray, shape (N, M)

    :param freq: Frequency in Hertz
    :type freq: array-like , shape (N, )
    :
    :return: Apparent resistivity in :math:`\Omega.m`
    :rtype: ndarray, shape (N, M)

    :example:
    >>> import numpy as np
    >>> z = np.array([2 + 1j *3 ])
    >>> f = np.array ([1014])
    >>> z2rhoa(z, f)
    ... array([[1623.73691735]])

    """

    z = np.array(z, dtype = 'complex' ) ;  freq = np.array(freq)

    if len(freq) != len(z):
        raise ValueError("frequency and tensor z must have the same length."
                         f"{len(freq)} & {len(z)} are given.")

    return np.abs(z)**2 / (2 * np.pi * freq[:, None] * mu0 )

def savitzky_golay1d (
        y: ArrayLike[DType[T]],
        window_size:int ,
        order: int,
        deriv: int =0,
        rate: int =1,
        mode: str ='same'
        )-> ArrayLike[DType[T]]:
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.

    The Savitzky-Golay filter removes high frequency noise from data. It has the
    advantage of preserving the original shape and features of the signal better
    than other types of filtering approaches, such as moving averages techniques.

    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    mode: str
         mode of the border prepending. Should be ``valid`` or ``same``.
         ``same`` is used for prepending or appending the first value of
         array for smoothing.Default is ``same``.
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly suited for
    smoothing noisy data. The main idea behind this approach is to make for
    each point a least-square fit with a polynomial of high order over a
    odd-sized window centered at the point.

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from watex.utils.exmath import savitzky_golay1d
    >>> t = np.linspace(-4, 4, 500)
    >>> y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    >>> ysg = savitzky_golay1d(y, window_size=31, order=4)
    >>> plt.plot(t, y, label='Noisy signal')
    >>> plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    >>> plt.plot(t, ysg, 'r', label='Filtered signal')
    >>> plt.legend()
    >>> plt.show()

    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    .. [3] https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter#Moving_average

    """

    try:
        window_size = abs(int(window_size))
        order = abs(int(order))
    except (ValueError, TypeError):
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)

    y = np.asarray(y)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.array(
        [[k**i for i in order_range]
         for k in range(-half_window, half_window+1)],
        dtype=float,
    )
    m = np.linalg.pinv(b)[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode=mode)
