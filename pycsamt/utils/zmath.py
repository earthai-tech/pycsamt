# -*- coding: utf-8 -*-
# Author: L. Kouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
helper functions for standard calculations.
"""

from __future__ import annotations 

import numpy as np
from numpy.typing import ArrayLike 
import math
import cmath

# Exception classes
class ZMathError(Exception): pass
class EmptyArrayError(ZMathError): pass
class InvalidSignificantFiguresError(ZMathError): pass
class PeriodListError(ZMathError): pass
class NearestIndexError(ZMathError): pass
class LogArrayError(ZMathError): pass
class MatrixInversionError(ZMathError): pass
class ImpedanceConversionError(ZMathError): pass
class DeterminantError(ZMathError): pass
class PolarToRectError(ZMathError): pass
class RectToPolarError(ZMathError): pass
class ZError2RPhiError(ZMathError): pass
class OldZError2RPhiError(ZMathError): pass
class RotationMatrixError(ZMathError): pass
class VectorRotationError(ZMathError): pass
class MatrixMultiplyError(ZMathError): pass
class ReorientDataError(ZMathError): pass


# Define constants
epsilon = 1e-9  # tolerance for small values
mu0 = 4e-7 * math.pi  # magnetic permeability in free space
_RAD_THR   = np.deg2rad(5.0)           # heuristics for auto-unit check
_DEG_SCALE = 180.0 / np.pi

__all__ = [
    'centre_point',
    'roundsf',
    'get_period_list',
    'nearest_index',
    'make_log_increasing_array',
    'invertmatrix_incl_errors',
    'rhophi2z',
    'compute_determinant_error',
    'propagate_error_polar2rect',
    'propagate_error_rect2polar',
    'z_error2r_phi_error',
    'old_z_error2r_phi_error',
    'rotatematrix_incl_errors',
    'rotatevector_incl_errors',
    'multiplymatrices_incl_errors',
    'reorient_data2D', 
    "rhophi_to_z",          # new name (was rhophi2z)
    "z_err_to_rphi_err",    # new name (was z_error2r_phi_error)
]

def centre_point(xarray, yarray):
    """
    Get center point of arrays of x and y values.
    """
    # Validate inputs
    if xarray.size == 0 or yarray.size == 0:
        raise EmptyArrayError(
            "Input arrays must be non-empty."
        )
    x_min = xarray.min()
    x_max = xarray.max()
    y_min = yarray.min()
    y_max = yarray.max()
    # Compute mid-points
    return (
        (x_max + x_min) / 2.0,
        (y_max + y_min) / 2.0
    )


def roundsf(number, sf):
    """
    Round a number to given significant figures.
    """
    # Validate significant figures
    try:
        sf_int = int(sf)
    except (ValueError, TypeError):
        raise InvalidSignificantFiguresError(
            f"sf must be int, got {sf!r}"
        )
    sf_int = max(sf_int, 1)
    # Handle zero
    if number == 0:
        return 0.0
    # Compute decimal places for rounding
    try:
        rounding = int(
            np.ceil(
                -np.log10(abs(number))
                + sf_int - 1
            )
        )
    except Exception as e:
        raise ZMathError(
            f"Error computing rounding: {e}"
        )
    return float(np.round(number, rounding))


def get_period_list(
    period_min,
    period_max,
    periods_per_decade,
    include_outside_range=True
):
    """
    Generate log-spaced list of periods.
    """
    # Validate bounds
    if period_min <= 0 or period_max <= 0:
        raise PeriodListError(
            "period_min and period_max must be > 0."
        )
    if period_min >= period_max:
        raise PeriodListError(
            "period_min must be < period_max."
        )
    # Compute log10 bounds
    try:
        log_min = math.log10(period_min)
        log_max = math.log10(period_max)
    except ValueError as e:
        raise PeriodListError(
            f"Invalid period values: {e}"
        )
    # Determine start log value
    if not log_min.is_integer():
        aligned_min = np.linspace(
            math.floor(log_min),
            math.ceil(log_min),
            periods_per_decade + 1
        )
        diffs_min = log_min - aligned_min
        idx = (
            np.where(diffs_min > 0)[0][-1]
            if include_outside_range
            else np.where(diffs_min < 0)[0][0]
        )
        start_log = aligned_min[idx]
    else:
        start_log = log_min
    # Determine stop log value
    if not log_max.is_integer():
        aligned_max = np.linspace(
            math.floor(log_max),
            math.ceil(log_max),
            periods_per_decade + 1
        )
        diffs_max = log_max - aligned_max
        idx = (
            np.where(diffs_max < 0)[0][0]
            if include_outside_range
            else np.where(diffs_max > 0)[0][-1]
        )
        stop_log = aligned_max[idx]
    else:
        stop_log = log_max
    # Number of periods to generate
    num = int(
        (stop_log - start_log)
        * periods_per_decade
        + 1
    )
    return np.logspace(
        start_log,
        stop_log,
        num
    )

def nearest_index(val, array):
    """
    Find index of nearest value in array.

    :param val: value to search for
    :param array: numpy array to search in
    :return: integer index of nearest element
    """
    arr = np.asarray(array)
    if arr.size == 0:
        raise NearestIndexError(
            "Input array must be non-empty"
        )
    # compute absolute differences
    diffs = np.abs(arr - val)
    # use argmin for nearest
    idx = int(diffs.argmin())
    return idx


def make_log_increasing_array(
    z1_layer,
    target_depth,
    n_layers,
    increment_factor=0.999
):
    """
    Create depth array with log-increasing cells.
    """
    # validate inputs
    if z1_layer <= 0 or target_depth <= 0:
        raise LogArrayError(
            "z1_layer and target_depth must be > 0"
        )
    if n_layers < 1 or not isinstance(
        n_layers, int
    ):
        raise LogArrayError(
            f"n_layers must be int >=1, got {n_layers}"
        )
    if not 0 < increment_factor < 1:
        raise LogArrayError(
            "increment_factor must be in (0,1)"
        )
    max_thick = target_depth
    attempt = 0
    max_attempts = int(1e6)
    while True:
        log_z = np.logspace(
            math.log10(z1_layer),
            math.log10(max_thick),
            num=n_layers
        )
        if log_z.sum() <= target_depth:
            return log_z
        max_thick *= increment_factor
        attempt += 1
        if attempt >= max_attempts:
            raise LogArrayError(
                "Unable to converge log-increasing array"
            )


def invertmatrix_incl_errors(
    inmatrix,
    inmatrix_err=None
):
    """
    Invert square matrix, propagating element errors.
    """
    if inmatrix is None:
        raise MatrixInversionError(
            "Input matrix must be provided"
        )
    mat = np.asarray(inmatrix)
    if mat.ndim != 2:
        raise MatrixInversionError(
            "Input must be 2D square matrix"
        )
    n, m = mat.shape
    if n != m:
        raise MatrixInversionError(
            "Matrix must be square"
        )
    # check error matrix shape
    err_mat = None
    if inmatrix_err is not None:
        err = np.asarray(inmatrix_err)
        if err.shape != mat.shape:
            raise MatrixInversionError(
                "Error matrix shape differs"
            )
        err_mat = err
    # determinant
    det = np.linalg.det(mat)
    if abs(det) < epsilon:
        raise MatrixInversionError(
            "Matrix is singular or nearly singular"
        )
    inv_mat = np.linalg.inv(mat)
    inv_err = None
    if err_mat is not None:
        inv_err = np.zeros_like(err_mat)
        # error propagation for 2x2 only
        if n != 2:
            raise MatrixInversionError(
                "Only 2x2 matrices supported for error"
            )
        for i in range(2):
            for j in range(2):
                total = 0.0
                for k in range(2):
                    for l in range(2):
                        total += abs(
                            -inv_mat[i, k]
                            * inv_mat[l, j]
                            * err_mat[k, l]
                        )
                inv_err[i, j] = total
    return inv_mat, inv_err


def rhophi2z(rho, phi, freq):
    """
    Convert Rho/Phi into complex Z (2x2 arrays).
    """
    rho_arr = np.asarray(rho)
    phi_arr = np.asarray(phi)
    if rho_arr.shape != (2, 2) or phi_arr.shape != (2, 2):
        raise ImpedanceConversionError(
            "rho and phi must be 2x2 arrays"
        )
    if not np.issubdtype(rho_arr.dtype, np.number) or \
       not np.issubdtype(phi_arr.dtype, np.number):
        raise ImpedanceConversionError(
            "rho and phi entries must be numeric"
        )
    Z = np.zeros((2, 2), dtype=complex)
    for i in range(2):
        for j in range(2):
            abs_z = math.sqrt(
                5.0 * freq * rho_arr[i, j]
            )
            angle = math.radians(phi_arr[i, j])
            Z[i, j] = cmath.rect(abs_z, angle)
    return Z


def compute_determinant_error(
    z_array,
    z_err_array,
    method='theoretical',
    repeats=1000
):
    """
    Compute error of det(Z) by theoretical or stochastic.
    Returns error and sqrt(error).
    """
    Z = np.asarray(z_array)
    err = np.asarray(z_err_array)
    if Z.shape != err.shape:
        raise DeterminantError(
            "z_array and z_err_array must match shape"
        )
    if method.lower().startswith('stoch'):
        sims = np.empty(
            (repeats, *Z.shape),
            dtype=complex
        )
        for r in range(repeats):
            noise = np.random.normal(
                loc=0.0,
                scale=err,
                size=Z.shape
            )
            sims[r] = Z + noise * (1 + 1j)
        dets = np.linalg.det(sims)
        error = np.std(dets, axis=0)
    else:
        # theoretical error for 2x2
        a, b = Z[0, 0], Z[0, 1]
        c, d = Z[1, 0], Z[1, 1]
        ea, eb = err[0, 0], err[0, 1]
        ec, ed = err[1, 0], err[1, 1]
        error = np.abs(
            ea * abs(d)
            + ed * abs(a)
            - eb * abs(c)
            - ec * abs(b)
        )
    return error, np.sqrt(error)


def propagate_error_polar2rect(
    r, r_error,
    phi, phi_error
):
    """
    Propagate errors from polar to cartesian coords.

    Returns x_err, y_err.
    """
    try:
        # compute corners of annulus sector
        corners = []
        for dr in (-r_error, r_error):
            for dphi in (-phi_error, phi_error):
                zc = cmath.rect(
                    r + dr,
                    math.radians(phi + dphi)
                )
                corners.append((zc.real, zc.imag))
        # outer boundary point
        zc = cmath.rect(
            r + r_error,
            math.radians(phi)
        )
        corners.append((zc.real, zc.imag))
    except Exception as e:
        raise PolarToRectError(
            f"Error computing corners: {e}"
        )
    xs, ys = zip(*corners)
    # nominal point
    z0 = cmath.rect(
        r,
        math.radians(phi)
    )
    x0, y0 = z0.real, z0.imag
    x_err = max(abs(x0 - xi) for xi in xs)
    y_err = max(abs(y0 - yi) for yi in ys)
    return x_err, y_err


def propagate_error_rect2polar(
    x, x_error,
    y, y_error
):
    """
    Propagate errors from cartesian to polar coords.

    Returns rho_err, phi_err.
    """
    # gather box corner points
    pts = [
        (x + x_error, y),
        (x - x_error, y),
        (x, y + y_error),
        (x, y - y_error),
        (x + x_error, y + y_error),
        (x + x_error, y - y_error),
        (x - x_error, y + y_error),
        (x - x_error, y - y_error)
    ]
    # check origin in box
    origin_in = (
        abs(x) <= x_error
        and abs(y) <= y_error
    )
    try:
        polars = [cmath.polar(complex(x0, y0)) for x0, y0 in pts]
    except Exception as e:
        raise RectToPolarError(
            f"Error computing polar values: {e}"
        )
    rhos = [p[0] for p in polars]
    phis = [math.degrees(p[1]) % 360 for p in polars]
    rho_err = 0.5 * (max(rhos) - min(rhos))
    phi_err = 0.5 * (max(phis) - min(phis))
    # handle wrap-around
    if max(phis) > 270 and min(phis) < 90:
        tmp1 = [p for p in phis if p < 90]
        tmp2 = [p for p in phis if p > 270]
        phi_err = 0.5 * ((max(tmp1) - min(tmp2)) % 360)
    if phi_err > 180:
        phi_err = (-phi_err) % 360
    if origin_in:
        rho_err = max(rhos)
        phi_err = 180.0
    return rho_err, phi_err


def z_error2r_phi_error(
    z_real, z_imag, error
):
    """
    Error from rectangular to polar coords.

    Returns resistivity rel err, phase abs err.
    """
    amp = np.abs(z_real + 1j * z_imag)
    with np.errstate(all='ignore'):
        rel_z_err = error / amp
    res_rel = 2.0 * rel_z_err
    if np.iterable(z_real):
        phi_err = np.degrees(np.arctan(rel_z_err))
        phi_err[res_rel > 1.0] = 90.0
    else:
        phi_err = (
            90.0 if res_rel > 1.0
            else np.degrees(np.arctan(rel_z_err))
        )
    return res_rel, phi_err


def old_z_error2r_phi_error(
    x, x_error,
    y, y_error
):
    """
    Legacy rect->polar error, with MT conventions.
    """
    pts = [
        (x + x_error, y),
        (x - x_error, y),
        (x, y + y_error),
        (x, y - y_error),
        (x - x_error, y - y_error),
        (x + x_error, y - y_error),
        (x + x_error, y + y_error),
        (x - x_error, y + y_error)
    ]
    # XXX TOD REVISE:
    origin_in = ( # noqa 
        abs(x) <= x_error
        and abs(y) <= y_error
    )
    polars = [cmath.polar(complex(u, v)) for u, v in pts]
    rhos = [p[0] for p in polars]
    phis = [math.degrees(p[1]) % 360 for p in polars] # noqa
    rho_err = 0.5 * (max(rhos) - min(rhos))
    rho = cmath.polar(complex(x, y))[0]
    rel_err = rho_err / rho if rho != 0 else 0.0
    if rel_err > 1.0:
        phi_err = 90.0
    else:
        phi_err = np.degrees(np.arcsin(rel_err))
    return rho_err, phi_err


def rotatematrix_incl_errors(
    inmatrix,
    angle,
    inmatrix_err=None
):
    """
    Rotate 2x2 matrix and propagate errors.
    Returns (rotated_matrix, err_matrix).
    """
    if inmatrix is None:
        raise RotationMatrixError(
            "Matrix must be provided"
        )
    M = np.asarray(inmatrix)
    if M.shape != (2, 2):
        raise RotationMatrixError(
            "inmatrix must be 2x2"
        )
    errM = None
    if inmatrix_err is not None:
        EM = np.asarray(inmatrix_err)
        if EM.shape != M.shape:
            raise RotationMatrixError(
                "Error matrix shape differs"
            )
        errM = np.real(EM)
    try:
        deg = angle % 360
    except Exception:
        raise RotationMatrixError(
            "Angle must be numeric degrees"
        )
    phi = math.radians(deg)
    c, s = math.cos(phi), math.sin(phi)
    R = np.array([
        [ c,  s],
        [-s,  c]
    ])
    Rt = np.linalg.inv(R)
    Mrot = R.dot(M).dot(Rt)
    ErrRot = None
    if errM is not None:
        e = errM
        ErrRot = np.zeros_like(e)
        ErrRot[0, 0] = math.sqrt(
            (c**2 * e[0,0])**2 +
            (c*s  * e[0,1])**2 +
            (c*s  * e[1,0])**2 +
            (s**2 * e[1,1])**2
        )
        ErrRot[0, 1] = math.sqrt(
            (c**2 * e[0,1])**2 +
            (c*s  * e[1,1])**2 +
            (c*s  * e[0,0])**2 +
            (s**2 * e[1,0])**2
        )
        ErrRot[1, 0] = math.sqrt(
            (c**2 * e[1,0])**2 +
            (c*s  * e[1,1])**2 +
            (c*s  * e[0,0])**2 +
            (s**2 * e[0,1])**2
        )
        ErrRot[1, 1] = math.sqrt(
            (c**2 * e[1,1])**2 +
            (c*s  * e[0,1])**2 +
            (c*s  * e[1,0])**2 +
            (s**2 * e[0,0])**2
        )
    return Mrot, ErrRot

def rotatevector_incl_errors(
    invector, angle, invector_err=None
):
    """
    Rotate 2D vector and propagate errors.
    Returns (vrot, verr).
    """
    if invector is None:
        raise VectorRotationError(
            "Vector must be provided"
        )
    v = np.asarray(invector)
    e = None
    if invector_err is not None:
        e = np.asarray(invector_err)
        if v.shape != e.shape:
            raise VectorRotationError(
                "Vector/error shape differs"
            )
    try:
        deg = angle % 360
    except Exception:
        raise VectorRotationError(
            "Angle must be numeric degrees"
        )
    phi = math.radians(deg)
    c, s = math.cos(phi), math.sin(phi)
    R = np.array([[c, s], [-s, c]])
    Rt = np.linalg.inv(R)
    if v.shape == (1,2):
        vrot = v.dot(Rt)
    else:
        vrot = R.dot(v)
    verr = None
    if e is not None:
        if e.shape == (1,2):
            verr = e.dot(np.abs(Rt))
        else:
            verr = np.abs(R).dot(e)
    return vrot, verr


def multiplymatrices_incl_errors(
    inmatrix1, inmatrix2,
    inmatrix1_err=None,
    inmatrix2_err=None
):
    """
    Multiply two 2x2 matrices and propagate
    element errors. Returns (product, error_matrix).
    """
    if inmatrix1 is None or inmatrix2 is None:
        raise MatrixMultiplyError(
            "Two matrices required"
        )
    A = np.asarray(inmatrix1)
    B = np.asarray(inmatrix2)
    if A.shape != (2,2) or B.shape != (2,2):
        raise MatrixMultiplyError(
            "Matrices must be 2x2"
        )
    prod = A.dot(B)
    if inmatrix1_err is None or inmatrix2_err is None:
        return prod, None
    EA = np.asarray(inmatrix1_err)
    EB = np.asarray(inmatrix2_err)
    if EA.shape != A.shape or EB.shape != B.shape:
        raise MatrixMultiplyError(
            "Error matrices shape differ"
        )
    var = np.zeros((2,2))
    var[0,0] = (
        (EA[0,0] * B[0,0])**2
        + (EA[0,1] * B[1,0])**2
        + (EB[0,0] * A[0,0])**2
        + (EB[1,0] * A[0,1])**2
    )
    var[0,1] = (
        (EA[0,0] * B[0,1])**2
        + (EA[0,1] * B[1,1])**2
        + (EB[0,1] * A[0,0])**2
        + (EB[1,1] * A[0,1])**2
    )
    var[1,0] = (
        (EA[1,0] * B[0,0])**2
        + (EA[1,1] * B[1,0])**2
        + (EB[0,0] * A[1,0])**2
        + (EB[1,0] * A[1,1])**2
    )
    var[1,1] = (
        (EA[1,0] * B[0,1])**2
        + (EA[1,1] * B[1,1])**2
        + (EB[0,1] * A[1,0])**2
        + (EB[1,1] * A[1,1])**2
    )
    return prod, np.sqrt(var)


def reorient_data2D(
    x_values, y_values,
    x_sensor_angle=0, y_sensor_angle=90
):
    """
    Re-orient time series data of a sensor
    pair to default (North/East) axes.

    Inputs:
      x_values, y_values: equal-length 1D arrays
      x_sensor_angle, y_sensor_angle: angles in degrees

    Returns:
      x_corr (North), y_corr (East)
    """
    x = np.asarray(x_values)
    y = np.asarray(y_values)
    # trim to common length
    if x.shape != y.shape:
        l = min(x.size, y.size)
        x = x[:l]
        y = y[:l]
    if x.ndim != 1:
        raise ReorientDataError(
            "Input values must be 1D arrays"
        )
    # stack into Nx2 array
    data = np.stack((x, y), axis=1)
    # compute sensor rotation matrix
    try:
        ax = math.radians(x_sensor_angle)
        ay = math.radians(y_sensor_angle)
    except Exception:
        raise ReorientDataError(
            "Angles must be numeric degrees"
        )
    R = np.array([
        [math.cos(ax),  math.sin(ax)],
        [math.cos(ay),  math.sin(ay)]
    ])
    # invert to default axes
    try:
        Rinv = np.linalg.inv(R)
    except Exception:
        raise ReorientDataError(
            "Sensor angles must define independent axes"
        )
    out = data.dot(Rinv.T)
    # split back to components
    return out[:,0], out[:,1]


def _as_array(x: ArrayLike, name: str) -> np.ndarray:
    """Convert *x* to 1-D float array with helpful error messages."""
    try:
        arr = np.asanyarray(x, dtype=float).ravel()
        if arr.size == 0:
            raise EmptyArrayError(f"{name!r} cannot be empty")
        return arr
    except ValueError as exc:                                  # pragma: no cover
        raise ZMathError(f"{name}: cannot convert to float") from exc

def rhophi_to_z(                       
        phase:       ArrayLike,
        freq:        ArrayLike,
        *,
        resistivity: ArrayLike | None = None,
        z_abs:       ArrayLike | None = None,
        deg:         bool             = False,
    ) -> tuple[np.ndarray, np.ndarray,
               np.ndarray, np.ndarray]:
    """
    Convert *(ρₐ, φ)* → complex **Z** or fill in missing quantity.

    Either *resistivity* **or** *z_abs* must be supplied.

    Parameters
    ----------
    phase, freq
        Per-sample phase (rad or deg, see *deg*) and transmit frequency
        (Hz).  Broadcastable to a common shape.
    resistivity
        Apparent resistivity (Ω·m).  Ignored if *z_abs* is given.
    z_abs
        |Z| (Volt ⋅ m⁻¹ ⋅ A⁻¹).  Overrides the *resistivity* path.
    deg
        ``True`` → *phase* is in **degrees** – it will be converted
        internally.

    Returns
    -------
    z_abs, z_real, z_imag, z_complex :
        Arrays with the broadcast common shape.

    Raises
    ------
    ImpedanceConversionError
        On shape mismatches or missing inputs.
    """
    ph  = _as_array(phase, "phase")
    frq = _as_array(freq,  "freq")

    if deg:                               # user forced degree input
        ph = np.deg2rad(ph)
    elif np.nanmax(np.abs(ph)) > np.pi and np.nanmin(
            np.abs(ph)) > _RAD_THR:       # heuristics → looks like deg
        ph = np.deg2rad(ph)

    # broadcast to a single shape
    try:
        ph, frq = np.broadcast_arrays(ph, frq)
    except ValueError:                                       
        raise ImpedanceConversionError("phase / freq not broadcastable")

    if z_abs is None:
        if resistivity is None:
            raise ImpedanceConversionError(
                "provide either `resistivity` or `z_abs`")
        rho = _as_array(resistivity, "resistivity")
        ph, rho = np.broadcast_arrays(ph, rho)
        z_abs = np.sqrt(mu0 * 2 * np.pi * frq * rho)          # ρ → |Z|
    else:
        z_abs = _as_array(z_abs, "z_abs")
        ph, z_abs = np.broadcast_arrays(ph, z_abs)

    z_real = z_abs * np.cos(ph)
    z_imag = z_abs * np.sin(ph)
    z_cplx = z_real + 1j * z_imag
    return z_abs, z_real, z_imag, z_cplx


def z_err_to_rphi_err(                
    z_real:  ArrayLike,
    z_imag:  ArrayLike,
    z_err:   ArrayLike,
    *,
    deg_out: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Propagate |Z| error to ρₐ *and* phase uncertainty.

    Parameters
    ----------
    z_real, z_imag
        Real / imaginary parts of the complex impedance.
    z_err
        Absolute error on |Z| (same units as |Z|).
    deg_out
        ``True`` → phase error returned in **degrees**
        (otherwise milliradians).

    Returns
    -------
    rho_rel_err, phase_err :
        *rho_rel_err* is the **relative** ρₐ error (unitless).  
        *phase_err* is the absolute phase uncertainty (deg or mrad).

    Notes
    -----
    When the amplitude relative error ≥ 1 the phase uncertainty is
    capped at 90 ° ( ≈ 1570 mrad ).
    """
    xr = _as_array(z_real, "z_real")
    xi = _as_array(z_imag, "z_imag")
    dz = _as_array(z_err,  "z_err")

    try:
        xr, xi, dz = np.broadcast_arrays(xr, xi, dz)
    except ValueError:                                        # pragma: no cover
        raise ImpedanceConversionError("z_real/imag/err shapes mismatch")

    z_amp        = np.hypot(xr, xi)
    z_rel_err    = dz / np.maximum(z_amp, epsilon)
    rho_rel_err  = 2.0 * z_rel_err

    # angle between |Z| vector and the tangent of its error circle
    phase_err = np.arctan(z_rel_err)
    phase_err[rho_rel_err >= 1.0] = np.pi / 2            # 90 deg cap
    if deg_out:
        phase_err *= _DEG_SCALE
    else:                                                # keep in mrad
        phase_err *= 1e3

    return rho_rel_err, phase_err


