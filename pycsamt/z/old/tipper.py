# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
Tipper (induction vector) utilities.

This module defines :class:`Tipper`, a light-weight container for the
complex tipper components :math:`(T_x, T_y)` per frequency, together
with convenient computations for amplitude/phase and the Parkinson
induction arrow magnitude and direction.

The internal storage uses shape ``(n_freq, 1, 2)`` where the last axis
is ordered as ``[Tx, Ty]``. Several common input shapes are accepted
and normalized at construction.

Notes
-----
- Errors, when provided, are **absolute** standard deviations with the
  same shape as the tipper array.
- Magnitudes and directions follow the Parkinson convention: the
  induction arrow points **towards** conductors, hence the minus sign
  in the direction computation.

Examples
--------
Create from vectors of Tx and Ty and compute amplitude/phase:

>>> import numpy as np
>>> from pycsamt.z.tipper import Tipper
>>> Tx = np.array([0.1+0.05j, 0.2+0.0j])
>>> Ty = np.array([-0.05+0.02j, 0.0+0.1j])
>>> T = np.stack([Tx, Ty], axis=-1)          # shape (n, 2)
>>> tip = Tipper(tipper_array=T, freq=[10., 1.])
>>> tip.amplitude.shape, tip.phase.shape
((2, 1, 2), (2, 1, 2))
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

from ..log.logger import get_logger
from ..exceptions import ZError  # reuse generic container error
from ..utils.zmath import (
    propagate_error_rect2polar,
    rotatevector_incl_errors,
)
from .base import BaseEM
logger = get_logger(__name__)


class Tipper(BaseEM):
    """
    Complex tipper container with basic derived quantities.

    Parameters
    ----------
    tipper_array : array-like, optional
        Complex tipper values. Accepted shapes:

        - ``(n_freq, 1, 2)`` (preferred)
        - ``(n_freq, 2)`` or ``(2,)`` (promoted)
        - ``(1, 2)`` (single-frequency, promoted)

        The last axis is ordered as ``[Tx, Ty]``.
    tipper_err_array : array-like, optional
        Absolute errors for the tipper, same shape as the normalized
        tipper array after promotion. If ``None``, uncertainties are
        not propagated.
    freq : array-like of float, optional
        Frequency vector in Hz, length ``n_freq``.

    Attributes
    ----------
    tipper : ndarray or None
        Complex tipper array of shape ``(n_freq, 1, 2)``.
    tipper_err : ndarray or None
        Absolute errors on the tipper with same shape.
    freq : ndarray or None
        Frequency vector (Hz).
    rotation_angle : float or ndarray
        Rotation history (degrees). ``0.0`` if not rotated; when
        multiple frequencies are present, stored as a vector.

    Derived (set by :meth:`compute_amp_phase` and
    :meth:`compute_mag_direction`)
    ---------------------------------------------------------------
    amplitude, phase : ndarray or None
        Amplitude :math:`|T|` and phase (deg) for each component.
    amplitude_err, phase_err : ndarray or None
        Propagated absolute errors (if available).
    mag_real, mag_imag : ndarray or None
        Parkinson arrow magnitudes for real and imaginary parts.
    angle_real, angle_imag : ndarray or None
        Parkinson arrow directions (deg) for real and imaginary parts.
    mag_err, angle_err : ndarray or None
        Crude uncertainty estimates for arrow magnitude and direction.

    Notes
    -----
    Setting :pyattr:`tipper`, :pyattr:`tipper_err` or :pyattr:`freq`
    does **not** automatically recompute derived quantities; call
    :meth:`compute_amp_phase` or :meth:`compute_mag_direction`
    explicitly to refresh them.
    """

    # ------------------------------------------------------------------ #
    # Construction
    # ------------------------------------------------------------------ #
    def __init__(
        self,
        tipper_array: Optional[np.ndarray] = None,
        tipper_err_array: Optional[np.ndarray] = None,
        freq: Optional[Sequence[float]] = None,
    ) -> None:
        self._tipper: Optional[np.ndarray] = None
        self._tipper_err: Optional[np.ndarray] = None
        self._freq: Optional[np.ndarray] = None

        self.rotation_angle: float | np.ndarray = 0.0

        # Derived quantities
        self._amplitude: Optional[np.ndarray] = None
        self._amplitude_err: Optional[np.ndarray] = None
        self._phase: Optional[np.ndarray] = None
        self._phase_err: Optional[np.ndarray] = None

        self._mag_real: Optional[np.ndarray] = None
        self._mag_imag: Optional[np.ndarray] = None
        self._angle_real: Optional[np.ndarray] = None
        self._angle_imag: Optional[np.ndarray] = None
        self._mag_err: Optional[np.ndarray] = None
        self._angle_err: Optional[np.ndarray] = None

        # Normalize and set inputs
        if tipper_array is not None:
            self.tipper = tipper_array
        if tipper_err_array is not None:
            self.tipper_err = tipper_err_array
        if freq is not None:
            self.freq = freq

        # Initialize rotation storage length
        if self._tipper is not None and isinstance(self.rotation_angle, float):
            self.rotation_angle = np.zeros(self._tipper.shape[0], dtype=float)

        # Optional eager computation
        if self._tipper is not None:
            self.compute_amp_phase()
            self.compute_mag_direction()

    # ------------------------------------------------------------------ #
    # Core arrays
    # ------------------------------------------------------------------ #
    @property
    def tipper(self) -> Optional[np.ndarray]:
        """Complex tipper array of shape ``(n_freq, 1, 2)``."""
        return self._tipper

    @tipper.setter
    def tipper(self, arr: np.ndarray) -> None:
        t = np.asarray(arr)

        # Accept common shapes and normalize to (n, 1, 2)
        if t.ndim == 1 and t.shape[0] == 2:
            t = t[None, None, :]
        elif t.ndim == 2:
            if t.shape == (1, 2):
                t = t[None, :, :]
            elif t.shape[1] == 2:
                t = t[:, None, :]
            else:
                raise ZError(
                    "Tipper shape (n, 2) or (1, 2) expected; got "
                    f"{t.shape!r}."
                )
        elif t.ndim == 3:
            if t.shape[1:] != (1, 2):
                raise ZError(
                    "Tipper 3-D shape must be (n_freq, 1, 2); got "
                    f"{t.shape!r}."
                )
        else:
            raise ZError(f"Unsupported tipper shape: {t.shape!r}.")

        self._tipper = t.astype(complex, copy=False)

        # Adjust rotation history storage
        if isinstance(self.rotation_angle, float):
            self.rotation_angle = np.zeros(self._tipper.shape[0], dtype=float)

    @property
    def tipper_err(self) -> Optional[np.ndarray]:
        """Absolute errors on :pyattr:`tipper` or ``None``."""
        return self._tipper_err

    @tipper_err.setter
    def tipper_err(self, arr: Optional[np.ndarray]) -> None:
        if arr is None:
            self._tipper_err = None
            return
        e = np.asarray(arr)
        # Apply the same normalization rules as for tipper
        if e.ndim == 1 and e.shape[0] == 2:
            e = e[None, None, :]
        elif e.ndim == 2:
            if e.shape == (1, 2):
                e = e[None, :, :]
            elif e.shape[1] == 2:
                e = e[:, None, :]
            else:
                raise ZError(
                    "Tipper error shape (n, 2) or (1, 2) expected; got "
                    f"{e.shape!r}."
                )
        elif e.ndim == 3:
            if e.shape[1:] != (1, 2):
                raise ZError(
                    "Tipper error 3-D shape must be (n_freq, 1, 2); got "
                    f"{e.shape!r}."
                )
        else:
            raise ZError(f"Unsupported tipper error shape: {e.shape!r}.")

        if self._tipper is not None and e.shape != self._tipper.shape:
            raise ZError(
                "'tipper_err' must match 'tipper' shape: "
                f"{e.shape!r} vs {self._tipper.shape!r}."
            )
        self._tipper = self._tipper  # no-op, for clarity
        self._tipper_err = e.astype(float, copy=False)

    @property
    def freq(self) -> Optional[np.ndarray]:
        """Frequency vector (Hz) or ``None``."""
        return self._freq

    @freq.setter
    def freq(self, f: Optional[Sequence[float]]) -> None:
        if f is None:
            self._freq = None
            return
        ff = np.asarray(f, dtype=float).ravel()
        if self._tipper is not None and ff.size != self._tipper.shape[0]:
            raise ZError(
                "Length of 'freq' must match tipper stack: "
                f"{ff.size} vs {self._tipper.shape[0]}."
            )
        if np.any(ff <= 0.0):
            raise ZError("Frequencies must be strictly positive.")
        self._freq = ff

    # ------------------------------------------------------------------ #
    # Amplitude and phase
    # ------------------------------------------------------------------ #
    def compute_amp_phase(self) -> None:
        """
        Compute amplitude :math:`|T|` and phase (degrees).

        When errors are available, propagate uncertainties using
        :func:`propagate_error_rect2polar` on each component.

        Sets
        ----
        amplitude, phase : ndarray
            Arrays of shape ``(n_freq, 1, 2)``.
        amplitude_err, phase_err : ndarray or None
            Arrays with propagated errors, or ``None`` if input
            errors are not available.

        Notes
        -----
        If :pyattr:`tipper` is missing, the method returns
        without raising and leaves derived fields unchanged.
        """
        if self._tipper is None:
            return

        T = self._tipper
        E = self._tipper_err

        self._amplitude = np.abs(T)
        self._phase = np.rad2deg(np.angle(T))

        if E is None:
            self._amplitude_err = None
            self._phase_err = None
            return

        amp_err = np.zeros_like(T, dtype=float)
        ph_err = np.zeros_like(T, dtype=float)

        # Loop over frequencies and the two components (Tx, Ty)
        for k in range(T.shape[0]):
            for j in range(2):
                r = T[k, 0, j].real
                im = T[k, 0, j].imag
                dr = E[k, 0, j]
                dim = E[k, 0, j]
                r_err, phi_err = propagate_error_rect2polar(r, dr, im, dim)
                amp_err[k, 0, j] = r_err
                ph_err[k, 0, j] = phi_err

        self._amplitude_err = amp_err
        self._phase_err = ph_err

    # ------------------------------------------------------------------ #
    # Magnitude and direction of induction arrows (Parkinson convention)
    # ------------------------------------------------------------------ #
    def compute_mag_direction(self) -> None:
        """
        Compute Parkinson induction arrow magnitudes and directions.

        The arrow is defined from tipper real/imag parts as:

        .. math::

           M_\\Re = \\sqrt{\\Re(T_x)^2 + \\Re(T_y)^2},\\quad
           \\theta_\\Re = \\operatorname{atan2}(-\\Re(T_y), -\\Re(T_x))

        and similarly for the imaginary parts (``M_Im``, ``θ_Im``).
        The minus sign makes arrows point **towards** conductors.

        Sets
        ----
        mag_real, mag_imag : ndarray
            Magnitudes for real/imag parts, shape ``(n_freq,)``.
        angle_real, angle_imag : ndarray
            Directions in degrees for real/imag parts, shape
            ``(n_freq,)``.
        mag_err, angle_err : ndarray or None
            Crude uncertainty proxies from component-wise errors if
            available; otherwise ``None``.
        """
        if self._tipper is None:
            return

        Tx = self._tipper[:, 0, 0]
        Ty = self._tipper[:, 0, 1]

        self._mag_real = np.sqrt(Tx.real ** 2 + Ty.real ** 2)
        self._mag_imag = np.sqrt(Tx.imag ** 2 + Ty.imag ** 2)

        # Parkinson convention: negative components in atan2
        self._angle_real = np.degrees(np.arctan2(-Ty.real, -Tx.real))
        self._angle_imag = np.degrees(np.arctan2(-Ty.imag, -Tx.imag))

        if self._tipper_err is None:
            self._mag_err = None
            self._angle_err = None
            return

        Ex = self._tipper_err[:, 0, 0]
        Ey = self._tipper_err[:, 0, 1]

        # Simple combination as a proxy (component-wise 2-norm)
        self._mag_err = np.sqrt(Ex ** 2 + Ey ** 2)

        # Angle uncertainty proxy (small-angle approx, degrees)
        # atan2(dy, dx) ~ dy/dx near origin; keep bounded
        with np.errstate(divide="ignore", invalid="ignore"):
            ang = np.degrees(np.arctan2(Ex, Ey))
        self._angle_err = np.mod(ang, 45.0)  # legacy cap at 45°

    # ------------------------------------------------------------------ #
    # Read-only accessors for derived quantities
    # ------------------------------------------------------------------ #
    @property
    def amplitude(self) -> Optional[np.ndarray]:
        """Amplitude :math:`|T|` (same shape as :pyattr:`tipper`)."""
        return self._amplitude

    @property
    def phase(self) -> Optional[np.ndarray]:
        """Phase of tipper (degrees), same shape as :pyattr:`tipper`."""
        return self._phase

    @property
    def amplitude_err(self) -> Optional[np.ndarray]:
        """Absolute error on amplitude, or ``None``."""
        return self._amplitude_err

    @property
    def phase_err(self) -> Optional[np.ndarray]:
        """Absolute error on phase (degrees), or ``None``."""
        return self._phase_err

    @property
    def mag_real(self) -> Optional[np.ndarray]:
        """Real-part induction arrow magnitude (Parkinson)."""
        return self._mag_real

    @property
    def mag_imag(self) -> Optional[np.ndarray]:
        """Imag-part induction arrow magnitude (Parkinson)."""
        return self._mag_imag

    @property
    def angle_real(self) -> Optional[np.ndarray]:
        """Real-part induction arrow direction (deg, Parkinson)."""
        return self._angle_real

    @property
    def angle_imag(self) -> Optional[np.ndarray]:
        """Imag-part induction arrow direction (deg, Parkinson)."""
        return self._angle_imag

    @property
    def mag_err(self) -> Optional[np.ndarray]:
        """Uncertainty proxy for arrow magnitude, or ``None``."""
        return self._mag_err

    @property
    def angle_err(self) -> Optional[np.ndarray]:
        """Uncertainty proxy for arrow direction (deg), or ``None``."""
        return self._angle_err


    def set_amp_phase(
        self,
        r_array: np.ndarray,
        phi_array: np.ndarray,
    ) -> None:
        """
        Set tipper from amplitude ``r`` and phase ``phi`` (degrees).

        This converts ``(r, phi)`` to complex values with
        :math:`T = r\\,e^{j\\,\\phi}`, updates
        :pyattr:`tipper`, and recomputes amplitude/phase and
        Parkinson arrows.

        Parameters
        ----------
        r_array : array-like
            Amplitudes. Accepted shapes are the same as for
            :pyattr:`tipper` and are normalized internally to
            ``(n_freq, 1, 2)``. Values must be real.
        phi_array : array-like
            Phases in **degrees**. Same shape rules as ``r_array``.
            Values must be real.

        Raises
        ------
        ZError
            If shapes are incompatible or inputs are not real.

        Notes
        -----
        - If :pyattr:`tipper` is already set, both arrays must match
          its shape (after normalization).
        - If :pyattr:`tipper` is ``None``, the provided arrays define
          the shape of the new tipper.
        """
        # Normalize shapes similarly to `tipper` setter
        def _normalize(arr: np.ndarray, name: str) -> np.ndarray:
            a = np.asarray(arr)
            if np.iscomplexobj(a):
                # amplitudes and phases must be real-valued
                if np.linalg.norm(np.imag(a)) != 0:
                    raise ZError(f'"{name}" must be real-valued.')
                a = np.real(a)
            if a.ndim == 1 and a.shape[0] == 2:
                a = a[None, None, :]
            elif a.ndim == 2:
                if a.shape == (1, 2):
                    a = a[None, :, :]
                elif a.shape[1] == 2:
                    a = a[:, None, :]
                else:
                    raise ZError(
                        f'{name} shape must be (n, 2) or (1, 2); got '
                        f'{a.shape!r}.'
                    )
            elif a.ndim == 3:
                if a.shape[1:] != (1, 2):
                    raise ZError(
                        f'{name} 3-D shape must be (n, 1, 2); got '
                        f'{a.shape!r}.'
                    )
            else:
                raise ZError(f'Unsupported {name} shape: {a.shape!r}.')
            return a.astype(float, copy=False)

        r = _normalize(r_array, "r")
        phi = _normalize(phi_array, "phi")

        # If tipper already exists, enforce same shape
        if self._tipper is not None and (r.shape != self._tipper.shape):
            raise ZError(
                '"r" must match current tipper shape: '
                f"{r.shape!r} vs {self._tipper.shape!r}."
            )
        if r.shape != phi.shape:
            raise ZError(
                'Shapes of "r" and "phi" must match: '
                f"{r.shape!r} vs {phi.shape!r}."
            )

        # Build complex tipper from polar (degrees → radians)
        tip_new = r * np.exp(1j * np.deg2rad(phi))

        # Assign and recompute derived quantities
        self.tipper = tip_new
        self.compute_amp_phase()
        self.compute_mag_direction()

    def set_mag_direction(
        self,
        mag_real: np.ndarray,
        ang_real: np.ndarray,
        mag_imag: np.ndarray,
        ang_imag: np.ndarray,
    ) -> None:
        """
        Set tipper from Parkinson arrow magnitude and direction.

        Given real/imag arrow magnitudes ``M`` and directions
        ``θ`` (degrees) defined by

        .. math::

           M = \\sqrt{T_x^2 + T_y^2},\\quad
           \\theta = \\operatorname{atan2}(-T_y, -T_x),

        this method reconstructs the tipper components via

        .. math::

           T_x = -M \\cos\\theta,\\quad
           T_y = -M \\sin\\theta,

        applied to the **real** and **imaginary** parts separately.

        Parameters
        ----------
        mag_real, mag_imag : array-like
            Arrow magnitudes per frequency for real/imag parts.
            Shape ``(n_freq,)`` or broadcastable to it.
        ang_real, ang_imag : array-like
            Arrow directions (degrees) per frequency for real/imag
            parts. Shape ``(n_freq,)`` or broadcastable.

        Raises
        ------
        ZError
            If :pyattr:`tipper` is not initialized or shapes are
            incompatible.

        Notes
        -----
        No uncertainty propagation is applied here. After setting,
        amplitude/phase and arrow quantities are recomputed.
        """
        if self._tipper is None:
            raise ZError(
                "Tipper must be initialized before setting from "
                "magnitude/direction."
            )

        n = self._tipper.shape[0]

        def _as_vec(x: np.ndarray, name: str) -> np.ndarray:
            a = np.asarray(x, dtype=float).ravel()
            if a.size not in (1, n):
                raise ZError(
                    f'"{name}" must be a scalar or length-{n} vector; '
                    f"got length {a.size}."
                )
            return a if a.size == n else np.full(n, a[0], dtype=float)

        Mr = _as_vec(mag_real, "mag_real")
        Mi = _as_vec(mag_imag, "mag_imag")
        th_r = np.deg2rad(_as_vec(ang_real, "ang_real"))
        th_i = np.deg2rad(_as_vec(ang_imag, "ang_imag"))

        # Reconstruct components following Parkinson convention
        Tx = self._tipper[:, 0, 0].copy()
        Ty = self._tipper[:, 0, 1].copy()

        Tx = (-Mr * np.cos(th_r)).astype(float) + 1j * (
            -Mi * np.cos(th_i)
        ).astype(float)
        Ty = (-Mr * np.sin(th_r)).astype(float) + 1j * (
            -Mi * np.sin(th_i)
        ).astype(float)

        self._tipper[:, 0, 0] = Tx
        self._tipper[:, 0, 1] = Ty

        # Refresh derived quantities
        self.compute_mag_direction()
        self.compute_amp_phase()

    # ------------------------------------------------------------------ #
    # Rotation
    # ------------------------------------------------------------------ #
    def rotate(self, alpha: float | Sequence[float]) -> None:
        """
        Rotate tipper(s) by angle(s) ``alpha`` (degrees, CW positive).

        Angles are referenced to geographic North (X→North, Y→East),
        positive clockwise. Accepts a single angle applied to all
        frequencies or one angle per frequency.

        Parameters
        ----------
        alpha : float or sequence of float
            Rotation angle(s) in degrees.

        Raises
        ------
        ZError
            If :pyattr:`tipper` is missing or the number of angles is
            invalid.

        Notes
        -----
        The following attributes are updated:

        - :pyattr:`tipper`
        - :pyattr:`tipper_err` (if present)
        - :pyattr:`rotation_angle`

        Derived quantities are recomputed afterwards.
        """
        if self._tipper is None:
            raise ZError("Tipper is not set; cannot rotate.")

        n = self._tipper.shape[0]

        # Normalize alpha to a length-n vector, modulo 360
        if np.isscalar(alpha) or (
            isinstance(alpha, (list, tuple)) and len(alpha) == 1
        ):
            ang = float(np.asarray(alpha).ravel()[0]) % 360.0
            alphas = np.full(n, ang, dtype=float)
        else:
            a = np.asarray(alpha, dtype=float).ravel()
            if a.size not in (1, n):
                raise ZError(
                    f"Expected 1 angle or {n} angles; got {a.size}."
                )
            alphas = (a % 360.0) if a.size == n else np.full(n, a[0] % 360.0)

        # Update rotation history
        if isinstance(self.rotation_angle, float):
            self.rotation_angle = np.zeros(n, dtype=float)
        self.rotation_angle = (self.rotation_angle + alphas) % 360.0

        # Rotate each frequency vector (shape (1, 2))
        T_rot = np.empty_like(self._tipper, dtype=complex)
        Terr_rot = (
            None if self._tipper_err is None else np.empty_like(self._tipper_err)
        )

        for k in range(n):
            ang = 0.0 if np.isnan(alphas[k]) else float(alphas[k])
            if self._tipper_err is None:
                T_rot[k], _ = rotatevector_incl_errors(self._tipper[k, :, :], ang)
            else:
                T_rot[k], Terr_rot[k] = rotatevector_incl_errors(
                    self._tipper[k, :, :], ang, self._tipper_err[k, :, :]
                )

        self.tipper = T_rot
        if Terr_rot is not None:
            self.tipper_err = Terr_rot

        # Recompute derived quantities
        self.compute_mag_direction()
        self.compute_amp_phase()
