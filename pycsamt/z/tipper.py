# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
Tipper
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np

from ..exceptions import ZError
from ..log.logger import get_logger
from ..utils.zmath import (
    propagate_error_rect2polar,
    rotatevector_incl_errors,
)
from .base import BaseEM

logger = get_logger(__name__)


class Tipper(BaseEM):
    r"""
    Induction tipper container with derived quantities.

    The class stores complex tipper components
    :math:`(T_x, T_y)` per frequency and provides
    convenience computations:

    * amplitude / phase of each component
    * Parkinson induction arrows (magnitude and
      direction) for real and imaginary parts
    * rotation of the tipper with optional error
      propagation

    Internally, data are kept with shape
    ``(n_freq, 1, 2)`` where the last axis is ordered
    as ``[T_x, T_y]``.  Several common input shapes are
    accepted and normalized on assignment.

    Parameters
    ----------
    tipper_array : array-like, optional
        Complex tipper values.  Accepted shapes:

        * ``(n_freq, 1, 2)``  (preferred)
        * ``(n_freq, 2)``     (promoted)
        * ``(2,)`` or ``(1, 2)``  (single row)

        The last axis is ordered as ``[T_x, T_y]``.
    tipper_err_array : array-like, optional
        Absolute errors with the same normalized shape
        as ``tipper_array``.  If ``None``, uncertainties
        are not propagated.
    freq : array-like of float, optional
        Frequency vector in Hz.  Length must be
        ``n_freq`` and values strictly positive.

    Attributes
    ----------
    tipper : ndarray or None
        Complex tipper, shape ``(n_freq, 1, 2)``.
    tipper_err : ndarray or None
        Absolute errors on :pyattr:`tipper`.
    freq : ndarray or None
        Frequency vector (Hz).
    rotation_angle : float or ndarray
        Cumulative rotation (deg, clockwise positive).
        Stored as a scalar or length ``n_freq`` vector.
    amplitude, phase : ndarray or None
        :math:`|T|` and phase (deg) for each component.
        Shape ``(n_freq, 1, 2)``.
    amplitude_err, phase_err : ndarray or None
        Propagated absolute amplitude error and absolute
        phase error (deg), or ``None``.
    mag_real, mag_imag : ndarray or None
        Parkinson arrow magnitudes built from the real
        and imaginary parts, shape ``(n_freq,)``.
    angle_real, angle_imag : ndarray or None
        Parkinson arrow directions (deg), real and imag,
        shape ``(n_freq,)``.
    mag_err, angle_err : ndarray or None
        Heuristic uncertainty proxies for arrow metrics,
        or ``None``.

    Notes
    -----
    * **Shape normalization.**  Inputs with shape
      ``(n_freq, 2)`` or ``(2,)`` are promoted to
      ``(n_freq, 1, 2)``.  Errors follow the same rules.

    * **Parkinson convention.**  Arrow directions use
      a minus sign so arrows point **towards**
      conductors:

      .. math::

         \theta = \operatorname{atan2}(-Y,\,-X)

      where :math:`(X, Y)` are the component pairs
      taken from either the real or the imaginary part.

    * **Rotation.**  Angles are clockwise-positive and
      referenced to geographic North (X→North, Y→East).

    Examples
    --------
    Minimal construction and amplitude / phase:

    >>> import numpy as np
    >>> from pycsamt.z.tipper import Tipper
    >>> T = np.array([[1 + 1j, 0 + 1j]])
    >>> tip = Tipper(tipper_array=T)  # (1, 2)
    >>> tip.amplitude.shape, tip.phase.shape
    ((1, 1, 2), (1, 1, 2))

    Parkinson arrows for a known direction:

    >>> th = np.deg2rad(30.0)
    >>> Tarr = np.zeros((1, 1, 2), complex)
    >>> Tarr[0, 0, 0] = -np.cos(th) + 0j
    >>> Tarr[0, 0, 1] = -np.sin(th) + 0j
    >>> tip = Tipper(tipper_array=Tarr)
    >>> float(tip.mag_real[0])  # doctest: +ELLIPSIS
    1.0
    >>> float(tip.angle_real[0])  # doctest: +ELLIPSIS
    30.0

    See Also
    --------
    pycsamt.utils.zmath.rotatevector_incl_errors :
        Rotation with error propagation.
    pycsamt.z.base.BaseEM :
        Shared container utilities.

    References
    ----------
    .. [1] Parkinson, W. D. (1959). *Directions of rapid
       geomagnetic fluctuations*. Geophys. J. Int.
    .. [2] Chave, A. D., and Jones, A. G. (2012). *The
       Magnetotelluric Method: Theory and Practice*.
       Cambridge Univ. Press.
    """

    def __init__(
        self,
        tipper_array: np.ndarray | None = None,
        tipper_err_array: np.ndarray | None = None,
        freq: Sequence[float] | None = None,
        *,
        name: str | None = None,
        **kw: Any,
    ) -> None:
        super().__init__(name=name, **kw)

        self._tipper: np.ndarray | None = None
        self._tipper_err: np.ndarray | None = None
        self._freq: np.ndarray | None = None

        self.rotation_angle: float | np.ndarray = 0.0

        self._amplitude: np.ndarray | None = None
        self._amplitude_err: np.ndarray | None = None
        self._phase: np.ndarray | None = None
        self._phase_err: np.ndarray | None = None

        self._mag_real: np.ndarray | None = None
        self._mag_imag: np.ndarray | None = None
        self._angle_real: np.ndarray | None = None
        self._angle_imag: np.ndarray | None = None
        self._mag_err: np.ndarray | None = None
        self._angle_err: np.ndarray | None = None

        if tipper_array is not None:
            self.tipper = tipper_array
        if tipper_err_array is not None:
            self.tipper_err = tipper_err_array
        if freq is not None:
            self.freq = freq

        if self._tipper is not None and isinstance(self.rotation_angle, float):
            self.rotation_angle = np.zeros(self._tipper.shape[0], dtype=float)

        if self._tipper is not None:
            self.compute_amp_phase()
            self.compute_mag_direction()

    @property
    def tipper(self) -> np.ndarray | None:
        return self._tipper

    @tipper.setter
    def tipper(self, arr: np.ndarray) -> None:
        t = np.asarray(arr)

        if t.ndim == 1 and t.shape[0] == 2:
            t = t[None, None, :]
        elif t.ndim == 2:
            if t.shape == (1, 2):
                t = t[None, :, :]
            elif t.shape[1] == 2:
                t = t[:, None, :]
            else:
                raise ZError(
                    f"Tipper shape (n, 2) or (1, 2) expected; got {t.shape!r}."
                )
        elif t.ndim == 3:
            if t.shape[1:] != (1, 2):
                raise ZError(
                    "Tipper 3-D shape must be (n_freq, 1, 2); " f"got {t.shape!r}."
                )
        else:
            raise ZError(f"Unsupported tipper shape: {t.shape!r}.")

        self._tipper = t.astype(complex, copy=False)

        if isinstance(self.rotation_angle, float):
            self.rotation_angle = np.zeros(self._tipper.shape[0], dtype=float)

    @property
    def tipper_err(self) -> np.ndarray | None:
        return self._tipper_err

    @tipper_err.setter
    def tipper_err(self, arr: np.ndarray | None) -> None:
        if arr is None:
            self._tipper_err = None
            return

        e = np.asarray(arr)

        if e.ndim == 1 and e.shape[0] == 2:
            e = e[None, None, :]
        elif e.ndim == 2:
            if e.shape == (1, 2):
                e = e[None, :, :]
            elif e.shape[1] == 2:
                e = e[:, None, :]
            else:
                raise ZError(
                    "Tipper error shape (n, 2) or (1, 2) " f"expected; got {e.shape!r}."
                )
        elif e.ndim == 3:
            if e.shape[1:] != (1, 2):
                raise ZError(
                    "Tipper error 3-D shape must be "
                    f"(n_freq, 1, 2); got {e.shape!r}."
                )
        else:
            raise ZError(f"Unsupported tipper error shape: {e.shape!r}.")

        if self._tipper is not None and e.shape != self._tipper.shape:
            raise ZError(
                "'tipper_err' must match 'tipper' shape: "
                f"{e.shape!r} vs {self._tipper.shape!r}."
            )

        self._tipper_err = e.astype(float, copy=False)

    @property
    def freq(self) -> np.ndarray | None:
        return self._freq

    @freq.setter
    def freq(self, f: Sequence[float] | None) -> None:
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

    def compute_amp_phase(self) -> None:
        r"""
        Compute :math:`|T|` and phase (deg) for each component.

        For every frequency and for both components
        :math:`(T_x, T_y)`, the amplitude is
        :math:`|T| = \sqrt{(\Re T)^2 + (\Im T)^2}` and the
        phase is the argument of the complex number in
        degrees.

        When :pyattr:`tipper_err` is present, uncertainties
        are propagated component-wise using
        :func:`pycsamt.utils.zmath.propagate_error_rect2polar`.

        Returns
        -------
        None
            Results are stored in the attributes
            :pyattr:`amplitude`, :pyattr:`phase`,
            and, when applicable,
            :pyattr:`amplitude_err`, :pyattr:`phase_err`.

        Notes
        -----
        If :pyattr:`tipper` is ``None``, the method exits
        silently and leaves derived fields unchanged.

        Examples
        --------
        >>> import numpy as np
        >>> from pycsamt.z.tipper import Tipper
        >>> T = np.array([[[1 + 1j, 0 + 1j]]])
        >>> tip = Tipper(tipper_array=T)
        >>> tip.compute_amp_phase()
        >>> tip.amplitude.shape, tip.phase.shape
        ((1, 1, 2), (1, 1, 2))
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

        for k in range(T.shape[0]):
            for j in range(2):
                r = T[k, 0, j].real
                im = T[k, 0, j].imag
                dr = E[k, 0, j]
                dim = E[k, 0, j]
                r_e, p_e = propagate_error_rect2polar(r, dr, im, dim)
                amp_err[k, 0, j] = r_e
                ph_err[k, 0, j] = p_e

        self._amplitude_err = amp_err
        self._phase_err = ph_err

    def compute_mag_direction(self) -> None:
        r"""
        Compute Parkinson arrow magnitudes and directions.

        Real-part arrow:

        .. math::

           M_{\Re} = \sqrt{(\Re T_x)^2 + (\Re T_y)^2}

        .. math::

           \theta_{\Re}
           = \operatorname{atan2}(-\Re T_y,\,-\Re T_x)

        Imag-part arrow is computed identically, using
        :math:`\Im T_x` and :math:`\Im T_y`.

        The minus sign makes arrows point **towards**
        conductors (Parkinson convention).

        Returns
        -------
        None
            Results are stored in
            :pyattr:`mag_real`, :pyattr:`mag_imag`,
            :pyattr:`angle_real`, :pyattr:`angle_imag`.
            If :pyattr:`tipper_err` is present, heuristic
            proxies :pyattr:`mag_err` and :pyattr:`angle_err`
            are also set.

        Notes
        -----
        Angle uncertainties are estimated via a bounded
        small-angle proxy and wrapped within a legacy
        45° cap.

        Examples
        --------
        >>> import numpy as np
        >>> from pycsamt.z.tipper import Tipper
        >>> th = np.deg2rad(60.0)
        >>> Tarr = np.zeros((1, 1, 2), complex)
        >>> Tarr[0, 0, 0] = 0.0 + 1j * (-np.cos(th))
        >>> Tarr[0, 0, 1] = 0.0 + 1j * (-np.sin(th))
        >>> tip = Tipper(tipper_array=Tarr)
        >>> float(tip.angle_imag[0])  # doctest: +ELLIPSIS
        60.0
        """

        if self._tipper is None:
            return

        Tx = self._tipper[:, 0, 0]
        Ty = self._tipper[:, 0, 1]

        self._mag_real = np.sqrt(Tx.real**2 + Ty.real**2)
        self._mag_imag = np.sqrt(Tx.imag**2 + Ty.imag**2)

        self._angle_real = np.degrees(np.arctan2(-Ty.real, -Tx.real))
        self._angle_imag = np.degrees(np.arctan2(-Ty.imag, -Tx.imag))

        if self._tipper_err is None:
            self._mag_err = None
            self._angle_err = None
            return

        Ex = self._tipper_err[:, 0, 0]
        Ey = self._tipper_err[:, 0, 1]

        self._mag_err = np.sqrt(Ex**2 + Ey**2)

        with np.errstate(divide="ignore", invalid="ignore"):
            ang = np.degrees(np.arctan2(Ex, Ey))
        self._angle_err = np.mod(ang, 45.0)

    @property
    def amplitude(self) -> np.ndarray | None:
        return self._amplitude

    @property
    def phase(self) -> np.ndarray | None:
        return self._phase

    @property
    def amplitude_err(self) -> np.ndarray | None:
        return self._amplitude_err

    @property
    def phase_err(self) -> np.ndarray | None:
        return self._phase_err

    @property
    def mag_real(self) -> np.ndarray | None:
        return self._mag_real

    @property
    def mag_imag(self) -> np.ndarray | None:
        return self._mag_imag

    @property
    def angle_real(self) -> np.ndarray | None:
        return self._angle_real

    @property
    def angle_imag(self) -> np.ndarray | None:
        return self._angle_imag

    @property
    def mag_err(self) -> np.ndarray | None:
        return self._mag_err

    @property
    def angle_err(self) -> np.ndarray | None:
        return self._angle_err

    def set_amp_phase(
        self,
        r_array: np.ndarray,
        phi_array: np.ndarray,
    ) -> None:
        r"""
        Set tipper from amplitude :math:`r` and phase :math:`\phi`.

        Converts the provided real arrays to complex tipper
        values via

        .. math::

           T = r \, e^{j \, \phi},

        where :math:`\phi` is given in **degrees**.  Shapes
        are normalized to ``(n_freq, 1, 2)`` following the
        same rules as :pyattr:`tipper`.

        Parameters
        ----------
        r_array : array-like
            Real amplitudes.  Accepted shapes are the same
            as for :pyattr:`tipper`.  Must be real-valued.
        phi_array : array-like
            Real phases in **degrees**.  Same shape rules as
            ``r_array``.

        Raises
        ------
        ZError
            If shapes are incompatible or inputs are not
            real after normalization.

        Returns
        -------
        None
            The method updates :pyattr:`tipper` and then
            recomputes amplitude / phase and arrow metrics.

        Notes
        -----
        If :pyattr:`tipper` already exists, both arrays must
        match its normalized shape.

        Examples
        --------
        >>> import numpy as np
        >>> from pycsamt.z.tipper import Tipper
        >>> r = np.ones((2, 1, 2))
        >>> phi = np.zeros_like(r)
        >>> phi[:, 0, 1] = 90.0
        >>> tip = Tipper()
        >>> tip.set_amp_phase(r, phi)
        >>> tip.tipper.shape
        (2, 1, 2)
        """

        def _normalize(a: np.ndarray, name: str) -> np.ndarray:
            x = np.asarray(a)
            if np.iscomplexobj(x):
                if np.linalg.norm(np.imag(x)) != 0:
                    raise ZError(f'"{name}" must be real-valued.')
                x = np.real(x)
            if x.ndim == 1 and x.shape[0] == 2:
                x = x[None, None, :]
            elif x.ndim == 2:
                if x.shape == (1, 2):
                    x = x[None, :, :]
                elif x.shape[1] == 2:
                    x = x[:, None, :]
                else:
                    raise ZError(
                        f"{name} shape must be (n, 2) or (1, 2); " f"got {x.shape!r}."
                    )
            elif x.ndim == 3:
                if x.shape[1:] != (1, 2):
                    raise ZError(
                        f"{name} 3-D shape must be (n, 1, 2); got {x.shape!r}."
                    )
            else:
                raise ZError(f"Unsupported {name} shape: {x.shape!r}.")
            return x.astype(float, copy=False)

        r = _normalize(r_array, "r")
        phi = _normalize(phi_array, "phi")

        if self._tipper is not None and (r.shape != self._tipper.shape):
            raise ZError(
                '"r" must match current tipper shape: '
                f"{r.shape!r} vs {self._tipper.shape!r}."
            )
        if r.shape != phi.shape:
            raise ZError(
                'Shapes of "r" and "phi" must match: ' f"{r.shape!r} vs {phi.shape!r}."
            )

        tip_new = r * np.exp(1j * np.deg2rad(phi))

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
        r"""
        Set tipper from Parkinson magnitudes and directions.

        For each frequency, reconstruct components per the
        Parkinson convention:

        .. math::

           T_x = -M \cos\theta,\quad
           T_y = -M \sin\theta,

        applied separately to the **real** and the
        **imaginary** parts.  Angles are provided in
        **degrees**.

        Parameters
        ----------
        mag_real, mag_imag : array-like
            Arrow magnitudes for the real and imaginary
            parts.  Scalar or length ``n_freq``.
        ang_real, ang_imag : array-like
            Arrow directions (deg) for real and imaginary
            parts.  Scalar or length ``n_freq``.

        Raises
        ------
        ZError
            If :pyattr:`tipper` is not initialized or if the
            supplied vectors are not scalar or length
            ``n_freq``.

        Returns
        -------
        None
            The method updates :pyattr:`tipper` and then
            recomputes arrow metrics and amplitude / phase.

        Examples
        --------
        >>> import numpy as np
        >>> from pycsamt.z.tipper import Tipper
        >>> tip = Tipper(tipper_array=np.zeros((1, 1, 2), complex))
        >>> tip.set_mag_direction(1.0, 0.0, 2.0, 90.0)
        >>> tuple(np.round(tip.tipper[0, 0].real, 6))
        (-1.0, -0.0)
        """

        if self._tipper is None:
            raise ZError(
                "Tipper must be initialized before setting from " "magnitude/direction."
            )

        n = self._tipper.shape[0]

        def _as_vec(x: np.ndarray, name: str) -> np.ndarray:
            a = np.asarray(x, dtype=float).ravel()
            if a.size not in (1, n):
                raise ZError(
                    f'"{name}" must be a scalar or length-{n} '
                    f"vector; got length {a.size}."
                )
            return a if a.size == n else np.full(n, a[0], dtype=float)

        Mr = _as_vec(mag_real, "mag_real")
        Mi = _as_vec(mag_imag, "mag_imag")
        th_r = np.deg2rad(_as_vec(ang_real, "ang_real"))
        th_i = np.deg2rad(_as_vec(ang_imag, "ang_imag"))

        Tx = self._tipper[:, 0, 0].copy()
        Ty = self._tipper[:, 0, 1].copy()

        Tx = (-Mr * np.cos(th_r)).astype(float) + 1j * (-Mi * np.cos(th_i)).astype(
            float
        )
        Ty = (-Mr * np.sin(th_r)).astype(float) + 1j * (-Mi * np.sin(th_i)).astype(
            float
        )

        self._tipper[:, 0, 0] = Tx
        self._tipper[:, 0, 1] = Ty

        self.compute_mag_direction()
        self.compute_amp_phase()

    def rotate(self, alpha: float | Sequence[float]) -> None:
        r"""
        Rotate tipper(s) clockwise by the given angle(s).

        Angles are referenced to geographic North
        (X→North, Y→East).  Positive angles are clockwise.

        Parameters
        ----------
        alpha : float or sequence of float
            Single angle applied to all frequencies, or a
            length ``n_freq`` sequence of angles (deg).

        Raises
        ------
        ZError
            If :pyattr:`tipper` is missing, or if the number
            of angles is not 1 or ``n_freq``.

        Returns
        -------
        None
            The method updates :pyattr:`tipper`,
            :pyattr:`tipper_err` (if present), and
            :pyattr:`rotation_angle`, then recomputes
            amplitude / phase and arrow metrics.

        Notes
        -----
        Error propagation uses
        :func:`pycsamt.utils.zmath.rotatevector_incl_errors`.
        Angles are reduced modulo 360° in the rotation
        history.

        Examples
        --------
        >>> import numpy as np
        >>> from pycsamt.z.tipper import Tipper
        >>> T = np.array([[[1.0 + 0.0j, 0.0 + 0.0j]]])
        >>> tip = Tipper(tipper_array=T)
        >>> tip.rotate(90.0)
        >>> np.allclose(tip.rotation_angle, 90.0)
        True
        """

        if self._tipper is None:
            raise ZError("Tipper is not set; cannot rotate.")

        n = self._tipper.shape[0]

        if np.isscalar(alpha) or (isinstance(alpha, (list, tuple)) and len(alpha) == 1):
            ang = float(np.asarray(alpha).ravel()[0]) % 360.0
            alphas = np.full(n, ang, dtype=float)
        else:
            a = np.asarray(alpha, dtype=float).ravel()
            if a.size not in (1, n):
                raise ZError(f"Expected 1 angle or {n} angles; got {a.size}.")
            alphas = a % 360.0 if a.size == n else np.full(n, a[0] % 360.0)

        if isinstance(self.rotation_angle, float):
            self.rotation_angle = np.zeros(n, dtype=float)
        self.rotation_angle = (self.rotation_angle + alphas) % 360.0

        T_rot = np.empty_like(self._tipper, dtype=complex)
        Terr_rot = None if self._tipper_err is None else np.empty_like(self._tipper_err)

        for k in range(n):
            ang = 0.0 if np.isnan(alphas[k]) else float(alphas[k])
            if self._tipper_err is None:
                T_rot[k], _ = rotatevector_incl_errors(self._tipper[k, :, :], ang)
            else:
                T_rot[k], Terr_rot[k] = rotatevector_incl_errors(
                    self._tipper[k, :, :],
                    ang,
                    self._tipper_err[k, :, :],
                )

        self.tipper = T_rot
        if Terr_rot is not None:
            self.tipper_err = Terr_rot

        self.compute_mag_direction()
        self.compute_amp_phase()
