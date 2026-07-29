# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
Impedance tensor (:class:`Z`) built on top of :class:`ResPhase`.

This module defines the high-level :class:`Z` class, which
**inherits** from :class:`~pycsamt.z.resphase.ResPhase`. It adds
conveniences around construction, input validation, storage of
the rotation history, and property setters for the impedance
tensor and the frequency vector.

"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from ..exceptions import ZError
from ..log.logger import get_logger
from ..utils.zmath import (
    invertmatrix_incl_errors,
    rotatematrix_incl_errors,
)
from .resphase import ResPhase

logger = get_logger(__name__)


__all__ = ["Z"]


class Z(ResPhase):
    r"""
    High-level impedance tensor container built on
    :class:`~pycsamt.z.resphase.ResPhase`.

    ``Z`` stores complex impedances with shape ``(n, 2, 2)`` and
    manages absolute uncertainties of the same shape. A single
    ``(2, 2)`` matrix is promoted to a one-row stack.

    Frequency is held as a 1-D vector ``(n,)`` in Hz. When either
    :pyattr:`z` or :pyattr:`freq` is set, apparent resistivity
    :math:`\rho` and phase :math:`\phi` are recomputed via the
    parent class.

    Parameters
    ----------
    z_array : ndarray, shape ``(n, 2, 2)`` or ``(2, 2)``, optional
        Complex impedance tensor(s). Real arrays are cast to
        complex.
    z_err_array : ndarray, optional
        Absolute per-component errors on :pyattr:`z`. Shapes as
        for ``z_array`` after any promotion.
    freq : array-like of float, optional
        Frequency vector in Hz. Must be 1-D and strictly
        positive. Length must equal the number of rows in
        :pyattr:`z`.
    name : str, optional
        Display name forwarded to
        :class:`~pycsamt.z.base.BaseEM`.
    meta : dict, optional
        Arbitrary metadata forwarded to ``BaseEM``.
    verbose : int, default 0
        Verbosity forwarded to ``BaseEM``.

    Attributes
    ----------
    z : ndarray or None
        Complex impedance stack, shape ``(n, 2, 2)``.
    z_err : ndarray or None
        Absolute errors on :pyattr:`z`.
    freq : ndarray or None
        Frequency vector in Hz.
    rotation_angle : float or ndarray
        Accumulated rotation(s) in degrees (CW positive). A
        scalar before data are set; a length-``n`` vector once
        data exist.

    Notes
    -----
    The index convention is

    - ``Zxx`` → ``[:, 0, 0]``
    - ``Zxy`` → ``[:, 0, 1]``
    - ``Zyx`` → ``[:, 1, 0]``
    - ``Zyy`` → ``[:, 1, 1]``

    Uncertainties are **absolute** (not relative) and are
    propagated to :math:`\rho` and :math:`\phi` when available.

    Examples
    --------
    Create from a single tensor and frequency::

        >>> import numpy as np
        >>> from pycsamt.z.z import Z
        >>> z = np.array([[0+0j, 1+1j],
        ...               [-1-1j, 0+0j]])
        >>> obj = Z(z_array=z, freq=[1.0])
        >>> obj.resistivity.shape
        (1, 2, 2)

    See Also
    --------
    pycsamt.z.resphase.ResPhase
        Implements the :math:`\rho` / :math:`\phi` mapping and
        error propagation.
    pycsamt.z.base.BaseEM
        Minimal mixin used for naming, metadata, and summary.

    References
    ----------
    .. [1] Chave, A. D., & Jones, A. G. (2012). *The
           Magnetotelluric Method: Theory and Practice*. CUP.
    .. [2] Simpson, F., & Bahr, K. (2005). *Practical
           Magnetotellurics*. CUP.
    """

    def __init__(
        self,
        z_array: np.ndarray | None = None,
        z_err_array: np.ndarray | None = None,
        freq: np.ndarray | None = None,
        *,
        name: str | None = None,
        meta: dict | None = None,
        verbose: int = 0,
    ) -> None:
        # Ensure BaseEM is initialized (logger, name, meta,
        # verbose). ResPhase does not call super().__init__.
        ResPhase.__init__(
            self,
            name=name,
            meta={} if meta is None else meta,
            verbose=verbose,
        )
        self._z: np.ndarray | None = None
        self._z_err: np.ndarray | None = None
        self._freq: np.ndarray | None = None

        # Either scalar or length-n vector (deg, CW positive).
        self.rotation_angle: float | np.ndarray = 0.0

        # Normalize and set initial arrays.
        if z_array is not None:
            self.z = z_array

        if z_err_array is not None:
            z_err = np.asarray(z_err_array)
            if z_err.ndim == 2 and z_err.shape == (2, 2):
                z_err = z_err[None, ...]
            self._check_tensor_shape(z_err, name="z_err")
            self._z_err = z_err.astype(float, copy=False)

        if freq is not None:
            self.freq = freq

        # If both Z and freq are present, derive rho/phi.
        if (self._z is not None) and (self._freq is not None):
            try:
                self.compute_resistivity_phase()
            except Exception as exc:  # pragma: no cover
                self.log.error(
                    "Failed to compute rho/phi at init: %s",
                    exc,
                )

    def __len__(self) -> int:
        """Return number of frequency samples (``n_freq``)."""
        return self.n_freq

    @classmethod
    def from_res_phase(
        cls,
        rho: np.ndarray,
        phi: np.ndarray,
        freq: np.ndarray,
        *,
        rho_err: np.ndarray | None = None,
        phi_err: np.ndarray | None = None,
        name: str | None = None,
        meta: dict | None = None,
        verbose: int = 0,
    ) -> Z:
        r"""
        Build a :class:`Z` instance from
        :math:`(\rho, \phi, f)` (with errors).

        This convenience constructor reconstructs complex
        impedances from apparent resistivity :math:`\rho`,
        phase :math:`\phi` (deg), and frequency, and attaches
        optional uncertainties. Internally it calls
        :meth:`~pycsamt.z.resphase.ResPhase.set_res_phase`.

        Parameters
        ----------
        rho : ndarray, shape ``(n, 2, 2)``
            Apparent resistivity (Ω·m).
        phi : ndarray, shape ``(n, 2, 2)``
            Phase (degrees).
        freq : ndarray, shape ``(n,)``
            Frequency in Hz. Must be strictly positive.
        rho_err : ndarray, optional
            Absolute error on :math:`\rho` (same shape).
        phi_err : ndarray, optional
            Absolute error on :math:`\phi` (deg, same shape).
        name, meta, verbose :
            Forwarded to the constructor.

        Returns
        -------
        Z
            A fully initialized impedance container.

        Notes
        -----
        The magnitude is recovered via

        .. math::

           |Z| = \sqrt{5 \, f \, \rho} \, ,

        followed by Euler reconstruction using :math:`\phi`.

        Examples
        --------
        >>> import numpy as np
        >>> rho = np.ones((2, 2, 2))
        >>> phi = np.zeros_like(rho)
        >>> f = np.array([10.0, 1.0])
        >>> Z.from_res_phase(rho, phi, f).z.shape
        (2, 2, 2)

        See Also
        --------
        pycsamt.z.resphase.ResPhase.set_res_phase
            Underlying implementation used here.
        """

        z = cls(name=name, meta=meta, verbose=verbose)
        z.set_res_phase(
            rho,
            phi,
            freq,
            res_err_array=rho_err,
            phase_err_array=phi_err,
        )
        return z

    @property
    def freq(self) -> np.ndarray | None:
        """Frequency vector in Hz, shape ``(n_freq,)`` or ``None``."""
        return self._freq

    @freq.setter
    def freq(self, freq_arr: np.ndarray | None) -> None:
        if freq_arr is None:
            self._freq = None
            return

        f = np.asarray(freq_arr, dtype=float)
        if f.ndim != 1:
            raise ZError("Frequency must be a 1-D array.")
        if np.any(f <= 0.0):
            raise ZError("Frequencies must be strictly positive.")

        if self._z is not None and f.size != self._z.shape[0]:
            raise ZError(
                "Length of 'freq' does not match Z: " f"{f.size} vs {self._z.shape[0]}."
            )

        self._freq = f

        if isinstance(self.rotation_angle, float) and self._z is not None:
            self.rotation_angle = np.repeat(
                self.rotation_angle,
                self._z.shape[0],
            )

        if self._z is not None:
            self.compute_resistivity_phase()

    @property
    def z(self) -> np.ndarray | None:
        """Complex Z, shape ``(n_freq, 2, 2)`` or ``None``."""
        return self._z

    @z.setter
    def z(self, z_array: np.ndarray) -> None:
        z = np.asarray(z_array)
        if z.ndim == 2 and z.shape == (2, 2):
            z = z[None, ...]
        self._check_tensor_shape(z, name="z")
        self._z = z.astype(complex, copy=False)

        if isinstance(self.rotation_angle, float):
            self.rotation_angle = np.zeros(
                self._z.shape[0],
                dtype=float,
            )

        if self._freq is not None:
            try:
                self.compute_resistivity_phase()
            except Exception as exc:  # pragma: no cover
                self.log.error(
                    "Failed to compute rho/phi after setting Z: %s",
                    exc,
                )

    @staticmethod
    def _check_tensor_shape(arr: np.ndarray, *, name: str) -> None:
        if arr.ndim != 3 or arr.shape[1:] != (2, 2):
            raise ZError(f"'{name}' must have shape (n_freq, 2, 2); got {arr.shape!r}.")

    @property
    def z_err(self) -> np.ndarray | None:
        """
        Absolute errors on complex Z, shape ``(n, 2, 2)`` or ``None``.
        """
        return self._z_err

    @z_err.setter
    def z_err(self, z_err_array: np.ndarray | None) -> None:
        if z_err_array is None:
            self._z_err = None
            if self._z is not None and self._freq is not None:
                self.compute_resistivity_phase()
            return

        zerr = np.asarray(z_err_array)
        if zerr.ndim == 2 and zerr.shape == (2, 2):
            zerr = zerr[None, ...]
        self._check_tensor_shape(zerr, name="z_err")

        if self._z is not None and zerr.shape != self._z.shape:
            raise ZError(
                "'z_err' shape must match 'z': " f"{zerr.shape!r} vs {self._z.shape!r}."
            )

        self._z_err = zerr.astype(float, copy=False)

        if (self._z is not None) and (self._freq is not None):
            self.compute_resistivity_phase()

    @property
    def inverse(self) -> np.ndarray:
        """
        Inverse tensor :math:`Z^{-1}` for each frequency.

        If ``z_err`` is present, errors are propagated
        internally, but only the inverted tensor is
        returned (legacy API).
        """
        if self._z is None:
            raise ZError("Z is not set; cannot invert.")

        inv = np.empty_like(self._z, dtype=complex)
        if self._z_err is None:
            for k in range(self._z.shape[0]):
                try:
                    inv[k] = np.linalg.inv(self._z[k])
                except np.linalg.LinAlgError:
                    raise ZError(f"Impedance tensor at index {k} is singular.")
        else:
            for k in range(self._z.shape[0]):
                try:
                    inv_k, _ = invertmatrix_incl_errors(
                        self._z[k],
                        self._z_err[k],
                    )
                    inv[k] = inv_k
                except Exception as exc:
                    raise ZError(
                        f"Failed to invert tensor at index {k}: {exc}"
                    ) from exc
        return inv

    def rotate(self, alpha: Sequence[float] | float) -> None:
        r"""
        Rotate :pyattr:`z` by angle(s) ``alpha`` (degrees, CW
        positive).

        The rotation is referenced to geographic axes (X→North,
        Y→East). A single angle is applied to all frequencies; a
        length-``n`` sequence applies element-wise angles. Errors
        are rotated consistently when present.

        Parameters
        ----------
        alpha : float or sequence of float
            Rotation angle(s) in degrees. If a sequence is given,
            its length must be 1 or equal to ``n``.

        Raises
        ------
        pycsamt.exceptions.ZError
            If :pyattr:`z` is missing or the number of angles is
            invalid.

        Notes
        -----
        The following are updated:

        - :pyattr:`z`
        - :pyattr:`z_err` (if present)
        - :pyattr:`rotation_angle`
        - Derived :math:`\rho` and :math:`\phi`

        Rotation uses
        :func:`~pycsamt.utils.zmath.rotatematrix_incl_errors`
        per frequency.

        Examples
        --------
        >>> import numpy as np
        >>> from pycsamt.z.z import Z
        >>> z = np.zeros((2, 2, 2), complex)
        >>> f = np.array([10.0, 1.0])
        >>> obj = Z(z_array=z, freq=f)
        >>> obj.rotate(30.0)  # single angle for all rows
        >>> obj.rotation_angle.shape
        (2,)
        """

        if self._z is None:
            raise ZError("Z is not set; cannot rotate.")

        n = self._z.shape[0]
        if np.isscalar(alpha) or (isinstance(alpha, (list, tuple)) and len(alpha) == 1):
            ang = float(np.asarray(alpha).ravel()[0]) % 360.0
            alphas = np.full(n, ang, dtype=float)
        else:
            a = np.asarray(alpha, dtype=float).ravel()
            if a.size not in (1, n):
                raise ZError(f"Expected 1 angle or {n} angles; got {a.size}.")
            alphas = (a % 360.0) if a.size == n else np.full(n, a[0] % 360.0)

        if isinstance(self.rotation_angle, float):
            self.rotation_angle = np.zeros(n, dtype=float)
        self.rotation_angle = (self.rotation_angle + alphas) % 360.0

        z_rot = np.empty_like(self._z, dtype=complex)
        zerr_rot = None if self._z_err is None else np.empty_like(self._z_err)

        for k in range(n):
            ang = 0.0 if np.isnan(alphas[k]) else float(alphas[k])
            if self._z_err is None:
                z_rot[k], _ = rotatematrix_incl_errors(
                    self._z[k],
                    ang,
                )
            else:
                z_rot[k], zerr_rot[k] = rotatematrix_incl_errors(
                    self._z[k],
                    ang,
                    self._z_err[k],
                )

        self.z = z_rot
        if zerr_rot is not None:
            self.z_err = zerr_rot

    def remove_static_shift(
        self,
        reduce_res_factor_x: float | Sequence[float] = 1.0,
        reduce_res_factor_y: float | Sequence[float] = 1.0,
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""
        Remove static shift using resistivity-scale correction factors.

        Assumes the observed tensor is :math:`Z = S \cdot Z_0`, with:

        .. math::

           S = \begin{bmatrix}
                 \sqrt{f_x} & 0 \\
                 0           & \sqrt{f_y}
               \end{bmatrix}

        where ``f_x`` and ``f_y`` are **resistivity-scale** factors.
        The corrected tensor is then:

        .. math::

           Z_0 = S^{-1} \cdot Z

        Parameters
        ----------
        reduce_res_factor_x : float or sequence, default=1.0
            Factor(s) for X components (resistivity scale). Either a
            single value or one per frequency.
        reduce_res_factor_y : float or sequence, default=1.0
            Factor(s) for Y components (resistivity scale). Either a
            single value or one per frequency.

        Returns
        -------
        static_shift : ndarray, shape (n_freq, 2, 2)
            Static-shift matrices applied at each frequency.
        z_corrected : ndarray, shape (n_freq, 2, 2)
            Corrected impedance tensor :math:`Z_0`.

        Raises
        ------
        ZError
            If Z is missing or factor lengths are invalid.

        Notes
        -----
        Factors are in **ρ scale**; the matrix entries use their
        square roots. No uncertainty propagation is applied here.
        """
        if self._z is None:
            raise ZError("Z is not set; cannot remove static shift.")

        n = self._z.shape[0]

        # Normalize factors to vectors of length n.
        def _normalize_factor(v: float | Sequence[float], name: str) -> np.ndarray:
            vv = np.asarray(v, dtype=float).ravel()
            if vv.size == 1:
                out = np.full(n, vv[0], dtype=float)
            elif vv.size == n:
                out = vv
            else:
                raise ZError(
                    f"'{name}' must be a scalar or length-{n} sequence; "
                    f"got length {vv.size}."
                )
            if np.any(out <= 0.0):
                raise ZError(f"'{name}' must contain positive values.")
            return out

        fx = _normalize_factor(reduce_res_factor_x, "reduce_res_factor_x")
        fy = _normalize_factor(reduce_res_factor_y, "reduce_res_factor_y")

        # Build per-frequency S and apply S^{-1} * Z.
        zcorr = np.empty_like(self._z, dtype=complex)
        S = np.zeros((n, 2, 2), dtype=float)
        S[:, 0, 0] = np.sqrt(fx)
        S[:, 1, 1] = np.sqrt(fy)

        Sinv = np.zeros_like(S)
        # Since S is diagonal, inverse is simple.
        Sinv[:, 0, 0] = 1.0 / S[:, 0, 0]
        Sinv[:, 1, 1] = 1.0 / S[:, 1, 1]

        # Apply Sinv @ Z for each frequency.
        for k in range(n):
            zcorr[k] = Sinv[k].dot(self._z[k])

        return S, zcorr

    def remove_distortion(
        self,
        distortion_tensor: np.ndarray,
        distortion_err_tensor: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
        r"""
        Remove galvanic distortion ``D`` from the observed impedance
        tensor :math:`Z` to obtain the undistorted tensor
        :math:`Z_0 = D^{-1} Z`.

        Uncertainty propagation is included (first-order, 1-norm). If
        either ``z_err`` or ``distortion_err_tensor`` is missing, zeros
        are assumed for the corresponding errors.

        Parameters
        ----------
        distortion_tensor : ndarray, shape (2, 2) or (n, 2, 2)
            Real 2×2 galvanic distortion matrix ``D``. If a stack is
            provided, only the first slice is used (distortion is assumed
            time-invariant).
        distortion_err_tensor : ndarray, shape (2, 2) or (n, 2, 2), \
                optional
            Absolute errors on ``D``. If omitted, zeros are assumed.
            If a stack is provided, only the first slice is used.

        Returns
        -------
        D : ndarray, shape (2, 2)
            The (real) distortion tensor used.
        Z0 : ndarray, shape (n_freq, 2, 2)
            Corrected impedance tensor :math:`Z_0 = D^{-1} Z`.
        Z0_err : ndarray or None, shape (n_freq, 2, 2)
            Propagated absolute errors on :math:`Z_0`. ``None`` if
            both input errors were ``None``.

        Raises
        ------
        ZError
            If ``Z`` is missing, distortion shapes are invalid, or the
            distortion matrix is singular.

        Notes
        -----
        Error propagation for a component
        :math:`(Z_0)_{ij} = \sum_k (D^{-1})_{ik} Z_{kj}` uses the
        simple 1-norm bound:

        .. math::

            \Delta (Z_0)_{ij} \approx
            \sum_k \\big(
              |\Delta (D^{-1})_{ik}| \cdot |Z_{kj}| +
              |(D^{-1})_{ik}| \\cdot |\Delta Z_{kj}|
            \big).

        Examples
        --------
        >>> D = np.array([[1.2, 0.5], [0.35, 2.1]])
        >>> D_used, Z0, Z0_err = zobj.remove_distortion(D)
        """
        if self._z is None:
            raise ZError("Z is not set; cannot remove distortion.")

        # --- Normalize distortion tensors to single 2×2 real arrays --------
        D = np.asarray(distortion_tensor, dtype=float)
        if D.ndim == 3:
            logger.info("Distortion provided as a stack; using the first slice.")
            D = D[0]
        if D.shape != (2, 2):
            raise ZError(
                "distortion_tensor must have shape (2, 2) or (n, 2, 2); "
                f"got {D.shape!r}."
            )

        if distortion_err_tensor is None:
            D_err = np.zeros((2, 2), dtype=float)
        else:
            D_err = np.asarray(distortion_err_tensor, dtype=float)
            if D_err.ndim == 3:
                logger.info("Distortion error provided as a stack; using first slice.")
                D_err = D_err[0]
            if D_err.shape != (2, 2):
                raise ZError(
                    "distortion_err_tensor must have shape (2, 2) or "
                    f"(n, 2, 2); got {D_err.shape!r}."
                )

        # --- Invert D with error propagation (errors of DI ignored later) ---
        try:
            DI, DI_err = invertmatrix_incl_errors(D, D_err)
        except np.linalg.LinAlgError as exc:
            raise ZError("Distortion tensor is singular; cannot invert.") from exc
        except Exception as exc:  # pragma: no cover
            raise ZError(f"Failed to invert distortion tensor: {exc}") from exc

        # --- Build corrected Z0 = DI @ Z, with propagated errors ------------
        n = self._z.shape[0]
        Z0 = np.empty_like(self._z, dtype=complex)

        # If no error info anywhere, we can return None for Z0_err.
        z_err_in = self._z_err
        if z_err_in is None and not np.any(D_err):
            Z0[:] = DI @ self._z.transpose(0, 2, 1)
            Z0 = Z0.transpose(0, 2, 1)  # back to (n, 2, 2)
            return D, Z0, None

        # Otherwise, prepare zero arrays for missing inputs.
        if z_err_in is None:
            z_err = np.zeros_like(self._z, dtype=float)
        else:
            z_err = np.asarray(z_err_in, dtype=float)
            if z_err.shape != self._z.shape:
                raise ZError(
                    f"'z_err' must match 'z' shape: {z_err.shape!r} vs "
                    f"{self._z.shape!r}."
                )

        Z0_err = np.zeros_like(self._z, dtype=float)

        # Compute per-frequency corrections with 1-norm error bound.
        for k in range(n):
            Zk = self._z[k]
            Z0[k] = DI @ Zk

            # error on each (i, j) uses k-sum over input 2×2:
            # sum_k |DI_err[i,k] * Z[k,j]| + |DI[i,k] * Z_err[k,j]|
            for i in range(2):
                for j in range(2):
                    term = (
                        np.abs(DI_err[i, 0] * Zk[0, j])
                        + np.abs(DI[i, 0] * z_err[k, 0, j])
                        + np.abs(DI_err[i, 1] * Zk[1, j])
                        + np.abs(DI[i, 1] * z_err[k, 1, j])
                    )
                    Z0_err[k, i, j] = term

        return D, Z0, Z0_err

    def _compute_det_variance(self) -> np.ndarray | None:
        r"""
        Approximate variance of :math:`\det(Z)` from Z and ``z_err``.

        Uses a simple central-difference style perturbation:

        .. math::

           \sigma_{\det} \approx
           \tfrac{1}{2}\,\big|\det(Z+\Delta Z) -
           \det(Z-\Delta Z)\big|

        and returns :math:`\sigma_{\det}^2` as a crude variance
        proxy. If no error is available the method returns ``None``.

        Returns
        -------
        ndarray or None
            Array of shape ``(n_freq,)`` with the approximate
            variance of the determinant, or ``None`` if no error
            information is available.
        """
        if self._z_err is None or self._z is None:
            return None

        det_plus = np.linalg.det(self._z + self._z_err)
        det_minus = np.linalg.det(self._z - self._z_err)
        sigma_det = 0.5 * np.abs(det_plus - det_minus)

        return sigma_det**2

    # 1-D / 2-D projections
    @property
    def only_1d(self) -> np.ndarray:
        r"""
        Return a **1-D**-like version of :math:`Z`.

        The diagonal entries are set to zero. The off-diagonal
        entries retain their complex signs, but their magnitudes are
        set to the mean of the **original off-diagonal magnitudes**
        at each frequency:

        .. code-block:: text

           |Z01|, |Z10| → m = 0.5 (|Z01| + |Z10|)
           Z01 := sign(Z01) * m
           Z10 := sign(Z10) * m

        Returns
        -------
        ndarray
            New tensor with shape ``(n_freq, 2, 2)``.
        """
        if self._z is None:
            raise ZError("Z is not set; cannot build 1-D projection.")

        z1d = self._z.copy()
        for k in range(z1d.shape[0]):
            z1d[k, 0, 0] = 0.0
            z1d[k, 1, 1] = 0.0
            z01 = z1d[k, 0, 1]
            z10 = z1d[k, 1, 0]
            mean_abs = 0.5 * (np.abs(z01) + np.abs(z10))
            # preserve complex sign (phase) of each entry
            z1d[k, 0, 1] = np.exp(1j * np.angle(z01)) * mean_abs
            z1d[k, 1, 0] = np.exp(1j * np.angle(z10)) * mean_abs
        return z1d

    @property
    def only_2d(self) -> np.ndarray:
        """
        Return a **2-D**-like version of :math:`Z`.

        The diagonal entries are set to zero; off-diagonal terms are
        kept unchanged.

        Returns
        -------
        ndarray
            New tensor with shape ``(n_freq, 2, 2)``.
        """
        if self._z is None:
            raise ZError("Z is not set; cannot build 2-D projection.")

        z2d = self._z.copy()
        z2d[:, 0, 0] = 0.0
        z2d[:, 1, 1] = 0.0
        return z2d

    # Basic algebraic invariants
    @property
    def trace(self) -> np.ndarray:
        r"""
        Trace of :math:`Z` at each frequency.

        Returns
        -------
        ndarray
            Array of shape ``(n_freq,)`` with
            :math:`\operatorname{tr}(Z)`.
        """
        if self._z is None:
            raise ZError("Z is not set; cannot compute trace.")
        return np.einsum("kii->k", self._z)

    @property
    def trace_err(self) -> np.ndarray | None:
        r"""
        Approximate error on :math:`\operatorname{tr}(Z)`.

        Returns
        -------
        ndarray or None
            If ``z_err`` is available, returns an array with
            shape ``(n_freq,)`` computed by summing the diagonal
            errors :math:`\Delta Z_{00} + \Delta Z_{11}`. Otherwise
            ``None``.
        """
        if self._z_err is None:
            return None
        return self._z_err[:, 0, 0] + self._z_err[:, 1, 1]

    @property
    def skew(self) -> np.ndarray:
        r"""
        Linear-algebra skew :math:`Z_{01} - Z_{10}`.

        .. note::
           This is *not* the MT skew used in dimensionality
           analysis; it is the simple matrix skew.

        Returns
        -------
        ndarray
            Array of shape ``(n_freq,)`` with
            :math:`Z_{01} - Z_{10}`.
        """
        if self._z is None:
            raise ZError("Z is not set; cannot compute skew.")
        return self._z[:, 0, 1] - self._z[:, 1, 0]

    @property
    def skew_err(self) -> np.ndarray | None:
        r"""
        Approximate error on the linear-algebra skew.

        Returns
        -------
        ndarray or None
            If ``z_err`` is available, returns an array with
            shape ``(n_freq,)`` computed as
            :math:`\Delta Z_{01} + \Delta Z_{10}` (1-norm bound).
            Otherwise ``None``.
        """
        if self._z_err is None:
            return None
        return self._z_err[:, 0, 1] + self._z_err[:, 1, 0]

    @property
    def det(self) -> np.ndarray:
        r"""
        Determinant of :math:`Z` at each frequency.

        Returns
        -------
        ndarray
            Array of shape ``(n_freq,)`` with :math:`\det(Z)`.
        """
        if self._z is None:
            raise ZError("Z is not set; cannot compute determinant.")
        return np.linalg.det(self._z)

    @property
    def det_err(self) -> np.ndarray | None:
        r"""
        Approximate error on :math:`\det(Z)`.

        Uses a central-difference style perturbation:

        .. math::

           \Delta \det(Z) \approx
           \tfrac{1}{2}\\,\big|\det(Z+\Delta Z) -
           \det(Z-\Delta Z)\big|.

        Returns
        -------
        ndarray or None
            Array of shape ``(n_freq,)`` with an error proxy, or
            ``None`` if ``z_err`` is not available.
        """
        if self._z_err is None or self._z is None:
            return None
        det_plus = np.linalg.det(self._z + self._z_err)
        det_minus = np.linalg.det(self._z - self._z_err)
        return 0.5 * np.abs(det_plus - det_minus)

    @property
    def norm(self) -> np.ndarray:
        r"""
        Frobenius norm :math:`\lVert Z \rVert_F` at each frequency.

        Returns
        -------
        ndarray
            Array of shape ``(n_freq,)`` with Frobenius norms.
        """
        if self._z is None:
            raise ZError("Z is not set; cannot compute norm.")
        return np.linalg.norm(self._z, axis=(1, 2))

    @property
    def norm_err(self) -> np.ndarray | None:
        """
        Approximate error on the Frobenius norm.

        Uses a first-order bound combining real and imaginary parts.

        Returns
        -------
        ndarray or None
            Array of shape ``(n_freq,)`` or ``None`` if ``z_err`` is
            not available.
        """
        if self._z_err is None or self._z is None:
            return None

        nrm = self.norm
        out = np.zeros_like(nrm, dtype=float)

        for k in range(self._z.shape[0]):
            zk = self._z[k]
            ek = self._z_err[k]
            # derivative of sqrt(sum x_i^2) ~ (1/||z||) * sum |e_i * x_i|
            # applied to real and imag parts
            rad = 0.0
            for i in range(2):
                for j in range(2):
                    rad += (ek[i, j] * np.real(zk[i, j])) ** 2
                    rad += (ek[i, j] * np.imag(zk[i, j])) ** 2
            out[k] = 0.0 if nrm[k] == 0.0 else np.sqrt(rad) / nrm[k]

        return out

    @property
    def invariants(self) -> dict[str, np.ndarray]:
        r"""
        Compute several algebraic invariants of :math:`Z`.

        Returns
        -------
        dict
            Dictionary with keys:

            - ``'z1'`` :math:`= (Z_{01} - Z_{10}) / 2`
            - ``'det'`` determinant of :math:`Z`
            - ``'det_real'`` determinant of :math:`\Re(Z)`
            - ``'det_imag'`` determinant of :math:`\Im(Z)`
            - ``'trace'`` trace of :math:`Z`
            - ``'skew'`` linear-algebra skew
            - ``'norm'`` Frobenius norm
            - ``'lambda_plus'``, ``'lambda_minus'`` :
              :math:`z_1 \pm \sqrt{z_1^2 / \det(Z)}`
            - ``'sigma_plus'``, ``'sigma_minus'`` :
              scalar combinations of norm and determinant magnitudes.

        Notes
        -----
        ``lambda_±`` are simple combinations frequently used in
        analytic manipulations of 2×2 matrices. ``sigma_±`` follow the
        historical form found in legacy code; they reduce to simple
        combinations of :math:`\lVert Z \rVert_F` and :math:`|\det|`.
        """
        if self._z is None:
            raise ZError("Z is not set; cannot compute invariants.")

        inv: dict[str, np.ndarray] = {}

        z1 = (self._z[:, 0, 1] - self._z[:, 1, 0]) / 2.0
        inv["z1"] = z1

        detz = self.det
        inv["det"] = detz

        inv["det_real"] = np.linalg.det(np.real(self._z))
        inv["det_imag"] = np.linalg.det(np.imag(self._z))

        inv["trace"] = self.trace
        inv["skew"] = self.skew
        nrm = self.norm
        inv["norm"] = nrm

        with np.errstate(divide="ignore", invalid="ignore"):
            sqrt_term = np.sqrt((z1 * z1) / detz)
        inv["lambda_plus"] = z1 + sqrt_term
        inv["lambda_minus"] = z1 - sqrt_term

        # Keep the historical form as in legacy implementation.
        inv["sigma_plus"] = (
            0.5 * nrm**2 + np.sqrt(0.25 * nrm**4) + np.abs(detz**2)
        )
        inv["sigma_minus"] = (
            0.5 * nrm**2 - np.sqrt(0.25 * nrm**4) + np.abs(detz**2)
        )

        return inv

    @property
    def z_xx(self) -> np.ndarray:
        if self._z is None:
            raise ZError("Z is not set.")
        return self._z[:, 0, 0]

    @property
    def z_xy(self) -> np.ndarray:
        if self._z is None:
            raise ZError("Z is not set.")
        return self._z[:, 0, 1]

    @property
    def z_yx(self) -> np.ndarray:
        if self._z is None:
            raise ZError("Z is not set.")
        return self._z[:, 1, 0]

    @property
    def z_yy(self) -> np.ndarray:
        if self._z is None:
            raise ZError("Z is not set.")
        return self._z[:, 1, 1]

    @property
    def z_err_xx(self) -> np.ndarray | None:
        return None if self._z_err is None else self._z_err[:, 0, 0]

    @property
    def z_err_xy(self) -> np.ndarray | None:
        return None if self._z_err is None else self._z_err[:, 0, 1]

    @property
    def z_err_yx(self) -> np.ndarray | None:
        return None if self._z_err is None else self._z_err[:, 1, 0]

    @property
    def z_err_yy(self) -> np.ndarray | None:
        return None if self._z_err is None else self._z_err[:, 1, 1]

    # Backward-compat alias
    remove_ss = remove_static_shift
