# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

from typing import Any

import numpy as np

from ..exceptions import (
    PhaseError,
    ResistivityError,
    ZError,
)
from ..log.logger import get_logger
from ..utils.zmath import (
    propagate_error_polar2rect,
    z_error2r_phi_error,
)
from .base import BaseEM

logger = get_logger(__name__)


class ResPhase(BaseEM):
    r"""
    Resistivity/phase container backed by complex **Z**.

    ``ResPhase`` computes apparent resistivity
    :math:`\rho` and phase :math:`\phi` from a complex
    impedance tensor **Z**, and supports the inverse map
    :math:`(\rho,\ \phi) \rightarrow \mathbf{Z}` with error
    propagation.  It is a light, independent container;
    :class:`~pycsamt.z.z.Z` inherits from it to add
    higher-level conveniences.

    Parameters
    ----------
    z_array : ndarray, shape (n_freq, 2, 2), optional
        Complex impedance tensor **Z**.  If omitted, call
        :meth:`compute_resistivity_phase` later.
    z_err_array : ndarray, shape (n_freq, 2, 2), optional
        Absolute per-component uncertainty on **Z**.  If omitted,
        uncertainties on :math:`\rho` and :math:`\phi` remain
        ``None``.
    freq : ndarray, shape (n_freq,), optional
        Frequency vector in Hz.  Must be 1-D, finite and strictly
        positive.
    **kwargs
        Forwarded to :class:`~pycsamt.z.base.BaseEM` (e.g., ``name``,
        ``meta``).

    Attributes
    ----------
    resistivity : ndarray, shape (n_freq, 2, 2)
        Apparent resistivity :math:`\rho` (Ω·m).  Set by
        :meth:`compute_resistivity_phase` or
        :meth:`set_res_phase`.
    phase : ndarray, shape (n_freq, 2, 2)
        Phase :math:`\phi` in degrees.  Set alongside
        :pyattr:`resistivity`.
    resistivity_err : ndarray or None, shape (n_freq, 2, 2)
        Absolute uncertainty on :math:`\rho` (Ω·m) or ``None`` if no
        **Z** errors were provided.
    phase_err : ndarray or None, shape (n_freq, 2, 2)
        Absolute phase uncertainty (deg) or ``None``.
    z : ndarray or None, shape (n_freq, 2, 2)
        Complex **Z**, when known (set by the inverse path).
    z_err : ndarray or None, shape (n_freq, 2, 2)
        Absolute uncertainty on **Z**, when propagated.
    freq : ndarray or None, shape (n_freq,)
        Frequency vector in Hz (1-D, finite, > 0).
    n_freq : int
        Inferred number of frequencies (from :pyattr:`freq` or the
        first dimension of known arrays).

    Notes
    -----
    **Forward path.**  :math:`\rho` and :math:`\phi` are computed
    from **Z** using

    .. math::

       \rho \;=\; 0.2\,\frac{|Z|^{2}}{f},
       \qquad
       \phi \;=\; \angle Z \quad (\text{in degrees}).

    **Error propagation (forward).**  If **Z** errors are available,
    the per-entry relative amplitude error is :math:`\Delta Z/|Z|`.
    We use
    :func:`~pycsamt.utils.zmath.z_error2r_phi_error`
    to map this to the :math:`\rho` relative error (×2) and to an
    absolute phase uncertainty (deg, capped at
    :math:`90^{\circ}`).

    **Inverse path.**  Given :math:`\rho` and :math:`\phi`, we
    recover :math:`|Z|` from

    .. math::

       |Z| \;=\; \sqrt{\,5\,f\,\rho\,},

    then build **Z** in Euler form.  When :math:`\rho` and
    :math:`\phi` errors are given, the :math:`|Z|` error follows

    .. math::

       \frac{d|Z|}{|Z|} \;=\; \tfrac{1}{2}\,\frac{d\rho}{\rho},

    and, together with the phase error, is converted to a single
    absolute **Z** error per component via
    :func:`~pycsamt.utils.zmath.propagate_error_polar2rect`.

    Examples
    --------
    Compute :math:`\rho` and :math:`\phi` from a stack of **Z**::

        >>> import numpy as np
        >>> from pycsamt.z.resphase import ResPhase
        >>> z = np.ones((2, 2, 2), complex)
        >>> f = np.array([10.0, 1.0])
        >>> rp = ResPhase()
        >>> rp.compute_resistivity_phase(z_array=z, freq=f)
        >>> rp.resistivity.shape, rp.phase.shape
        ((2, 2, 2), (2, 2, 2))

    Reconstruct **Z** from :math:`\rho,\ \phi` (no uncertainties)::

        >>> rho = (0.2 / f)[:, None, None] * np.ones((2, 2, 2))
        >>> phi = np.zeros_like(rho)
        >>> rp.set_res_phase(rho, phi, f)
        >>> rp._z.shape
        (2, 2, 2)

    See Also
    --------
    pycsamt.z.z.Z
        High-level impedance container built on ``ResPhase``.
    pycsamt.utils.zmath.z_error2r_phi_error
        Forward error mapping for :math:`\rho` and :math:`\phi`.
    pycsamt.utils.zmath.propagate_error_polar2rect
        Polar → rectangular error propagation for **Z**.

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
        **kwargs: Any,
    ) -> None:
        super().__init__(**kwargs)

        self._z: np.ndarray | None = None
        self._z_err: np.ndarray | None = None

        self._resistivity: np.ndarray | None = None
        self._phase: np.ndarray | None = None

        self._resistivity_err: np.ndarray | None = None
        self._phase_err: np.ndarray | None = None

        if z_array is not None:
            self._z = np.asarray(z_array, dtype=complex)
        if z_err_array is not None:
            self._z_err = np.asarray(z_err_array, dtype=float)
        if freq is not None:
            self.freq = freq

    @property
    def resistivity(self) -> np.ndarray:
        if self._resistivity is None:
            raise ResistivityError("resistivity not computed/attached")
        return self._resistivity

    @resistivity.setter
    def resistivity(self, res_array: np.ndarray) -> None:
        self._resistivity = np.asarray(res_array, dtype=float)

    @property
    def resistivity_err(self) -> np.ndarray | None:
        return self._resistivity_err

    @resistivity_err.setter
    def resistivity_err(self, res_err_array: np.ndarray | None) -> None:
        if res_err_array is None:
            self._resistivity_err = None
        else:
            self._resistivity_err = np.asarray(res_err_array, dtype=float)

    @property
    def phase(self) -> np.ndarray:
        if self._phase is None:
            raise PhaseError("phase not computed/attached")
        return self._phase

    @phase.setter
    def phase(self, phase_array: np.ndarray) -> None:
        self._phase = np.asarray(phase_array, dtype=float)

    @property
    def phase_err(self) -> np.ndarray | None:
        return self._phase_err

    @phase_err.setter
    def phase_err(self, phase_err_array: np.ndarray | None) -> None:
        if phase_err_array is None:
            self._phase_err = None
        else:
            self._phase_err = np.asarray(phase_err_array, dtype=float)

    def compute_resistivity_phase(
        self,
        z_array: np.ndarray | None = None,
        z_err_array: np.ndarray | None = None,
        freq: np.ndarray | None = None,
    ) -> None:
        r"""
        Compute :math:`rho` and :math:`phi`(and their errors)
        from complex **Z**.

        Any provided inputs override the instance state.  On success,
        :pyattr:`resistivity`, :pyattr:`phase`, and, when applicable,
        :pyattr:`resistivity_err` and :pyattr:`phase_err` are set.

        Parameters
        ----------
        z_array : ndarray, shape (n_freq, 2, 2), optional
            Complex impedance tensor **Z**.  If given, it replaces the
            internal value used for the computation.
        z_err_array : ndarray, shape (n_freq, 2, 2), optional
            Absolute **Z** error.  If omitted, ρ and φ uncertainties are
            set to ``None``.
        freq : ndarray, shape (n_freq,), optional
            Frequency vector in Hz.  Must be 1-D, finite and > 0.

        Returns
        -------
        None
            Results are stored on the instance.

        Raises
        ------
        ZError
            If **Z** is missing, shapes are inconsistent, values are not
            finite, or frequencies are not strictly positive.

        Notes
        -----
        We use :math:`\rho = 0.2\,|Z|^2 / f` and
        :math:`\phi = \angle Z` (in degrees).

        If ``z_err_array`` is given, per-entry uncertainties are
        computed via
        :func:`~pycsamt.utils.zmath.z_error2r_phi_error`.  The
        resistivity error is **absolute** (Ω·m).  The phase error is
        **absolute** (deg) and is capped at :math:`90^\circ`.

        Examples
        --------
        >>> import numpy as np
        >>> z = np.ones((2, 2, 2), complex)
        >>> f = np.array([10.0, 1.0])
        >>> rp = ResPhase()
        >>> rp.compute_resistivity_phase(z_array=z, freq=f)
        >>> rp.resistivity.shape
        (2, 2, 2)
        """

        if z_array is not None:
            self._z = np.asarray(z_array, dtype=complex)
        if z_err_array is not None:
            self._z_err = np.asarray(z_err_array, dtype=float)
        if freq is not None:
            self.freq = freq

        if self._z is None or self.freq is None:
            raise ZError("missing Z and/or 'freq' to compute ρ and φ")

        z = np.asarray(self._z, dtype=complex)
        f = np.asarray(self.freq, dtype=float)

        if z.ndim != 3 or z.shape[1:] != (2, 2):
            raise ZError(f"Z must have shape (n,2,2); got {z.shape!r}")
        if f.ndim != 1 or f.shape[0] != z.shape[0]:
            raise ZError("freq must be 1-D with length equal to Z.shape[0]")
        if np.any(f <= 0.0) or not np.all(np.isfinite(f)):
            raise ZError("freq must be finite and > 0")
        abs_z2 = np.abs(z) ** 2
        self._resistivity = 0.2 * abs_z2 / f[:, None, None]
        self._phase = np.degrees(np.angle(z))

        if self._z_err is None:
            self._resistivity_err = None
            self._phase_err = None
            return

        z_err = np.asarray(self._z_err, dtype=float)
        if z_err.shape != z.shape:
            raise ZError(
                "Z error must have same shape as Z; got "
                f"{z_err.shape!r} vs {z.shape!r}"
            )
        valid_z = np.isfinite(z.real) & np.isfinite(z.imag)
        valid_zerr = np.isfinite(z_err) & (z_err >= 0)
        if np.any(~valid_zerr & valid_z):
            raise ZError("z_err must be finite and non-negative")

        rho_e = np.full_like(self._resistivity, np.nan, dtype=float)
        phi_e = np.full_like(self._phase, np.nan, dtype=float)

        for k in range(f.size):
            for i in range(2):
                for j in range(2):
                    if not (valid_z[k, i, j] and valid_zerr[k, i, j]):
                        continue
                    re = z[k, i, j].real
                    im = z[k, i, j].imag
                    dz = z_err[k, i, j]
                    r_rel, ph_err = z_error2r_phi_error(re, im, dz)
                    val = self._resistivity[k, i, j]
                    rho_e[k, i, j] = (
                        0.0
                        if (val == 0.0) or (not np.isfinite(r_rel))
                        else val * r_rel
                    )
                    phi_e[k, i, j] = ph_err

        self._resistivity_err = rho_e
        self._phase_err = phi_e

    # ---------------------------------
    # Inverse: (ρ, φ) → Z (+ Z errs)
    # ---------------------------------
    def set_res_phase(
        self,
        res_array: np.ndarray,
        phase_array: np.ndarray,
        freq: np.ndarray,
        res_err_array: np.ndarray | None = None,
        phase_err_array: np.ndarray | None = None,
    ) -> None:
        r"""
        Attach :math:`rho` and :math:`phi` (with optional errors)
        and reconstruct **Z**.

        This inverse path accepts apparent resistivity (ρ) and phase (φ)
        at each frequency, reconstructs |Z| via ``|Z| = sqrt(5 f ρ)``,
        and builds the complex tensor **Z**.  If both ρ and φ errors are
        supplied, a per-entry absolute **Z** uncertainty is propagated in
        polar coordinates and converted to rectangular form.

        Parameters
        ----------
        res_array : ndarray, shape (n_freq, 2, 2)
            Apparent resistivity (Ω·m).  Must be real and finite.
        phase_array : ndarray, shape (n_freq, 2, 2)
            Phase in degrees.  Must be real and finite.
        freq : ndarray, shape (n_freq,)
            Frequency in Hz (1-D, finite, strictly positive).
        res_err_array : ndarray, shape (n_freq, 2, 2), optional
            Absolute error on ρ (Ω·m).  If omitted, **Z** error is left
            ``None``.
        phase_err_array : ndarray, shape (n_freq, 2, 2), optional
            Absolute phase error in degrees.

        Returns
        -------
        None
            Results are stored on the instance (:pyattr:`_z`,
            :pyattr:`_z_err`, :pyattr:`resistivity`, :pyattr:`phase`).

        Raises
        ------
        ResistivityError
            If ρ contains complex values.
        PhaseError
            If φ contains complex values.
        ZError
            If shapes are inconsistent, values are non-finite, or
            frequencies are not strictly positive.

        Notes
        -----
        The relationship between :math:`|Z|` and :math:`\rho` implies

        ..math::

            \frac{d|Z|}{|Z|} = \tfrac{1}{2}\,\frac{d\rho}{\rho}.

        When error arrays are supplied, |Z| error follows the above and
        is combined with phase error by
        :func:`~pycsamt.utils.zmath.propagate_error_polar2rect` to yield
        a single absolute **Z** error per component.

        Examples
        --------
        >>> import numpy as np
        >>> f = np.array([10.0, 1.0])
        >>> rho = (0.2 / f)[:, None, None] * np.ones((2, 2, 2))
        >>> phi = np.zeros_like(rho)
        >>> rp = ResPhase()
        >>> rp.set_res_phase(rho, phi, f)
        >>> rp._z.shape
        (2, 2, 2)
        """

        rho = np.asarray(res_array, dtype=float)
        phi = np.asarray(phase_array, dtype=float)
        f = np.asarray(freq, dtype=float)

        if np.iscomplexobj(res_array):
            raise ResistivityError("ρ must be real-valued")
        if np.iscomplexobj(phase_array):
            raise PhaseError("phase must be real-valued")

        if f.ndim != 1 or f.size != rho.shape[0]:
            raise ZError("freq must be 1-D with length equal to ρ.shape[0]")
        if np.any(f <= 0.0) or not np.all(np.isfinite(f)):
            raise ZError("freq must be finite and > 0")
        if not np.all(np.isfinite(rho)) or not np.all(np.isfinite(phi)):
            raise ZError("ρ and φ must be finite")

        self._resistivity = rho
        self._phase = phi
        self.freq = f

        self._resistivity_err = (
            None
            if res_err_array is None
            else np.asarray(res_err_array, dtype=float)
        )
        self._phase_err = (
            None
            if phase_err_array is None
            else np.asarray(phase_err_array, dtype=float)
        )

        abs_z = np.sqrt(5.0 * f[:, None, None] * rho)
        self._z = abs_z * np.exp(1j * np.radians(phi))

        self._z_err = None

        if self._resistivity_err is None or self._phase_err is None:
            return

        rho_err = np.asarray(self._resistivity_err, dtype=float)
        phi_err = np.asarray(self._phase_err, dtype=float)

        if rho_err.shape != rho.shape or phi_err.shape != phi.shape:
            raise ZError(
                "error arrays must match shapes of ρ and φ: "
                f"{rho_err.shape!r} vs {rho.shape!r}, "
                f"{phi_err.shape!r} vs {phi.shape!r}"
            )
        if not np.all(np.isfinite(rho_err)) or not np.all(
            np.isfinite(phi_err)
        ):
            raise ZError("ρ_err and φ_err must be finite")
        if np.any(rho_err < 0) or np.any(phi_err < 0):
            raise ZError("ρ_err and φ_err must be non-negative")

        with np.errstate(divide="ignore", invalid="ignore"):
            rel_rho = np.where(rho != 0.0, rho_err / rho, 0.0)
        abs_z_err = 0.5 * abs_z * rel_rho

        z_err = np.zeros_like(self._z, dtype=float)
        for k in range(f.size):
            for i in range(2):
                for j in range(2):
                    ex, ey = propagate_error_polar2rect(
                        abs_z[k, i, j],
                        abs_z_err[k, i, j],
                        phi[k, i, j],
                        phi_err[k, i, j],
                    )
                    z_err[k, i, j] = float(np.hypot(ex, ey))
        self._z_err = z_err

    # ---------------------------
    # Convenience component views
    # ---------------------------
    def _need_rho(self) -> None:
        if self._resistivity is None:
            raise ResistivityError("resistivity not computed/attached")

    def _need_phi(self) -> None:
        if self._phase is None:
            raise PhaseError("phase not computed/attached")

    @property
    def res_xx(self) -> np.ndarray:
        self._need_rho()
        return self._resistivity[:, 0, 0]

    @property
    def res_xy(self) -> np.ndarray:
        self._need_rho()
        return self._resistivity[:, 0, 1]

    @property
    def res_yx(self) -> np.ndarray:
        self._need_rho()
        return self._resistivity[:, 1, 0]

    @property
    def res_yy(self) -> np.ndarray:
        self._need_rho()
        return self._resistivity[:, 1, 1]

    @property
    def phase_xx(self) -> np.ndarray:
        self._need_phi()
        return self._phase[:, 0, 0]

    @property
    def phase_xy(self) -> np.ndarray:
        self._need_phi()
        return self._phase[:, 0, 1]

    @property
    def phase_yx(self) -> np.ndarray:
        self._need_phi()
        return self._phase[:, 1, 0]

    @property
    def phase_yy(self) -> np.ndarray:
        self._need_phi()
        return self._phase[:, 1, 1]

    @property
    def res_err_xx(self) -> np.ndarray | None:
        if self._resistivity_err is None:
            return None
        return self._resistivity_err[:, 0, 0]

    @property
    def res_err_xy(self) -> np.ndarray | None:
        if self._resistivity_err is None:
            return None
        return self._resistivity_err[:, 0, 1]

    @property
    def res_err_yx(self) -> np.ndarray | None:
        if self._resistivity_err is None:
            return None
        return self._resistivity_err[:, 1, 0]

    @property
    def res_err_yy(self) -> np.ndarray | None:
        if self._resistivity_err is None:
            return None
        return self._resistivity_err[:, 1, 1]

    @property
    def phase_err_xx(self) -> np.ndarray | None:
        if self._phase_err is None:
            return None
        return self._phase_err[:, 0, 0]

    @property
    def phase_err_xy(self) -> np.ndarray | None:
        if self._phase_err is None:
            return None
        return self._phase_err[:, 0, 1]

    @property
    def phase_err_yx(self) -> np.ndarray | None:
        if self._phase_err is None:
            return None
        return self._phase_err[:, 1, 0]

    @property
    def phase_err_yy(self) -> np.ndarray | None:
        if self._phase_err is None:
            return None
        return self._phase_err[:, 1, 1]

    # -------------------------
    # Determinant-based metrics
    # -------------------------
    @property
    def _zdet(self) -> np.ndarray:
        if self._z is None:
            raise ZError("Z is not set")
        return np.array([np.linalg.det(zz) ** 0.5 for zz in self._z])

    @property
    def _zdet_var(self) -> np.ndarray:
        if self._z_err is not None:
            return np.array(
                [abs(np.linalg.det(zzv)) ** 0.5 for zzv in self._z_err]
            )
        return np.ones_like(self._zdet, dtype=float)

    @property
    def phase_det(self) -> np.ndarray:
        zd = self._zdet
        return np.degrees(np.angle(zd))

    @property
    def phase_det_err(self) -> np.ndarray:
        zd = self._zdet
        zv = self._zdet_var
        with np.errstate(divide="ignore", invalid="ignore"):
            out = np.arcsin(np.clip(zv / np.abs(zd), -1.0, 1.0))
        return np.degrees(out)

    @property
    def res_det(self) -> np.ndarray:
        zd = self._zdet
        f = np.asarray(self.freq, dtype=float)
        return 0.2 * (np.abs(zd) ** 2) / f

    @property
    def res_det_err(self) -> np.ndarray:
        zd = self._zdet
        zv = self._zdet_var
        f = np.asarray(self.freq, dtype=float)
        return 0.2 * (np.abs(zd + zv) ** 2) / f - self.res_det


# Backward-compat aliases
ResPhase.compute_rho_phi = ResPhase.compute_resistivity_phase
ResPhase.set_rho_phi = ResPhase.set_res_phase
