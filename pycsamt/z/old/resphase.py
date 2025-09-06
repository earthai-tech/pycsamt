# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0


"""
Resistivity/phase container built from complex impedances Z.

This module defines :class:`ResPhase`, a light container that
derives apparent resistivity (ρ) and phase (φ) from the complex
impedance tensor Z (and optionally their uncertainties), and
supports the inverse operation (ρ, φ) → Z with error propagation.

The class is intentionally minimal and independent; the higher
level :class:`Z` class should **inherit** from :class:`ResPhase`
and extend behavior as needed.
"""

from __future__ import annotations

from typing import Any, Optional

import numpy as np

from ..exceptions import ZError, PhaseError, ResistivityError
from ..utils.zmath import (
    propagate_error_polar2rect,
    z_error2r_phi_error,
)
from ..log.logger import get_logger
from .base import BaseEM 

logger = get_logger(__name__)


class ResPhase (BaseEM):
    """
    Container for apparent resistivity (ρ) and phase (φ) with Z link.

    Parameters
    ----------
    z_array : ndarray, shape (n_freq, 2, 2), optional
        Complex impedance tensor values Z(ω).
    z_err_array : ndarray, shape (n_freq, 2, 2), optional
        Per-component uncertainty on Z (same shape as ``z_array``).
        Values are expected as **absolute** errors on complex Z.
    freq : ndarray, shape (n_freq,), optional
        Frequency vector in Hz, one value per row in Z.
    **kwargs
        Additional attributes attached verbatim to the instance.

    Notes
    -----
    Apparent resistivity is computed as::

        ρ = 0.2 * |Z|^2 / f

    which is equivalent to :math:`|Z| = \\sqrt{5 f \\rho}` used for
    the inverse transform when reconstructing Z from (ρ, φ).
    """

    # --------------------------
    # Construction & attributes
    # --------------------------
    def __init__(
        self,
        z_array: Optional[np.ndarray] = None,
        z_err_array: Optional[np.ndarray] = None,
        freq: Optional[np.ndarray] = None,
        **kwargs: Any,
    ) -> None:
        self._z: Optional[np.ndarray] = z_array
        self._z_err: Optional[np.ndarray] = z_err_array

        self._resistivity: Optional[np.ndarray] = None
        self._phase: Optional[np.ndarray] = None

        self._resistivity_err: Optional[np.ndarray] = None
        self._phase_err: Optional[np.ndarray] = None

        self._freq = None if freq is None else np.asarray(freq, dtype=float)

        # Attach any extra metadata passed in kwargs.
        for key, val in kwargs.items():
            setattr(self, key, val)

    # ----------
    # Properties
    # ----------
    @property
    def freq(self) -> Optional[np.ndarray]:
        return self._freq
    
    @freq.setter
    def freq(self, f: Optional[np.ndarray]) -> None:
        if f is None:
            self._freq = None
            return
        f = np.asarray(f, dtype=float)
        if f.ndim != 1 or np.any(f <= 0.0):
            raise ZError("Frequency must be 1-D and strictly positive.")
        self._freq = f
    
    @property
    def resistivity(self) -> Optional[np.ndarray]:
        """Apparent resistivity array (ρ), shape (n_freq, 2, 2)."""
        if self._resistivity is None:
            raise ResistivityError("Resistivity not computed/attached yet.")
        return self._resistivity

    @resistivity.setter
    def resistivity(self, res_array: np.ndarray) -> None:
        self._resistivity = res_array

    @property
    def resistivity_err(self) -> Optional[np.ndarray]:
        """Uncertainty on ρ, shape (n_freq, 2, 2)."""
        return self._resistivity_err

    @resistivity_err.setter
    def resistivity_err(self, res_err_array: Optional[np.ndarray]) -> None:
        self._resistivity_err = res_err_array

    @property
    def phase(self) -> np.ndarray:
        """Phase array (degrees), shape (n_freq, 2, 2)."""
        return self._phase

    @phase.setter
    def phase(self, phase_array: np.ndarray) -> None:
        self._phase = phase_array

    @property
    def phase_err(self) -> Optional[np.ndarray]:
        """Uncertainty on phase (degrees), shape (n_freq, 2, 2)."""
        return self._phase_err

    @phase_err.setter
    def phase_err(self, phase_err_array: Optional[np.ndarray]) -> None:
        self._phase_err = phase_err_array

    # -----------------------------
    # Forward: Z → (ρ, φ) (+ errors)
    # -----------------------------
    def compute_resistivity_phase(
        self,
        z_array: Optional[np.ndarray] = None,
        z_err_array: Optional[np.ndarray] = None,
        freq: Optional[np.ndarray] = None,
    ) -> None:
        """
        Compute (ρ, φ) and their errors from Z (and Z errors).

        Parameters
        ----------
        z_array : ndarray, optional
            Complex impedance tensor, shape ``(n_freq, 2, 2)``.
            If provided, overrides the internal Z.
        z_err_array : ndarray, optional
            Absolute error on Z, same shape. If omitted, error
            outputs are set to zero arrays.
        freq : ndarray, optional
            Frequency vector in Hz. If provided, overrides the
            internal frequency.

        Raises
        ------
        ZError
            If Z or frequency is missing or has inconsistent shape.
        """
        if z_array is not None:
            self._z = z_array
        if z_err_array is not None:
            self._z_err = z_err_array
        if freq is not None:
            self.freq = freq

        if self._z is None or self.freq is None:
            raise ZError("Missing Z and/or 'freq' to compute ρ and φ.")

        z = np.asarray(self._z, dtype=complex)
        f = np.asarray(self.freq, dtype=float)

        if z.ndim != 3 or z.shape[1:] != (2, 2):
            raise ZError(
                "Z must have shape (n_freq, 2, 2); got "
                f"{z.shape!r}."
            )
        if f.ndim != 1 or f.shape[0] != z.shape[0]:
            raise ZError(
                "Frequency must be 1-D with length equal to "
                "Z.shape[0]."
            )
        if np.any(f <= 0.0):
            raise ZError("Frequencies must be strictly positive.")

        # ρ = 0.2 * |Z|^2 / f  (broadcast f along 2×2 plane)
        abs_z2 = np.abs(z) ** 2
        self._resistivity = 0.2 * abs_z2 / f[:, None, None]

        # φ in degrees
        self._phase = np.degrees(np.angle(z))

        # Initialize errors as zeros (float) by default.
        self._resistivity_err = np.zeros_like(self._resistivity, dtype=float)
        self._phase_err = np.zeros_like(self._phase, dtype=float)

        # If no Z error provided, we are done.
        if self._z_err is None:
            return

        z_err = np.asarray(self._z_err, dtype=float)
        if z_err.shape != z.shape:
            raise ZError(
                "Z error must have same shape as Z; got "
                f"{z_err.shape!r} vs {z.shape!r}."
            )

        # Per-frequency/component error propagation using helper that
        # converts |Z| error → relative ρ error and phase error.
        for k in range(f.size):
            for i in range(2):
                for j in range(2):
                    re = z[k, i, j].real
                    im = z[k, i, j].imag
                    dz = z_err[k, i, j]

                    r_rel, phi_err = z_error2r_phi_error(re, im, dz)
                    # ρ_err is relative → scale ρ by r_rel
                    self._resistivity_err[k, i, j] = (
                        self._resistivity[k, i, j] * r_rel
                    )
                    self._phase_err[k, i, j] = phi_err

    # ---------------------------------
    # Inverse: (ρ, φ) → Z (+ Z errors)
    # ---------------------------------
    def set_res_phase(
        self,
        res_array: np.ndarray,
        phase_array: np.ndarray,
        freq: np.ndarray,
        res_err_array: Optional[np.ndarray] = None,
        phase_err_array: Optional[np.ndarray] = None,
    ) -> None:
        """
        Attach (ρ, φ) (and errors) and reconstruct Z with errors.

        Parameters
        ----------
        res_array : ndarray, shape (n_freq, 2, 2)
            Apparent resistivity (Ohm·m).
        phase_array : ndarray, shape (n_freq, 2, 2)
            Phase (degrees).
        freq : ndarray, shape (n_freq,)
            Frequency vector (Hz).
        res_err_array : ndarray, optional
            Uncertainty on ρ (same shape).
        phase_err_array : ndarray, optional
            Uncertainty on phase (degrees, same shape).

        Raises
        ------
        ResistivityError
            If ρ contains non-real values.
        PhaseError
            If phase contains non-real values.
        ZError
            If frequency has invalid shape or non-positive values.
        """
        rho = np.asarray(res_array, dtype=float)
        phi = np.asarray(phase_array, dtype=float)
        f = np.asarray(freq, dtype=float)

        if np.iscomplexobj(res_array):
            raise ResistivityError("Resistivity array must be real.")
        if np.iscomplexobj(phase_array):
            raise PhaseError("Phase array must be real.")

        if f.ndim != 1 or f.size != rho.shape[0]:
            raise ZError(
                "Frequency must be 1-D with length equal to "
                "ρ.shape[0]."
            )
        if np.any(f <= 0.0):
            raise ZError("Frequencies must be strictly positive.")

        self._resistivity = rho
        self._phase = phi
        self.freq = f
        self._resistivity_err = (
            None if res_err_array is None 
            else np.asarray(res_err_array, dtype=float)
        )
        self._phase_err = (
            None if phase_err_array is None 
            else np.asarray(phase_err_array, dtype=float)
        )

        # |Z| from ρ via |Z| = sqrt(5 f ρ)
        abs_z = np.sqrt(5.0 * f[:, None, None] * rho)
        self._z = abs_z * np.exp(1j * np.radians(phi))

        # Z error from (ρ_err, φ_err) using error propagation in polar
        # coordinates. Initialize with zeros; fill when both errors exist.
        self._z_err = np.zeros_like(self._z, dtype=float)

        if self._resistivity_err is None or self._phase_err is None:
            # Nothing else to propagate.
            return

        rho_err = np.asarray(self._resistivity_err, dtype=float)
        phi_err = np.asarray(self._phase_err, dtype=float)

        if rho_err.shape != rho.shape or phi_err.shape != phi.shape:
            raise ZError(
                "Error arrays must match shapes of ρ and φ: "
                f"{rho_err.shape!r} vs {rho.shape!r}, "
                f"{phi_err.shape!r} vs {phi.shape!r}."
            )

        # Relationship between abs(Z) and ρ implies:
        #   abs(Z) ∝ sqrt(ρ)  →  d|Z|/|Z| = 0.5 * dρ/ρ
        with np.errstate(divide="ignore", invalid="ignore"):
            rel_rho = np.where(rho != 0.0, rho_err / rho, 0.0)
        abs_z_err = 0.5 * abs_z * rel_rho  # absolute error on |Z|

        # Propagate polar (|Z|, φ) errors to rectangular Z error.
        for k in range(f.size):
            for i in range(2):
                for j in range(2):
                    self._z_err[k, i, j] = max(
                        propagate_error_polar2rect(
                            abs_z[k, i, j],
                            abs_z_err[k, i, j],
                            phi[k, i, j],
                            phi_err[k, i, j],
                        )
                    )

    # ---------------------------
    # Convenience component views
    # ---------------------------
    @property
    def res_xx(self) -> np.ndarray:
        return self._resistivity[:, 0, 0]

    @property
    def res_xy(self) -> np.ndarray:
        return self._resistivity[:, 0, 1]

    @property
    def res_yx(self) -> np.ndarray:
        return self._resistivity[:, 1, 0]

    @property
    def res_yy(self) -> np.ndarray:
        return self._resistivity[:, 1, 1]

    @property
    def phase_xx(self) -> np.ndarray:
        return self._phase[:, 0, 0]

    @property
    def phase_xy(self) -> np.ndarray:
        return self._phase[:, 0, 1]

    @property
    def phase_yx(self) -> np.ndarray:
        return self._phase[:, 1, 0]

    @property
    def phase_yy(self) -> np.ndarray:
        return self._phase[:, 1, 1]

    @property
    def res_err_xx(self) -> Optional[np.ndarray]:
        return ( 
            None if self._resistivity_err is None 
            else self._resistivity_err[:, 0, 0]
            )

    @property
    def res_err_xy(self) -> Optional[np.ndarray]:
        return ( 
            None if self._resistivity_err is None 
            else self._resistivity_err[:, 0, 1]
            )

    @property
    def res_err_yx(self) -> Optional[np.ndarray]:
        return ( 
            None if self._resistivity_err is None 
            else self._resistivity_err[:, 1, 0]
            )

    @property
    def res_err_yy(self) -> Optional[np.ndarray]:
        return ( 
            None if self._resistivity_err is None 
            else self._resistivity_err[:, 1, 1]
            )

    @property
    def phase_err_xx(self) -> Optional[np.ndarray]:
        return ( 
            None if self._phase_err is None 
            else self._phase_err[:, 0, 0]
            )

    @property
    def phase_err_xy(self) -> Optional[np.ndarray]:
        return ( 
            None if self._phase_err is None 
            else self._phase_err[:, 0, 1]
            )

    @property
    def phase_err_yx(self) -> Optional[np.ndarray]:
        return ( 
            None if self._phase_err is None 
            else self._phase_err[:, 1, 0]
            )

    @property
    def phase_err_yy(self) -> Optional[np.ndarray]:
        return ( 
            None if self._phase_err is None else 
            self._phase_err[:, 1, 1]
            )

    # -------------------------
    # Determinant-based metrics
    # -------------------------
    @property
    def _zdet(self) -> np.ndarray:
        """
        Square-root determinant of Z for each frequency.

        Returns
        -------
        ndarray, shape (n_freq,)
            Complex values :math:`\\sqrt{\\det(Z)}`.
        """
        if self._z is None:
            raise ZError("Z is not set.")
        return np.array([np.linalg.det(zz) ** 0.5 for zz in self._z])

    @property
    def _zdet_var(self) -> np.ndarray:
        """
        Proxy uncertainty on :math:`\\sqrt{\\det(Z)}`.

        Returns
        -------
        ndarray, shape (n_freq,)
            If ``Z_err`` is available, returns
            :math:`|\\sqrt{\\det(Z_{err})}|`. Otherwise, returns a
            vector of ones (same length as ``_zdet``).
        """
        if self._z_err is not None:
            return np.array(
                [abs(np.linalg.det(zzv)) ** 0.5 for zzv in self._z_err]
            )
        return np.ones_like(self._zdet, dtype=float)

    @property
    def phase_det(self) -> np.ndarray:
        """
        Phase of :math:`\\sqrt{\\det(Z)}` in degrees.
        """
        zd = self._zdet
        return np.degrees(np.angle(zd))

    @property
    def phase_det_err(self) -> np.ndarray:
        """
        Approximate error on determinant phase (degrees).
        """
        zd = self._zdet
        zv = self._zdet_var
        with np.errstate(divide="ignore", invalid="ignore"):
            out = np.arcsin(np.clip(zv / np.abs(zd), -1.0, 1.0))
        return np.degrees(out)

    @property
    def res_det(self) -> np.ndarray:
        """
        Apparent resistivity from :math:`\\sqrt{\\det(Z)}`.
        """
        zd = self._zdet
        f = np.asarray(self.freq, dtype=float)
        return 0.2 * (np.abs(zd) ** 2) / f

    @property
    def res_det_err(self) -> np.ndarray:
        """
        Approximate uncertainty on ``res_det``.
        """
        zd = self._zdet
        zv = self._zdet_var
        f = np.asarray(self.freq, dtype=float)
        return 0.2 * (np.abs(zd + zv) ** 2) / f - self.res_det


# -----------------------
# Backward-compat aliases
# -----------------------
# Short method names (optional convenience).
ResPhase.compute_rho_phi = ResPhase.compute_resistivity_phase
ResPhase.set_rho_phi = ResPhase.set_res_phase
