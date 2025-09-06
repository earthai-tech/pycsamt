# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

from typing import Any, Optional

import numpy as np

from ..exceptions import (
    ZError,
    PhaseError,
    ResistivityError,
)
from ..utils.zmath import (
    propagate_error_polar2rect,
    z_error2r_phi_error,
)
from ..log.logger import get_logger
from .base import BaseEM

logger = get_logger(__name__)


class ResPhase(BaseEM):
    def __init__(
        self,
        z_array: Optional[np.ndarray] = None,
        z_err_array: Optional[np.ndarray] = None,
        freq: Optional[np.ndarray] = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(**kwargs)
        self._z: Optional[np.ndarray] = None
        self._z_err: Optional[np.ndarray] = None

        self._resistivity: Optional[np.ndarray] = None
        self._phase: Optional[np.ndarray] = None

        self._resistivity_err: Optional[np.ndarray] = None
        self._phase_err: Optional[np.ndarray] = None

        if z_array is not None:
            self._z = np.asarray(z_array, dtype=complex)
        if z_err_array is not None:
            self._z_err = np.asarray(z_err_array, dtype=float)
        if freq is not None:
            self.freq = freq

    # -----------------------
    # Core properties (guard)
    # -----------------------
    @property
    def resistivity(self) -> np.ndarray:
        if self._resistivity is None:
            raise ResistivityError(
                "resistivity not computed/attached"
            )
        return self._resistivity

    @resistivity.setter
    def resistivity(self, res_array: np.ndarray) -> None:
        self._resistivity = np.asarray(res_array, dtype=float)

    @property
    def resistivity_err(self) -> Optional[np.ndarray]:
        return self._resistivity_err

    @resistivity_err.setter
    def resistivity_err(
        self, res_err_array: Optional[np.ndarray]
    ) -> None:
        if res_err_array is None:
            self._resistivity_err = None
        else:
            self._resistivity_err = np.asarray(
                res_err_array, dtype=float
            )

    @property
    def phase(self) -> np.ndarray:
        if self._phase is None:
            raise PhaseError("phase not computed/attached")
        return self._phase

    @phase.setter
    def phase(self, phase_array: np.ndarray) -> None:
        self._phase = np.asarray(phase_array, dtype=float)

    @property
    def phase_err(self) -> Optional[np.ndarray]:
        return self._phase_err

    @phase_err.setter
    def phase_err(
        self, phase_err_array: Optional[np.ndarray]
    ) -> None:
        if phase_err_array is None:
            self._phase_err = None
        else:
            self._phase_err = np.asarray(
                phase_err_array, dtype=float
            )

    # -----------------------------
    # Forward: Z → (ρ, φ) (+ errs)
    # -----------------------------
    def compute_resistivity_phase(
        self,
        z_array: Optional[np.ndarray] = None,
        z_err_array: Optional[np.ndarray] = None,
        freq: Optional[np.ndarray] = None,
    ) -> None:
        if z_array is not None:
            self._z = np.asarray(z_array, dtype=complex)
        if z_err_array is not None:
            self._z_err = np.asarray(z_err_array, dtype=float)
        if freq is not None:
            self.freq = freq

        if self._z is None or self.freq is None:
            raise ZError(
                "missing Z and/or 'freq' to compute ρ and φ"
            )

        z = np.asarray(self._z, dtype=complex)
        f = np.asarray(self.freq, dtype=float)

        if z.ndim != 3 or z.shape[1:] != (2, 2):
            raise ZError(
                "Z must have shape (n,2,2); got "
                f"{z.shape!r}"
            )
        if f.ndim != 1 or f.shape[0] != z.shape[0]:
            raise ZError(
                "freq must be 1-D with length equal to Z.shape[0]"
            )
        if np.any(f <= 0.0) or not np.all(np.isfinite(f)):
            raise ZError("freq must be finite and > 0")
        if not np.all(np.isfinite(z.real + z.imag)):
            raise ZError("Z must be finite")

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
        if not np.all(np.isfinite(z_err)) or np.any(z_err < 0):
            raise ZError("z_err must be finite and non-negative")

        rho_e = np.zeros_like(self._resistivity, dtype=float)
        phi_e = np.zeros_like(self._phase, dtype=float)

        for k in range(f.size):
            for i in range(2):
                for j in range(2):
                    re = z[k, i, j].real
                    im = z[k, i, j].imag
                    dz = z_err[k, i, j]
                    r_rel, ph_err = z_error2r_phi_error(
                        re, im, dz
                    )
                    rho_e[k, i, j] = (
                        self._resistivity[k, i, j] * r_rel
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
        res_err_array: Optional[np.ndarray] = None,
        phase_err_array: Optional[np.ndarray] = None,
    ) -> None:
        rho = np.asarray(res_array, dtype=float)
        phi = np.asarray(phase_array, dtype=float)
        f = np.asarray(freq, dtype=float)

        if np.iscomplexobj(res_array):
            raise ResistivityError("ρ must be real-valued")
        if np.iscomplexobj(phase_array):
            raise PhaseError("phase must be real-valued")

        if f.ndim != 1 or f.size != rho.shape[0]:
            raise ZError(
                "freq must be 1-D with length equal to ρ.shape[0]"
            )
        if np.any(f <= 0.0) or not np.all(np.isfinite(f)):
            raise ZError("freq must be finite and > 0")
        if not np.all(np.isfinite(rho)) or not np.all(
            np.isfinite(phi)
        ):
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
        if (not np.all(np.isfinite(rho_err))
                or not np.all(np.isfinite(phi_err))):
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
            raise ResistivityError(
                "resistivity not computed/attached"
            )

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
    def res_err_xx(self) -> Optional[np.ndarray]:
        if self._resistivity_err is None:
            return None
        return self._resistivity_err[:, 0, 0]

    @property
    def res_err_xy(self) -> Optional[np.ndarray]:
        if self._resistivity_err is None:
            return None
        return self._resistivity_err[:, 0, 1]

    @property
    def res_err_yx(self) -> Optional[np.ndarray]:
        if self._resistivity_err is None:
            return None
        return self._resistivity_err[:, 1, 0]

    @property
    def res_err_yy(self) -> Optional[np.ndarray]:
        if self._resistivity_err is None:
            return None
        return self._resistivity_err[:, 1, 1]

    @property
    def phase_err_xx(self) -> Optional[np.ndarray]:
        if self._phase_err is None:
            return None
        return self._phase_err[:, 0, 0]

    @property
    def phase_err_xy(self) -> Optional[np.ndarray]:
        if self._phase_err is None:
            return None
        return self._phase_err[:, 0, 1]

    @property
    def phase_err_yx(self) -> Optional[np.ndarray]:
        if self._phase_err is None:
            return None
        return self._phase_err[:, 1, 0]

    @property
    def phase_err_yy(self) -> Optional[np.ndarray]:
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
        return np.array(
            [np.linalg.det(zz) ** 0.5 for zz in self._z]
        )

    @property
    def _zdet_var(self) -> np.ndarray:
        if self._z_err is not None:
            return np.array(
                [
                    abs(np.linalg.det(zzv)) ** 0.5
                    for zzv in self._z_err
                ]
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
            out = np.arcsin(
                np.clip(zv / np.abs(zd), -1.0, 1.0)
            )
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
