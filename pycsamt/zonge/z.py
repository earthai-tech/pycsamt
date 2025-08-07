# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
Lightweight container for the complex impedance tensor *Z*.

Only the strict tensor bookkeeping lives here; heavy I/O and
higher-level helpers (e.g. rotations, quality metrics) are bolted on
later in the workflow.
"""
from __future__ import annotations

from typing import Sequence, Dict, List, Any
import numpy as np

from .base   import AVGComponentBase
from .utils  import chunk_by_frequency, number_stations
from ..exceptions import ZDataError, ZError

from .tensor import ImpedanceTensor, TensorFactory 

from ..utils.zmath   import z_error_to_rho_phi  
from ..constants import MU_0

__all__ = ["Z"]

class Z(AVGComponentBase):
    """
    Impedance-tensor stack.

    Parameters
    ----------
    z_array : array_like[complex], shape ``(n_freq, 2, 2)`` or
              ``(n_freq, n_station, 2, 2)``
        Complex impedance tensors.  A single ``(2, 2)`` matrix is promoted
        to a stack with ``n_freq = 1``.
    freq    : array_like[float], optional
        Discrete transmit frequencies (Hz).  **Must** match
        ``n_freq``.
    station_ids : Sequence[str] | None, default *None*
        Friendly identifiers (“S00”, “S01”, …).  When omitted they are
        auto-generated in the *kind-2* fashion.
    """

    def __init__(
        self,
        z_array:       Sequence[Any] | np.ndarray | None = None,
        *,
        freq:          Sequence[float] | np.ndarray | None = None,
        station_ids:   Sequence[str] | None = None,
        ) -> None:

        super().__init__()
        self._z: np.ndarray | None        = None
        self.freq: np.ndarray | None      = None
        self.station_ids: List[str]       = list(station_ids) if station_ids \
                                                else []
        self.loc: Dict[str, np.ndarray]   = {}

        if z_array is not None:
            self.read(z_array, freq=freq, station_ids=station_ids)

    def read(
        self,
        z_array:     Sequence[Any] | np.ndarray,
        *,
        freq:        Sequence[float] | np.ndarray | None = None,
        station_ids: Sequence[str]           | None = None,
        ) -> None:
        """
        Load a tensor stack and derive the *loc* mapping.

        Raises
        ------
        ZDataError
            When dimensions between *z_array*, *freq* and *station_ids*
            are inconsistent.
        """
        arr = np.asarray(z_array, dtype=complex)
        if arr.ndim == 2:                       # promote a single matrix
            arr = arr[None, ...]               # → (1, 2, 2)
        if arr.ndim == 3:                       # add station axis
            arr = arr[:, None, ...]            # → (n_freq, 1, 2, 2)
        if arr.shape[-2:] != (2, 2):
            raise ZDataError("Last two dims must be (2, 2)")

        n_freq, n_station = arr.shape[:2]

        # ­--- frequencies
        if freq is None:
            raise ZDataError("`freq` must be supplied the first time")
        freq = np.asarray(freq, dtype=float).ravel()
        if freq.size != n_freq:
            raise ZDataError("len(freq) ≠ n_freq")

        # ­-- station identifiers
        if station_ids is None:
            station_ids, _ = number_stations(n_station, n_freq)
        if len(station_ids) != n_station:
            raise ZDataError("station-id count mismatch")

        # commit
        self._z          = arr
        self.freq        = freq
        self.station_ids = list(station_ids)

        # build fast lookup dict:  (n_freq, 2, 2) per station
        chunks = chunk_by_frequency(arr, n_freq)  # returns list(n_station)
        self.loc = {sid: slab for sid, slab in zip(self.station_ids, chunks)}

   
    def write(self) -> List[str]:
        """Return a flat *repr* mainly for unit-test round-trips."""
        if self._z is None:
            return []
        flat = self._z.reshape(self._z.shape[0], self._z.shape[1], 4)
        return [",".join(f"{v.real:+.6e}{v.imag:+.6e}j" for v in row.ravel())
                for row in flat]


    def __getitem__(self, station: str) -> np.ndarray:
        """Quick access: ``Z['S05']`` -> ``(n_freq, 2, 2)`` tensor stack."""
        return self.loc[station]

    def __len__(self) -> int:
        return 0 if self._z is None else self._z.shape[1]  # station count

    def __repr__(self) -> str:                              # noqa: D401
        if self._z is None:
            return "Z(empty)"
        n_freq, n_station = self._z.shape[:2]
        return f"Z(n_freq={n_freq}, n_station={n_station})"


    @property
    def mag(self) -> np.ndarray:
        """Absolute value ``|Z|`` – shape *(n_freq, n_station, 2, 2)*."""
        if self._z is None:
            raise ZDataError("Z array not loaded yet")
        return np.abs(self._z)

    @property
    def rho(self) -> np.ndarray:
        r"""Apparent resistivity,  
        :math:`\rho_a = \dfrac{|Z|^2}{\mu_0\,\omega}`.

        Returns
        -------
        ndarray
            Same shape as :pyattr:`mag`.
        """
        if self._z is None or self.freq is None:
            raise ZDataError("Z/Freq not initialised")
        if not hasattr(self, "_rho_cache"):
            # broadcasting: (n_freq, 1, 1, 1) against full tensor
            omega               = 2.0 * np.pi * self.freq[:, None, None, None]
            self._rho_cache     = self.mag ** 2 / (MU_0 * omega)
        return self._rho_cache                                  # type: ignore[attr-defined]

    @property
    def phase(self) -> np.ndarray:
        """Phase angle in **degrees** – same shape as :pyattr:`mag`."""
        if self._z is None:
            raise ZDataError("Z array not loaded yet")
        if not hasattr(self, "_phase_cache"):
            self._phase_cache = np.degrees(np.angle(self._z))
        return self._phase_cache                                 # type: ignore[attr-defined]

    def get_station_view(self, station: str) -> Dict[str, np.ndarray]:
        """Return every derived slice for a single *station*."""
        slab = self[station]                 # ← uses __getitem__
        idx  = self.station_ids.index(station)
        return {
            "z"    : slab,
            "mag"  : self.mag[:, idx, ...],
            "rho"  : self.rho[:, idx, ...],
            "phase": self.phase[:, idx, ...],
        }

    # Uncertainty propagation (σ|Z| → σρ, σφ)
    @property
    def z_err(self) -> np.ndarray | None:
        """Absolute σ|Z| values – same shape as :pyattr:`_z`."""
        return getattr(self, "_z_err", None)
    
    @z_err.setter
    def z_err(self, arr: Sequence[Any] | np.ndarray) -> None:
        self._z_err = np.asarray(arr, dtype=float).reshape(self._z.shape)
    
    @property
    def rho_err(self) -> np.ndarray:
        """Relative ρ-error (Δρ/ρ).  Lazily computed from *z_err*."""
        if self.z_err is None:
            raise ZError("Set `z_err` before requesting `rho_err`")
        if not hasattr(self, "_rho_err_cache"):
            rel, _phi = z_error_to_rho_phi(self._z.real,
                                           self._z.imag,
                                           self.z_err)
            self._rho_err_cache = rel.reshape(self._z.shape)
        return self._rho_err_cache                           # type: ignore[attr-defined]
    
    @property
    def phase_err(self) -> np.ndarray:
        """Absolute phase error in **degrees** –  same shape as :pyattr:`phase`."""
        if self.z_err is None:
            raise ZError("Set `z_err` before requesting `phase_err`")
        if not hasattr(self, "_phi_err_cache"):
            _rel, phi = z_error_to_rho_phi(self._z.real,
                                           self._z.imag,
                                           self.z_err)
            self._phi_err_cache = phi.reshape(self._z.shape)
        return self._phi_err_cache                           # type: ignore[attr-defined]

    # Quick-and-dirty QC masks
    def bad_mask(
        self,
        max_rho_err: float = 0.2,
        max_phi_err: float = 5.0
        ) -> np.ndarray:
        """
        Boolean mask of unreliable data points.
    
        A point is flagged *True* whenever
        ``rho_err > max_rho_err`` **OR** ``phase_err > max_phi_err``.
        """
        return (self.rho_err > max_rho_err) | (self.phase_err > max_phi_err)
    
    # I guess since we are dealing with only Z class we, need on tensor z and zerror if exist 
    # and not rho and phase. 
    # therefore we need to mofify the function method TensorFactory.build. 
    # because We wILL ADD thsis two classes on Resistivy class, Phase Class also 
    # so the TensorFactory.build should be more flexible and versatile. 

    def as_tensor(self) -> "ImpedanceTensor":
        """Return the current Z-stack as a 3-D
        ``(n_freq, 2, 2)`` :class:`~pycsamt.zonge.tensor.ImpedanceTensor`.

        The object bundles *z*, |Z|, ρₐ, φ and – if available –
        the propagated one-sigma errors.
        """
        if self._z is None:
            raise ZDataError("Nothing loaded – call `read()` first")

        return TensorFactory.build(
            z      = {"xx": self._z[:, :, 0, 0],
                      "xy": self._z[:, :, 0, 1],
                      "yx": self._z[:, :, 1, 0],
                      "yy": self._z[:, :, 1, 1]},

            z_err  = None if self.z_err is None else {
                      "xx": self.z_err[:, :, 0, 0],
                      "xy": self.z_err[:, :, 0, 1],
                      "yx": self.z_err[:, :, 1, 0],
                      "yy": self.z_err[:, :, 1, 1]},
        )


    @classmethod
    def from_tensor(
        cls,
        tensor: "ImpedanceTensor",
        *,
        freq: Sequence[float],
        station_ids: Sequence[str] | None = None
    ) -> "Z":
        """Instantiate a :class:`Z` from an existing
        :class:`~pycsamt.zonge.factory.ImpedanceTensor`.

        Notes
        -----
        * ``freq`` must have the same length as ``tensor.z`` along
          the first axis.
        * If *station_ids* is omitted a default “S00, S01, …” sequence
          is generated.
        """
        z_obj = cls()
        z_obj.read(
            tensor.z,
            freq        = freq,
            station_ids = station_ids,
        )
        # optional uncertainties
        if tensor.z_err is not None:
            z_obj.z_err = tensor.z_err
        return z_obj
