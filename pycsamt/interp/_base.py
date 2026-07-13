# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""ResistivityModel — method-agnostic 2-D EM model container.

Any inversion output (Occam2D, ModEM, AI-based) can be adapted into a
:class:`ResistivityModel` so the rest of :mod:`pycsamt.interp` works
without knowing which code produced the model.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from pycsamt.models.occam2d.results import InversionResult

__all__ = ["ResistivityModel"]


@dataclass
class ResistivityModel:
    """Unified 2-D resistivity model container.

    All spatial arrays use metres; depth is positive downward.

    Parameters
    ----------
    x_centers : ndarray (n_x,)
        Horizontal cell-centre positions along the profile, metres.
    z_centers : ndarray (n_z,)
        Depth cell-centre positions, metres (positive downward).
    rho_2d : ndarray (n_z, n_x)
        :math:`\\log_{10}(\\rho / \\Omega\\mathrm{m})` on the cell grid.
    station_x : ndarray (n_sta,)
        Station positions along the profile, metres.
    station_names : list of str
        Station labels, same length as *station_x*.
    method : str
        Source inversion method tag (``'occam2d'``, ``'modem'``,
        ``'ai'``, ``'generic'``).
    rms : float
        Final RMS misfit from the inversion; ``nan`` if unknown.
    """

    x_centers: np.ndarray
    z_centers: np.ndarray
    rho_2d: np.ndarray
    station_x: np.ndarray = field(default_factory=lambda: np.array([]))
    station_names: list[str] = field(default_factory=list)
    method: str = "generic"
    rms: float = float("nan")

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------

    @property
    def n_x(self) -> int:
        return int(self.rho_2d.shape[1])

    @property
    def n_z(self) -> int:
        return int(self.rho_2d.shape[0])

    @property
    def depth_max(self) -> float:
        return float(self.z_centers[-1]) if len(self.z_centers) else 0.0

    @property
    def profile_length(self) -> float:
        if len(self.x_centers) < 2:
            return 0.0
        return float(self.x_centers[-1] - self.x_centers[0])

    # ------------------------------------------------------------------
    # Column access
    # ------------------------------------------------------------------

    def column_nearest(self, x: float) -> np.ndarray:
        """Return the log10-rho column (n_z,) nearest to position *x*."""
        ix = int(np.argmin(np.abs(self.x_centers - x)))
        return self.rho_2d[:, ix].copy()

    def station_column(self, station_name: str) -> np.ndarray:
        """Return the column nearest to the named station."""
        if station_name not in self.station_names:
            raise KeyError(f"Unknown station: {station_name!r}")
        x = self.station_x[self.station_names.index(station_name)]
        return self.column_nearest(x)

    # ------------------------------------------------------------------
    # Adapters
    # ------------------------------------------------------------------

    @classmethod
    def from_occam2d(cls, result: InversionResult) -> ResistivityModel:
        """Build from an :class:`~pycsamt.models.occam2d.InversionResult`.

        Parameters
        ----------
        result : InversionResult
            A loaded Occam2D post-inversion result.

        Raises
        ------
        ValueError
            If *result* does not contain a valid ``rho_2d`` grid.
        """
        if result.rho_2d is None or result.mesh is None:
            raise ValueError(
                "InversionResult has no rho_2d or mesh — ensure the "
                "workdir contains mesh, model, and iter files."
            )
        mesh = result.mesh
        x_c = (mesh.x_nodes[:-1] + mesh.x_nodes[1:]) / 2.0
        z_c = (mesh.z_nodes[:-1] + mesh.z_nodes[1:]) / 2.0

        sta_x: np.ndarray = np.array([])
        sta_names: list[str] = []
        if result.data is not None:
            sta_x = np.asarray(result.data.offsets, dtype=float)
            sta_names = list(result.data.sites)

        return cls(
            x_centers=x_c,
            z_centers=z_c,
            rho_2d=result.rho_2d.copy(),
            station_x=sta_x,
            station_names=sta_names,
            method="occam2d",
            rms=float(result.final_rms),
        )

    @classmethod
    def from_array(
        cls,
        rho_2d: np.ndarray,
        x_centers: np.ndarray,
        z_centers: np.ndarray,
        *,
        station_x: np.ndarray | None = None,
        station_names: list[str] | None = None,
        method: str = "generic",
        rms: float = float("nan"),
    ) -> ResistivityModel:
        """Build directly from numpy arrays.

        Parameters
        ----------
        rho_2d : ndarray (n_z, n_x)
            :math:`\\log_{10}(\\rho)`.
        x_centers, z_centers : ndarray
            Cell-centre coordinates, metres.
        station_x : ndarray, optional
            Station positions.  Defaults to ``x_centers``.
        station_names : list of str, optional
            Station labels.
        method, rms
            Metadata.
        """
        rho_2d = np.asarray(rho_2d, dtype=float)
        x_c = np.asarray(x_centers, dtype=float)
        z_c = np.asarray(z_centers, dtype=float)
        sta_x = (
            np.asarray(station_x, dtype=float)
            if station_x is not None
            else x_c
        )
        if station_names is None:
            station_names = [f"S{i:03d}" for i in range(len(sta_x))]
        return cls(
            x_centers=x_c,
            z_centers=z_c,
            rho_2d=rho_2d,
            station_x=sta_x,
            station_names=list(station_names),
            method=method,
            rms=float(rms),
        )
