# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""ResistivityModel — method-agnostic 2-D EM model container.

Any inversion output (Occam2D, ModEM, AI-based) can be adapted into a
:class:`ResistivityModel` so the rest of :mod:`pycsamt.interp` works
without knowing which code produced the model.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np

from ..api.property import MetadataMixin, PyCSAMTObject

if TYPE_CHECKING:
    from pycsamt.models.occam2d.results import InversionResult

__all__ = ["ResistivityModel"]


@dataclass(repr=False)
class ResistivityModel(PyCSAMTObject, MetadataMixin):
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
    metadata : dict
        Free-form provenance (coordinate reference, backend version,
        original file paths, and similar). Empty unless populated
        explicitly or via :meth:`~pycsamt.api.property.MetadataMixin.update_metadata`.
    """

    x_centers: np.ndarray
    z_centers: np.ndarray
    rho_2d: np.ndarray
    station_x: np.ndarray = field(default_factory=lambda: np.array([]))
    station_names: list[str] = field(default_factory=list)
    method: str = "generic"
    rms: float = float("nan")
    metadata: dict[str, Any] = field(default_factory=dict)

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

    def clip_to_stations(self, *, pad_m: float = 150.0) -> ResistivityModel:
        """Return a copy restricted to the station-carrying core.

        A model built through :meth:`from_occam2d` spans Occam2D's own
        wide (often multi-kilometre) boundary-padding columns on both
        sides of the stations -- real, but rarely useful for display or
        for anything that treats ``x_centers`` as "the profile". This
        drops everything outside ``[station_x.min() - pad_m,
        station_x.max() + pad_m]``.

        Parameters
        ----------
        pad_m : float, default 150.0
            Margin kept on each side of the outermost stations, metres.

        Raises
        ------
        ValueError
            If ``station_x`` is empty (nothing to clip around).

        Examples
        --------
        >>> import numpy as np
        >>> from pycsamt.interp import ResistivityModel
        >>> x = np.linspace(-5000.0, 5000.0, 50)
        >>> z = np.array([10.0, 20.0])
        >>> rho = np.zeros((2, 50))
        >>> m = ResistivityModel.from_array(rho, x, z, station_x=np.array([0.0, 100.0]))
        >>> clipped = m.clip_to_stations(pad_m=150.0)
        >>> clipped.x_centers.min() >= -150.0, clipped.x_centers.max() <= 250.0
        (True, True)
        """
        if self.station_x.size == 0:
            raise ValueError(
                "clip_to_stations needs a non-empty station_x to clip around."
            )
        lo = float(self.station_x.min()) - pad_m
        hi = float(self.station_x.max()) + pad_m
        keep = (self.x_centers >= lo) & (self.x_centers <= hi)
        return ResistivityModel.from_array(
            self.rho_2d[:, keep], self.x_centers[keep], self.z_centers,
            station_x=self.station_x, station_names=self.station_names,
            method=self.method, rms=self.rms,
        )

    def clip_to_depth(self, max_depth_m: float) -> ResistivityModel:
        """Return a copy restricted to ``z_centers <= max_depth_m``.

        Solvers that need a numerically distant lower boundary (Occam2D,
        ModEM, MARE2DEM) commonly extend their mesh to tens of
        kilometres depth even when the geological zone of interest is a
        couple of kilometres at most. Plotting the full mesh then
        squeezes that zone into a sliver of the figure. This keeps only
        the shallow cells actually worth displaying.

        Parameters
        ----------
        max_depth_m : float
            Deepest cell centre to keep, metres.

        Raises
        ------
        ValueError
            If no cell centre is at or above ``max_depth_m``.

        Examples
        --------
        >>> import numpy as np
        >>> from pycsamt.interp import ResistivityModel
        >>> x = np.array([0.0, 100.0])
        >>> z = np.array([10.0, 500.0, 9000.0])
        >>> rho = np.zeros((3, 2))
        >>> m = ResistivityModel.from_array(rho, x, z)
        >>> m.clip_to_depth(1000.0).z_centers
        array([ 10., 500.])
        """
        keep = self.z_centers <= max_depth_m
        if not np.any(keep):
            raise ValueError(
                f"No cell centre at or above max_depth_m={max_depth_m} "
                f"(z_centers spans {self.z_centers.min()}-{self.z_centers.max()})."
            )
        return ResistivityModel.from_array(
            self.rho_2d[keep, :], self.x_centers, self.z_centers[keep],
            station_x=self.station_x, station_names=self.station_names,
            method=self.method, rms=self.rms,
        )

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

        Notes
        -----
        ``mesh.x_nodes`` is in the mesh's own coordinate frame, which
        starts at the outer edge of :meth:`~pycsamt.models.occam2d.OccamMesh.from_data`'s
        geometrically-expanding horizontal padding -- not at the first
        station. ``result.data.offsets`` (``station_x`` below), by
        contrast, is real station chainage. Left uncorrected, every
        station's nearest ``x_centers`` column collapses onto the same
        outermost padding cell (``column_nearest``/``station_column``
        silently return the same, meaningless column for every station).
        Padding is built identically on both sides
        (``from_data``'s ``left_pad``/``right_pad`` share one ``pad``
        list), so the station-carrying core is exactly centred in the
        mesh; that symmetry is what lets the shift below be recovered
        from ``x_nodes`` and ``data.offsets`` alone, without assuming a
        specific padding cell count.
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

        if sta_x.size >= 2:
            total_width = float(mesh.x_nodes[-1] - mesh.x_nodes[0])
            core_width = float(sta_x.max() - sta_x.min())
            pad_each_side = (total_width - core_width) / 2.0
            shift = float(sta_x.min()) - (float(mesh.x_nodes[0]) + pad_each_side)
            x_c = x_c + shift

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

    @classmethod
    def from_any(
        cls,
        model: Any,
        *,
        station_x: np.ndarray | None = None,
        station_names: list[str] | None = None,
        unit: str = "m",
    ) -> ResistivityModel:
        """Build from *any* supported 2-D inversion result.

        A single entry point on top of :func:`pycsamt.topo.section.extract_grid`
        (the same "accept anything" adapter :func:`~pycsamt.topo.plot_topo_section`
        uses for topography draping), so callers that just want a
        :class:`ResistivityModel` — for :mod:`pycsamt.interp.export`'s
        Surfer writers, calibration, classification, or anything else in
        this package — do not need to know which of Occam2D, ModEM,
        the backend-neutral :class:`~pycsamt.inversion.results.InversionResult`,
        or an AI agent result actually produced it.

        Parameters
        ----------
        model : object
            Already a :class:`ResistivityModel` (returned as-is), a raw
            ``(x_centers, z_centers, rho_2d)`` triple, a backend-neutral
            :class:`~pycsamt.inversion.results.InversionResult`, a
            native Occam2D or 2-D ModEM ``InversionResult``, or an AI
            agent result exposing ``pred_rho``.
        station_x : ndarray, optional
            Overrides the station positions carried by *model*.
        station_names : list of str, optional
            Overrides the station labels carried by *model*.
        unit : {"m", "km"}
            Unit of *model*'s own coordinates when *model* is a raw
            array triple. Ignored for input forms that carry their own
            unit; the returned model is always in metres.

        Raises
        ------
        TypeError
            If *model* is not a recognised input form.
        NotImplementedError
            For a native MARE2DEM result (unstructured triangular mesh
            — regrid onto a regular grid first).
        ValueError
            For a native 3-D ModEM result (extract a 2-D cut first).

        Examples
        --------
        >>> import numpy as np
        >>> from pycsamt.interp import ResistivityModel
        >>> x = np.array([0.0, 100.0, 200.0])
        >>> z = np.array([10.0, 50.0])
        >>> rho_log10 = np.array([[2.0, 2.1, 2.2], [2.5, 2.6, 2.7]])
        >>> model = ResistivityModel.from_any((x, z, rho_log10))
        >>> model.method, model.n_x, model.n_z
        ('array', 3, 2)
        """
        if isinstance(model, cls):
            return model

        from ..topo.section import extract_grid

        info = extract_grid(
            model, station_x=station_x, station_names=station_names, unit=unit
        )
        to_m = 1000.0 if info.unit == "km" else 1.0
        return cls.from_array(
            info.rho_log10,
            info.x_centers * to_m,
            info.z_centers * to_m,
            station_x=info.station_x * to_m,
            station_names=info.station_names,
            method=info.method,
            rms=info.rms,
        )
