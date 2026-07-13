# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Estimate the area of interest for a MARE2DEM survey.

Port of ``m2d_estimateAreaOfInterest.m``.

Returns the recommended y-range and z-range (depth, positive down)
for plotting or mesh generation given the survey geometry.  The rules
mirror the MATLAB logic:

* Single station → 2 × 2 km square centred on the station.
* Profile (dy >> dz) → 10% lateral buffer, depth = 75% of profile width.
* Square geometry → equal buffer on all sides.
"""

from __future__ import annotations

import numpy as np

from ..iotools.emdata import EMDataFile

__all__ = ["estimate_area_of_interest"]


def estimate_area_of_interest(
    em: EMDataFile,
) -> tuple[np.ndarray | None, np.ndarray | None]:
    """Estimate the y and z limits for display or mesh generation.

    Port of ``m2d_estimateAreaOfInterest.m``.

    Parameters
    ----------
    em : EMDataFile
        Survey data file supplying MT, CSEM, and DC receiver /
        transmitter positions.

    Returns
    -------
    ylim : numpy.ndarray, shape (2,) or None
        ``[y_min, y_max]`` recommended profile extent in metres.
        ``None`` when no survey geometry is available.
    zlim : numpy.ndarray, shape (2,) or None
        ``[z_min, z_max]`` recommended depth extent in metres.
        ``None`` when no survey geometry is available.

    Examples
    --------
    >>> from pycsamt.models.mare2dem import read_emdata
    >>> from pycsamt.models.mare2dem.geom.area_of_interest import estimate_area_of_interest
    >>> em = read_emdata("survey.emdata")
    >>> ylim, zlim = estimate_area_of_interest(em)
    >>> ylim
    array([-5500.,  5500.])
    """
    points: list[np.ndarray] = []

    if em.csem is not None:
        if len(em.csem.receivers):
            points.append(em.csem.receivers[:, 1:3])  # y, z
        if len(em.csem.transmitters):
            points.append(em.csem.transmitters[:, 1:3])

    if em.mt is not None and len(em.mt.receivers):
        points.append(em.mt.receivers[:, 1:3])

    if em.dc is not None:
        if len(em.dc.tx_electrodes):
            points.append(em.dc.tx_electrodes[:, 1:3])
        if len(em.dc.rx_electrodes):
            points.append(em.dc.rx_electrodes[:, 1:3])

    if not points:
        return None, None

    survey = np.vstack(points)
    ymin, zmin = survey.min(axis=0)
    ymax, zmax = survey.max(axis=0)
    dy = ymax - ymin
    dz = zmax - zmin

    if dy == 0 and dz == 0:
        # single station — 2 × 2 km box
        return (
            np.array([ymin - 1000.0, ymin + 1000.0]),
            np.array([zmin - 1000.0, zmin + 1000.0]),
        )

    if dy > 2 * dz:
        # elongated profile along y
        ylim = np.array([ymin - dy / 10, ymax + dy / 10])
        zlim = np.array([zmin - dy / 10, zmax + dy * 0.75])

        # CSEM-only: limit max depth to max Tx–Rx range
        if em.mt is None or not len(em.mt.receivers):
            if (
                em.csem is not None
                and len(em.csem.receivers)
                and len(em.csem.transmitters)
                and em.data is not None
                and len(em.data)
            ):
                csem_mask = em.data[:, 0] < 100
                if csem_mask.any():
                    itx = em.data[csem_mask, 2].astype(int) - 1
                    irx = em.data[csem_mask, 3].astype(int) - 1
                    n_tx = len(em.csem.transmitters)
                    n_rx = len(em.csem.receivers)
                    valid = (
                        (itx >= 0) & (itx < n_tx) & (irx >= 0) & (irx < n_rx)
                    )
                    if valid.any():
                        r = np.abs(
                            em.csem.receivers[irx[valid], 1]
                            - em.csem.transmitters[itx[valid], 1]
                        ).max()
                        zlim = np.array([zmin - r / 10, zmax + r])
    else:
        # roughly square geometry
        dr = max(dy, dz)
        if dy > dz:
            ylim = np.array([ymin - dr / 10, ymax + dr / 10])
            dd = dy - dz
            zlim = np.array(
                [zmin - dd / 2 - dr / 10, zmax + dd / 2 + dr / 10]
            )
        else:
            zlim = np.array([zmin - dr / 10, zmax + dr / 10])
            dd = dz - dy
            ylim = np.array(
                [ymin - dd / 2 - dr / 10, ymax + dd / 2 + dr / 10]
            )

    return ylim, zlim
