# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Generate a plausible synthetic elevation profile for display.

For a synthetic survey with no real field topography, this replaces the ad
hoc ``base + A*sin(x/L) + B*cos(x/L)`` one-liner already duplicated across
several documentation-generation scripts with one tested function. It has
no bearing on any forward solve or inversion -- see
:func:`~pycsamt.topo.section.plot_topo_section` for draping a resulting
profile onto a 2-D section, and note that most 2-D forward solvers
(:class:`~pycsamt.forward.em2d.MT2DForward` in particular) require
receivers at ``z=0`` regardless of what this function returns.
"""

from __future__ import annotations

import numpy as np

__all__ = ["synthetic_elevation_profile"]


def synthetic_elevation_profile(
    chainage_m,
    *,
    base_m: float = 100.0,
    amplitude_m: float = 30.0,
    period_m: float = 900.0,
    secondary_period_m: float = 260.0,
    phase_m: float = 0.0,
) -> np.ndarray:
    """Return a smooth, deterministic elevation profile, metres.

    The shape is two summed sinusoids -- a dominant one at *period_m* and
    a secondary one (0.4x the amplitude) at *secondary_period_m* -- giving
    a single broad rise or fall with smaller-scale undulation rather than
    a perfectly regular wave.

    Parameters
    ----------
    chainage_m : array_like
        Along-profile position(s), metres.
    base_m : float, default 100.0
        Mean elevation, metres.
    amplitude_m : float, default 30.0
        Amplitude of the dominant sinusoid, metres.
    period_m : float, default 900.0
        Period of the dominant sinusoid, metres.
    secondary_period_m : float, default 260.0
        Period of the secondary (0.4x amplitude) sinusoid, metres.
    phase_m : float, default 0.0
        Along-profile shift applied before evaluating the profile --
        two calls with different *phase_m* (everything else equal) give
        related but genuinely different curves, e.g. for two nearby
        survey lines without duplicating this function.

    Returns
    -------
    ndarray
        Elevation, metres, same shape as *chainage_m*.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.topo.synthetic import synthetic_elevation_profile
    >>> chainage = np.linspace(0, 2400, 5)
    >>> np.round(synthetic_elevation_profile(chainage), 1)
    array([112. , 110.5, 128. , 136.9, 101.9])
    """
    x = np.asarray(chainage_m, dtype=float) + phase_m
    return (
        base_m
        + amplitude_m * np.sin(x / period_m)
        + 0.4 * amplitude_m * np.cos(x / secondary_period_m)
    )
