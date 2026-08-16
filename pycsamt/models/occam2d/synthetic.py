# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Package a real forward-modelled response for :class:`~pycsamt.models.occam2d.InputBuilder`.

There is no public, physically-measured data source behind a synthetic
inversion experiment -- only a forward-modelled
:class:`~pycsamt.forward.em2d.ForwardResponse2D`. :class:`~pycsamt.site.Site`
has no array-based constructor (only ``Site(edi: EDIFile)``, wrapping a real
or in-memory-synthesized EDI file), so building one for pure array data would
mean populating an ``EDIFile``'s ``HEAD``/``Z`` sections by hand for no
benefit when there is no EDI file to round-trip.

:class:`~pycsamt.models.occam2d.data.OccamData.from_edi` -- the method
:class:`~pycsamt.models.occam2d.builder.InputBuilder` actually calls -- reads
through ``_normalise_source``, which accepts any iterable of objects
exposing ``name``, ``coords``, ``freq``, ``rho``, ``phase``, ``rho_err``,
and ``phase_err``, not only a real ``Site``. :class:`SyntheticSite` is the
minimal object that contract actually needs, and :func:`sites_from_response`
builds one per station directly from a real
:class:`~pycsamt.forward.em2d.ForwardResponse2D`.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np

__all__ = ["SyntheticSite", "sites_from_response"]

_M_PER_DEG = 111_195.0


class SyntheticSite:
    """Minimal object accepted by ``OccamData.from_edi`` (via
    :class:`~pycsamt.models.occam2d.builder.InputBuilder`).

    Not a real :class:`pycsamt.site.Site` -- there is no EDI file behind
    forward-modelled array data to justify building one. ``coords`` fakes
    a longitude so :func:`~pycsamt.models.occam2d.data.OccamData`'s
    lat/lon-to-metres offset recovery lands on the real along-profile
    chainage that was passed in.

    Parameters
    ----------
    name : str
        Station label.
    x_m : float
        Along-profile position, metres.
    freq : array_like, shape (n_freq,)
        Frequencies, Hz.
    rho : ndarray, shape (n_freq, 2, 2)
        Apparent resistivity, Ω·m. Only ``[:, 0, 1]`` (TE / Zxy) and
        ``[:, 1, 0]`` (TM / Zyx) are read.
    phase : ndarray, shape (n_freq, 2, 2)
        Phase, degrees, raw (not pre-shifted -- ``OccamData.from_edi``
        applies the conventional TM ``+180`` degree normalisation into
        the first quadrant itself).
    rho_err, phase_err : ndarray, shape (n_freq, 2, 2)
        Per-datum errors, same units as *rho*/*phase*.
    """

    def __init__(self, name, x_m, freq, rho, phase, rho_err, phase_err):
        self.name = name
        self.coords = (0.0, float(x_m) / _M_PER_DEG, 0.0)
        self.freq = np.asarray(freq, dtype=float)
        self.rho = np.asarray(rho, dtype=float)
        self.phase = np.asarray(phase, dtype=float)
        self.rho_err = np.asarray(rho_err, dtype=float)
        self.phase_err = np.asarray(phase_err, dtype=float)


def sites_from_response(
    resp,
    station_x: Sequence[float],
    station_names: Sequence[str],
    *,
    rho_err_frac: float = 0.05,
    phase_err_deg: float = 1.5,
) -> list[SyntheticSite]:
    """Build one :class:`SyntheticSite` per station from a real
    :class:`~pycsamt.forward.em2d.ForwardResponse2D`.

    Parameters
    ----------
    resp : ForwardResponse2D
        Real forward-modelled response, e.g. from
        ``MT2DForward(freqs, grid).run()``. ``resp.rho_a_te``/``rho_a_tm``
        and ``resp.phase_te``/``phase_tm`` must have shape
        ``(n_freq, n_stations)``, matching ``station_x``'s length.
    station_x : sequence of float
        Along-profile position of each station, metres, same order and
        length as *resp*'s station axis.
    station_names : sequence of str
        Station labels, same length as *station_x*.
    rho_err_frac : float, default 0.05
        Assigned relative resistivity error (5%).
    phase_err_deg : float, default 1.5
        Assigned absolute phase error, degrees.

    Returns
    -------
    list of SyntheticSite

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.forward.em2d import MT2DForward
    >>> from pycsamt.forward.grid2d import Grid2D
    >>> from pycsamt.models.occam2d.synthetic import sites_from_response
    >>> grid = Grid2D.halfspace(rho=100.0, nx=20, x_max=2000.0, n_stations=5)
    >>> resp = MT2DForward(np.array([100.0, 10.0]), grid, verbose=False).run()
    >>> sites = sites_from_response(resp, grid.x_stations, ["S0", "S1", "S2", "S3", "S4"])
    >>> len(sites), sites[0].name
    (5, 'S0')
    >>> sites[0].rho.shape
    (2, 2, 2)
    """
    n_freq = len(resp.freqs)
    sites = []
    for i, (x, name) in enumerate(zip(station_x, station_names)):
        rho = np.zeros((n_freq, 2, 2))
        phase = np.zeros((n_freq, 2, 2))
        rho_err = np.zeros((n_freq, 2, 2))
        phase_err = np.zeros((n_freq, 2, 2))
        rho[:, 0, 1] = resp.rho_a_te[:, i]
        phase[:, 0, 1] = resp.phase_te[:, i]
        rho[:, 1, 0] = resp.rho_a_tm[:, i]
        phase[:, 1, 0] = resp.phase_tm[:, i]
        rho_err[:, 0, 1] = rho[:, 0, 1] * rho_err_frac
        rho_err[:, 1, 0] = rho[:, 1, 0] * rho_err_frac
        phase_err[:, 0, 1] = phase_err_deg
        phase_err[:, 1, 0] = phase_err_deg
        sites.append(SyntheticSite(name, x, resp.freqs, rho, phase, rho_err, phase_err))
    return sites
