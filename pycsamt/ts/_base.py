# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Multichannel magnetotelluric time-series container.

:class:`TSData` is the common in-memory form produced by every
reader in :mod:`pycsamt.ts.readers`.  It stores one 1-D float
array per channel (``NaN`` marks missing samples), a single
sampling interval, and the acquisition metadata needed later to
synthesize EDI header sections (station, coordinates, azimuths,
dipole lengths, units).

It intentionally mirrors the spirit of
:class:`pycsamt.seg.time_series.TimeSeries` (the EDI ``>TSERIES``
container) but is designed for *field* recordings: fixed dt for
all channels, missing-data awareness, and rich site metadata.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..log.logger import get_logger
from ..z.base import BaseEM

logger = get_logger(__name__)

__all__ = ["TSData", "CHANNEL_ORDER"]

#: Canonical channel ordering used across :mod:`pycsamt.ts`
#: (matches the SEG EDI recommendation Hx, Hy, Hz, Ex, Ey and
#: the LiMS/EMSLAB multiplex order).  Remote channels follow.
CHANNEL_ORDER: Tuple[str, ...] = (
    "HX", "HY", "HZ", "EX", "EY", "RHX", "RHY",
)


def _norm_cid(cid: str) -> str:
    """Normalize a channel id to upper-case letters/digits."""
    return "".join(
        c for c in str(cid).upper() if c.isalnum()
    )


class TSData(BaseEM):
    r"""
    Container for multichannel MT field time series.

    Parameters
    ----------
    data : dict of str to ndarray, optional
        Mapping ``channel -> samples`` (1-D float).  Missing
        samples are ``NaN``.  Channel keys are normalized to
        upper case (``"hx"`` → ``"HX"``).
    dt : float, optional
        Sampling interval in seconds (common to all channels).
    name : str, optional
        Display name (defaults to ``station`` when set).
    verbose : int, default 0
        Verbosity forwarded to :class:`~pycsamt.z.base.BaseEM`.
    **meta
        Free-form metadata.  Recognized keys are promoted to
        attributes (``station``, ``lat``, ``lon``, ``elev``,
        ``declination``, ``coordsys``, ``start``, ``stop``);
        anything else lands in :attr:`meta`.

    Attributes
    ----------
    ids : list of str
        Channel identifiers in insertion order.
    data : dict of str to ndarray
        Per-channel samples, ``NaN`` = missing.
    dt : float or None
        Sampling interval (s).
    station : str or None
        Station/site identifier (EDI ``DATAID``).
    lat, lon, elev : float or None
        Geographic coordinates (decimal degrees, meters).
    declination : float or None
        Geomagnetic declination (degrees).
    coordsys : str or None
        Azimuth reference (e.g. ``"MAGNETIC NORTH"``).
    start, stop : str or None
        Acquisition start/end timestamps
        (``YYYY-MM-DD hh:mm:ss``).
    azim : dict of str to float
        Channel azimuths in degrees.
    units : dict of str to str
        Channel units (``"nT"``, ``"mV/km"``).
    gain : dict of str to float
        Channel gains as recorded.
    baseline : dict of str to float
        Baseline (absolute-field offset) per channel if known.
    dipole : dict of str to float
        Electric dipole lengths in meters (keys ``EX``/``EY``).
    sensor : dict of str to str
        Sensor identifiers per channel.
    missing : float or None
        Original missing-data sentinel of the source file
        (samples are already converted to ``NaN`` on read).
    meta : dict
        Everything else (source path, comments, filters, ...).

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.ts import TSData
    >>> ts = TSData(
    ...     data={"HX": np.zeros(8), "EX": np.zeros(8)},
    ...     dt=5.0, station="kap103",
    ... )
    >>> ts.n_chan, ts.n_samples
    (2, 8)

    See Also
    --------
    pycsamt.ts.readers.read_ts
        Build a :class:`TSData` from a field file.
    pycsamt.seg.time_series.TimeSeries
        EDI ``>TSERIES`` container this class can convert to.
    """

    def __init__(
        self,
        data: Optional[Dict[str, np.ndarray]] = None,
        dt: Optional[float] = None,
        *,
        name: Optional[str] = None,
        verbose: int = 0,
        **meta,
    ) -> None:
        station = meta.pop("station", None)
        super().__init__(
            name=name or station, verbose=verbose
        )
        self.ids: List[str] = []
        self.data: Dict[str, np.ndarray] = {}
        self.dt: Optional[float] = (
            float(dt) if dt is not None else None
        )

        self.station: Optional[str] = station
        self.lat: Optional[float] = meta.pop("lat", None)
        self.lon: Optional[float] = meta.pop("lon", None)
        self.elev: Optional[float] = meta.pop("elev", None)
        self.declination: Optional[float] = meta.pop(
            "declination", None
        )
        self.coordsys: Optional[str] = meta.pop(
            "coordsys", None
        )
        self.start: Optional[str] = meta.pop("start", None)
        self.stop: Optional[str] = meta.pop("stop", None)

        self.azim: Dict[str, float] = dict(
            meta.pop("azim", {}) or {}
        )
        self.units: Dict[str, str] = dict(
            meta.pop("units", {}) or {}
        )
        self.gain: Dict[str, float] = dict(
            meta.pop("gain", {}) or {}
        )
        self.baseline: Dict[str, float] = dict(
            meta.pop("baseline", {}) or {}
        )
        self.dipole: Dict[str, float] = dict(
            meta.pop("dipole", {}) or {}
        )
        self.sensor: Dict[str, str] = dict(
            meta.pop("sensor", {}) or {}
        )
        self.missing: Optional[float] = meta.pop(
            "missing", None
        )
        self.meta: Dict[str, object] = dict(meta)

        for cid, arr in (data or {}).items():
            self.add_channel(cid, arr)

    # ------------------------------------------------ basic API
    def add_channel(
        self, cid: str, samples: Sequence[float]
    ) -> "TSData":
        """Register (or replace) a channel array."""
        key = _norm_cid(cid)
        arr = np.asarray(samples, dtype=float).ravel()
        if key not in self.data:
            self.ids.append(key)
        self.data[key] = arr
        return self

    def get(self, cid: str) -> np.ndarray:
        """Return the sample array for channel ``cid``."""
        return self.data[_norm_cid(cid)]

    def channels(self) -> List[str]:
        """Ordered channel identifiers."""
        return list(self.ids)

    @property
    def n_chan(self) -> int:
        return len(self.ids)

    @property
    def n_samples(self) -> int:
        """Length of the longest channel (0 when empty)."""
        if not self.data:
            return 0
        return max(a.size for a in self.data.values())

    @property
    def duration(self) -> Optional[float]:
        """Record length in seconds (``n_samples * dt``)."""
        if self.dt is None or self.n_samples == 0:
            return None
        return float(self.n_samples) * float(self.dt)

    def time(self, cid: Optional[str] = None) -> np.ndarray:
        """Relative time vector (s) for ``cid`` (or longest)."""
        n = (
            self.get(cid).size if cid is not None
            else self.n_samples
        )
        dt = 1.0 if self.dt is None else float(self.dt)
        return np.arange(n, dtype=float) * dt

    def matrix(
        self,
        ids: Optional[Sequence[str]] = None,
    ) -> Tuple[np.ndarray, List[str]]:
        """
        Return samples as a 2-D array ``(n_samples, n_chan)``.

        Channels shorter than the longest one are right-padded
        with ``NaN``.  Returns ``(M, order)``.
        """
        order = (
            [_norm_cid(c) for c in ids]
            if ids is not None else list(self.ids)
        )
        n = max(self.data[c].size for c in order)
        M = np.full((n, len(order)), np.nan, dtype=float)
        for j, c in enumerate(order):
            v = self.data[c]
            M[: v.size, j] = v
        return M, order

    def missing_fraction(self, cid: str) -> float:
        """Fraction of ``NaN`` samples in channel ``cid``."""
        x = self.get(cid)
        if x.size == 0:
            return 1.0
        return float(np.isnan(x).mean())

    def slice(self, start: int, stop: int) -> "TSData":
        """Return a sample-index slice as a new :class:`TSData`."""
        out = self.copy_meta()
        for cid in self.ids:
            out.add_channel(cid, self.data[cid][start:stop])
        return out

    def copy_meta(self) -> "TSData":
        """New empty :class:`TSData` carrying the same metadata."""
        out = TSData(
            dt=self.dt,
            name=self.name,
            verbose=self.verbose,
            station=self.station,
        )
        for attr in (
            "lat", "lon", "elev", "declination", "coordsys",
            "start", "stop", "missing",
        ):
            setattr(out, attr, getattr(self, attr))
        for attr in (
            "azim", "units", "gain", "baseline", "dipole",
            "sensor", "meta",
        ):
            setattr(out, attr, dict(getattr(self, attr)))
        return out

    # ------------------------------------- bridge to seg module
    def to_seg_timeseries(
        self,
        *,
        max_samples: Optional[int] = None,
        step: int = 1,
    ):
        """
        Convert to :class:`pycsamt.seg.time_series.TimeSeries`
        so the record can be embedded in an EDI file as
        ``>=TSERIESSECT`` / ``>TSERIES`` blocks.

        Parameters
        ----------
        max_samples : int, optional
            Keep at most this many samples per channel (head of
            the record after ``step`` decimation).  Guards
            against writing multi-hundred-MB EDIs by accident.
        step : int, default 1
            Plain stride decimation applied first (no filter —
            for previews, not for further processing).
        """
        from ..seg.time_series import TimeSeries

        ts = TimeSeries(name=self.name, verbose=self.verbose)
        dt = (
            1.0 if self.dt is None else float(self.dt) * step
        )
        ts._sect_dt = dt
        ts.sectid = self.station or self.name or "TS"
        for cid in self.ids:
            x = self.data[cid][::step]
            if max_samples is not None:
                x = x[: int(max_samples)]
            # EDI numeric blocks cannot hold NaN: use EMPTY
            x = np.where(np.isfinite(x), x, 1.0e32)
            ts.ids.append(cid)
            ts.data[cid] = x
            ts.dt_map[cid] = dt
            ts.npts_map[cid] = int(x.size)
        return ts

    # ------------------------------------------------- dunders
    def __contains__(self, cid: str) -> bool:
        return _norm_cid(cid) in self.data

    def __len__(self) -> int:
        return self.n_samples

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        chan = ",".join(self.ids)
        return (
            f"TSData(station={self.station!r}, chan=[{chan}], "
            f"n={self.n_samples}, dt={self.dt})"
        )
