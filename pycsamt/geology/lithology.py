# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Lithology — resistivity-to-geology classification for EM methods.

The built-in :class:`RockDatabase` covers the resistivity ranges
encountered in MT, AMT, CSAMT, and CSEM surveys. Its entries live in
:mod:`pycsamt.geology.rock_library` (not in this module, so the table can
grow independently of the classification logic here), drawing on:

* Palacky (1988) — resistivity of geological formations
* Telford, Geldart & Sheriff (1990) — Applied Geophysics
* Keller & Frischknecht (1966) — Electrical Methods in Geophysical
  Prospecting
* Slichter & Telkes (1942) — electrical properties of minerals
* Updated ranges for fracture zones, aquifers, and economic targets

:class:`StratigraphicLog` packages a per-station depth profile with
lithology assignments and can be exported to CSV, LAS 2.0, or Oasis
Montaj XYZ format via :mod:`pycsamt.interp.export`.

Example
-------
>>> db = RockDatabase.default()
>>> db.classify(250.0)
RockEntry(name='Granite (weathered)', rho_min=50, rho_max=2000, ...)
>>> log = StratigraphicLog.from_column(
...     "S17", x=1050.0, z_centers=z, rho_log10=col, db=db
... )
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Union

import numpy as np

from ..api.property import MetadataMixin, PyCSAMTObject
from .rock_library import BUILTIN_ROCKS

if TYPE_CHECKING:
    from .rock_providers import RockPropertyProvider

__all__ = ["RockEntry", "RockDatabase", "Layer", "StratigraphicLog"]

PathLike = Union[str, Path]


# ---------------------------------------------------------------------------
# RockEntry
# ---------------------------------------------------------------------------


@dataclass
class RockEntry(PyCSAMTObject):
    """A single entry in the rock physics database.

    Parameters
    ----------
    name : str
        Geological unit / lithology name.
    rho_min, rho_max : float
        Resistivity range in Ω·m (linear scale).
    color : str
        Hex colour code for plotting (e.g. ``'#28B463'``).
    description : str
        Optional free-text note.
    code : int
        Integer code used in LAS exports.
    source : str
        Optional literature citation for this range (e.g.
        ``'Palacky (1988)'``). Empty when unspecified, as for entries
        loaded from a CSV without a ``source`` column.
    """

    name: str
    rho_min: float
    rho_max: float
    color: str = "#AAAAAA"
    description: str = ""
    code: int = 0
    source: str = ""

    @property
    def rho_mid(self) -> float:
        return float(np.sqrt(self.rho_min * self.rho_max))

    @property
    def log_rho_mid(self) -> float:
        return float(np.log10(self.rho_mid))

    def contains(self, rho_ohm_m: float) -> bool:
        return self.rho_min <= rho_ohm_m <= self.rho_max


# ---------------------------------------------------------------------------
# RockDatabase
# ---------------------------------------------------------------------------


class RockDatabase(PyCSAMTObject, MetadataMixin):
    """Extensible rock physics database for EM resistivity interpretation.

    Parameters
    ----------
    entries : list of RockEntry
        The database entries.  ``default()`` returns the built-in set.

    Example
    -------
    >>> db = RockDatabase.default()
    >>> db.classify(180.0).name
    'Granite (weathered)'
    """

    def __init__(
        self,
        entries: list[RockEntry],
        *,
        metadata: dict | None = None,
    ) -> None:
        self._entries: list[RockEntry] = list(entries)
        # Pre-compute log-midpoints for fast nearest-match
        self._log_mids = np.array([e.log_rho_mid for e in self._entries])
        self.metadata: dict = dict(metadata) if metadata else {}

    @property
    def entries(self) -> tuple[RockEntry, ...]:
        """Read-only view of the database entries, in insertion order."""
        return tuple(self._entries)

    # ------------------------------------------------------------------
    # Factory
    # ------------------------------------------------------------------

    @classmethod
    def default(cls) -> RockDatabase:
        """Return a database pre-loaded with the built-in rock entries.

        Entries come from :data:`pycsamt.geology.rock_library.BUILTIN_ROCKS`;
        see that module for the literature the ranges are drawn from.
        """
        entries = [
            RockEntry(
                name=r[0],
                rho_min=r[1],
                rho_max=r[2],
                color=r[3],
                description=r[4],
                code=r[5],
                source=r[6] if len(r) > 6 else "",
            )
            for r in BUILTIN_ROCKS
        ]
        return cls(entries, metadata={"origin": "default"})

    @classmethod
    def from_csv(cls, path: PathLike) -> RockDatabase:
        """Load from a CSV file.

        Required columns: ``name, rho_min, rho_max``
        Optional columns: ``color, description, code, source``
        """
        p = Path(path)
        entries: list[RockEntry] = []
        with p.open(newline="") as fh:
            reader = csv.DictReader(fh)
            for i, row in enumerate(reader):
                entries.append(
                    RockEntry(
                        name=row.get("name", f"Rock_{i}").strip(),
                        rho_min=float(row["rho_min"]),
                        rho_max=float(row["rho_max"]),
                        color=row.get("color", "#AAAAAA").strip(),
                        description=row.get("description", "").strip(),
                        code=int(row["code"]) if "code" in row else i + 1,
                        source=row.get("source", "").strip(),
                    )
                )
        return cls(
            entries,
            metadata={"origin": "csv", "path": str(p)},
        )

    @classmethod
    def from_provider(cls, provider: RockPropertyProvider) -> RockDatabase:
        """Build a database from any :class:`RockPropertyProvider`.

        This is the generic entry point behind :meth:`from_url`; call it
        directly when *provider* is something other than a plain URL
        fetch, such as a provider pre-configured with authentication or a
        provider composing several sources.
        """
        entries, metadata = provider.fetch()
        return cls(list(entries), metadata=dict(metadata))

    @classmethod
    def from_url(
        cls,
        url: str,
        *,
        cache_dir: PathLike | None = None,
        ttl_seconds: float = 86400.0,
        timeout: float = 10.0,
        force: bool = False,
        fallback: bool = True,
    ) -> RockDatabase:
        """Fetch a database from *url*, with local caching and fallback.

        *url* must serve a JSON array of objects using the same fields as
        :class:`RockEntry` (``name``, ``rho_min``, ``rho_max``, and
        optionally ``color``, ``description``, ``code``, ``source``). No
        public rock-resistivity service is bundled or assumed by default
        -- point this at a project- or organisation-controlled endpoint
        (an internal API, a shared JSON file on a file server, and so on).

        Parameters
        ----------
        url : str
            Location to fetch the JSON entry list from.
        cache_dir : path-like, optional
            Override the local cache location. Defaults to
            ``$PYCSAMT_ROCKDB_CACHE`` or ``~/.pycsamt/rock_db``.
        ttl_seconds : float
            Reuse a cached response younger than this many seconds
            instead of re-fetching. Default one day.
        timeout : float
            Network timeout in seconds for the fetch itself.
        force : bool
            Re-fetch even if a fresh cache entry exists.
        fallback : bool
            If the fetch fails (network, timeout, malformed response) and
            no usable cache entry exists, fall back to
            :meth:`default` instead of raising. Set to ``False`` to
            surface the failure instead.

        See Also
        --------
        pycsamt.geology.rock_providers.RemoteRockPropertyProvider
            The provider implementing the fetch/cache/fallback policy
            used here.
        """
        from .rock_providers import RemoteRockPropertyProvider

        provider = RemoteRockPropertyProvider(
            url,
            cache_dir=cache_dir,
            ttl_seconds=ttl_seconds,
            timeout=timeout,
            force=force,
            fallback=fallback,
        )
        return cls.from_provider(provider)

    # ------------------------------------------------------------------
    # Classification
    # ------------------------------------------------------------------

    def classify(
        self,
        rho_ohm_m: float,
        method: str = "nearest",
    ) -> RockEntry:
        """Return the best-matching rock entry for *rho_ohm_m*.

        Parameters
        ----------
        rho_ohm_m : float
            Resistivity in Ω·m (linear).
        method : {'nearest', 'overlap'}
            ``'nearest'``: log-distance to midpoint.
            ``'overlap'``: first entry whose range brackets *rho_ohm_m*.
        """
        if np.isnan(rho_ohm_m) or rho_ohm_m <= 0:
            return self._entries[0]

        if method == "overlap":
            for e in self._entries:
                if e.contains(rho_ohm_m):
                    return e

        # nearest in log-space
        log_rho = np.log10(max(rho_ohm_m, 1e-6))
        idx = int(np.argmin(np.abs(self._log_mids - log_rho)))
        return self._entries[idx]

    def classify_column(
        self,
        rho_log10: np.ndarray,
    ) -> list[RockEntry]:
        """Classify every cell in a log10-rho depth column."""
        return [self.classify(float(10.0**v)) for v in rho_log10]

    def __len__(self) -> int:
        return len(self._entries)

    def __repr__(self) -> str:
        return f"RockDatabase({len(self._entries)} entries)"


# ---------------------------------------------------------------------------
# Layer  /  StratigraphicLog
# ---------------------------------------------------------------------------


@dataclass
class Layer(PyCSAMTObject):
    """One geological unit in a pseudo-stratigraphic log.

    Parameters
    ----------
    top, bottom : float
        Depth in metres (positive downward).
    rho_log10 : float
        Representative :math:`\\log_{10}(\\rho)` of the layer.
    lithology : str
        Rock name from :class:`RockDatabase`.
    color : str
        Hex colour for plotting.
    confidence : float
        Fraction of depth cells whose DB classification matches the
        reported lithology (0 – 1).
    """

    top: float
    bottom: float
    rho_log10: float
    lithology: str
    color: str = "#AAAAAA"
    confidence: float = 1.0

    @property
    def thickness(self) -> float:
        return self.bottom - self.top

    @property
    def rho_ohm_m(self) -> float:
        return float(10.0**self.rho_log10)


class StratigraphicLog(PyCSAMTObject):
    """Per-station pseudo-stratigraphic depth profile.

    Constructed from a 1-D log10-rho column and a :class:`RockDatabase`,
    it merges adjacent cells that share the same lithology into discrete
    :class:`Layer` objects.

    Parameters
    ----------
    station_name : str
    station_x : float
    z_centers : ndarray (n_z,)
        Depth cell centres, metres.
    rho_log10 : ndarray (n_z,)
        :math:`\\log_{10}(\\rho)` for each depth cell.
    layers : list of Layer
        Merged geological units (assembled by :meth:`from_column`).
    """

    def __init__(
        self,
        station_name: str,
        station_x: float,
        z_centers: np.ndarray,
        rho_log10: np.ndarray,
        layers: list[Layer],
    ) -> None:
        self.station_name = station_name
        self.station_x = float(station_x)
        self.z_centers = np.asarray(z_centers, dtype=float)
        self.rho_log10 = np.asarray(rho_log10, dtype=float)
        self.layers = layers

    # ------------------------------------------------------------------

    @classmethod
    def from_column(
        cls,
        station_name: str,
        x: float,
        z_centers: np.ndarray,
        rho_log10: np.ndarray,
        db: RockDatabase | None = None,
        *,
        merge_tolerance: float = 0.2,
    ) -> StratigraphicLog:
        """Build a log from a 1-D resistivity column.

        Parameters
        ----------
        station_name : str
        x : float
            Station position, metres.
        z_centers : ndarray (n_z,)
        rho_log10 : ndarray (n_z,)
        db : RockDatabase, optional
            Defaults to :meth:`RockDatabase.default`.
        merge_tolerance : float
            Log10-rho difference threshold for merging adjacent cells
            into one layer (default 0.2 decade).
        """
        if db is None:
            db = RockDatabase.default()

        z = np.asarray(z_centers, dtype=float)
        rho = np.asarray(rho_log10, dtype=float)

        entries = db.classify_column(rho)

        dz = np.diff(z)
        half_dz = np.append(dz / 2, dz[-1] / 2) if len(dz) else np.array([1.0])

        layers: list[Layer] = []
        i = 0
        while i < len(z):
            e0 = entries[i]
            j = i + 1
            while j < len(z):
                if entries[j].name != e0.name:
                    break
                if abs(rho[j] - rho[i]) > merge_tolerance:
                    break
                j += 1
            top = float(z[i] - half_dz[i])
            bottom = float(z[j - 1] + half_dz[j - 1])
            rho_rep = float(np.nanmean(rho[i:j]))
            n_match = sum(1 for k in range(i, j) if entries[k].name == e0.name)
            conf = n_match / max(j - i, 1)
            layers.append(
                Layer(
                    top=max(top, 0.0),
                    bottom=bottom,
                    rho_log10=rho_rep,
                    lithology=e0.name,
                    color=e0.color,
                    confidence=conf,
                )
            )
            i = j

        return cls(station_name, x, z, rho, layers)

    # ------------------------------------------------------------------
    # Export helpers
    # ------------------------------------------------------------------

    def to_dataframe(self):
        """Return layers as a :class:`pandas.DataFrame`."""
        try:
            import pandas as pd
        except ImportError as exc:
            raise ImportError(
                "pandas required for StratigraphicLog.to_dataframe"
            ) from exc
        rows = [
            {
                "station": self.station_name,
                "x": self.station_x,
                "top": ly.top,
                "bottom": ly.bottom,
                "thickness": ly.thickness,
                "rho_log10": ly.rho_log10,
                "rho_ohm_m": ly.rho_ohm_m,
                "lithology": ly.lithology,
                "color": ly.color,
                "confidence": ly.confidence,
            }
            for ly in self.layers
        ]
        return pd.DataFrame(rows)

    def to_dict(self) -> dict:
        return {
            "station_name": self.station_name,
            "station_x": self.station_x,
            "layers": [
                {
                    "top": ly.top,
                    "bottom": ly.bottom,
                    "rho_log10": ly.rho_log10,
                    "lithology": ly.lithology,
                }
                for ly in self.layers
            ],
        }

    def __repr__(self) -> str:
        return (
            f"StratigraphicLog({self.station_name!r}, x={self.station_x:.1f} m, "
            f"{len(self.layers)} layers)"
        )
