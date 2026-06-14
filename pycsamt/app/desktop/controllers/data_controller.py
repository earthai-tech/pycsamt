# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
DataController — load, cache, and filter survey data.

No Qt imports: usable from the Dash web app and from unit tests
without a running QApplication.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd


class DataController:
    """
    Loads EDI / AVG / J files into a Sites collection and exposes
    a filtered pandas DataFrame for the StationModel.

    Parameters
    ----------
    progress_callback : callable(int), optional
        Called with 0–100 during loading so a progress bar can update.
    """

    #: Columns exposed to the StationModel
    STATION_COLUMNS = ["ID", "Line", "Latitude", "Longitude", "Elevation", "N_freq", "Tipper"]

    def __init__(self, progress_callback=None) -> None:
        self._progress_cb = progress_callback
        self._sites = None
        self._df: Optional[pd.DataFrame] = None
        self._station_to_line: dict[str, str] = {}

    # ── Loading ───────────────────────────────────────────────────────

    def load(
        self,
        paths: List[str | Path],
        path_to_line: Optional[Dict[str, str]] = None,
    ) -> object:
        """
        Load EDI (or other supported) files and return a ``Sites`` object.

        Parameters
        ----------
        paths : list of str or Path
            File paths to load.
        path_to_line : dict {path_str: line_name}, optional
            Maps each EDI file path to its survey-line name.  When supplied,
            a ``Line`` column is populated in :attr:`dataframe`.
        """
        from pycsamt.seg.edi import EDIFile
        from pycsamt.site.base import Sites

        paths = [Path(p) for p in paths]

        # Build stem → line lookup (EDIFile.station == Path.stem typically)
        self._station_to_line = {
            Path(p).stem: line
            for p, line in (path_to_line or {}).items()
        }

        edis = []
        total = len(paths)
        for i, p in enumerate(paths):
            edis.append(EDIFile(p))
            if self._progress_cb is not None:
                self._progress_cb(int((i + 1) / total * 90))

        self._sites = Sites(edis)
        self._df = self._build_dataframe()

        if self._progress_cb is not None:
            self._progress_cb(100)

        return self._sites

    # ── DataFrame ─────────────────────────────────────────────────────

    def _build_dataframe(self) -> pd.DataFrame:
        """Build a summary DataFrame from the loaded Sites."""
        if self._sites is None:
            return pd.DataFrame(columns=self.STATION_COLUMNS)

        rows: List[Dict] = []
        for edi in self._sites.as_list():
            # sites.as_list() returns EDIFile objects; Sites.get() returns a
            # proper Site whose .summary() yields a dict with lat/lon/name/etc.
            try:
                site = self._sites.get(edi.station)
                s = site.summary()
            except Exception:
                # Fallback: read metadata directly from the EDIFile header
                try:
                    head = edi.get_section("head")
                    loc = head.Location
                except Exception:
                    loc = None
                s = {
                    "name":   edi.station,
                    "lat":    loc.lat  if loc is not None else float("nan"),
                    "lon":    loc.lon  if loc is not None else float("nan"),
                    "elev":   loc.elev if loc is not None else float("nan"),
                    "nfreq":  edi.n_freq,
                    "tipper": edi.has_tipper,
                }
            station_name = s.get("name", "")
            rows.append({
                "ID":        station_name,
                "Line":      self._station_to_line.get(station_name, "—"),
                "Latitude":  s.get("lat", float("nan")),
                "Longitude": s.get("lon", float("nan")),
                "Elevation": s.get("elev", float("nan")),
                "N_freq":    s.get("nfreq", 0),
                "Tipper":    bool(s.get("tipper", False)),
            })

        df = pd.DataFrame(rows, columns=self.STATION_COLUMNS)
        df["Latitude"]  = pd.to_numeric(df["Latitude"],  errors="coerce")
        df["Longitude"] = pd.to_numeric(df["Longitude"], errors="coerce")
        df["Elevation"] = pd.to_numeric(df["Elevation"], errors="coerce")
        return df

    # ── Filtered views ────────────────────────────────────────────────

    @property
    def dataframe(self) -> pd.DataFrame:
        """Full summary DataFrame (one row per station)."""
        if self._df is None:
            return pd.DataFrame(columns=self.STATION_COLUMNS)
        return self._df.copy()

    @property
    def sites(self):
        return self._sites

    def filter_by_ids(self, ids: List[str]) -> pd.DataFrame:
        """Return a sub-DataFrame for the given station IDs."""
        df = self.dataframe
        return df[df["ID"].isin(ids)].reset_index(drop=True)

    def station_coords(self) -> pd.DataFrame:
        """Return only ID / Latitude / Longitude columns — used by MapPanel."""
        return self.dataframe[["ID", "Latitude", "Longitude"]].dropna()

    def is_loaded(self) -> bool:
        return self._sites is not None
