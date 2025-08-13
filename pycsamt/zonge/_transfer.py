# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

import math
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import ( 
    Any, 
    Dict, 
    Optional, 
    Sequence, 
    Union
)

import numpy as np
import pandas as pd

from ..log.logger import get_logger 
from .utils import load_avg, to_xarray, AvgDataError

logger = get_logger(__name__)

__all__= ['LegacyAVGBase']

@dataclass
class LegacyAVGBase:
    """
    Base transformer to convert **legacy** Zonge ``*.avg`` tables
    (AMTAVG/MTEdit style) into a **modern** xarray dataset.

    The class:
      • normalizes legacy columns → canonical, lower-case names
      • guarantees minimalist coordinates (station, freq, comp)
      • injects modern fields that may be missing in legacy data
      • derives conservative quantities (e.g., |Z| in ohms)
      • composes survey attrs with safe placeholders

    Subclass this base to customize any hook:
      - ``_normalize_columns``
      - ``_ensure_coords``
      - ``_inject_placeholders``
      - ``_derive_safe_quantities``
      - ``_build_attrs``

    Notes on units
    --------------
    Legacy files are not always explicit about units.  This base
    assumes the **modern CSAVGW** conventions when deriving |Z|:
      - ``E.mag`` in **nV/(m·A)** → converted to **V/(m·A)**
      - ``B.mag`` in **pT/A**     → converted to **T/A**
      - ``|Z| = E/B`` in **ohms**
    If only ``rho`` and ``freq`` exist, we use:
      - ``|Z| = sqrt(μ0 * ω * ρa)``  (SI coherent)
    We **do not** synthesize ``zmag`` (historic "km/sec") to avoid
    mixing conventions.  If ``Z.mag`` exists, it is preserved.

    Resulting Dataset
    -----------------
    The dataset is shaped as:
        ``station × freq × comp``
    and contains all numeric legacy variables plus modern stubs
    (weights, errors, rho_sc, etc.), ready for downstream tools.
    """

    # Magnetic permeability (H/m) in SI.
    mu0: float = 4.0e-7 * math.pi

    # When both routes to |Z| exist, prefer E/B over ρa–ω.
    prefer_emh_over_rho: bool = True

    # Column normalization map (legacy → canonical).
    col_map: Dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Populate a robust default column map."""
        if self.col_map:
            return
        # Keep this compact but comprehensive.  The goal is to
        # tolerate typical legacy headers without being brittle.
        self.col_map = {
            # coordinates
            "Station": "station",
            "Stn": "station",
            "Freq": "freq",
            "Freq.": "freq",
            "Comp": "comp",
            # Tx / fields and phases
            "Amps": "amps",
            "Tx.Amp": "amps",
            "Emag": "emag",
            "E.mag": "emag",
            "Ephz": "ephz",
            "E.phz": "ephz",
            "Hmag": "hmag",
            "Hphz": "hphz",
            "B.mag": "hmag",
            "B.phz": "hphz",
            # impedance, rho, phase
            "Resistivity": "rho",
            "ARes.mag": "rho",
            "Phase": "phase",
            "Z.phz": "phase",
            "Z.mag": "zmag",
            "|Z|": "zabs",
            # errors / weights (legacy + modern labels)
            "%Emag": "e.%err",
            "sEphz": "e.perr",
            "%Hmag": "h.%err",
            "sHphz": "h.perr",
            "%Rho": "rho.%err",
            "sPhz": "z.perr",
            "E.%err": "e.%err",
            "E.perr": "e.perr",
            "B.%err": "h.%err",
            "B.perr": "h.perr",
            "Z.%err": "z.%err",
            "Z.perr": "z.perr",
            "ARes.%err": "rho.%err",
            "Z.mwgt": "z.mwgt",
            "Z.pwgt": "z.pwgt",
            "E.wgt": "e.wgt",
            "H.wgt": "h.wgt",
            # static-corrected rho
            "SRes": "rho_sc",
            "TMARES": "rho_sc",
            "TMARES/SRES": "rho_sc",
            # legacy skip flag
            "skp": "skp",
        }

    # ---------- public API -----------------------------------

    def from_path(
        self,
        path: Union[str, Path],
        *,
        utm_zone: Optional[int] = None,
    ):
        """
        Read a legacy ``*.avg`` file and return a dataset.

        Parameters
        ----------
        path
            Filesystem path to the legacy AVG file.
        utm_zone
            Optional UTM zone override passed to ``load_avg``.

        Returns
        -------
        xarray.Dataset
            Multi-dimensional dataset (station × freq × comp).
        """
        df, meta = load_avg(Path(path), utm_zone=utm_zone)
        return self.from_dataframe(df, meta=meta)

    def from_dataframe(
        self,
        df: pd.DataFrame,
        *,
        meta: Optional[Dict[str, Any]] = None,
        coords: Sequence[str] = ("station", "freq", "comp"),
        data_vars: Optional[Sequence[str]] = None,
    ):
        """
        Transform an already-parsed legacy table into a dataset.

        Parameters
        ----------
        df
            Tidy legacy data as a ``pandas.DataFrame``.
        meta
            Optional metadata dict (e.g., from ``load_avg``).
        coords
            Coordinate columns to use for gridding.
        data_vars
            Optional explicit list of data variables to include.
            If ``None``, all numeric columns not used as coords
            are included.

        Returns
        -------
        xarray.Dataset
            Dataset with conservative placeholders in place.
        """
        meta = dict(meta or {})
        df = self._normalize_columns(df)
        df = self._ensure_coords(df)
        df = self._apply_legacy_skip(df)
        df = self._inject_placeholders(df)
        df = self._derive_safe_quantities(df, meta)

        # Build dataset attrs once data is stabilized.
        attrs = self._build_attrs(df, meta)

        # Choose data variables if not provided.
        if data_vars is None:
            num = df.select_dtypes(include=[np.number]).columns
            data_vars = [c for c in num if c not in coords]

        ds = to_xarray(
            df,
            coords=coords,
            data_vars=data_vars,
            attrs=attrs,
        )
        return ds

    # ---------- hooks / internals -----------------------------

    def _normalize_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Rename legacy / mixed-case columns into canonical names
        and trim whitespace.  No type coercion here.
        """
        if df.empty:
            raise AvgDataError("Empty legacy table.")
        mapping = {c: self.col_map.get(c, c) for c in df.columns}
        out = df.rename(columns=mapping).copy()
        out.columns = [str(c).strip() for c in out.columns]
        return out

    def _ensure_coords(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Guarantee coordinate columns exist and have reasonable
        numeric types.  Insert defaults when absent.
        """
        out = df.copy()

        # Component often missing in kind-1; default to ExHy.
        if "comp" not in out.columns:
            out["comp"] = "ExHy"

        # Coerce station and freq to numeric (keep NaNs).
        for col in ("station", "freq"):
            if col in out.columns:
                out[col] = pd.to_numeric(out[col], errors="coerce")

        # If station looks integer-grid, round it for stability.
        if "station" in out.columns:
            s = out["station"]
            if s.notna().any():
                r = s.round().astype("Int64")
                if (s.dropna() - r.dropna()).abs().max() < 1e-6:
                    out["station"] = r.astype(float)

        return out

    def _apply_legacy_skip(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Convert legacy ``skp`` (0=drop, 1=skip, 2=good) into the
        modern **weight flags**.  If weights already exist, we
        keep them and only add missing ones.
        """
        out = df.copy()

        if "skp" in out.columns:
            skp = pd.to_numeric(out["skp"], errors="coerce")
            # Good → 1.0, else 0.0
            w = (skp == 2).astype(float)
            for col in ("z.mwgt", "z.pwgt"):
                if col not in out.columns:
                    out[col] = w
            # Also publish a convenience boolean "use"
            if "use" not in out.columns:
                out["use"] = w.astype(bool)

        return out

    def _inject_placeholders(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Add modern columns with safe default placeholders so the
        rest of the pipeline can rely on their presence.
        """
        out = df.copy()

        # Weights: default to "use" (1.0) if missing.
        for w in ("z.mwgt", "z.pwgt", "e.wgt", "h.wgt"):
            if w not in out.columns:
                out[w] = 1.0

        # Error columns expected by modern tools.
        for c in (
            "e.%err", "e.perr", "h.%err", "h.perr",
            "z.%err", "z.perr", "rho.%err",
        ):
            if c not in out.columns:
                out[c] = np.nan

        # Static-corrected rho placeholder.
        if "rho_sc" not in out.columns:
            out["rho_sc"] = np.nan

        # Optional coherence and GDP/acq fields.
        for c in ("coh", "gdp_blk", "gdp_chn", "gdp_time"):
            if c not in out.columns:
                out[c] = np.nan

        # Provide '|Z|' slot; fill later when we can.
        if "zabs" not in out.columns:
            out["zabs"] = np.nan

        return out

    def _derive_safe_quantities(
        self,
        df: pd.DataFrame,
        meta: Dict[str, Any],
    ) -> pd.DataFrame:
        """
        Compute conservative derived values only when units are
        unambiguous.  Priority for |Z| derivation:
          1) E/B route (ohms), if both present
          2) ρa–ω route (ohms), if rho and freq present
        If both routes produce values, the ``prefer_emh_over_rho``
        flag decides which to keep.
        """
        out = df.copy()

        # --- Route 1: |Z| from E and B in SI
        z_from_emh = pd.Series(np.nan, index=out.index, dtype=float)
        if {"emag", "hmag"}.issubset(out.columns):
            e_si = pd.to_numeric(out["emag"], errors="coerce") * 1e-9
            b_si = pd.to_numeric(out["hmag"], errors="coerce") * 1e-12
            with np.errstate(divide="ignore", invalid="ignore"):
                z = e_si / b_si
            m = np.isfinite(e_si) & np.isfinite(b_si) & (b_si != 0)
            z_from_emh.loc[m] = z[m]

        # --- Route 2: |Z| from rho and freq in SI
        z_from_rho = pd.Series(np.nan, index=out.index, dtype=float)
        if {"rho", "freq"}.issubset(out.columns):
            rho = pd.to_numeric(out["rho"], errors="coerce")
            f = pd.to_numeric(out["freq"], errors="coerce")
            omega = 2.0 * math.pi * f
            with np.errstate(invalid="ignore"):
                z2 = np.sqrt(self.mu0 * omega * rho)
            z_from_rho = z2

        # Merge routes into zabs without overwriting good values.
        if self.prefer_emh_over_rho:
            primary, secondary = z_from_emh, z_from_rho
        else:
            primary, secondary = z_from_rho, z_from_emh

        # Fill existing zabs only when missing.
        zabs = out["zabs"].copy()
        p_mask = np.isfinite(primary)
        s_mask = np.isfinite(secondary) & zabs.isna() & ~p_mask
        zabs.loc[p_mask] = primary.loc[p_mask]
        zabs.loc[s_mask] = secondary.loc[s_mask]
        out["zabs"] = zabs

        # Do **not** synthesize "zmag" due to unit ambiguity.
        return out

    def _build_attrs(
        self,
        df: pd.DataFrame,
        meta: Dict[str, Any],
    ) -> Dict[str, Any]:
        """
        Compose dataset attributes from file meta and data
        introspection.  Adds conservative defaults for missing
        survey descriptors so downstream consumers have a stable
        header to read.
        """
        attrs: Dict[str, Any] = {}

        # 1) Copy existing meta (omit heavy per-block details).
        for k, v in meta.items():
            if k != "blocks":
                attrs[k] = v

        # 2) Survey descriptors: fill if not provided.
        attrs.setdefault("Survey.Type", "CSAMT")
        attrs.setdefault("Survey.Array", "Scalar")

        # CSAVGW-style units (safe defaults for legacy files).
        attrs.setdefault("Unit.Length", "m")
        attrs.setdefault("Unit.E", "nV/Am")
        attrs.setdefault("Unit.B", "pT/A")
        attrs.setdefault("Unit.Phase", "mrad")

        # 3) Station bounds and increment if grid-like.
        if "station" in df.columns:
            st = df["station"].dropna().unique()
            st.sort()
            if st.size:
                attrs.setdefault("Stn.Left", float(st.min()))
                attrs.setdefault("Stn.Right", float(st.max()))
                attrs.setdefault("Stn.Beg", float(st.min()))
                attrs.setdefault("Stn.GdpBeg", float(st.min()))
                if st.size > 1:
                    diffs = np.diff(st)
                    inc = np.median(diffs)
                    close = np.allclose(
                        diffs, inc, rtol=0.0, atol=1e-6
                    )
                    if close:
                        attrs.setdefault("Stn.Inc", float(inc))
                        attrs.setdefault("Stn.GdpInc", float(inc))

        # 4) Minimal provenance tag.
        now = datetime.now(timezone.utc).strftime(
            "%Y-%m-%dT%H:%M:%SZ"
        )
        hist = attrs.get("history", "")
        tag = f"legacy→xarray:{now}"
        attrs["history"] = f"{hist} | {tag}" if hist else tag

        return attrs
    
    # ---------- public ergonomic adapters --------------------
    def to_xarray(
        self,
        obj: Union[pd.DataFrame, str, Path],
        *,
        meta: Optional[Dict[str, Any]] = None,
        coords: Sequence[str] = ("station", "freq", "comp"),
        data_vars: Optional[Sequence[str]] = None,
        utm_zone: Optional[int] = None,
    ):
        """
        Convert a legacy table (DataFrame) **or** a path to a legacy
        ``*.avg`` file into a modern xarray.Dataset.

        - If ``obj`` is a path/str → reads the file via ``from_path``.
        - If ``obj`` is a DataFrame → forwards to ``from_dataframe``.

        Returns
        -------
        xarray.Dataset
        """
        if isinstance(obj, (str, Path)):
            # Path case: let from_path do the reading.
            return self.from_path(obj, utm_zone=utm_zone)

        if isinstance(obj, pd.DataFrame):
            # DataFrame case: transform in-memory table.
            return self.from_dataframe(
                obj,
                meta=meta,
                coords=coords,
                data_vars=data_vars,
            )

        raise TypeError(
            "to_xarray expects a pandas.DataFrame or a file path; "
            f"got {type(obj)!r}"
        )

    def transform(
        self,
        df: pd.DataFrame,
        *,
        meta: Optional[Dict[str, Any]] = None,
        coords: Sequence[str] = ("station", "freq", "comp"),
        data_vars: Optional[Sequence[str]] = None,
    ):
        """
        Compatibility alias that mirrors scikit-learn-style
        ``transform``. Simply forwards to ``from_dataframe``.
        """
        return self.from_dataframe(
            df,
            meta=meta,
            coords=coords,
            data_vars=data_vars,
        )

    def __call__(self, obj, **kwargs):
        """
        Callable shortcut: ``LegacyAVGBase()(obj, ...)`` behaves like
        ``to_xarray(obj, ...)``.
        """
        return self.to_xarray(obj, **kwargs)
