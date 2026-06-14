# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
stratagem.qc
============

Quality-control and frequency-filtering classes for Stratagem AMT surveys.

:class:`QualityController`
    Station-level QC: builds a per-station report (SNR, fraction of good
    frequencies, phase-tensor skew) and flags stations that fall below
    configurable thresholds.  Optionally enriches the report with
    hardware-level stack counts from a :class:`~pycsamt.stratagem.io.StratagemRawReader`.

:class:`FrequencyFilter`
    Frequency-level editing: removes incoherent, low-SNR, or out-of-band
    frequency bins.  When a :class:`~pycsamt.stratagem.io.StratagemRawReader`
    is supplied, hardware-measured zero-stack rows are also masked before
    any statistical criteria are applied.

Both classes delegate their core algorithms to :mod:`pycsamt.emtools.qc`,
:mod:`pycsamt.emtools.frequency`, and
:mod:`pycsamt.emtools.remove_noise` — they add only the Stratagem-specific
wiring (hardware mask alignment, station-index mapping, result persistence).
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd

from ..api.property import MetadataMixin, PyCSAMTObject
from ..emtools._core import _iter_items
from ..emtools.frequency import select_band
from ..emtools.qc import build_qc_table, qc_flags
from ..emtools.remove_noise import mask_incoherent_freqs
from ..exceptions import NotFittedError

__all__ = ["FrequencyFilter", "QualityController"]


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _extract_edis(sites_or_list) -> list:
    """Extract EDIFile objects from a Sites wrapper or a plain list."""
    result = []
    for item in _iter_items(sites_or_list):
        edi = getattr(item, "edi", None)
        if edi is not None and getattr(edi, "Z", None) is not None:
            result.append(edi)
        else:
            result.append(item)
    return result


def _align_hardware_mask(
    raw_freqs: np.ndarray,
    raw_mask: np.ndarray,      # shape (n_raw_freqs,) for one station
    edi_freqs: np.ndarray,
) -> np.ndarray:
    """Map a hardware SNR mask from the raw frequency grid to an EDI grid.

    Uses nearest-neighbour matching in linear frequency space.  A raw mask
    value of False propagates to all EDI frequencies that map to that raw bin.

    Parameters
    ----------
    raw_freqs : ndarray, shape (n_raw,)
    raw_mask  : ndarray of bool, shape (n_raw,)
    edi_freqs : ndarray, shape (n_edi,)

    Returns
    -------
    ndarray of bool, shape (n_edi,)
        True where the nearest raw frequency bin had valid data.
    """
    if raw_freqs.size == 0:
        return np.ones(edi_freqs.size, dtype=bool)

    edi_mask = np.ones(edi_freqs.size, dtype=bool)
    for k, ef in enumerate(edi_freqs):
        nearest = int(np.argmin(np.abs(raw_freqs - ef)))
        edi_mask[k] = bool(raw_mask[nearest])
    return edi_mask


def _apply_freq_mask_to_edi(edi, keep: np.ndarray) -> None:
    """Zero out impedance tensor rows where *keep* is False (in-place).

    Parameters
    ----------
    edi : EDIFile
    keep : ndarray of bool, shape (n_freqs,)
        True = keep this frequency row; False = set to NaN.
    """
    z = getattr(edi.Z, "z", None)
    if z is None:
        return
    mask = ~keep
    if not np.any(mask):
        return
    z2 = z.copy()
    z2[mask] = np.nan
    edi.Z.z = z2

    z_err = getattr(edi.Z, "z_err", None)
    if z_err is not None:
        ze2 = z_err.copy()
        ze2[mask] = np.nan
        try:
            edi.Z.z_err = ze2
        except Exception:
            # z_err setter may trigger rho/phase recomputation on NaN Z;
            # bypass via internal attribute when that fails
            try:
                edi.Z.__dict__["_z_err"] = ze2
            except Exception:
                pass

    tip = getattr(edi.Tip, "tipper", None)
    if tip is not None and tip.shape[0] == z.shape[0]:
        t2 = tip.copy()
        t2[mask] = np.nan
        try:
            edi.Tip.tipper = t2
        except Exception:
            pass


# ---------------------------------------------------------------------------
# QualityController
# ---------------------------------------------------------------------------

class QualityController(PyCSAMTObject, MetadataMixin):
    """Station-level quality-control report for Stratagem AMT surveys.

    Wraps :func:`~pycsamt.emtools.qc.build_qc_table` and
    :func:`~pycsamt.emtools.qc.qc_flags` with optional hardware-level
    enrichment from a :class:`~pycsamt.stratagem.io.StratagemRawReader`.

    Parameters
    ----------
    min_frac_ok : float, default 0.6
        Minimum fraction of valid (non-NaN) impedance rows; stations below
        this are flagged ``low_coverage``.
    min_snr_med : float, default 2.0
        Minimum median SNR; stations below this are flagged ``low_snr``.
    max_skew_med : float, default 6.0
        Maximum median absolute phase-tensor skew angle (°); stations
        exceeding this are flagged ``high_skew``.
    include_skew : bool, default True
        Include phase-tensor skew in the report.  Requires a valid
        impedance tensor.
    verbose : int, default 0

    Attributes
    ----------
    report_ : pandas.DataFrame
        Per-station QC metrics.  Columns: ``station``, ``n_freq``,
        ``n_ok``, ``frac_ok``, ``snr_med``, ``pmin``, ``pmax``,
        and (when ``include_skew=True``) ``skew_med``, ``skew_iqr``.
        When a ``StratagemRawReader`` is supplied to :meth:`fit`, three
        additional columns are appended: ``hw_freqs``,
        ``hw_usable_freqs``, ``hw_coverage``.
    flags_ : pandas.DataFrame
        Per-station flag strings in the ``flags`` column.

    Examples
    --------
    >>> from pycsamt.stratagem import EDIBatch, CoordinateInjector
    >>> from pycsamt.stratagem.qc import QualityController
    >>> batch = EDIBatch("2/2EDI").fit()
    >>> inj = CoordinateInjector(epsg=32649).fit(batch, "2.csv")
    >>> qc = QualityController().fit(inj.edi_objects_)
    >>> qc.report_.head()
    >>> qc.summary()
    """

    __repr_fields__ = (
        "min_frac_ok", "min_snr_med", "max_skew_med", "n_stations_",
    )

    def __init__(
        self,
        *,
        min_frac_ok: float = 0.6,
        min_snr_med: float = 2.0,
        max_skew_med: float = 6.0,
        include_skew: bool = True,
        verbose: int = 0,
    ) -> None:
        self.min_frac_ok = min_frac_ok
        self.min_snr_med = min_snr_med
        self.max_skew_med = max_skew_med
        self.include_skew = include_skew
        self.verbose = verbose

    # ------------------------------------------------------------------
    def fit(
        self,
        edi_objects: list,
        raw_reader: "Optional[StratagemRawReader]" = None,
    ) -> "QualityController":
        """Build the QC report.

        Parameters
        ----------
        edi_objects : list of EDIFile
            Stations to assess.  Typically from
            :attr:`~pycsamt.stratagem.gis_correct.CoordinateInjector.edi_objects_`
            or :attr:`~pycsamt.stratagem.io.EDIBatch.edi_objects_`.
        raw_reader : StratagemRawReader, optional
            When supplied, hardware stack counts and SNR masks are joined
            into :attr:`report_` as extra columns ``hw_freqs``,
            ``hw_usable_freqs``, and ``hw_coverage``.

        Returns
        -------
        self
        """
        self.n_stations_ = len(edi_objects)

        self.report_ = build_qc_table(
            edi_objects,
            include_skew=self.include_skew,
            api=False,
            verbose=self.verbose,
        )

        self.flags_ = qc_flags(
            edi_objects,
            min_frac_ok=self.min_frac_ok,
            min_snr_med=self.min_snr_med,
            max_skew_med=self.max_skew_med,
            verbose=self.verbose,
        )

        # ── optional hardware enrichment ──────────────────────────────
        if raw_reader is not None and hasattr(raw_reader, "snr_mask_"):
            self._enrich_with_hardware(edi_objects, raw_reader)

        if self.verbose:
            n_flagged = int((self.flags_["flags"] != "").sum())
            print(
                f"[QualityController] {len(self.report_)} stations assessed, "
                f"{n_flagged} flagged"
            )

        return self

    def _enrich_with_hardware(self, edi_objects: list, raw_reader) -> None:
        """Append hardware coverage columns to report_."""
        # Use station-number matching rather than index alignment.
        edi_to_raw = raw_reader.match_to_edis(edi_objects)
        hw_freqs, hw_usable, hw_cov = [], [], []

        for i in range(len(edi_objects)):
            if i in edi_to_raw:
                raw_idx = edi_to_raw[i]
                mask = raw_reader.snr_mask_[raw_idx]
                hw_freqs.append(int(raw_reader.n_freqs_))
                hw_usable.append(int(mask.sum()))
                hw_cov.append(float(mask.mean()))
            else:
                hw_freqs.append(None)
                hw_usable.append(None)
                hw_cov.append(None)

        # align on station column if lengths match
        if len(hw_freqs) == len(self.report_):
            self.report_ = self.report_.copy()
            self.report_["hw_freqs"] = hw_freqs
            self.report_["hw_usable_freqs"] = hw_usable
            self.report_["hw_coverage"] = hw_cov

    # ------------------------------------------------------------------
    def summary(self) -> str:
        """Return a compact text summary of the QC results."""
        if not hasattr(self, "report_"):
            raise NotFittedError("Call fit() first.")

        r = self.report_
        f = self.flags_
        n_total = len(r)
        n_flagged = int((f["flags"] != "").sum())

        lines = [
            f"QualityController: {n_total} stations",
            f"  flagged  : {n_flagged} ({100 * n_flagged / max(1, n_total):.0f}%)",
        ]
        if not r.empty:
            lines += [
                f"  frac_ok  : {r['frac_ok'].mean():.2f} mean",
                f"  snr_med  : {r['snr_med'].median():.1f} median",
            ]
            if "skew_med" in r.columns:
                lines.append(f"  skew_med : {r['skew_med'].median():.1f}° median")

        # flag breakdown
        all_flags: Dict[str, int] = {}
        for row_flags in f["flags"].dropna():
            for flag in str(row_flags).split(","):
                flag = flag.strip()
                if flag:
                    all_flags[flag] = all_flags.get(flag, 0) + 1
        if all_flags:
            lines.append("  flag breakdown:")
            for k, v in sorted(all_flags.items(), key=lambda x: -x[1]):
                lines.append(f"    {k}: {v}")

        return "\n".join(lines)

    def flagged_stations(self) -> List[str]:
        """Return station names with at least one QC flag.

        Returns
        -------
        list of str
        """
        if not hasattr(self, "flags_"):
            raise NotFittedError("Call fit() first.")
        mask = self.flags_["flags"].astype(str).str.strip() != ""
        return self.flags_.loc[mask, "station"].tolist()


# ---------------------------------------------------------------------------
# FrequencyFilter
# ---------------------------------------------------------------------------

class FrequencyFilter(PyCSAMTObject):
    """Remove bad frequency bins from Stratagem AMT data.

    Combines three filtering strategies that are applied in order:

    1. **Hardware mask** (optional) — zero-stack rows from raw Stratagem
       files are masked before any statistical analysis.  Requires a
       fitted :class:`~pycsamt.stratagem.io.StratagemRawReader`.
    2. **Band selection** — frequencies outside ``[fmin, fmax]`` are dropped.
    3. **Incoherent-frequency mask** — frequencies that fail the SNR
       threshold across more than ``(1 - min_frac)`` of stations are masked.

    All masking is performed in-place on the ``EDIFile.Z.z`` arrays of
    the supplied objects.  Use ``copy=True`` in :meth:`fit` to avoid
    mutating the originals.

    Parameters
    ----------
    fmin : float, optional
        Lower frequency bound (Hz).  Default: no lower bound.
    fmax : float, optional
        Upper frequency bound (Hz).  Default: no upper bound.
    snr_thresh : float, default 2.5
        Per-station SNR threshold for incoherent-frequency masking.
    min_frac : float, default 0.4
        Minimum fraction of stations that must pass ``snr_thresh`` for
        a frequency to be retained.
    use_hardware_mask : bool, default True
        When a ``raw_reader`` is given to :meth:`fit`, apply the
        hardware SNR mask.
    verbose : int, default 0

    Attributes
    ----------
    edi_objects_ : list of EDIFile
        Filtered EDI objects (in-place modified unless ``copy=True``).
    n_masked_hw_ : int
        Number of (station, frequency) pairs masked by hardware SNR.
    n_masked_stat_ : int
        Number masked by the statistical incoherence criterion.
    n_dropped_band_ : int
        Number of frequency rows removed by band selection.

    Examples
    --------
    >>> filt = FrequencyFilter(fmin=10.0, fmax=10000.0)
    >>> filt.fit(inj.edi_objects_, raw_reader=rdr)
    FrequencyFilter(fmin=10.0, fmax=10000.0, ...)
    >>> paths = filt.out("2/2EDIF")
    """

    __repr_fields__ = (
        "fmin", "fmax", "snr_thresh", "n_masked_hw_", "n_dropped_band_",
    )

    def __init__(
        self,
        *,
        fmin: Optional[float] = None,
        fmax: Optional[float] = None,
        snr_thresh: float = 2.5,
        min_frac: float = 0.4,
        use_hardware_mask: bool = True,
        verbose: int = 0,
    ) -> None:
        self.fmin = fmin
        self.fmax = fmax
        self.snr_thresh = snr_thresh
        self.min_frac = min_frac
        self.use_hardware_mask = use_hardware_mask
        self.verbose = verbose

    # ------------------------------------------------------------------
    def fit(
        self,
        edi_objects: list,
        raw_reader: "Optional[StratagemRawReader]" = None,
        *,
        copy: bool = False,
    ) -> "FrequencyFilter":
        """Apply frequency filters.

        Parameters
        ----------
        edi_objects : list of EDIFile
        raw_reader : StratagemRawReader, optional
            Provides hardware SNR masks aligned to station order.
        copy : bool, default False
            When True, deep-copies the Z data of each EDIFile before
            masking so the originals are not mutated.

        Returns
        -------
        self
        """
        if copy:
            edi_objects = [deepcopy(e) for e in edi_objects]

        self.n_masked_hw_ = 0
        self.n_dropped_band_ = 0
        self.n_masked_stat_ = 0

        # ── 1. hardware mask ──────────────────────────────────────────
        if (
            self.use_hardware_mask
            and raw_reader is not None
            and hasattr(raw_reader, "snr_mask_")
        ):
            raw_freqs = raw_reader.freqs_
            # Use station-number matching so that raw[i] aligns with the
            # correct EDI even when sequences start at different offsets
            # (e.g. raw stations 1-87, EDI files starting from station 2).
            edi_to_raw = raw_reader.match_to_edis(edi_objects)

            for j, edi in enumerate(edi_objects):
                if j not in edi_to_raw:
                    continue
                raw_idx = edi_to_raw[j]
                edi_freqs = getattr(edi.Z, "freq", None)
                if edi_freqs is None or raw_freqs.size == 0:
                    continue

                hw_keep = _align_hardware_mask(
                    raw_freqs, raw_reader.snr_mask_[raw_idx], edi_freqs
                )
                n_bad = int((~hw_keep).sum())
                if n_bad:
                    _apply_freq_mask_to_edi(edi, hw_keep)
                    self.n_masked_hw_ += n_bad

        # ── 2. band selection ─────────────────────────────────────────
        if self.fmin is not None or self.fmax is not None:
            before_counts = [
                int(np.sum(np.isfinite(getattr(e.Z, "z", np.array([])).ravel())))
                for e in edi_objects
            ]
            filtered = select_band(
                edi_objects,
                fmin=self.fmin,
                fmax=self.fmax,
                inplace=True,
                verbose=0,
            )
            after_counts = [
                int(np.sum(np.isfinite(getattr(e.Z, "z", np.array([])).ravel())))
                for e in edi_objects
            ]
            self.n_dropped_band_ = sum(
                max(0, b - a) for b, a in zip(before_counts, after_counts)
            )

        # ── 3. incoherent-frequency mask ──────────────────────────────
        before_nan = sum(
            int(np.sum(~np.isfinite(getattr(e.Z, "z", np.array([])).ravel())))
            for e in edi_objects
        )
        mask_incoherent_freqs(
            edi_objects,
            snr_thresh=self.snr_thresh,
            min_frac=self.min_frac,
            inplace=True,
            verbose=0,
        )
        after_nan = sum(
            int(np.sum(~np.isfinite(getattr(e.Z, "z", np.array([])).ravel())))
            for e in edi_objects
        )
        self.n_masked_stat_ = max(0, after_nan - before_nan)

        self.edi_objects_ = edi_objects

        if self.verbose:
            print(
                f"[FrequencyFilter] hw={self.n_masked_hw_} masked, "
                f"band_drop={self.n_dropped_band_}, "
                f"incoherent={self.n_masked_stat_}"
            )

        return self

    # ------------------------------------------------------------------
    def out(
        self,
        savepath: "Union[str, Path, None]" = None,
        *,
        overwrite: bool = False,
    ) -> "Union[List, List[Path]]":
        """Write filtered EDI files to disk or return objects.

        Parameters
        ----------
        savepath : path-like, optional
            Output directory.  When ``None``, returns the list of
            filtered :class:`~pycsamt.seg.edi.EDIFile` objects instead
            of writing to disk.
        overwrite : bool, default False

        Returns
        -------
        list of EDIFile (when savepath is None) or list of Path
        """
        if not hasattr(self, "edi_objects_"):
            raise NotFittedError("Call fit() first.")

        if savepath is None:
            return self.edi_objects_

        out_dir = Path(savepath).expanduser().resolve()
        out_dir.mkdir(parents=True, exist_ok=True)
        written: List[Path] = []

        for edi in self.edi_objects_:
            fname = (
                edi.path.name if getattr(edi, "path", None) is not None
                else f"{edi.station or 'station'}.edi"
            )
            out_path = out_dir / fname
            if out_path.exists() and not overwrite:
                written.append(out_path)
                continue
            try:
                edi.write(new_edifn=fname, savepath=str(out_dir))
                written.append(out_path)
            except Exception as exc:
                if self.verbose:
                    print(f"[FrequencyFilter] write failed {fname}: {exc}")

        if self.verbose:
            print(f"[FrequencyFilter] wrote {len(written)} files → {out_dir}")

        return written
