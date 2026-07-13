from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

# re-use package editors when it saves code
from ..api.station import PYCSAMT_STATION_RENDERING
from ..api.view import (
    PYCSAMT_API_VIEW,
    maybe_wrap_frame,
    wrap_result,
)
from ._core import (
    _apply_each,
    _get_t_block,
    _get_z_block,
    _iter_items,
    _name,
    ensure_sites,
)

_BACKWARD_SINCE = "2.0.0"
_BACKWARD_REMOVE = "2.17.0"


@dataclass
class FrequencyEditResult:
    """Container returned by confidence-based frequency editing."""

    sites: Any
    report: Any
    decisions: Any
    mode: str
    method: str
    ci_hi: float
    ci_lo: float
    reject: str
    interpolation: str

    @property
    def n_dropped(self) -> int:
        """Total number of dropped station-frequency rows."""
        if getattr(self.decisions, "empty", True):
            return 0
        return int((self.decisions["action"] == "dropped").sum())

    @property
    def n_masked(self) -> int:
        """Total number of masked station-frequency rows."""
        if getattr(self.decisions, "empty", True):
            return 0
        return int((self.decisions["action"] == "masked").sum())

    @property
    def n_recovered(self) -> int:
        """Total number of recovered station-frequency rows."""
        if getattr(self.decisions, "empty", True):
            return 0
        return int((self.decisions["action"] == "recovered").sum())

    def summary(self) -> str:
        """Return a compact text summary of the edit result."""
        return (
            "FrequencyEditResult("
            f"mode={self.mode!r}, method={self.method!r}, "
            f"dropped={self.n_dropped}, masked={self.n_masked}, "
            f"recovered={self.n_recovered})"
        )

    def __repr__(self) -> str:  # noqa: D105
        return self.summary()


# ------------------------------- helpers -------------------------------- #


def _sorted_unique(fr: np.ndarray) -> np.ndarray:
    fr = np.asarray(fr, dtype=float)
    fr = fr[np.isfinite(fr)]
    if fr.size == 0:
        return fr
    fr = np.unique(fr)
    fr = fr[fr > 0.0]
    return np.sort(fr)


def _nearest_idx(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    # for each y, nearest index in x
    idx = np.searchsorted(x, y)
    idx = np.clip(idx, 1, x.size - 1)
    left = x[idx - 1]
    right = x[idx]
    pick_left = (y - left) <= (right - y)
    idx[pick_left] -= 1
    return idx


def _interp_complex(
    x: np.ndarray,
    y: np.ndarray,
    xnew: np.ndarray,
    *,
    method: str = "nearest",
) -> np.ndarray:
    if y.ndim == 1:
        y = y[:, None]
    r = y.real
    im = y.imag
    if method == "nearest":
        idx = _nearest_idx(x, xnew)
        rr = r[idx]
        ii = im[idx]
    else:
        rr = np.vstack(
            [np.interp(xnew, x, r[:, j]) for j in range(r.shape[1])]
        ).T
        ii = np.vstack(
            [np.interp(xnew, x, im[:, j]) for j in range(im.shape[1])]
        ).T
    out = rr + 1j * ii
    return out.squeeze()


def _interp_rows_by_freq(
    values: np.ndarray,
    fr: np.ndarray,
    fill: np.ndarray,
    good: np.ndarray,
    *,
    method: str = "linear",
) -> np.ndarray:
    """Recover selected rows by interpolating finite trusted rows."""
    out_dtype = (
        np.result_type(values.dtype, np.complex128)
        if np.iscomplexobj(values)
        else values.dtype
    )
    out = values.astype(out_dtype, copy=True)
    fr = np.asarray(fr, dtype=float)
    fill = np.asarray(fill, dtype=bool)
    good = np.asarray(good, dtype=bool)
    if values.shape[0] != fr.size:
        return out
    if not np.any(fill):
        return out
    x_good = np.log10(np.maximum(fr[good], 1e-24))
    x_fill = np.log10(np.maximum(fr[fill], 1e-24))
    if x_good.size < 2 or x_fill.size == 0:
        out[fill] = np.nan
        return out
    flat = values.reshape(values.shape[0], -1)
    flat_out = out.reshape(out.shape[0], -1)
    for j in range(flat.shape[1]):
        y = flat[:, j]
        valid = good & np.isfinite(y)
        if np.count_nonzero(valid) < 2:
            flat_out[fill, j] = np.nan
            continue
        interp = _interp_complex(
            np.log10(np.maximum(fr[valid], 1e-24)),
            y[valid],
            x_fill,
            method=method,
        )
        flat_out[fill, j] = (
            interp if np.iscomplexobj(flat_out) else interp.real
        )
    return flat_out.reshape(values.shape)


def _confidence_decision_table(
    sites: Any,
    *,
    method: str,
    ci_hi: float,
    ci_lo: float,
    weights: dict[str, float] | None,
    recursive: bool,
    on_dup: str,
    strict: bool,
    verbose: int,
):
    from .qc import frequency_confidence_table

    return frequency_confidence_table(
        sites,
        method=method,
        weights=weights,
        ci_hi=ci_hi,
        ci_lo=ci_lo,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )


def _apply_station_rendering(
    ax: Any,
    positions: Sequence[float],
    labels: Sequence[Any],
    *,
    station_label_step: int | None,
    station_preset: str,
    station_style: Any | None,
) -> None:
    """Apply the package station-axis API to frequency-edit plots."""
    import copy

    from pycsamt.api.station import PYCSAMT_STATION_RENDERING

    style = station_style or PYCSAMT_STATION_RENDERING.style_for(
        station_preset,
    )
    style = copy.copy(style)
    style.side = "top"
    style.max_labels = max(int(style.max_labels), len(labels))
    if station_label_step is None:
        style.every = 1
    else:
        style.every = int(station_label_step)
    style.apply(
        ax,
        positions,
        labels,
        xlim=(-0.5, len(labels) - 0.5),
    )


def _station_confidence_vectors(
    table: Any,
    station: str,
    fr: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return confidence and matching flags for station frequencies."""
    confidence = np.full(fr.size, np.nan, dtype=float)
    flags = np.full(fr.size, "", dtype=object)
    if table is None or getattr(table, "empty", True):
        return confidence, flags
    sub = table[table["station"].astype(str) == str(station)]
    if sub.empty:
        return confidence, flags
    ftab = sub["frequency_hz"].to_numpy(dtype=float)
    ctab = sub["confidence"].to_numpy(dtype=float)
    flagtab = sub["flags"].astype(str).to_numpy()
    for i, freq in enumerate(fr):
        if not np.isfinite(freq):
            continue
        idx = int(np.nanargmin(np.abs(ftab - freq)))
        if np.isclose(ftab[idx], freq, rtol=1e-6, atol=1e-12):
            confidence[i] = ctab[idx]
            flags[i] = flagtab[idx]
    return confidence, flags


def _apply_row_mask_to_block(
    obj: Any, fields: Sequence[str], keep: np.ndarray, fr: np.ndarray
) -> None:
    """Apply a row-keep mask to tensor/tipper arrays and their frequency."""
    if _set_masked_strict_block(obj, fields, keep, fr):
        return
    for field in fields:
        value = getattr(obj, field, None)
        if isinstance(value, np.ndarray) and value.shape[0] == fr.size:
            _set_array_field(obj, field, value[keep])
    _set_block_freq(obj, fr[keep])


def _set_bad_rows_to_nan(
    obj: Any, fields: Sequence[str], bad: np.ndarray, fr: np.ndarray
) -> None:
    """Set selected tensor/tipper rows to NaN."""
    for field in fields:
        value = getattr(obj, field, None)
        if isinstance(value, np.ndarray) and value.shape[0] == fr.size:
            new = value.copy()
            new[bad] = np.nan
            _set_array_field(obj, field, new)


def _set_array_field(obj: Any, field: str, value: np.ndarray) -> None:
    """Set an array field, tolerating NaN-containing values.

    ``Z``-like containers with ``compute_resistivity_phase`` accept
    NaN-containing ``z``/``z_err`` writes fine (the setter applies the
    array first, then logs internally if the derived resistivity/phase
    recompute can't handle the NaNs) — masking and recovery both rely on
    writing NaN rows here, so this must not refuse the write.
    """
    try:
        setattr(obj, field, value)
    except Exception:
        pass


def _set_masked_strict_block(
    obj: Any, fields: Sequence[str], keep: np.ndarray, fr: np.ndarray
) -> bool:
    """Atomically subset strict Z-like containers."""
    if not hasattr(obj, "compute_resistivity_phase"):
        return False
    if "z" not in fields or not hasattr(obj, "_z"):
        return False
    z = getattr(obj, "z", None)
    if not isinstance(z, np.ndarray) or z.shape[0] != fr.size:
        return False
    z_new = z[keep]
    if not np.isfinite(z_new).all():
        return False
    try:
        obj._freq = np.asarray(fr[keep], dtype=float)
        obj._z = z_new.astype(complex, copy=False)
        z_err = getattr(obj, "z_err", None)
        if (
            "z_err" in fields
            and isinstance(z_err, np.ndarray)
            and z_err.shape[0] == fr.size
        ):
            obj._z_err = z_err[keep].astype(float, copy=False)
        obj.compute_resistivity_phase()
    except Exception:
        return False
    return True


def _set_block_freq(obj: Any, fr: np.ndarray) -> None:
    try:
        obj.freq = fr
    except Exception:
        pass


def _regrid_z(
    z: np.ndarray,
    fr: np.ndarray,
    fnew: np.ndarray,
    *,
    method: str,
) -> np.ndarray:
    # z: (n,2,2) complex
    n = fnew.size
    out = np.empty((n, 2, 2), dtype=complex)
    x = np.log10(fr)
    xn = np.log10(fnew)
    for a in range(2):
        for b in range(2):
            y = z[:, a, b]
            out[:, a, b] = _interp_complex(x, y, xn, method=method)
    return out


def _regrid_t(
    t: np.ndarray,
    fr: np.ndarray,
    fnew: np.ndarray,
    *,
    method: str,
) -> np.ndarray:
    # t: (n,2) complex
    n = fnew.size
    out = np.empty((n, 2), dtype=complex)
    x = np.log10(fr)
    xn = np.log10(fnew)
    for k in range(2):
        y = t[:, k]
        out[:, k] = _interp_complex(x, y, xn, method=method)
    return out


def _apply_grid_one(
    ed: Any,
    fnew: np.ndarray,
    *,
    method: str,
) -> None:
    Z, z, frz = _get_z_block(ed)
    if Z is not None:
        z2 = _regrid_z(z, frz, fnew, method=method)
        try:
            Z.z = z2
        except Exception:
            pass
        z_err = getattr(Z, "z_err", None)
        if isinstance(z_err, np.ndarray) and z_err.shape[0] == frz.size:
            try:
                Z.z_err = _regrid_z(z_err, frz, fnew, method=method).real
            except Exception:
                try:
                    Z._z_err = _regrid_z(z_err, frz, fnew, method=method).real
                except Exception:
                    pass
        _set_block_freq(Z, fnew)
    T, t, frt = _get_t_block(ed)
    if T is not None:
        t2 = _regrid_t(t, frt, fnew, method=method)
        try:
            T.tipper = t2
        except Exception:
            pass
        tipper_err = getattr(T, "tipper_err", None)
        if (
            isinstance(tipper_err, np.ndarray)
            and tipper_err.shape[0] == frt.size
        ):
            try:
                T.tipper_err = _regrid_t(
                    tipper_err, frt, fnew, method=method
                ).real
            except Exception:
                try:
                    T._tipper_err = _regrid_t(
                        tipper_err, frt, fnew, method=method
                    ).real
                except Exception:
                    pass
        _set_block_freq(T, fnew)


def _valid_freq_of(ed: Any) -> np.ndarray | None:
    Z, z, fr = _get_z_block(ed)
    if isinstance(fr, np.ndarray) and fr.size:
        return _sorted_unique(fr)
    T, t, fr = _get_t_block(ed)
    if isinstance(fr, np.ndarray) and fr.size:
        return _sorted_unique(fr)
    return None


def _union_grid(SitesObj) -> np.ndarray:
    frs: list[np.ndarray] = []
    for ed in _iter_items(SitesObj):
        fr = _valid_freq_of(ed)
        if isinstance(fr, np.ndarray) and fr.size:
            frs.append(fr)
    if not frs:
        return np.array([], dtype=float)
    g = np.unique(np.concatenate(frs))
    return _sorted_unique(g)


def _intersect_grid(SitesObj) -> np.ndarray:
    frs: list[np.ndarray] = []
    for ed in _iter_items(SitesObj):
        fr = _valid_freq_of(ed)
        if isinstance(fr, np.ndarray) and fr.size:
            frs.append(fr)
    if not frs:
        return np.array([], dtype=float)
    g = frs[0]
    for fr in frs[1:]:
        g = np.intersect1d(g, fr)
    return _sorted_unique(g)


# ------------------------------ public API ------------------------------- #


def select_band(
    sites: Any,
    *,
    fmin: float | None = None,
    fmax: float | None = None,
    band_hz: Sequence[float] | None = None,
    keep: Sequence[float] | None = None,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    # Resolve band_hz alias → fmin / fmax.
    # band_hz=(lo, hi) is a convenience form; canonical fmin/fmax always win.
    if band_hz is not None:
        _lo, _hi = band_hz
        if fmin is None:
            fmin = _lo
        if fmax is None:
            fmax = _hi
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    keep_values = None if keep is None else np.asarray(keep, dtype=float)

    def _keep_mask(fr: np.ndarray) -> np.ndarray:
        fr = np.asarray(fr, dtype=float)
        if keep_values is not None:
            mask = np.zeros(fr.size, dtype=bool)
            for value in keep_values[np.isfinite(keep_values)]:
                mask |= np.isclose(fr, value, rtol=1e-6, atol=1e-12)
            return mask
        mask = np.ones(fr.size, dtype=bool)
        if fmin is not None:
            mask &= fr >= float(fmin)
        if fmax is not None:
            mask &= fr <= float(fmax)
        return mask

    def _one(Si):
        ed = next(_iter_items(Si))
        Z, z, fr = _get_z_block(ed)
        if Z is not None:
            _apply_row_mask_to_block(Z, ("z", "z_err"), _keep_mask(fr), fr)
        T, t, fr = _get_t_block(ed)
        if T is not None:
            _apply_row_mask_to_block(
                T, ("tipper", "tipper_err"), _keep_mask(fr), fr
            )
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


def drop_duplicates(
    sites: Any,
    *,
    tol: float = 1e-10,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    def _one(Si):
        ed = next(_iter_items(Si))
        for objname in ("Z", "z", "Tipper", "tipper"):
            O = getattr(ed, objname, None)
            if O is None:
                continue
            fr = getattr(O, "freq", None)
            if not isinstance(fr, np.ndarray):
                continue
            fr = np.asarray(fr, dtype=float)
            if fr.size == 0:
                continue
            frs = _sorted_unique(fr)
            if frs.size == fr.size and np.allclose(fr, frs):
                continue
            # nearest index map
            idx = _nearest_idx(frs, fr)
            # keep first occurrence only
            seen = set()
            keep_idx = []
            for i, j in enumerate(idx):
                if abs(fr[i] - frs[j]) <= tol and j not in seen:
                    keep_idx.append(i)
                    seen.add(j)
            keep_idx = np.array(keep_idx, dtype=int)
            for fld in ("z", "z_err", "tipper", "tipper_err"):
                val = getattr(O, fld, None)
                if isinstance(val, np.ndarray) and val.shape[0] == fr.size:
                    setattr(O, fld, val[keep_idx])
            _set_block_freq(O, fr[keep_idx])
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


def drop_low_confidence_frequencies(
    sites: Any,
    *,
    method: str = "composite",
    threshold: float = 0.50,
    weights: dict[str, float] | None = None,
    also: str = "both",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    """Drop rows whose frequency confidence is below ``threshold``.

    The confidence scores are computed with
    :func:`pycsamt.emtools.qc.frequency_confidence_table`.  The operation is
    station-aware: each station keeps or drops its own bad frequency rows.
    A new :class:`~pycsamt.site.base.Sites` object is returned unless
    ``inplace=True``.
    """
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    table = _confidence_decision_table(
        S,
        method=method,
        ci_hi=max(float(threshold), 0.50),
        ci_lo=float(threshold),
        weights=weights,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )

    def _one(Si):
        ed = next(_iter_items(Si))
        station = _name(ed, 0)
        if also in ("z", "both"):
            Z, z, fr = _get_z_block(ed)
            if Z is not None:
                conf, _ = _station_confidence_vectors(table, station, fr)
                keep = ~np.isfinite(conf) | (conf >= float(threshold))
                _apply_row_mask_to_block(Z, ("z", "z_err"), keep, fr)
        if also in ("tipper", "both"):
            T, t, fr = _get_t_block(ed)
            if T is not None:
                conf, _ = _station_confidence_vectors(table, station, fr)
                keep = ~np.isfinite(conf) | (conf >= float(threshold))
                _apply_row_mask_to_block(
                    T,
                    ("tipper", "tipper_err"),
                    keep,
                    fr,
                )
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


def mask_low_confidence_frequencies(
    sites: Any,
    *,
    method: str = "composite",
    threshold: float = 0.50,
    weights: dict[str, float] | None = None,
    also: str = "both",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    """Set low-confidence frequency rows to NaN without changing the grid."""
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    table = _confidence_decision_table(
        S,
        method=method,
        ci_hi=max(float(threshold), 0.50),
        ci_lo=float(threshold),
        weights=weights,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )

    def _one(Si):
        ed = next(_iter_items(Si))
        station = _name(ed, 0)
        if also in ("z", "both"):
            Z, z, fr = _get_z_block(ed)
            if Z is not None:
                conf, _ = _station_confidence_vectors(table, station, fr)
                bad = np.isfinite(conf) & (conf < float(threshold))
                _set_bad_rows_to_nan(Z, ("z", "z_err"), bad, fr)
        if also in ("tipper", "both"):
            T, t, fr = _get_t_block(ed)
            if T is not None:
                conf, _ = _station_confidence_vectors(table, station, fr)
                bad = np.isfinite(conf) & (conf < float(threshold))
                _set_bad_rows_to_nan(T, ("tipper", "tipper_err"), bad, fr)
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


def recover_low_confidence_frequencies(
    sites: Any,
    *,
    method: str = "composite",
    ci_hi: float = 0.90,
    ci_lo: float = 0.50,
    weights: dict[str, float] | None = None,
    interpolation: str = "linear",
    reject: str = "mask",
    also: str = "both",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    """Recover recoverable frequency rows using trusted neighboring rows.

    Rows with confidence in ``[ci_lo, ci_hi)`` are treated as recoverable and
    are interpolated in log-frequency from rows with confidence
    ``>= ci_hi``.  Rows below ``ci_lo`` are considered rejected and are either
    masked, dropped, or kept depending on ``reject``.
    """
    reject = str(reject).lower()
    if reject not in {"mask", "drop", "keep"}:
        msg = "reject must be one of 'mask', 'drop', or 'keep'."
        raise ValueError(msg)
    if interpolation not in {"linear", "nearest"}:
        msg = "interpolation must be 'linear' or 'nearest'."
        raise ValueError(msg)
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    table = _confidence_decision_table(
        S,
        method=method,
        ci_hi=ci_hi,
        ci_lo=ci_lo,
        weights=weights,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )

    def _recover_block(obj, fields, fr, conf):
        trusted = np.isfinite(conf) & (conf >= float(ci_hi))
        recover = (
            np.isfinite(conf) & (conf >= float(ci_lo)) & (conf < float(ci_hi))
        )
        reject_mask = np.isfinite(conf) & (conf < float(ci_lo))
        if recover.any():
            for field in fields:
                value = getattr(obj, field, None)
                if (
                    isinstance(value, np.ndarray)
                    and value.shape[0] == fr.size
                ):
                    recovered = _interp_rows_by_freq(
                        value,
                        fr,
                        recover,
                        trusted,
                        method=interpolation,
                    )
                    _set_array_field(obj, field, recovered)
        if reject == "mask" and reject_mask.any():
            _set_bad_rows_to_nan(obj, fields, reject_mask, fr)
        elif reject == "drop" and reject_mask.any():
            keep = ~reject_mask
            _apply_row_mask_to_block(obj, fields, keep, fr)

    def _one(Si):
        ed = next(_iter_items(Si))
        station = _name(ed, 0)
        if also in ("z", "both"):
            Z, z, fr = _get_z_block(ed)
            if Z is not None:
                conf, _ = _station_confidence_vectors(table, station, fr)
                _recover_block(Z, ("z", "z_err"), fr, conf)
        if also in ("tipper", "both"):
            T, t, fr = _get_t_block(ed)
            if T is not None:
                conf, _ = _station_confidence_vectors(table, station, fr)
                _recover_block(T, ("tipper", "tipper_err"), fr, conf)
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


def edit_frequencies_by_confidence(
    sites: Any,
    *,
    mode: str = "recover",
    before_sites: Any | None = None,
    method: str = "composite",
    threshold: float = 0.50,
    ci_hi: float = 0.90,
    ci_lo: float = 0.50,
    weights: dict[str, float] | None = None,
    interpolation: str = "linear",
    reject: str = "drop",
    also: str = "both",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    api: bool | None = None,
) -> Any:
    """Edit frequency rows and return diagnostics in one workflow.

    This is the high-level confidence-editing entry point.  It applies one
    of the low-level edit strategies and immediately computes a station
    report and a station-frequency decision table.  Use ``before_sites`` when
    ``sites`` is already an in-memory object and a reliable before/after
    comparison is required, because lower-level site editors can mutate the
    wrapped impedance objects while constructing the edited return value.

    Parameters
    ----------
    sites : path-like, EDI-like, Sites, or sequence
        Input data to edit.  Path-like inputs are normally safe to use
        directly because they can be loaded independently by the package.
        In-memory objects should be paired with ``before_sites`` when the
        report must preserve an untouched baseline.
    mode : {'recover', 'drop', 'mask'}, default 'recover'
        Frequency-editing strategy.  ``'recover'`` interpolates recoverable
        rows in log-frequency and handles rejected rows according to
        ``reject``.  ``'drop'`` removes rows below ``threshold``.
        ``'mask'`` keeps the frequency grid but replaces low-confidence
        tensor rows by missing values when the container allows it.
    before_sites : optional
        Independent baseline used only for reporting and decision tracking.
        If omitted, ``sites`` is used as the baseline.
    method : str, default 'composite'
        Confidence metric passed to
        :func:`pycsamt.emtools.qc.frequency_confidence_table`.
    threshold : float, default 0.50
        Confidence threshold used by ``mode='drop'`` and ``mode='mask'``.
    ci_hi, ci_lo : float, default 0.90 and 0.50
        High-confidence and low-confidence limits used by ``mode='recover'``
        and by the diagnostic report.
    weights : dict or None, default None
        Optional confidence-metric weights.
    interpolation : {'linear', 'nearest'}, default 'linear'
        Interpolation strategy for recoverable rows in ``mode='recover'``.
    reject : {'drop', 'mask', 'keep'}, default 'drop'
        Handling of rows below ``ci_lo`` in ``mode='recover'``.
    also : {'z', 'tipper', 'both'}, default 'both'
        Data blocks edited when present.
    inplace : bool, default False
        Forwarded to the low-level edit function.
    recursive, on_dup, strict, verbose
        Site-loading options forwarded to :func:`ensure_sites`.

    Returns
    -------
    FrequencyEditResult
        Edited sites together with station-level and station-frequency
        diagnostics.
    """
    mode = str(mode).lower()
    if mode not in {"recover", "drop", "mask"}:
        msg = "mode must be one of 'recover', 'drop', or 'mask'."
        raise ValueError(msg)
    baseline = before_sites if before_sites is not None else sites
    if mode == "recover":
        edited = recover_low_confidence_frequencies(
            sites,
            method=method,
            ci_hi=ci_hi,
            ci_lo=ci_lo,
            weights=weights,
            interpolation=interpolation,
            reject=reject,
            also=also,
            inplace=inplace,
            recursive=recursive,
            on_dup=on_dup,
            strict=strict,
            verbose=verbose,
        )
    elif mode == "drop":
        edited = drop_low_confidence_frequencies(
            sites,
            method=method,
            threshold=threshold,
            weights=weights,
            also=also,
            inplace=inplace,
            recursive=recursive,
            on_dup=on_dup,
            strict=strict,
            verbose=verbose,
        )
    else:
        edited = mask_low_confidence_frequencies(
            sites,
            method=method,
            threshold=threshold,
            weights=weights,
            also=also,
            inplace=inplace,
            recursive=recursive,
            on_dup=on_dup,
            strict=strict,
            verbose=verbose,
        )
    report = frequency_edit_report(
        baseline,
        edited,
        method=method,
        ci_hi=ci_hi,
        ci_lo=ci_lo,
        weights=weights,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    decisions = frequency_edit_decision_table(
        baseline,
        edited,
        method=method,
        ci_hi=ci_hi,
        ci_lo=ci_lo,
        weights=weights,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    result = FrequencyEditResult(
        sites=edited,
        report=report,
        decisions=decisions,
        mode=mode,
        method=method,
        ci_hi=float(ci_hi),
        ci_lo=float(ci_lo),
        reject=str(reject),
        interpolation=str(interpolation),
    )
    if api is False:
        return result
    if api is None and not PYCSAMT_API_VIEW.enabled():
        return result
    return wrap_result(
        {
            "sites": result.sites,
            "report": result.report,
            "decisions": result.decisions,
            "mode": result.mode,
            "method": result.method,
            "ci_hi": result.ci_hi,
            "ci_lo": result.ci_lo,
            "reject": result.reject,
            "interpolation": result.interpolation,
            "n_dropped": result.n_dropped,
            "n_masked": result.n_masked,
            "n_recovered": result.n_recovered,
        },
        name="frequency_edit",
        kind="emtools.frequency.edit",
        meta={
            "mode": result.mode,
            "method": result.method,
            "ci_hi": result.ci_hi,
            "ci_lo": result.ci_lo,
            "reject": result.reject,
            "interpolation": result.interpolation,
        },
    )


def frequency_edit_report(
    before_sites: Any,
    after_sites: Any,
    *,
    method: str = "composite",
    ci_hi: float = 0.90,
    ci_lo: float = 0.50,
    weights: dict[str, float] | None = None,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    api: bool | None = None,
):
    """Summarize station-level changes after frequency editing.

    The report compares the native frequency rows and finite tensor rows
    before and after an edit such as dropping, masking, or recovery.  It also
    carries the median confidence from
    :func:`pycsamt.emtools.qc.frequency_confidence_table`.
    """
    import pandas as pd

    before = ensure_sites(
        before_sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    after = ensure_sites(
        after_sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    tb_before = _confidence_decision_table(
        before,
        method=method,
        ci_hi=ci_hi,
        ci_lo=ci_lo,
        weights=weights,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )
    tb_after = _confidence_decision_table(
        after,
        method=method,
        ci_hi=ci_hi,
        ci_lo=ci_lo,
        weights=weights,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )

    def _counts(sites_obj):
        rows = []
        for i, ed in enumerate(_iter_items(sites_obj)):
            station = _name(ed, i)
            Z, z, fr = _get_z_block(ed)
            if Z is None:
                rows.append(
                    dict(
                        station=station,
                        n_freq=0,
                        n_finite=0,
                        frac_finite=np.nan,
                    )
                )
                continue
            finite_rows = np.isfinite(z.reshape(z.shape[0], -1)).all(axis=1)
            rows.append(
                dict(
                    station=station,
                    n_freq=int(fr.size),
                    n_finite=int(np.count_nonzero(finite_rows)),
                    frac_finite=float(
                        np.count_nonzero(finite_rows) / max(1, fr.size),
                    ),
                )
            )
        return pd.DataFrame.from_records(rows)

    def _conf_summary(table, suffix):
        if table.empty:
            return pd.DataFrame(
                columns=[
                    "station",
                    f"confidence_median_{suffix}",
                    f"safe_fraction_{suffix}",
                    f"recoverable_fraction_{suffix}",
                    f"reject_fraction_{suffix}",
                ]
            )
        grouped = table.assign(
            safe=table["confidence"] >= ci_hi,
            recoverable=(table["confidence"] >= ci_lo)
            & (table["confidence"] < ci_hi),
            reject=table["confidence"] < ci_lo,
        ).groupby("station")
        return grouped.agg(
            **{
                f"confidence_median_{suffix}": ("confidence", "median"),
                f"safe_fraction_{suffix}": ("safe", "mean"),
                f"recoverable_fraction_{suffix}": ("recoverable", "mean"),
                f"reject_fraction_{suffix}": ("reject", "mean"),
            }
        ).reset_index()

    left = _counts(before).rename(
        columns={
            "n_freq": "n_freq_before",
            "n_finite": "n_finite_before",
            "frac_finite": "frac_finite_before",
        }
    )
    right = _counts(after).rename(
        columns={
            "n_freq": "n_freq_after",
            "n_finite": "n_finite_after",
            "frac_finite": "frac_finite_after",
        }
    )
    out = left.merge(right, on="station", how="outer")
    out = out.merge(
        _conf_summary(tb_before, "before"), on="station", how="left"
    )
    out = out.merge(
        _conf_summary(tb_after, "after"), on="station", how="left"
    )
    out["n_dropped"] = out["n_freq_before"] - out["n_freq_after"]
    out["n_masked_or_unfinite"] = (
        out["n_finite_before"] - out["n_finite_after"]
    )
    out["confidence_delta"] = (
        out["confidence_median_after"] - out["confidence_median_before"]
    )
    df = out.sort_values("station").reset_index(drop=True)

    return maybe_wrap_frame(
        df,
        api=api,
        name="frequency_edit_report",
        kind="emtools.frequency.edit_report",
        source=before_sites,
        description="Station-level changes after frequency editing.",
    )


def frequency_edit_decision_table(
    before_sites: Any,
    after_sites: Any,
    *,
    method: str = "composite",
    ci_hi: float = 0.90,
    ci_lo: float = 0.50,
    weights: dict[str, float] | None = None,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    api: bool | None = None,
):
    """Return one row per original station-frequency edit decision."""
    import pandas as pd

    before = ensure_sites(
        before_sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    after = ensure_sites(
        after_sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    conf = _confidence_decision_table(
        before,
        method=method,
        ci_hi=ci_hi,
        ci_lo=ci_lo,
        weights=weights,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )
    after_map = {}
    for i, ed in enumerate(_iter_items(after)):
        station = _name(ed, i)
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue
        finite = np.isfinite(z.reshape(z.shape[0], -1)).all(axis=1)
        after_map[station] = (fr, finite, z)

    rows = []
    for i, ed in enumerate(_iter_items(before)):
        station = _name(ed, i)
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue
        finite_before = np.isfinite(z.reshape(z.shape[0], -1)).all(axis=1)
        after_fr, after_finite, after_z = after_map.get(
            station,
            (
                np.asarray([], dtype=float),
                np.asarray([], dtype=bool),
                np.asarray([], dtype=complex),
            ),
        )
        cvec, flagvec = _station_confidence_vectors(conf, station, fr)
        for j, freq in enumerate(fr):
            present_after = False
            finite_after = False
            changed_after = False
            if after_fr.size:
                idx = int(np.nanargmin(np.abs(after_fr - freq)))
                if np.isclose(after_fr[idx], freq, rtol=1e-6, atol=1e-12):
                    present_after = True
                    finite_after = bool(after_finite[idx])
                    if after_z.ndim == z.ndim and after_z.shape[0] > idx:
                        changed_after = not np.allclose(
                            z[j],
                            after_z[idx],
                            rtol=1e-7,
                            atol=1e-12,
                            equal_nan=True,
                        )
            if not present_after:
                action = "dropped"
            elif finite_before[j] and not finite_after:
                action = "masked"
            elif (not finite_before[j]) and finite_after:
                action = "recovered"
            elif (
                np.isfinite(cvec[j])
                and float(ci_lo) <= cvec[j] < float(ci_hi)
                and changed_after
            ):
                action = "recovered"
            else:
                action = "kept"
            rows.append(
                dict(
                    station=station,
                    frequency_hz=float(freq),
                    period_s=float(1.0 / freq) if freq else np.nan,
                    log10_period=(
                        float(np.log10(1.0 / freq)) if freq > 0 else np.nan
                    ),
                    confidence=float(cvec[j]),
                    flags=str(flagvec[j]),
                    finite_before=bool(finite_before[j]),
                    present_after=bool(present_after),
                    finite_after=bool(finite_after),
                    action=action,
                )
            )
    df = pd.DataFrame.from_records(rows)

    return maybe_wrap_frame(
        df,
        api=api,
        name="frequency_edit_decision_table",
        kind="emtools.frequency.edit_decisions",
        source=before_sites,
        description="Per-frequency edit decisions for each station.",
    )


def plot_frequency_edit_summary(
    before_sites: Any,
    after_sites: Any,
    *,
    method: str = "composite",
    ci_hi: float = 0.90,
    ci_lo: float = 0.50,
    figsize: tuple[float, float] = (9.0, 4.0),
    station_label_step: int | None = 1,
    station_preset: str = "pseudosection",
    station_style: Any | None = None,
    ax: Any | None = None,
):
    """Plot station-level before/after frequency-edit summary."""
    import matplotlib.pyplot as plt

    from pycsamt.api.style import PYCSAMT_STYLE

    report = frequency_edit_report(
        before_sites,
        after_sites,
        method=method,
        ci_hi=ci_hi,
        ci_lo=ci_lo,
    )
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if report.empty:
        ax.text(0.5, 0.5, "no stations", ha="center", va="center")
        return ax
    x = np.arange(len(report))
    before_kw = PYCSAMT_STYLE.correction.before.plot_kwargs()
    after_kw = PYCSAMT_STYLE.correction.after.plot_kwargs()
    ax.plot(
        x,
        report["n_freq_before"].to_numpy(dtype=float),
        **before_kw,
    )
    ax.plot(
        x,
        report["n_freq_after"].to_numpy(dtype=float),
        **after_kw,
    )
    dropped = report["n_dropped"].fillna(0).to_numpy(dtype=float)
    masked = report["n_masked_or_unfinite"].fillna(0).to_numpy(dtype=float)
    ax.bar(x, dropped, color="#d62728", alpha=0.22, label="dropped")
    ax.bar(
        x,
        masked,
        bottom=dropped,
        color="#ff99c8",
        alpha=0.28,
        label="masked/unfinite delta",
    )
    station_names = report["station"].astype(str).tolist()
    _apply_station_rendering(
        ax,
        x,
        station_names,
        station_label_step=station_label_step,
        station_preset=station_preset,
        station_style=station_style,
    )
    ax.set_ylabel("Frequency rows")
    ax.set_title("Frequency edit summary", fontsize=10)
    ax.grid(True, ls=":", alpha=0.35)
    ax.legend(fontsize=8)
    return ax


def plot_frequency_edit_decisions(
    before_sites: Any,
    after_sites: Any,
    *,
    method: str = "composite",
    ci_hi: float = 0.90,
    ci_lo: float = 0.50,
    figsize: tuple[float, float] = (10.0, 5.0),
    station_label_step: int | None = 1,
    station_preset: str = "pseudosection",
    station_style: Any | None = None,
    ax: Any | None = None,
):
    """Plot dropped, masked, recovered, and kept frequency decisions."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap

    table = frequency_edit_decision_table(
        before_sites,
        after_sites,
        method=method,
        ci_hi=ci_hi,
        ci_lo=ci_lo,
    )
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if table.empty:
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return ax
    action_codes = {
        "dropped": 0,
        "masked": 1,
        "recovered": 2,
        "kept": 3,
    }
    stations = table["station"].drop_duplicates().astype(str).tolist()
    periods = np.sort(table["log10_period"].dropna().unique())
    matrix = np.full((periods.size, len(stations)), np.nan, dtype=float)
    for j, station in enumerate(stations):
        sub = table[table["station"].astype(str) == station]
        lookup = {
            float(row.log10_period): action_codes.get(str(row.action), np.nan)
            for _, row in sub.iterrows()
            if np.isfinite(row.log10_period)
        }
        for i, period in enumerate(periods):
            matrix[i, j] = lookup.get(float(period), np.nan)
    cmap = ListedColormap(["#8b0026", "#ff99c8", "#20b455", "#d8d8d8"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)
    im = ax.imshow(
        matrix,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
        cmap=cmap,
        norm=norm,
        extent=(-0.5, len(stations) - 0.5, periods.min(), periods.max()),
    )
    _apply_station_rendering(
        ax,
        np.arange(len(stations)),
        stations,
        station_label_step=station_label_step,
        station_preset=station_preset,
        station_style=station_style,
    )
    ax.set_ylabel(r"$\log_{10}T$ (s)")
    ax.set_title("Frequency edit decisions", fontsize=10)
    cbar = ax.figure.colorbar(im, ax=ax, pad=0.015, ticks=[0, 1, 2, 3])
    cbar.ax.set_yticklabels(["dropped", "masked", "recovered", "kept"])
    return ax


def regrid_to(
    sites: Any,
    target_freq: Sequence[float],
    *,
    method: str = "nearest",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    fnew = _sorted_unique(np.asarray(target_freq, dtype=float))

    def _one(Si):
        ed = next(_iter_items(Si))
        if fnew.size == 0:
            return Si
        _apply_grid_one(ed, fnew, method=method)
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


def regrid_logspace(
    sites: Any,
    *,
    fmin: float | None = None,
    fmax: float | None = None,
    band_hz: Sequence[float] | None = None,
    per_decade: int = 6,
    n_per_decade: int | None = None,
    method: str = "nearest",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    # Resolve band_hz alias → fmin / fmax (canonical wins).
    if band_hz is not None:
        _lo, _hi = band_hz
        if fmin is None:
            fmin = _lo
        if fmax is None:
            fmax = _hi
    # Resolve n_per_decade alias → per_decade.
    # n_per_decade only takes effect when per_decade is still at its default (6).
    if n_per_decade is not None and per_decade == 6:
        per_decade = int(n_per_decade)
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    # infer band from union
    U = _union_grid(S)
    if U.size == 0:
        return S
    lo = fmin if fmin is not None else U.min()
    hi = fmax if fmax is not None else U.max()
    lo, hi = float(lo), float(hi)
    if lo <= 0 or hi <= 0 or hi <= lo:
        return S
    ndec = np.log10(hi) - np.log10(lo)
    npts = max(2, int(np.ceil(ndec * per_decade)))
    fnew = np.logspace(np.log10(lo), np.log10(hi), npts)
    return regrid_to(
        S,
        fnew,
        method=method,
        inplace=inplace,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )


def decimate_step(
    sites: Any,
    *,
    step: int = 2,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    def _one(Si):
        ed = next(_iter_items(Si))
        Z, z, frz = _get_z_block(ed)
        if Z is not None:
            idx = np.arange(0, frz.size, step, dtype=int)
            Z.z = z[idx]
            _set_block_freq(Z, frz[idx])
        T, t, frt = _get_t_block(ed)
        if T is not None:
            idx = np.arange(0, frt.size, step, dtype=int)
            T.tipper = t[idx]
            _set_block_freq(T, frt[idx])
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


def smooth_mavg(
    sites: Any,
    *,
    k: int = 3,
    window: int | None = None,
    on: str = "both",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    # Resolve window alias → k.
    # window only takes effect when k is still at its default (3).
    if window is not None and k == 3:
        k = int(window)
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    k = max(1, int(k))
    if k == 1:
        return S

    def _roll(y: np.ndarray) -> np.ndarray:
        if y.ndim == 1:
            y = y[:, None]
        out = np.empty_like(y)
        for j in range(y.shape[1]):
            v = y[:, j].astype(complex)
            re = np.convolve(v.real, np.ones(k) / k, "same")
            im = np.convolve(v.imag, np.ones(k) / k, "same")
            out[:, j] = re + 1j * im
        return out.squeeze()

    def _one(Si):
        ed = next(_iter_items(Si))
        if on in ("z", "both"):
            Z, z, fr = _get_z_block(ed)
            if Z is not None:
                z2 = z.copy()
                for a in range(2):
                    for b in range(2):
                        z2[:, a, b] = _roll(z[:, a, b])
                Z.z = z2
        if on in ("tipper", "both"):
            T, t, fr = _get_t_block(ed)
            if T is not None:
                t2 = t.copy()
                for j in range(2):
                    t2[:, j] = _roll(t[:, j])
                T.tipper = t2
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


def align_grid(
    sites: Any,
    *,
    mode: str = "union",
    ref_station: str | None = None,
    method: str = "nearest",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if mode == "union":
        fnew = _union_grid(S)
    elif mode == "intersection":
        fnew = _intersect_grid(S)
    elif mode == "ref":
        fnew = None
        for ed in _iter_items(S):
            nm = getattr(ed, "station", None) or getattr(ed, "name", None)
            if nm == ref_station:
                fr = _valid_freq_of(ed)
                if isinstance(fr, np.ndarray):
                    fnew = fr
                break
        if fnew is None:
            return S
    else:
        raise ValueError("mode must be union|intersection|ref")
    if fnew.size == 0:
        return S
    return regrid_to(
        S,
        fnew,
        method=method,
        inplace=inplace,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )


# ------------------------- coverage + quality ---------------------------- #

import matplotlib.pyplot as plt


def _rel_err_block(z: np.ndarray, ze: np.ndarray | None) -> np.ndarray:
    if ze is None:
        return np.zeros(z.shape[0], dtype=float)
    # avg relative err across off-diagonals (fallback: all)
    zx = z[:, 0, 1]
    zy = z[:, 1, 0]
    ex = ze[:, 0, 1] if ze.ndim == 3 else None
    ey = ze[:, 1, 0] if ze.ndim == 3 else None
    if ex is not None and ey is not None:
        rx = np.abs(ex) / (np.abs(zx) + 1e-12)
        ry = np.abs(ey) / (np.abs(zy) + 1e-12)
        r = 0.5 * (rx + ry)
        return np.asarray(r, dtype=float)
    # fallback to full-tensor mean
    num = np.linalg.norm(ze.reshape(ze.shape[0], -1), axis=1)
    den = np.linalg.norm(z.reshape(z.shape[0], -1), axis=1) + 1e-12
    return num / den


def plot_coverage_quality_heatmap(
    sites: Any,
    *,
    axis: str = "period",
    figsize: tuple[float, float] = (7.5, 4.5),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    # assemble presence & quality per site on union grid
    labs: list[str] = []
    frs: list[np.ndarray] = []
    quals: list[np.ndarray] = []
    for i, ed in enumerate(_iter_items(S)):
        nm = getattr(ed, "station", None) or getattr(ed, "name", None)
        nm = nm or f"site_{i}"
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue
        ze = getattr(Z, "z_err", None)
        r = _rel_err_block(z, ze)
        q = 1.0 / (1.0 + r)
        frs.append(_sorted_unique(fr))
        # re-align q to unique set by averaging duplicates
        uf = frs[-1]
        if uf.size != fr.size:
            idx = _nearest_idx(uf, fr)
            acc = np.zeros(uf.size, dtype=float)
            cnt = np.zeros(uf.size, dtype=float)
            for k, j in enumerate(idx):
                acc[j] += q[k]
                cnt[j] += 1.0
            q = acc / (cnt + 1e-12)
        else:
            q = q
        quals.append(q)
        labs.append(nm)
    grid = _union_grid(S)
    if grid.size == 0 or not labs:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return ax
    # project site vectors onto union grid
    M = np.zeros((len(labs), grid.size), dtype=float)
    for i, (fr, q) in enumerate(zip(frs, quals)):
        idx = _nearest_idx(grid, fr)
        M[i, idx] = q
    if axis == "period":
        ylab = "period (s)"
        y_values = 1.0 / np.maximum(grid, 1e-24)
        row_order = np.argsort(y_values)
        y_values = y_values[row_order]
        image = M.T[row_order]
    else:
        ylab = "freq (Hz)"
        y_values = grid
        image = M.T
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(
        image,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
        vmin=0.0,
        vmax=1.0,
    )
    ax.set_ylabel(ylab)
    PYCSAMT_STATION_RENDERING.apply(
        ax,
        np.arange(len(labs)),
        labs,
        preset="pseudosection",
        xlim=(-0.5, len(labs) - 0.5),
    )
    yt = np.linspace(0, image.shape[0] - 1, num=min(8, image.shape[0]))
    yv = np.linspace(y_values.min(), y_values.max(), num=yt.size)
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{v:.2g}" for v in yv])
    if axis == "period" and not ax.yaxis_inverted():
        ax.invert_yaxis()
    cb = plt.colorbar(im, ax=ax)
    cb.set_label("quality (1/(1+rel_err))")
    return ax


# ------------------------- apparent depth image -------------------------- #


def _rho_det_from_z(z: np.ndarray, fr: np.ndarray) -> np.ndarray:
    # ρa_xy/yx = 0.2 * |Z|^2 / f (ohm·m); det via geom. mean
    zx = z[:, 0, 1]
    zy = z[:, 1, 0]
    rx = 0.2 * (np.abs(zx) ** 2) / (fr + 1e-24)
    ry = 0.2 * (np.abs(zy) ** 2) / (fr + 1e-24)
    rdet = np.sqrt(rx * ry)
    return rdet


def _depth_from_rho_f(rho: np.ndarray, fr: np.ndarray) -> np.ndarray:
    # δ ≈ 503 * sqrt(ρ / f)  (meters)
    with np.errstate(invalid="ignore"):
        d = 503.0 * np.sqrt(np.maximum(rho, 0.0) / (fr + 1e-24))
    return d


def plot_apparent_depth_psection(
    sites: Any,
    *,
    axis_y: str = "period",
    agg: str = "median",
    figsize: tuple[float, float] = (7.5, 4.5),
    log_color: bool = True,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    # gather rows: station, freq, period, depth
    rows: list[tuple[str, float, float, float]] = []
    for i, ed in enumerate(_iter_items(S)):
        nm = getattr(ed, "station", None) or getattr(ed, "name", None)
        nm = nm or f"site_{i}"
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue
        rho = _rho_det_from_z(z, fr)
        dep = _depth_from_rho_f(rho, fr)
        per = 1.0 / fr
        for f, p, d in zip(fr, per, dep):
            if not np.isfinite(d):
                continue
            rows.append((nm, f, p, d))
    if not rows:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return ax
    # pivot to station × axis_y
    dtype = [
        ("station", "U64"),
        ("freq", "f8"),
        ("period", "f8"),
        ("depth", "f8"),
    ]
    arr = np.array(rows, dtype=dtype)
    import pandas as pd

    df = pd.DataFrame(arr)
    ykey = "period" if axis_y == "period" else "freq"
    piv = df.pivot_table(
        index=ykey,
        columns="station",
        values="depth",
        aggfunc=agg,
    )
    piv = piv.sort_index()
    if axis_y == "period":
        piv = piv.sort_index(ascending=True)
    # piv has shape (n_periods, n_stations); do NOT transpose so that
    # imshow maps x → stations and y → periods.
    Zm = piv.to_numpy(dtype=float)
    if log_color:
        Zplot = np.log10(np.maximum(Zm, 1e-3))
        cblab = "log10 depth (m)"
    else:
        Zplot = Zm
        cblab = "depth (m)"
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(
        Zplot,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
    )
    ax.set_ylabel("Period (s)" if axis_y == "period" else "Frequency (Hz)")
    PYCSAMT_STATION_RENDERING.apply(
        ax,
        np.arange(Zplot.shape[1], dtype=float),
        list(piv.columns),
        preset="pseudosection",
        xlim=(-0.5, Zplot.shape[1] - 0.5),
    )
    yt = np.linspace(
        0, Zplot.shape[0] - 1, num=min(8, Zplot.shape[0])
    )  # shape[0] = n_periods
    yv = np.linspace(
        piv.index.min(), piv.index.max(), num=min(8, len(piv.index))
    )
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{v:.3g}" for v in yv])
    if axis_y == "period" and not ax.yaxis_inverted():
        ax.invert_yaxis()
    cb = plt.colorbar(im, ax=ax)
    cb.set_label(cblab)
    return ax


# ------------------------- band-collapsed microplots --------------------- #


def _det_phase_from_z(z: np.ndarray) -> np.ndarray:
    detz = z[:, 0, 0] * z[:, 1, 1] - z[:, 0, 1] * z[:, 1, 0]
    return np.degrees(np.angle(detz))


def _tip_amp(t: np.ndarray | None) -> np.ndarray | None:
    if t is None:
        return None
    return np.sqrt(np.abs(t[:, 0]) ** 2 + np.abs(t[:, 1]) ** 2)


def _bands_auto(all_period: np.ndarray, n: int) -> np.ndarray:
    lo = float(np.nanmin(all_period))
    hi = float(np.nanmax(all_period))
    lo = max(lo, 1e-6)
    edges = np.logspace(np.log10(lo), np.log10(hi), n + 1)
    return edges


def plot_band_microstrips(
    sites: Any,
    *,
    bands: Sequence[tuple[float, float]] | None = None,
    n_bands: int = 6,
    figsize: tuple[float, float] = (9.0, 6.0),
    marker_size: float = 16.0,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    # collect per-site arrays
    sites_data = []
    allP = []
    for i, ed in enumerate(_iter_items(S)):
        nm = getattr(ed, "station", None) or getattr(ed, "name", None)
        nm = nm or f"site_{i}"
        Z, z, fr = _get_z_block(ed)
        T, t, frt = _get_t_block(ed)
        if Z is None:
            continue
        per = 1.0 / fr
        rho = _rho_det_from_z(z, fr)
        lgr = np.log10(np.maximum(rho, 1e-12))
        ph = _det_phase_from_z(z)
        ta = _tip_amp(t) if T is not None else None
        sites_data.append((nm, per, lgr, ph, ta))
        allP.append(per)
    if not sites_data:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return ax
    allP = np.concatenate(allP)
    if bands is None:
        edges = _bands_auto(allP, n_bands)
        bands = list(zip(edges[:-1], edges[1:]))
    # aggregate per band
    rows = []
    for nm, per, lgr, ph, ta in sites_data:
        for bi, (lo, hi) in enumerate(bands):
            m = (per >= lo) & (per < hi)
            if not np.any(m):
                rows.append((nm, bi, np.nan, np.nan, np.nan))
                continue
            r = float(np.nanmedian(lgr[m]))
            p = float(np.nanmedian(ph[m]))
            if ta is not None and ta.size == per.size:
                tt = float(np.nanmedian(ta[m]))
            else:
                tt = np.nan
            rows.append((nm, bi, r, p, tt))
    # normalize per metric for color
    import pandas as pd

    df = pd.DataFrame(
        rows, columns=["station", "band", "logrho", "phi", "tamp"]
    )

    # global min/max (robust) for color normalization
    def _rng(s: pd.Series) -> tuple[float, float]:
        v = s.to_numpy()
        v = v[np.isfinite(v)]
        if v.size == 0:
            return 0.0, 1.0
        return (np.nanpercentile(v, 5), np.nanpercentile(v, 95))

    r0, r1 = _rng(df["logrho"])
    p0, p1 = _rng(df["phi"])
    t0, t1 = _rng(df["tamp"])

    def _norm(x, a, b):
        return np.clip((x - a) / (b - a + 1e-12), 0.0, 1.0)

    df["cr"] = _norm(df["logrho"], r0, r1)
    df["cp"] = _norm(df["phi"], p0, p1)
    df["ct"] = _norm(df["tamp"], t0, t1)
    # plot: 3 dots per band (circle, square, triangle)
    st = list(df["station"].unique())
    st_map = {s: i for i, s in enumerate(st)}
    y0 = np.arange(len(st), dtype=float)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    for _, row in df.iterrows():
        yi = st_map[row["station"]]
        xi = float(row["band"])
        # offsets for 3 metrics
        yR = yi - 0.22
        yP = yi + 0.00
        yT = yi + 0.22
        ax.scatter(
            [xi],
            [yR],
            s=marker_size,
            c=[[row["cr"]]],
            cmap="viridis",
            marker="o",
            vmin=0.0,
            vmax=1.0,
        )
        ax.scatter(
            [xi],
            [yP],
            s=marker_size,
            c=[[row["cp"]]],
            cmap="viridis",
            marker="s",
            vmin=0.0,
            vmax=1.0,
        )
        ax.scatter(
            [xi],
            [yT],
            s=marker_size,
            c=[[row["ct"]]],
            cmap="viridis",
            marker="^",
            vmin=0.0,
            vmax=1.0,
        )
    ax.set_xlim(-0.5, len(bands) - 0.5)
    ax.set_ylim(-0.6, len(st) - 0.4)
    ax.set_xlabel("band")
    ax.set_ylabel("station")
    ax.set_yticks(y0)
    ax.set_yticklabels(st)
    # simple mini-legend using text markers
    ax.text(1.01, 0.85, "o log10 ρ_det", transform=ax.transAxes)
    ax.text(1.01, 0.70, "s φ_det", transform=ax.transAxes)
    ax.text(1.01, 0.55, "^ |T|", transform=ax.transAxes)
    return ax
