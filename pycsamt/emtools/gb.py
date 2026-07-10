from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from ..api.view import maybe_wrap_frame
from ._core import (
    _apply_each,
    _get_z_block,
    _iter_items,
    _name,
    ensure_sites,
)


def _gb_items(sites: Any, resolved: Any) -> list[Any]:
    """Return resolved sites, or direct EDI-like items if coercion was empty."""
    items = list(_iter_items(resolved))
    if items:
        return items
    if isinstance(sites, (str, bytes)):
        return items
    direct = list(_iter_items(sites))
    return [ed for ed in direct if _get_z_block(ed)[0] is not None]


@dataclass
class GroomBaileyResult:
    """Container returned by :func:`groom_bailey_decomposition`."""

    sites: Any
    table: pd.DataFrame
    applied: bool
    method: str

    @property
    def n_station(self) -> int:
        """Number of stations with fitted distortion parameters."""
        return int(len(self.table))

    def summary(self) -> str:
        """Return a compact text summary."""
        med = (
            float(np.nanmedian(self.table["rms_fit"]))
            if not self.table.empty and "rms_fit" in self.table
            else np.nan
        )
        return (
            "GroomBaileyResult("
            f"stations={self.n_station}, applied={self.applied}, "
            f"median_rms={med:.4g})"
        )

    def __repr__(self) -> str:  # noqa: D105
        return self.summary()


def _rotmat(deg: float) -> np.ndarray:
    th = np.radians(float(deg))
    c, s = np.cos(th), np.sin(th)
    return np.array([[c, s], [-s, c]], dtype=float)


def _rotate_tensor(z: np.ndarray, deg: float) -> np.ndarray:
    R = _rotmat(deg)
    Rt = R.T
    return R[None, :, :] @ z @ Rt[None, :, :]


def _select_band(
    z: np.ndarray,
    fr: np.ndarray,
    band: tuple[float, float] | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    fr = np.asarray(fr, dtype=float).ravel()
    n = min(z.shape[0], fr.size)
    z = np.asarray(z[:n], dtype=np.complex128)
    fr = fr[:n]
    mask = np.isfinite(fr) & (fr > 0.0)
    mask &= np.isfinite(z).all(axis=(1, 2))
    if band is not None:
        per = 1.0 / fr
        lo, hi = float(band[0]), float(band[1])
        mask &= (per >= lo) & (per <= hi)
    return z[mask], fr[mask], mask


def _anti_diagonal_from_distortion(
    z: np.ndarray,
    D: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Best anti-diagonal regional tensor for fixed real distortion D."""
    c00, c01 = D[0, 0], D[0, 1]
    c10, c11 = D[1, 0], D[1, 1]
    den_u = max(c00 * c00 + c10 * c10, 1e-24)
    den_v = max(c01 * c01 + c11 * c11, 1e-24)
    u = (c00 * z[:, 0, 1] + c10 * z[:, 1, 1]) / den_u
    v = (c01 * z[:, 0, 0] + c11 * z[:, 1, 0]) / den_v
    z2 = np.zeros_like(z)
    z2[:, 0, 1] = u
    z2[:, 1, 0] = v
    return z2, u, v


def _solve_real_row(
    obs0: np.ndarray,
    obs1: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    rows = []
    rhs = []
    for a, b, uu, vv, w in zip(obs0, obs1, u, v, weights):
        sw = float(np.sqrt(max(w, 0.0)))
        vals = (
            (np.array([0.0, np.real(vv)]), np.real(a)),
            (np.array([0.0, np.imag(vv)]), np.imag(a)),
            (np.array([np.real(uu), 0.0]), np.real(b)),
            (np.array([np.imag(uu), 0.0]), np.imag(b)),
        )
        for row, val in vals:
            rows.append(sw * row)
            rhs.append(sw * val)
    A = np.asarray(rows, dtype=float)
    y = np.asarray(rhs, dtype=float)
    if A.size == 0 or np.linalg.matrix_rank(A) < 2:
        return np.array([1.0, 0.0], dtype=float)
    sol, *_ = np.linalg.lstsq(A, y, rcond=None)
    return sol.astype(float)


def _normalise_distortion(
    D: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    det = float(np.linalg.det(D))
    if not np.isfinite(det) or abs(det) < 1e-12:
        return D, u, v
    scale = np.sqrt(abs(det))
    if not np.isfinite(scale) or scale <= 0:
        return D, u, v
    return D / scale, u * scale, v * scale


def _fit_gb_distortion(
    z: np.ndarray,
    *,
    max_iter: int,
    tol: float,
    robust: bool,
) -> dict[str, Any]:
    D = np.eye(2, dtype=float)
    weights = np.ones(z.shape[0], dtype=float)
    last = np.inf
    z2 = np.zeros_like(z)
    for _ in range(max(1, int(max_iter))):
        z2, u, v = _anti_diagonal_from_distortion(z, D)
        row0 = _solve_real_row(z[:, 0, 0], z[:, 0, 1], u, v, weights)
        row1 = _solve_real_row(z[:, 1, 0], z[:, 1, 1], u, v, weights)
        D = np.vstack([row0, row1])
        D, u, v = _normalise_distortion(D, u, v)
        z2 = np.zeros_like(z)
        z2[:, 0, 1] = u
        z2[:, 1, 0] = v
        model = D[None, :, :] @ z2
        resid = z - model
        denom = np.sqrt(np.nanmean(np.abs(z) ** 2)) + 1e-24
        rms = float(np.sqrt(np.nanmean(np.abs(resid) ** 2)) / denom)
        if abs(last - rms) <= float(tol):
            break
        last = rms
        if robust:
            r = np.sqrt(np.nanmean(np.abs(resid) ** 2, axis=(1, 2)))
            med = np.nanmedian(r)
            mad = 1.4826 * np.nanmedian(np.abs(r - med)) + 1e-24
            c = 1.345 * mad
            weights = np.ones_like(r)
            high = r > c
            weights[high] = c / np.maximum(r[high], 1e-24)
            weights = np.clip(weights, 0.05, 1.0)
    try:
        corrected = np.linalg.inv(D)[None, :, :] @ z
    except np.linalg.LinAlgError:
        corrected = np.full_like(z, np.nan)
    return dict(D=D, regional=z2, corrected=corrected, rms_fit=rms)


def _gb_parameters(D: np.ndarray) -> dict[str, float]:
    gain = float(np.sqrt(abs(np.linalg.det(D))))
    Dn = D / gain if np.isfinite(gain) and gain > 0 else D.copy()
    twist_rad = np.arctan2(
        Dn[0, 1] - Dn[1, 0],
        Dn[0, 0] + Dn[1, 1],
    )
    twist_deg = float(np.degrees(twist_rad))
    M = _rotmat(-twist_deg) @ Dn
    denom = float(M[0, 0] + M[1, 1])
    if abs(denom) < 1e-12:
        shear = np.nan
        anisotropy = np.nan
    else:
        shear = float((M[0, 1] + M[1, 0]) / denom)
        anisotropy = float((M[0, 0] - M[1, 1]) / denom)
    shear = float(np.clip(shear, -0.99, 0.99)) if np.isfinite(shear) else shear
    anisotropy = (
        float(np.clip(anisotropy, -0.99, 0.99))
        if np.isfinite(anisotropy) else anisotropy
    )
    return dict(
        gain=gain,
        twist_deg=twist_deg,
        shear=shear,
        shear_angle_deg=float(np.degrees(np.arctan(shear)))
        if np.isfinite(shear) else np.nan,
        anisotropy=anisotropy,
    )


def _diag_ratio(z: np.ndarray) -> float:
    diag = np.sqrt(np.abs(z[:, 0, 0]) ** 2 + np.abs(z[:, 1, 1]) ** 2)
    off = np.sqrt(np.abs(z[:, 0, 1]) ** 2 + np.abs(z[:, 1, 0]) ** 2)
    val = diag / np.maximum(off, 1e-24)
    return float(np.nanmedian(val)) if np.isfinite(val).any() else np.nan


def groom_bailey_table(
    sites: Any,
    *,
    band: tuple[float, float] | None = None,
    rotate_deg: float | None = None,
    min_freq: int = 4,
    max_iter: int = 30,
    tol: float = 1e-6,
    robust: bool = True,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    api: bool | None = None,
) -> Any:
    r"""Estimate Groom-Bailey-style galvanic distortion parameters.

    The fitted model is

    .. math:: Z_\mathrm{obs}(f) \approx D\,Z_{2D}(f),

    where ``D`` is a real, frequency-independent 2x2 distortion matrix
    and ``Z_2D`` is anti-diagonal at each frequency.  The fitted matrix is
    decomposed into gain, twist, shear, and anisotropy-style parameters.
    """
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    rows: list[dict[str, Any]] = []
    for i, ed in enumerate(_gb_items(sites, S)):
        station = _name(ed, i)
        Z, z, fr = _get_z_block(ed)
        if Z is None or z is None or fr is None:
            continue
        z_work = np.asarray(z, dtype=np.complex128)
        if rotate_deg is not None:
            z_work = _rotate_tensor(z_work, float(rotate_deg))
        z_fit, fr_fit, _ = _select_band(z_work, fr, band)
        if z_fit.shape[0] < int(min_freq):
            rows.append(dict(
                station=station,
                status="insufficient_frequencies",
                n_freq=int(z_fit.shape[0]),
            ))
            continue
        fit = _fit_gb_distortion(
            z_fit,
            max_iter=max_iter,
            tol=tol,
            robust=robust,
        )
        D = fit["D"]
        params = _gb_parameters(D)
        corrected = fit["corrected"]
        rows.append(dict(
            station=station,
            status="ok",
            n_freq=int(z_fit.shape[0]),
            period_min_s=float(np.nanmin(1.0 / fr_fit)),
            period_max_s=float(np.nanmax(1.0 / fr_fit)),
            rotate_deg=(
                float(rotate_deg) if rotate_deg is not None else np.nan
            ),
            distortion_xx=float(D[0, 0]),
            distortion_xy=float(D[0, 1]),
            distortion_yx=float(D[1, 0]),
            distortion_yy=float(D[1, 1]),
            gain=params["gain"],
            twist_deg=params["twist_deg"],
            shear=params["shear"],
            shear_angle_deg=params["shear_angle_deg"],
            anisotropy=params["anisotropy"],
            rms_fit=float(fit["rms_fit"]),
            diagonal_ratio_before=_diag_ratio(z_fit),
            diagonal_ratio_after=_diag_ratio(corrected),
            robust=bool(robust),
            method="gb_real_distortion_2d",
        ))
    df = pd.DataFrame.from_records(rows)

    return maybe_wrap_frame(
        df,
        api=api,
        name="groom_bailey_table",
        kind="emtools.gb.table",
        source=sites,
        description="Groom-Bailey-style galvanic distortion parameters.",
    )


def _distortion_map(table: pd.DataFrame) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    if table is None or table.empty:
        return out
    for _, row in table.iterrows():
        if str(row.get("status", "ok")) != "ok":
            continue
        D = np.array(
            [
                [row.get("distortion_xx", np.nan),
                 row.get("distortion_xy", np.nan)],
                [row.get("distortion_yx", np.nan),
                 row.get("distortion_yy", np.nan)],
            ],
            dtype=float,
        )
        if np.isfinite(D).all() and abs(np.linalg.det(D)) > 1e-12:
            out[str(row["station"])] = D
    return out


def apply_groom_bailey(
    sites: Any,
    table: pd.DataFrame | None = None,
    *,
    band: tuple[float, float] | None = None,
    rotate_deg: float | None = None,
    min_freq: int = 4,
    max_iter: int = 30,
    tol: float = 1e-6,
    robust: bool = True,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Any:
    """Remove fitted Groom-Bailey galvanic distortion from impedance tensors."""
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    if table is None:
        table = groom_bailey_table(
            S,
            band=band,
            rotate_deg=rotate_deg,
            min_freq=min_freq,
            max_iter=max_iter,
            tol=tol,
            robust=robust,
            recursive=False,
            on_dup=on_dup,
            strict=False,
            verbose=verbose,
        )
    dmap = _distortion_map(table)

    def _one(Si):
        ed = next(_iter_items(Si))
        station = _name(ed, 0)
        D = dmap.get(str(station))
        Z, z, _ = _get_z_block(ed)
        if D is None or Z is None or z is None:
            return Si
        try:
            DI = np.linalg.inv(D)
        except np.linalg.LinAlgError:
            return Si
        z2 = DI[None, :, :] @ np.asarray(z, dtype=np.complex128)
        ze = getattr(Z, "z_err", None)
        try:
            Z.z = z2
            if ze is not None:
                Z.z_err = np.abs(DI)[None, :, :] @ np.asarray(ze, dtype=float)
        except Exception:
            pass
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


def groom_bailey_decomposition(
    sites: Any,
    *,
    apply: bool = False,
    band: tuple[float, float] | None = None,
    rotate_deg: float | None = None,
    min_freq: int = 4,
    max_iter: int = 30,
    tol: float = 1e-6,
    robust: bool = True,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> GroomBaileyResult:
    """Estimate, and optionally apply, Groom-Bailey distortion correction."""
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    table = groom_bailey_table(
        S,
        band=band,
        rotate_deg=rotate_deg,
        min_freq=min_freq,
        max_iter=max_iter,
        tol=tol,
        robust=robust,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )
    out_sites = S
    if apply:
        out_sites = apply_groom_bailey(
            S,
            table=table,
            inplace=inplace,
            recursive=False,
            on_dup=on_dup,
            strict=False,
            verbose=verbose,
        )
    return GroomBaileyResult(
        sites=out_sites,
        table=table.copy() if hasattr(table, "copy") else pd.DataFrame(table),
        applied=bool(apply),
        method="gb_real_distortion_2d",
    )


__all__ = [
    "GroomBaileyResult",
    "groom_bailey_table",
    "apply_groom_bailey",
    "groom_bailey_decomposition",
]
