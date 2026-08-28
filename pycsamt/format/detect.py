# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.format.detect
=====================

Smart source detection for the PCSF / PCSM conversion pipeline.

:func:`detect_source` looks at a path — a single file *or* a whole solver
working directory — and reports, without loading any heavy array data,
what it is and how it should be converted:

===============  ==========================================================
``category``     meaning
===============  ==========================================================
``"pcsf"``       an existing ``.pcsf`` (HDF5) container
``"pcsm"``       an existing ``.pcsm`` / ``.pcsm.gz`` ASCII sibling
``"solver"``     a classical-solver result (``backend`` is one of
                 ``occam2d`` / ``modem`` / ``mare2dem``)
``"ai_arrays"``  a plain array bundle (``.npz`` / ``.npy``) written by an
                 AI/DL inversion — routed through
                 :mod:`pycsamt.format.adapters.generic`
``"unknown"``    nothing recognisable
===============  ==========================================================

The result also carries a ``target_geometry`` (``grid2d`` / ``grid3d`` /
``mesh_unstructured``) for ``solver`` / ``ai_arrays`` sources, the peeked
``geometry`` for ``pcsf`` / ``pcsm`` sources, and a ``hints`` dict with
any concrete paths / array keys the converter will need (e.g. the
MARE2DEM ``.poly`` PSLG, or the resolved resistivity key inside an
``.npz``).

This module is deliberately dependency-light: it never imports a solver
package, ``h5py``, or ``triangle`` — only :func:`pycsamt.format.peek_kind`
(pure text/attr read) and, for ``.npz`` inspection, ``numpy``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

__all__ = [
    "SourceKind",
    "detect_source",
    "describe_source",
    "AI_RESISTIVITY_KEYS",
    "AI_LOG10_RESISTIVITY_KEYS",
    "AI_LN_RESISTIVITY_KEYS",
]

# ---------------------------------------------------------------------------
# Extension groups
# ---------------------------------------------------------------------------

_PCSF_SUFFIXES = (".pcsf",)
_PCSM_SUFFIXES = (".pcsm", ".pcsm.gz")
_AI_ARRAY_SUFFIXES = (".npz", ".npy")

# ---------------------------------------------------------------------------
# Solver fingerprints (case-insensitive filename / suffix signatures)
# ---------------------------------------------------------------------------

_OCCAM_NAME_SIGNATURES = frozenset(
    {
        "occamdatafile.dat",
        "occam2dmesh",
        "occam2dmodel",
        "startup",
        "occamstartup",
        "occam2dstartup",
        "occamdata.dat",
        "iter.dat",
    }
)
_OCCAM_SUFFIX_SIGNATURES = frozenset({".iter"})

_MODEM_NAME_SIGNATURES = frozenset(
    {
        "modemdata.dat",
        "modem.inv",
        "modular_nlcg.log",
        "modem.cov",
    }
)
_MODEM_SUFFIX_SIGNATURES = frozenset({".rho", ".rms", ".inv"})

_MARE_NAME_SIGNATURES = frozenset({"mare2dem.settings"})
_MARE_SUFFIX_SIGNATURES = frozenset({".poly", ".resistivity", ".emdata"})

_BACKEND_GEOMETRY = {
    "occam2d": "grid2d",
    "modem": "grid3d",
    "mare2dem": "mesh_unstructured",
}

# ---------------------------------------------------------------------------
# AI array-key aliases (lower-cased; matched case-insensitively)
# ---------------------------------------------------------------------------

AI_RESISTIVITY_KEYS = (
    "resistivity",
    "rho",
    "rho_ohm_m",
    "res",
    "resistivity_ohm_m",
    "model",
    "m",
    "prediction",
    "pred",
    "y_pred",
    "output",
)
AI_LOG10_RESISTIVITY_KEYS = (
    "log10_resistivity",
    "log10_rho",
    "log10rho",
    "log_resistivity",
    "log_rho",
    "logrho",
    "log10_model",
)
AI_LN_RESISTIVITY_KEYS = (
    "ln_resistivity",
    "ln_rho",
    "lnrho",
    "log_e_rho",
)

_X_KEYS = ("x", "x_centers", "x_center", "xc", "cell_x", "easting", "x_m")
_Y_KEYS = ("y", "y_centers", "y_center", "yc", "cell_y", "northing", "y_m")
_Z_KEYS = ("z", "z_centers", "z_center", "zc", "cell_z", "depth", "z_m")
_XN_KEYS = ("x_nodes", "x_edges", "xn", "x_grid")
_YN_KEYS = ("y_nodes", "y_edges", "yn", "y_grid")
_ZN_KEYS = ("z_nodes", "z_edges", "zn", "z_grid")
_NODES_KEYS = ("nodes", "vertices", "points", "xy", "node_xy", "coords")
_CONN_KEYS = (
    "connectivity",
    "triangles",
    "tris",
    "elements",
    "cells",
    "simplices",
    "faces",
)
_REGION_KEYS = ("region_ids", "regions", "region", "labels", "cell_region")
_UNCERT_KEYS = ("uncertainty", "std", "sigma", "stddev", "uncert")
_SENS_KEYS = ("sensitivity", "jacobian_diag", "sens")


# ---------------------------------------------------------------------------
# Result type
# ---------------------------------------------------------------------------


@dataclass
class SourceKind:
    """Outcome of :func:`detect_source`.

    Attributes
    ----------
    category : str
        ``"pcsf"`` | ``"pcsm"`` | ``"solver"`` | ``"ai_arrays"`` |
        ``"unknown"``.
    path : pathlib.Path
        The path that was probed (a file, or a solver working directory).
    is_dir : bool
        Whether *path* is a directory.
    backend : str or None
        ``"occam2d"`` | ``"modem"`` | ``"mare2dem"`` for ``solver``
        sources, else ``None``.
    geometry : str or None
        For ``pcsf`` / ``pcsm``: the peeked geometry kind.
    target_geometry : str or None
        For ``solver`` / ``ai_arrays``: the PCSF geometry the conversion
        will produce.
    confidence : str
        ``"high"`` | ``"medium"`` | ``"low"``.
    detail : str
        Human-readable one-line explanation.
    hints : dict
        Extra paths / keys the converter needs (see module docstring).
    """

    category: str
    path: Path
    is_dir: bool
    backend: str | None = None
    geometry: str | None = None
    target_geometry: str | None = None
    confidence: str = "high"
    detail: str = ""
    hints: dict[str, Any] = field(default_factory=dict)

    @property
    def convertible(self) -> bool:
        """True when this source can be fed to ``pycsamt format convert``."""
        return self.category in {"pcsf", "pcsm", "solver", "ai_arrays"}

    def to_dict(self) -> dict[str, Any]:
        """JSON-friendly plain-dict view."""
        return {
            "category": self.category,
            "path": str(self.path),
            "is_dir": self.is_dir,
            "backend": self.backend,
            "geometry": self.geometry,
            "target_geometry": self.target_geometry,
            "confidence": self.confidence,
            "detail": self.detail,
            "hints": {k: _jsonify(v) for k, v in self.hints.items()},
        }


def _jsonify(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (list, tuple)):
        return [_jsonify(v) for v in value]
    if isinstance(value, dict):
        return {k: _jsonify(v) for k, v in value.items()}
    return value


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def detect_source(
    path: str | Path,
    *,
    solver_hint: str | None = None,
) -> SourceKind:
    """Detect what *path* is and how it should convert to PCSF/PCSM.

    Parameters
    ----------
    path : path-like
        A ``.pcsf`` / ``.pcsm`` / ``.npz`` / ``.npy`` file, a
        solver-specific file (``.iter``, ``.rho``, ``.poly``,
        ``.resistivity``, ...), or a solver working directory.
    solver_hint : {"occam2d", "modem", "mare2dem"}, optional
        Force the solver backend instead of fingerprinting it. Ignored
        for PCSF/PCSM/array sources.

    Returns
    -------
    SourceKind

    Raises
    ------
    FileNotFoundError
        If *path* does not exist.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"No such file or directory: {p}")

    if solver_hint is not None:
        solver_hint = solver_hint.lower()
        if solver_hint not in _BACKEND_GEOMETRY:
            raise ValueError(
                f"Unknown solver_hint {solver_hint!r}; expected one of "
                f"{sorted(_BACKEND_GEOMETRY)}."
            )

    if p.is_dir():
        return _detect_dir(p, solver_hint)
    return _detect_file(p, solver_hint)


# ---------------------------------------------------------------------------
# File detection
# ---------------------------------------------------------------------------


def _lname(p: Path) -> str:
    return p.name.lower()


def _has_pcsm_suffix(name: str) -> bool:
    return name.endswith(".pcsm") or name.endswith(".pcsm.gz")


def _detect_file(p: Path, solver_hint: str | None) -> SourceKind:
    name = _lname(p)
    suffix = p.suffix.lower()

    # -- PCSF / PCSM ------------------------------------------------------
    if name.endswith(".pcsf"):
        return _peek_pcsf_or_pcsm(p, "pcsf")
    if _has_pcsm_suffix(name):
        return _peek_pcsf_or_pcsm(p, "pcsm")

    # -- AI array bundles ----------------------------------------------
    if suffix in _AI_ARRAY_SUFFIXES:
        return _detect_ai_arrays(p)

    # -- Solver-specific single files -> classify via the parent dir ---
    if solver_hint is not None or suffix in (
        _OCCAM_SUFFIX_SIGNATURES
        | _MODEM_SUFFIX_SIGNATURES
        | _MARE_SUFFIX_SIGNATURES
        | {".dat", ".settings", ".log"}
    ):
        backend = solver_hint or _backend_from_suffix(suffix)
        parent = _detect_dir(p.parent, solver_hint or backend)
        if parent.category == "solver":
            parent.detail = (
                f"{parent.detail} (pointed at {p.name}; using its folder "
                f"{p.parent})"
            )
            parent.hints.setdefault("picked_file", p)
            if backend == "mare2dem" and suffix == ".poly":
                parent.hints["poly"] = p
            return parent
        if backend is not None:
            return SourceKind(
                category="solver",
                path=p.parent,
                is_dir=True,
                backend=backend,
                target_geometry=_BACKEND_GEOMETRY[backend],
                confidence="low",
                detail=(
                    f"{p.name} looks like a {backend} file but its folder "
                    f"has no corroborating signature."
                ),
                hints={"picked_file": p}
                | ({"poly": p} if suffix == ".poly" else {}),
            )

    return SourceKind(
        category="unknown",
        path=p,
        is_dir=False,
        confidence="low",
        detail=f"Unrecognised file type: {p.suffix or p.name!r}.",
    )


def _backend_from_suffix(suffix: str) -> str | None:
    if suffix in _OCCAM_SUFFIX_SIGNATURES:
        return "occam2d"
    if suffix in _MODEM_SUFFIX_SIGNATURES:
        return "modem"
    if suffix in _MARE_SUFFIX_SIGNATURES:
        return "mare2dem"
    return None


def _peek_pcsf_or_pcsm(p: Path, category: str) -> SourceKind:
    from .text import peek_kind

    label = category.upper()
    try:
        kind = peek_kind(p)
        conf = "high"
        detail = f"{label} container, geometry={kind}."
    except Exception as exc:  # noqa: BLE001 - a corrupt file must still classify
        kind = None
        conf = "low"
        detail = f"{label} file, but its header could not be read: {exc}"
    return SourceKind(
        category=category,
        path=p,
        is_dir=False,
        geometry=kind,
        confidence=conf,
        detail=detail,
    )


# ---------------------------------------------------------------------------
# Directory detection
# ---------------------------------------------------------------------------


def _detect_dir(d: Path, solver_hint: str | None) -> SourceKind:
    try:
        entries = [x for x in d.iterdir() if x.is_file()]
    except OSError as exc:
        return SourceKind(
            category="unknown",
            path=d,
            is_dir=True,
            confidence="low",
            detail=f"Cannot list directory: {exc}",
        )

    names = {x.name.lower() for x in entries}
    suffixes = {x.suffix.lower() for x in entries}

    if solver_hint is not None:
        return _solver_dir_result(d, solver_hint, entries, forced=True)

    matches: list[str] = []
    if names & _MARE_NAME_SIGNATURES or (
        ".poly" in suffixes and ".resistivity" in suffixes
    ):
        matches.append("mare2dem")
    elif ".resistivity" in suffixes and ".emdata" in suffixes:
        matches.append("mare2dem")

    if names & _MODEM_NAME_SIGNATURES or (
        ".rho" in suffixes and ".dat" in suffixes
    ):
        matches.append("modem")

    if names & _OCCAM_NAME_SIGNATURES or (suffixes & _OCCAM_SUFFIX_SIGNATURES):
        matches.append("occam2d")

    if not matches:
        # A folder that simply holds one PCSF/PCSM file.
        pcs = [
            x
            for x in entries
            if x.name.lower().endswith((".pcsf", ".pcsm", ".pcsm.gz"))
        ]
        if len(pcs) == 1:
            cat = "pcsm" if _has_pcsm_suffix(pcs[0].name.lower()) else "pcsf"
            return _peek_pcsf_or_pcsm(pcs[0], cat)
        return SourceKind(
            category="unknown",
            path=d,
            is_dir=True,
            confidence="low",
            detail=(
                "No Occam2D / ModEM / MARE2DEM signature found in "
                f"{d}. Pass --solver to force one."
            ),
        )

    if len(matches) > 1:
        backend = matches[0]
        res = _solver_dir_result(d, backend, entries, forced=False)
        res.confidence = "medium"
        res.detail = (
            f"Multiple solver signatures found ({', '.join(matches)}); "
            f"assuming {backend}. Pass --solver to override."
        )
        return res

    return _solver_dir_result(d, matches[0], entries, forced=False)


def _solver_dir_result(
    d: Path,
    backend: str,
    entries: list[Path],
    *,
    forced: bool,
) -> SourceKind:
    hints: dict[str, Any] = {}
    if backend == "mare2dem":
        polys = sorted(x for x in entries if x.suffix.lower() == ".poly")
        if polys:
            hints["poly"] = polys[0]
        res_files = sorted(
            x for x in entries if x.suffix.lower() == ".resistivity"
        )
        if res_files:
            hints["resistivity_files"] = res_files
    elif backend == "occam2d":
        iters = sorted(x for x in entries if x.suffix.lower() == ".iter")
        if iters:
            hints["iteration_files"] = iters
    elif backend == "modem":
        rhos = sorted(x for x in entries if x.suffix.lower() == ".rho")
        if rhos:
            hints["rho_files"] = rhos

    detail = (
        f"Forced solver backend: {backend}."
        if forced
        else f"Detected {backend} working directory."
    )
    return SourceKind(
        category="solver",
        path=d,
        is_dir=True,
        backend=backend,
        target_geometry=_BACKEND_GEOMETRY[backend],
        confidence="high" if not forced else "high",
        detail=detail,
        hints=hints,
    )


# ---------------------------------------------------------------------------
# AI array-bundle detection
# ---------------------------------------------------------------------------


def _match_key(
    available: dict[str, str], aliases: tuple[str, ...]
) -> str | None:
    """Return the first real key whose lower-case form is in *aliases*."""
    for alias in aliases:
        if alias in available:
            return available[alias]
    return None


def _detect_ai_arrays(p: Path) -> SourceKind:
    try:
        import numpy as np
    except Exception as exc:  # pragma: no cover - numpy is a hard dep
        return SourceKind(
            category="unknown",
            path=p,
            is_dir=False,
            confidence="low",
            detail=f"numpy unavailable, cannot inspect {p.name}: {exc}",
        )

    if p.suffix.lower() == ".npy":
        try:
            arr = np.load(p, allow_pickle=False, mmap_mode="r")
        except Exception as exc:  # noqa: BLE001
            return SourceKind(
                category="unknown",
                path=p,
                is_dir=False,
                confidence="low",
                detail=f"Could not read {p.name}: {exc}",
            )
        ndim = arr.ndim
        geom = "grid3d" if ndim == 3 else "grid2d" if ndim == 2 else None
        if geom is None:
            return SourceKind(
                category="unknown",
                path=p,
                is_dir=False,
                confidence="low",
                detail=(
                    f"{p.name} holds a {ndim}-D array; expected a 2-D "
                    "(grid2d) or 3-D (grid3d) resistivity volume."
                ),
            )
        return SourceKind(
            category="ai_arrays",
            path=p,
            is_dir=False,
            target_geometry=geom,
            confidence="low",
            detail=(
                f"Bare .npy {ndim}-D array assumed to be linear "
                f"resistivity on a unit-spaced {geom}. Prefer an .npz "
                "carrying explicit coordinate arrays."
            ),
            hints={
                "resistivity_key": None,
                "encoding": "linear",
                "synthetic_coords": True,
            },
        )

    # -- .npz ---------------------------------------------------------
    try:
        with np.load(p, allow_pickle=False) as npz:
            files = list(npz.files)
            shapes = {k: tuple(npz[k].shape) for k in files}
    except Exception as exc:  # noqa: BLE001
        return SourceKind(
            category="unknown",
            path=p,
            is_dir=False,
            confidence="low",
            detail=f"Could not read {p.name} as a .npz bundle: {exc}",
        )

    lower = {k.lower(): k for k in files}

    rho_key = _match_key(lower, AI_RESISTIVITY_KEYS)
    encoding = "linear"
    if rho_key is None:
        rho_key = _match_key(lower, AI_LOG10_RESISTIVITY_KEYS)
        if rho_key is not None:
            encoding = "log10"
    if rho_key is None:
        rho_key = _match_key(lower, AI_LN_RESISTIVITY_KEYS)
        if rho_key is not None:
            encoding = "ln"

    if rho_key is None:
        return SourceKind(
            category="unknown",
            path=p,
            is_dir=False,
            confidence="low",
            detail=(
                f"{p.name} has no recognisable resistivity array. Keys: "
                f"{', '.join(files)}. Expected one of "
                f"{', '.join(AI_RESISTIVITY_KEYS[:5])}, ..."
            ),
        )

    nodes_key = _match_key(lower, _NODES_KEYS)
    conn_key = _match_key(lower, _CONN_KEYS)
    x_key = _match_key(lower, _X_KEYS)
    y_key = _match_key(lower, _Y_KEYS)
    z_key = _match_key(lower, _Z_KEYS)

    rho_ndim = len(shapes[rho_key])

    hints: dict[str, Any] = {
        "resistivity_key": rho_key,
        "encoding": encoding,
        "keys": files,
        "shapes": {k: list(v) for k, v in shapes.items()},
        "uncertainty_key": _match_key(lower, _UNCERT_KEYS),
        "sensitivity_key": _match_key(lower, _SENS_KEYS),
    }

    if nodes_key and conn_key:
        hints.update(
            nodes_key=nodes_key,
            connectivity_key=conn_key,
            region_ids_key=_match_key(lower, _REGION_KEYS),
        )
        return SourceKind(
            category="ai_arrays",
            path=p,
            is_dir=False,
            target_geometry="mesh_unstructured",
            confidence="high",
            detail=(
                f"AI array bundle: mesh ({nodes_key}+{conn_key}), "
                f"resistivity={rho_key} ({encoding})."
            ),
            hints=hints,
        )

    if (x_key and y_key and z_key) or rho_ndim == 3:
        hints.update(
            x_key=x_key,
            y_key=y_key,
            z_key=z_key,
            x_nodes_key=_match_key(lower, _XN_KEYS),
            y_nodes_key=_match_key(lower, _YN_KEYS),
            z_nodes_key=_match_key(lower, _ZN_KEYS),
        )
        have_coords = bool(x_key and y_key and z_key)
        hints["synthetic_coords"] = not have_coords
        return SourceKind(
            category="ai_arrays",
            path=p,
            is_dir=False,
            target_geometry="grid3d",
            confidence="high" if have_coords else "medium",
            detail=(
                f"AI array bundle: grid3d, resistivity={rho_key} "
                f"({encoding}), "
                + (
                    f"coords={x_key}/{y_key}/{z_key}."
                    if have_coords
                    else "no coordinate arrays — unit spacing assumed."
                )
            ),
            hints=hints,
        )

    if (x_key and z_key) or rho_ndim == 2:
        hints.update(
            x_key=x_key,
            z_key=z_key,
            x_nodes_key=_match_key(lower, _XN_KEYS),
            z_nodes_key=_match_key(lower, _ZN_KEYS),
        )
        have_coords = bool(x_key and z_key)
        hints["synthetic_coords"] = not have_coords
        return SourceKind(
            category="ai_arrays",
            path=p,
            is_dir=False,
            target_geometry="grid2d",
            confidence="high" if have_coords else "medium",
            detail=(
                f"AI array bundle: grid2d, resistivity={rho_key} "
                f"({encoding}), "
                + (
                    f"coords={x_key}/{z_key}."
                    if have_coords
                    else "no coordinate arrays — unit spacing assumed."
                )
            ),
            hints=hints,
        )

    return SourceKind(
        category="unknown",
        path=p,
        is_dir=False,
        confidence="low",
        detail=(
            f"{p.name}: found resistivity ({rho_key}, {rho_ndim}-D) but "
            "could not infer a geometry (no coord arrays, no mesh keys)."
        ),
    )


# ---------------------------------------------------------------------------
# Pretty description
# ---------------------------------------------------------------------------


def describe_source(sk: SourceKind) -> str:
    """Return a short multi-line human summary of *sk*."""
    lines = [
        f"path        : {sk.path}",
        f"category    : {sk.category}",
    ]
    if sk.backend:
        lines.append(f"backend     : {sk.backend}")
    if sk.geometry:
        lines.append(f"geometry    : {sk.geometry}")
    if sk.target_geometry:
        lines.append(f"-> geometry : {sk.target_geometry}")
    lines.append(f"confidence  : {sk.confidence}")
    lines.append(f"detail      : {sk.detail}")
    if sk.hints:
        keyhints = ", ".join(
            f"{k}={_jsonify(v)}"
            for k, v in sk.hints.items()
            if v is not None and k not in {"keys", "shapes"}
        )
        if keyhints:
            lines.append(f"hints       : {keyhints}")
    return "\n".join(lines)
