# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""Adapt pyCSAMT transfer functions to validated Occam1D soundings.

This module is the boundary between heterogeneous site/EDI representations
and :class:`Occam1DData`. It performs component selection, canonical pyCSAMT
impedance conversion, uncertainty propagation, frequency filtering, and
provenance capture. Native text serialization remains in :mod:`.data`.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from os import PathLike
from pathlib import Path
from typing import Any

import numpy as np

from ...z.utils import rho_phi_from_z
from .config import Occam1DConfig
from .data import Occam1DData
from .schema import MODE_OPTIONS

__all__ = ["extract_sounding", "normalise_sites"]

_COMPONENTS = {
    "xy": (0, 1),
    "te": (0, 1),
    "yx": (1, 0),
    "tm": (1, 0),
}


def _first_attr(obj: Any, *names: str):
    """Return the first present, non-``None`` attribute."""
    if obj is None:
        return None
    for name in names:
        value = getattr(obj, name, None)
        if value is not None:
            return value
    return None


def _site_arrays(site):
    """Discover transfer-function arrays across supported site dialects."""
    container = getattr(site, "edi", site)
    z_container = getattr(container, "Z", None)
    sources = (site, container, z_container)

    def discover(*names):
        for source in sources:
            value = _first_attr(source, *names)
            if value is not None:
                return value
        return None

    return {
        "frequency": discover("freq", "frequency"),
        "impedance": discover("z", "impedance"),
        "z_error": discover(
            "z_err", "z_error", "impedance_err", "_z_err"
        ),
        "resistivity": discover("rho", "resistivity", "res"),
        "phase": discover("phase", "phi", "phase_deg"),
        "rho_error": discover(
            "rho_err", "resistivity_err", "_rho_err"
        ),
        "phase_error": discover(
            "phase_err", "phi_err", "_phase_err"
        ),
    }


def _frequency_vector(values) -> np.ndarray:
    """Validate the source frequency vector before tensor conversion."""
    if values is None:
        raise ValueError("The selected site has no frequency vector.")
    try:
        frequency = np.array(values, dtype=float, copy=True)
    except (TypeError, ValueError) as error:
        raise TypeError("site frequency must contain real numbers.") from error
    if frequency.ndim != 1 or not frequency.size:
        raise ValueError(
            "site frequency must be a non-empty one-dimensional vector."
        )
    return frequency


def _tensor(
    array,
    name,
    n_frequency,
    *,
    complex_values=False,
    allow_vector=False,
):
    """Normalize a transfer-function quantity to ``(n, 2, 2)``."""
    if array is None:
        return None
    dtype = complex if complex_values else float
    try:
        values = np.array(array, dtype=dtype, copy=True)
    except (TypeError, ValueError) as error:
        raise TypeError(f"{name} must contain numeric values.") from error
    if allow_vector and values.shape == (n_frequency,):
        return values
    if values.shape == (n_frequency, 4):
        values = values.reshape(n_frequency, 2, 2)
    if values.shape == (2, 2) and n_frequency == 1:
        values = values[np.newaxis, ...]
    if values.shape != (n_frequency, 2, 2):
        raise ValueError(
            f"{name} must have shape (n_frequencies, 2, 2) or "
            "(n_frequencies, 4); got "
            f"{values.shape!r}."
        )
    return values


def _mask_missing_offdiagonal(impedance):
    r"""Replace sentinel-filled off-diagonal impedance rows with NaN.

    EDI readers commonly represent a missing transfer-function estimate
    (the file's own ``EMPTY`` sentinel, e.g. ``1e32``) as an exact
    ``0+0j`` rather than NaN. A genuine :math:`Z_{xy}` or :math:`Z_{yx}`
    estimate from real data is never exactly zero to floating-point
    precision, so an exact zero is an unambiguous "no estimate here"
    signal, not a physically tiny impedance. Left unmasked, it silently
    produces a spurious zero determinant apparent resistivity while its
    phase (:func:`numpy.angle` of ``0+0j`` is ``0.0``, which is finite)
    passes the usual "is this observation finite" check -- exactly the
    combination :class:`~.data.Occam1DData` rejects as an invalid,
    non-positive resistivity. Masking to NaN here lets the existing
    finite-value filtering in :func:`extract_sounding` drop the row
    instead, both here and for a single off-diagonal component.
    """
    missing = (impedance[:, 0, 1] == 0) | (impedance[:, 1, 0] == 0)
    if not np.any(missing):
        return impedance
    masked = impedance.copy()
    masked[missing, 0, 1] = complex(np.nan, np.nan)
    masked[missing, 1, 0] = complex(np.nan, np.nan)
    return masked


def _canonical_from_impedance(z, z_error, frequency):
    """Use pyCSAMT's canonical impedance conversion and error propagation."""
    try:
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            return rho_phi_from_z(z, frequency, z_error)
    except Exception as error:
        raise ValueError(
            "Could not convert the site's impedance tensor to apparent "
            "resistivity and phase."
        ) from error


def _relative_rho_error(absolute, resistivity):
    """Convert absolute resistivity error to a relative fraction."""
    if absolute is None:
        return None
    return np.divide(
        np.abs(absolute),
        np.abs(resistivity),
        out=np.full(np.shape(resistivity), np.nan, dtype=float),
        where=np.isfinite(resistivity) & (np.abs(resistivity) > 0),
    )


def _determinant_values(z, z_error, frequency):
    r"""Return determinant response and first-order uncertainties.

    The determinant impedance follows the pyCSAMT/MT convention
    :math:`Z_d=\sqrt{-Z_{xy}Z_{yx}}`. Independent off-diagonal relative
    errors are propagated in quadrature through the square root.
    """
    zxy = z[:, 0, 1]
    zyx = z[:, 1, 0]
    determinant = np.sqrt(-zxy * zyx)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        rho = np.abs(determinant) ** 2 * 0.2 / frequency
    phase = np.angle(determinant, deg=True)
    if z_error is None:
        return rho, phase, None, None
    rel_xy = np.divide(
        np.abs(z_error[:, 0, 1]),
        np.abs(zxy),
        out=np.full(zxy.shape, np.nan, dtype=float),
        where=np.abs(zxy) > 0,
    )
    rel_yx = np.divide(
        np.abs(z_error[:, 1, 0]),
        np.abs(zyx),
        out=np.full(zyx.shape, np.nan, dtype=float),
        where=np.abs(zyx) > 0,
    )
    relative_z = 0.5 * np.sqrt(rel_xy**2 + rel_yx**2)
    rho_error = 2.0 * relative_z
    phase_error = np.degrees(np.arctan(relative_z))
    return rho, phase, rho_error, phase_error


def _floor_errors(values, floor, observed):
    """Apply a positive floor only where an observation is usable."""
    result = np.full(np.shape(observed), np.nan, dtype=float)
    mask = np.isfinite(observed)
    if values is None:
        result[mask] = floor
        return result
    errors = np.abs(np.asarray(values, dtype=float))
    valid = mask & np.isfinite(errors)
    result[valid] = np.maximum(errors[valid], floor)
    result[mask & ~valid] = floor
    return result


def _component_values(arrays, frequency, row, column):
    """Resolve component observations, preferring supplied physical arrays."""
    z = arrays["impedance"]
    z_error = arrays["z_error"]
    calculated = (None, None, None, None)
    if z is not None:
        calculated = _canonical_from_impedance(z, z_error, frequency)
    (
        calculated_rho,
        calculated_phase,
        calculated_rho_error,
        calc_phase_error,
    ) = calculated

    rho_tensor = arrays["resistivity"]
    phase_tensor = arrays["phase"]

    def component(values):
        if values is None:
            return None
        return values if values.ndim == 1 else values[:, row, column]

    rho = (
        component(rho_tensor)
        if rho_tensor is not None
        else None if calculated_rho is None else calculated_rho[:, row, column]
    )
    phase = (
        component(phase_tensor)
        if phase_tensor is not None
        else None
        if calculated_phase is None
        else calculated_phase[:, row, column]
    )
    if rho is None or phase is None:
        raise ValueError(
            "The selected component lacks resistivity/phase observations and "
            "cannot be derived from impedance."
        )

    supplied_rho_error = arrays["rho_error"]
    supplied_phase_error = arrays["phase_error"]
    if supplied_rho_error is not None:
        rho_error = _relative_rho_error(
            component(supplied_rho_error), rho
        )
        rho_source = "supplied absolute resistivity error"
    elif calculated_rho_error is not None:
        rho_error = _relative_rho_error(
            calculated_rho_error[:, row, column],
            calculated_rho[:, row, column],
        )
        rho_source = "propagated impedance error"
    else:
        rho_error = None
        rho_source = "configured floor"
    if supplied_phase_error is not None:
        phase_error = np.abs(component(supplied_phase_error))
        phase_source = "supplied phase error"
    elif calc_phase_error is not None:
        phase_error = np.abs(calc_phase_error[:, row, column])
        phase_source = "propagated impedance error"
    else:
        phase_error = None
        phase_source = "configured floor"
    return rho, phase, rho_error, phase_error, rho_source, phase_source


def extract_sounding(site, config=None, *, mode=None) -> Occam1DData:
    r"""Extract one component or determinant sounding from a site object.

    Parameters
    ----------
    site : object
        A :class:`pycsamt.site.Site`, EDI-like object, or compatible object
        exposing frequency and either impedance or resistivity/phase tensors.
    config : Occam1DConfig, optional
        Validated frequency bounds, component, and uncertainty floors. A fresh
        default configuration is used when omitted.
    mode : {"xy", "yx", "te", "tm", "determinant"}, optional
        Case-insensitive override of ``config.mode``. ``te`` aliases ``xy``;
        ``tm`` aliases ``yx``.

    Returns
    -------
    Occam1DData
        Frequency-descending sounding with physical units, aligned positive
        uncertainties, and extraction provenance in ``metadata``.

    Raises
    ------
    TypeError
        If ``config`` or numerical source values have unsupported types.
    ValueError
        If tensors are malformed, the component is unavailable, frequencies
        are duplicated, or filtering leaves no usable observations.

    Notes
    -----
    Supplied physical resistivity/phase values take precedence over values
    derived from impedance. Supplied ``rho_error`` is interpreted according
    to the pyCSAMT convention as absolute error in ohm metres and converted to
    the relative fraction required by :class:`Occam1DData`.

    Determinant phase is mapped to the acute passive-earth branch and YX/TM
    phase modulo 180 degrees.  This is essential for the determinant:
    ``numpy.sqrt`` returns the principal square-root branch and EDI impedance
    sign conventions can represent the same passive response with a negative
    phase.  The 1-D forward kernel uses the equivalent ``[0, 90]`` branch.
    Missing/invalid uncertainties receive the configured floor only where the
    corresponding observation is finite.
    """
    if config is not None and not isinstance(config, Occam1DConfig):
        raise TypeError("config must be an Occam1DConfig instance or None.")
    cfg = config or Occam1DConfig()
    cfg.validate()
    if mode is not None and not isinstance(mode, str):
        raise TypeError("mode must be a string or None.")
    selected_mode = (mode or cfg.mode).strip().lower()
    if selected_mode not in MODE_OPTIONS:
        choices = ", ".join(sorted(MODE_OPTIONS))
        raise ValueError(f"mode must be one of: {choices}.")

    raw = _site_arrays(site)
    frequency = _frequency_vector(raw["frequency"])
    n_frequency = frequency.size
    arrays = {
        "impedance": _tensor(
            raw["impedance"], "impedance", n_frequency, complex_values=True
        ),
        "z_error": _tensor(raw["z_error"], "z_error", n_frequency),
        "resistivity": _tensor(
            raw["resistivity"],
            "resistivity",
            n_frequency,
            allow_vector=True,
        ),
        "phase": _tensor(
            raw["phase"], "phase", n_frequency, allow_vector=True
        ),
        "rho_error": _tensor(
            raw["rho_error"],
            "rho_error",
            n_frequency,
            allow_vector=True,
        ),
        "phase_error": _tensor(
            raw["phase_error"],
            "phase_error",
            n_frequency,
            allow_vector=True,
        ),
    }
    if arrays["impedance"] is not None:
        arrays["impedance"] = _mask_missing_offdiagonal(arrays["impedance"])
    if arrays["z_error"] is not None and arrays["impedance"] is None:
        raise ValueError("z_error was supplied without an impedance tensor.")

    if selected_mode == "determinant":
        if arrays["impedance"] is None:
            raise ValueError(
                "Determinant mode requires the complex impedance tensor."
            )
        rho, phase, rho_error, phase_error = _determinant_values(
            arrays["impedance"], arrays["z_error"], frequency
        )
        error_source = (
            "propagated impedance error"
            if arrays["z_error"] is not None
            else "configured floor"
        )
        rho_source = phase_source = error_source
        phase = np.abs(phase)
    else:
        row, column = _COMPONENTS[selected_mode]
        (
            rho,
            phase,
            rho_error,
            phase_error,
            rho_source,
            phase_source,
        ) = _component_values(arrays, frequency, row, column)
        if selected_mode in {"yx", "tm"}:
            phase = np.mod(phase, 180.0)

    rho = np.asarray(rho, dtype=float)
    phase = np.asarray(phase, dtype=float)
    rho_error = _floor_errors(rho_error, cfg.error_floor_rho, rho)
    phase_error = _floor_errors(
        phase_error, cfg.error_floor_phase, phase
    )
    usable = (
        (np.isfinite(rho) & (rho > 0)) | np.isfinite(phase)
    )
    selected = np.isfinite(frequency) & (frequency > 0) & usable
    if cfg.freq_min is not None:
        selected &= frequency >= cfg.freq_min
    if cfg.freq_max is not None:
        selected &= frequency <= cfg.freq_max
    if not np.any(selected):
        raise ValueError(
            "Frequency and observation filtering removed every datum."
        )
    chosen_frequency = frequency[selected]
    if np.unique(chosen_frequency).size != chosen_frequency.size:
        raise ValueError(
            "Selected site frequencies must be unique for Occam1D."
        )
    order = np.argsort(chosen_frequency, kind="stable")[::-1]
    name = _first_attr(site, "name", "station", "dataid") or "site"
    metadata = {
        "source_class": type(site).__qualname__,
        "source_module": type(site).__module__,
        "source_frequency_count": n_frequency,
        "selected_frequency_count": int(np.count_nonzero(selected)),
        "resistivity_error_source": rho_source,
        "phase_error_source": phase_source,
        "frequency_order": "descending",
    }
    return Occam1DData(
        chosen_frequency[order],
        rho[selected][order],
        phase[selected][order],
        rho_error[selected][order],
        phase_error[selected][order],
        mode=selected_mode,
        station=str(name),
        metadata=metadata,
    )


def normalise_sites(source) -> list:
    """Normalize supported single-site and collection inputs to a list.

    Parameters
    ----------
    source : object, iterable, mapping, or path-like
        A site, EDI-like object, :class:`pycsamt.site.Sites`, mapping of names
        to sites, iterable of sites, or source accepted by ``Sites.from_any``.

    Returns
    -------
    list
        Materialized site objects in source order.

    Raises
    ------
    TypeError
        If ``source`` is ``None`` or cannot be normalized.
    ValueError
        If a path-like source cannot be loaded as site data.
    """
    if source is None:
        raise TypeError("source cannot be None.")
    if isinstance(source, (str, bytes, PathLike)):
        path = Path(source).expanduser()
        if path.is_file() and path.suffix.casefold() == ".edi":
            from ...seg.edi import EDIFile

            try:
                return [EDIFile(path)]
            except Exception as error:
                raise ValueError(
                    f"Could not read EDI site data from {source!r}."
                ) from error
        from ...site import Sites

        try:
            return list(Sites.from_any(source))
        except Exception as error:
            raise ValueError(
                f"Could not load site data from {source!r}."
            ) from error
    if hasattr(source, "freq") or hasattr(source, "edi"):
        return [source]
    if hasattr(source, "Z") and not hasattr(source, "_items"):
        from ...site.base import Site

        try:
            return [Site(source)]
        except Exception as error:
            raise ValueError(
                "Could not adapt the EDI-like source to a Site object."
            ) from error
    if isinstance(source, Mapping):
        return list(source.values())
    if hasattr(source, "_items") or hasattr(source, "edic"):
        return list(source)
    if isinstance(source, Iterable):
        return list(source)
    raise TypeError(
        "source must be a site, site collection, iterable, mapping, or path."
    )
