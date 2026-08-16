"""Data-driven ("smart") QC plot gating for method-aware presets.

A plain :class:`~pycsamt.pipeline.StepSpec` QC hook is called as ``fn(sites)``
only (see :meth:`~pycsamt.pipeline.Step.generate_qc_plots`) — there is no
channel to pass extra parameters through. The wrappers in this module use
that single-argument contract to decide, from the data itself, whether a
plot is actually meaningful before calling it:

* Phase-tensor-ellipse plotting only makes sense once both off-diagonal
  impedance components (Zxy and Zyx) are populated on at least one station
  — a single-component TE/TM-only CSAMT line cannot produce a meaningful
  ellipse.
* Tipper plots only make sense once a real (non-empty) vertical-field
  channel is present on at least one station — most AMT surveys carry none.

Each wrapper returns ``None`` when its precondition isn't met.
:meth:`~pycsamt.pipeline.Step.generate_qc_plots` already treats a ``None``
return as "skip this plot, no figure saved, no error" — the same mechanism
every other QC hook relies on — so no changes to the step-execution engine
were needed to support this.

The detection helpers below (``_any_multicomponent``, ``_any_tipper``) reuse
the exact low-level idioms already used elsewhere in ``emtools`` for the
same checks (``_get_z_block``/direct Zxy·Zyx presence in
``fieldzone.py``/``source_effects.py``; ``_get_t_block`` in
``strike.py``/``inspect.py``), rather than the higher-level
``SiteMixin.has_component``/``quality_flags`` API, so behaviour matches what
``correct_near_field`` and ``plot_strike_analysis`` already consider
"present" on the exact same data.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from ..emtools._core import _get_t_block, _get_z_block, _iter_items, ensure_sites

__all__ = [
    "phase_tensor_smart",
    "induction_multiperiod_map_smart",
    "induction_section_smart",
    "response_tipper_smart",
    "tipper_components_smart",
    "tipper_hodograms_smart",
]


def _any_multicomponent(sites: Any) -> bool:
    """True if any station has both Zxy and Zyx populated (finite, nonzero)."""
    S = ensure_sites(sites, strict=False)
    for ed in _iter_items(S):
        _Z, z, _fr = _get_z_block(ed)
        if z is None:
            continue
        zxy = np.abs(z[:, 0, 1])
        zyx = np.abs(z[:, 1, 0])
        if np.any((zxy > 0) & np.isfinite(zxy) & (zyx > 0) & np.isfinite(zyx)):
            return True
    return False


def _any_tipper(sites: Any) -> bool:
    """True if any station carries a real (non-empty) tipper channel."""
    S = ensure_sites(sites, strict=False)
    for ed in _iter_items(S):
        _T, t, fr = _get_t_block(ed)
        if t is None or fr is None:
            continue
        return True
    return False


def _first_of(result: Any) -> Any:
    """Unwrap a ``(Figure, axes_array)``-style return down to the figure."""
    if isinstance(result, tuple) and result:
        return result[0]
    return result


# ---------------------------------------------------------------------------
# Phase tensor — gated on multi-component data
# ---------------------------------------------------------------------------


def phase_tensor_smart(sites: Any) -> Any:
    """``plot_phase_tensor_psection`` if the data is multi-component, else None.

    A single-component TE/TM-only CSAMT line has only one off-diagonal Z
    component populated, so a phase-tensor ellipse cannot be constructed
    meaningfully — skip rather than draw a degenerate figure.
    """
    if not _any_multicomponent(sites):
        return None
    from ..emtools.tensor import plot_phase_tensor_psection

    return plot_phase_tensor_psection(sites)


# ---------------------------------------------------------------------------
# Tipper-dependent plots — gated on tipper presence
# ---------------------------------------------------------------------------


def induction_multiperiod_map_smart(sites: Any) -> Any:
    """``plot_induction_multiperiod_map`` if tipper is present, else None."""
    if not _any_tipper(sites):
        return None
    from ..emtools.tf import plot_induction_multiperiod_map

    return _first_of(plot_induction_multiperiod_map(sites))


def induction_section_smart(sites: Any) -> Any:
    """``plot_induction_section`` if tipper is present, else None."""
    if not _any_tipper(sites):
        return None
    from ..emtools.tf import plot_induction_section

    return plot_induction_section(sites)


def response_tipper_smart(sites: Any) -> Any:
    """``plot_response_tipper`` if tipper is present, else None."""
    if not _any_tipper(sites):
        return None
    from ..emtools.plot import plot_response_tipper

    return plot_response_tipper(sites)


def tipper_components_smart(sites: Any) -> Any:
    """``plot_tipper_components`` if tipper is present, else None."""
    if not _any_tipper(sites):
        return None
    from ..emtools.inspect import plot_tipper_components

    return plot_tipper_components(sites)


def tipper_hodograms_smart(sites: Any) -> Any:
    """``plot_tipper_hodograms`` if tipper is present, else None."""
    if not _any_tipper(sites):
        return None
    from ..emtools.tf import plot_tipper_hodograms

    return plot_tipper_hodograms(sites)
