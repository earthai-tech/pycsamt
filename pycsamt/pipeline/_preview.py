"""Random-station raw-vs-processed preview for method-aware presets.

:func:`pick_preview_stations` deterministically samples a handful of station
names from a ``Sites`` collection; :func:`plot_raw_preview` and
:func:`plot_processed_preview` reuse the existing
:func:`pycsamt.emtools.plot.plot_raw_sites_1d` (which already supports both
raw and processed rendering via its ``raw`` flag — no new plotting code was
needed) to draw the same stations before and after processing.

Both plot functions are registered as QC hooks on the ``preview`` category
steps (``PRE001``/``PRE002``) — placed first and last in a preset
respectively — but are also plain, directly callable functions for anyone
who wants a custom station count or explicit station list outside the
pipeline, since a QC hook is always invoked as ``fn(sites)`` with no way to
forward extra keyword arguments (see :mod:`pycsamt.pipeline._smart_qc` for
the same constraint).

Both steps re-derive the station selection independently with the same
default seed, so if an intermediate processing step drops one of the
previewed stations, the "after" plot's reselection (drawn from the smaller
surviving set, same seed) can end up choosing a different station than the
"before" plot. This is a deliberate trade-off: solving it exactly would
require carrying state across step boundaries, which the pipeline's
otherwise-stateless step-execution model does not support.
"""

from __future__ import annotations

import random
from typing import Any

from ..emtools._core import _name, ensure_sites, _iter_items
from ..emtools.plot import plot_raw_sites_1d

__all__ = ["pick_preview_stations", "plot_raw_preview", "plot_processed_preview"]


def pick_preview_stations(
    sites: Any, *, n: int = 3, random_state: int = 0
) -> list[str]:
    """Deterministically pick up to *n* station names from *sites*.

    Uses the same station-name resolution as every other ``emtools``
    plotting function (:func:`pycsamt.emtools._core._name`), so the
    returned names match what ``plot_raw_sites_1d(sites, stations=...)``
    expects.
    """
    S = ensure_sites(sites, strict=False)
    names = [_name(ed, i) for i, ed in enumerate(_iter_items(S))]
    n = min(max(0, n), len(names))
    return sorted(random.Random(random_state).sample(names, n))


def plot_raw_preview(sites: Any, *, n: int = 3, random_state: int = 0) -> Any:
    """Raw 1-D rho/phase panels for a deterministic random station subset."""
    stations = pick_preview_stations(sites, n=n, random_state=random_state)
    return plot_raw_sites_1d(sites, stations=stations, raw=True)


def plot_processed_preview(sites: Any, *, n: int = 3, random_state: int = 0) -> Any:
    """Processed 1-D rho/phase panels for the same deterministic station subset."""
    stations = pick_preview_stations(sites, n=n, random_state=random_state)
    return plot_raw_sites_1d(sites, stations=stations, raw=False)
