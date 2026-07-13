# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.metadata.geology
========================

Geological formation catalog for MT / AMT / TEM modelling.

The catalog centralises the ``GEOLOGY_PRIORS`` that used to live
in ``forward/synthetic.py`` and extends it with a richer set of
geological scenarios.  Every scenario is described by:

* resistivity range (Ω·m),
* depth range (m),
* n_layers hint for layered-Earth generation,
* free-text description and associated rock types.

Quick start
-----------
::

    from pycsamt.metadata.geology import CATALOG, Formation

    # Access a named scenario
    f = CATALOG.get("sedimentary")
    print(f.rho_mid)         # geometric-mean resistivity
    print(f.to_prior())      # dict compatible with LayeredModel.from_geology()

    # Fuzzy look-up by resistivity
    matches = CATALOG.lookup_by_resistivity(100.0, n=3)

    # All names
    print(CATALOG.names())

    # Add a custom formation
    CATALOG.add(Formation(
        name="custom_ore",
        resistivity_range=(0.01, 10.0),
        depth_range=(100, 500),
        n_layers_range=(2, 4),
        description="Massive sulphide ore body",
        rock_types=["massive sulphide", "graphite"],
    ))

Compatibility with forward.synthetic
-------------------------------------
:func:`geology_prior` returns a dict in the exact format expected by
``LayeredModel.from_geology(name)``::

    from pycsamt.metadata.geology import geology_prior
    prior = geology_prior("sedimentary")   # → used by LayeredModel
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

from ..api.view import maybe_wrap_frame

__all__ = [
    "Formation",
    "GeologyCatalog",
    "CATALOG",
    "geology_prior",
]


# ---------------------------------------------------------------------------
# Formation
# ---------------------------------------------------------------------------


@dataclass
class Formation:
    """A named geological scenario for layered-Earth modelling.

    Parameters
    ----------
    name : str
        Unique identifier (lower-case, e.g. ``"sedimentary"``).
    resistivity_range : tuple[float, float]
        Expected resistivity interval (Ω·m) as *(rho_min, rho_max)*.
    depth_range : tuple[float, float]
        Typical investigation depth range (m) as *(depth_min, depth_max)*.
    n_layers_range : tuple[int, int]
        Plausible number of layers for synthetic generation,
        *(n_min, n_max)*.
    description : str
        Free-text geological description.
    rock_types : list[str]
        Rock type names that correspond to entries in
        :data:`~pycsamt.metadata.rocks.GEO_ROCK_RESISTIVITY`.
    extra : dict
        Arbitrary additional parameters preserved on round-trips.

    Examples
    --------
    ::

        f = Formation(
            name="hydrothermal",
            resistivity_range=(1.0, 1e3),
            depth_range=(200, 3000),
            n_layers_range=(3, 6),
            description="Hydrothermal alteration zone",
        )
        print(f.rho_mid, f.log_rho_range)
    """

    name: str
    resistivity_range: tuple[float, float]
    depth_range: tuple[float, float]
    n_layers_range: tuple[int, int]
    description: str = ""
    rock_types: list[str] = field(default_factory=list)
    extra: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        rlo, rhi = self.resistivity_range
        if rlo <= 0 or rhi <= 0:
            raise ValueError("resistivity_range values must be > 0")
        if rlo > rhi:
            self.resistivity_range = (rhi, rlo)
        dlo, dhi = self.depth_range
        if dlo > dhi:
            self.depth_range = (dhi, dlo)
        nlo, nhi = self.n_layers_range
        if nlo > nhi:
            self.n_layers_range = (nhi, nlo)

    # ------------------------------------------------------------------
    # Computed properties
    # ------------------------------------------------------------------

    @property
    def log_rho_range(self) -> tuple[float, float]:
        """Log10-resistivity interval ``(log_rho_min, log_rho_max)``."""
        lo, hi = self.resistivity_range
        return (math.log10(lo), math.log10(hi))

    @property
    def rho_mid(self) -> float:
        """Geometric-mean resistivity (Ω·m)."""
        lo, hi = self.resistivity_range
        return math.sqrt(lo * hi)

    @property
    def depth_mid(self) -> float:
        """Arithmetic-mean depth (m)."""
        return (self.depth_range[0] + self.depth_range[1]) / 2.0

    @property
    def n_layers_mid(self) -> int:
        """Midpoint of the n_layers range, rounded."""
        return round((self.n_layers_range[0] + self.n_layers_range[1]) / 2.0)

    # ------------------------------------------------------------------
    # Compatibility helpers
    # ------------------------------------------------------------------

    def to_prior(self) -> dict[str, Any]:
        """Return a dict compatible with ``LayeredModel.from_geology(name)``.

        The returned dict has keys ``n_layers``, ``log_rho_range``,
        ``depth_max_range``, and ``description`` — exactly the format
        expected by ``GEOLOGY_PRIORS`` in ``forward/synthetic.py``.
        """
        return {
            "n_layers": self.n_layers_range,
            "log_rho_range": self.log_rho_range,
            "depth_max_range": self.depth_range,
            "description": self.description,
        }

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "resistivity_range": list(self.resistivity_range),
            "depth_range": list(self.depth_range),
            "n_layers_range": list(self.n_layers_range),
            "description": self.description,
            "rock_types": list(self.rock_types),
            **self.extra,
        }

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> Formation:
        known = {
            "name",
            "resistivity_range",
            "depth_range",
            "n_layers_range",
            "description",
            "rock_types",
        }
        extra = {k: v for k, v in d.items() if k not in known}
        return cls(
            name=d["name"],
            resistivity_range=tuple(d["resistivity_range"]),
            depth_range=tuple(d["depth_range"]),
            n_layers_range=tuple(d.get("n_layers_range", (3, 6))),
            description=d.get("description", ""),
            rock_types=list(d.get("rock_types", [])),
            extra=extra,
        )

    @classmethod
    def from_prior(cls, name: str, prior: dict[str, Any]) -> Formation:
        """Build a :class:`Formation` from a ``GEOLOGY_PRIORS`` entry."""
        log_lo, log_hi = prior["log_rho_range"]
        depth_lo, depth_hi = prior["depth_max_range"]
        n_lo, n_hi = prior["n_layers"]
        return cls(
            name=name,
            resistivity_range=(10.0**log_lo, 10.0**log_hi),
            depth_range=(depth_lo, depth_hi),
            n_layers_range=(n_lo, n_hi),
            description=prior.get("description", ""),
        )

    def __repr__(self) -> str:
        return (
            f"Formation({self.name!r}  "
            f"ρ=[{self.resistivity_range[0]:.2g}, {self.resistivity_range[1]:.2g}] Ω·m  "
            f"d=[{self.depth_range[0]:.0f}, {self.depth_range[1]:.0f}] m)"
        )


# ---------------------------------------------------------------------------
# Built-in formation definitions
# ---------------------------------------------------------------------------

_BUILTIN_FORMATIONS: list[dict[str, Any]] = [
    # ── Sedimentary / basin ─────────────────────────────────────────────
    dict(
        name="sedimentary",
        resistivity_range=(3.0, 3162.0),
        depth_range=(500, 3000),
        n_layers_range=(3, 7),
        description="Alternating clay/shale and sand/carbonate layers",
        rock_types=[
            "clay",
            "shale",
            "sedimentary rock",
            "dolomite/limestone",
        ],
    ),
    dict(
        name="evaporite",
        resistivity_range=(1e3, 1e6),
        depth_range=(200, 3000),
        n_layers_range=(2, 5),
        description="Salt/anhydrite evaporite basin",
        rock_types=["hard rock"],
    ),
    dict(
        name="coastal",
        resistivity_range=(0.5, 500.0),
        depth_range=(10, 500),
        n_layers_range=(3, 6),
        description="Coastal/delta mixing zone: fresh water over saline sediments",
        rock_types=["fresh water", "salt water", "clay", "gravel/sand"],
    ),
    # ── Crystalline / hard rock ─────────────────────────────────────────
    dict(
        name="crystalline",
        resistivity_range=(100.0, 31623.0),
        depth_range=(5000, 30000),
        n_layers_range=(3, 6),
        description="Resistive upper to conductive lower crust",
        rock_types=["igneous rock", "metamorphic rock", "hard rock"],
    ),
    dict(
        name="basement",
        resistivity_range=(1e3, 1e5),
        depth_range=(1000, 20000),
        n_layers_range=(2, 4),
        description="Deep resistive basement / craton",
        rock_types=["igneous rock", "metamorphic rock"],
    ),
    dict(
        name="volcanic",
        resistivity_range=(10.0, 1e4),
        depth_range=(200, 5000),
        n_layers_range=(3, 7),
        description="Volcanic pile: alternating lavas and pyroclastics",
        rock_types=["igneous rock", "hard rock"],
    ),
    # ── Geothermal / hydrothermal ────────────────────────────────────────
    dict(
        name="geothermal",
        resistivity_range=(2.0, 1e4),
        depth_range=(500, 5000),
        n_layers_range=(3, 5),
        description="Resistive cap over conductive geothermal reservoir",
        rock_types=["clay", "shale", "igneous rock"],
    ),
    dict(
        name="hydrothermal",
        resistivity_range=(1.0, 1e3),
        depth_range=(200, 3000),
        n_layers_range=(3, 6),
        description="Hydrothermal alteration zone with clay cap",
        rock_types=["clay", "saprolite", "shale"],
    ),
    # ── Marine / offshore ────────────────────────────────────────────────
    dict(
        name="marine",
        resistivity_range=(0.3, 1e3),
        depth_range=(100, 2000),
        n_layers_range=(3, 6),
        description="Seawater over possible HC reservoir (CSEM context)",
        rock_types=["sea water", "sedimentary rock", "gravel/sand"],
    ),
    # ── Permafrost / Arctic ──────────────────────────────────────────────
    dict(
        name="permafrost",
        resistivity_range=(10.0, 31623.0),
        depth_range=(50, 500),
        n_layers_range=(3, 5),
        description="Frozen resistive layer over conductive unfrozen sediments",
        rock_types=["permafrost", "tills"],
    ),
    # ── Ore-bearing / mining ─────────────────────────────────────────────
    dict(
        name="mineralized",
        resistivity_range=(0.01, 1e3),
        depth_range=(50, 1000),
        n_layers_range=(3, 6),
        description="Conductive ore-bearing zone embedded in resistive host rock",
        rock_types=[
            "massive sulphide",
            "ore minerals",
            "graphite",
            "igneous rock",
        ],
    ),
    dict(
        name="porphyry",
        resistivity_range=(10.0, 1e4),
        depth_range=(200, 2000),
        n_layers_range=(3, 6),
        description="Porphyry copper system: clay-alteration cap over "
        "mineralised core",
        rock_types=["clay", "igneous rock", "massive sulphide"],
    ),
    # ── Regolith / laterite ──────────────────────────────────────────────
    dict(
        name="laterite",
        resistivity_range=(10.0, 5000.0),
        depth_range=(10, 200),
        n_layers_range=(3, 6),
        description="Deeply weathered tropical profile: laterite / saprolite / "
        "bedrock",
        rock_types=["saprolite", "clay", "tills", "igneous rock"],
    ),
]


# ---------------------------------------------------------------------------
# GeologyCatalog
# ---------------------------------------------------------------------------


class GeologyCatalog:
    """Registry of geological formations for layered-Earth modelling.

    Parameters
    ----------
    formations : list of Formation, optional
        Initial formations.  When omitted the built-in catalog is loaded.

    Examples
    --------
    ::

        from pycsamt.metadata.geology import CATALOG

        # list all names
        print(CATALOG.names())

        # look up by name
        f = CATALOG.get("sedimentary")

        # find closest by resistivity
        matches = CATALOG.lookup_by_resistivity(50.0, n=3)

        # add a custom formation
        CATALOG.add(Formation(
            name="custom",
            resistivity_range=(5.0, 200.0),
            depth_range=(100, 500),
            n_layers_range=(3, 5),
        ))

        # export to a pandas DataFrame
        df = CATALOG.to_dataframe()
    """

    def __init__(
        self,
        formations: list[Formation] | None = None,
    ) -> None:
        self._store: dict[str, Formation] = {}
        if formations is not None:
            for f in formations:
                self._store[f.name] = f
        else:
            self._load_builtins()

    # ------------------------------------------------------------------
    # Mutation
    # ------------------------------------------------------------------

    def add(self, formation: Formation) -> None:
        """Register *formation*, overwriting any existing entry with the same name."""
        self._store[formation.name] = formation

    def remove(self, name: str) -> None:
        """Remove formation by name (raises KeyError when absent)."""
        del self._store[name]

    def reset(self) -> None:
        """Restore the built-in catalog, discarding custom additions."""
        self._store.clear()
        self._load_builtins()

    # ------------------------------------------------------------------
    # Queries
    # ------------------------------------------------------------------

    def get(self, name: str) -> Formation:
        """Return the :class:`Formation` for *name*.

        Raises
        ------
        KeyError
            When *name* is not in the catalog.
        """
        key = name.lower()
        if key not in self._store:
            available = sorted(self._store.keys())
            raise KeyError(
                f"Formation {name!r} not found.  Available: {available}"
            )
        return self._store[key]

    def names(self) -> list[str]:
        """Return a sorted list of all registered formation names."""
        return sorted(self._store.keys())

    def __contains__(self, name: str) -> bool:
        return name.lower() in self._store

    def __len__(self) -> int:
        return len(self._store)

    def __iter__(self):
        return iter(self._store.values())

    def lookup_by_resistivity(
        self,
        rho: float,
        n: int = 3,
    ) -> list[Formation]:
        """Return the *n* formations whose resistivity range best covers *rho*.

        Matching priority:

        1. Formations whose interval contains *rho* (sorted by narrowest range).
        2. Formations nearest to *rho* in log10 space.

        Parameters
        ----------
        rho : float
            Query resistivity in Ω·m.
        n : int
            Maximum number of results to return.
        """
        containing = [
            f
            for f in self._store.values()
            if f.resistivity_range[0] <= rho <= f.resistivity_range[1]
        ]
        if containing:
            return sorted(
                containing,
                key=lambda f: math.log10(
                    f.resistivity_range[1] / f.resistivity_range[0]
                ),
            )[:n]

        # Fall back: distance in log space
        log_rho = math.log10(max(rho, 1e-9))
        return sorted(
            self._store.values(),
            key=lambda f: min(
                abs(log_rho - math.log10(max(f.resistivity_range[0], 1e-9))),
                abs(log_rho - math.log10(max(f.resistivity_range[1], 1e-9))),
            ),
        )[:n]

    def lookup_by_depth(
        self,
        depth_m: float,
        n: int = 3,
    ) -> list[Formation]:
        """Return formations whose depth range contains *depth_m* (metres)."""
        containing = [
            f
            for f in self._store.values()
            if f.depth_range[0] <= depth_m <= f.depth_range[1]
        ]
        if containing:
            return sorted(
                containing,
                key=lambda f: f.depth_range[1] - f.depth_range[0],
            )[:n]
        return sorted(
            self._store.values(),
            key=lambda f: min(
                abs(depth_m - f.depth_range[0]),
                abs(depth_m - f.depth_range[1]),
            ),
        )[:n]

    def lookup_by_rock_type(self, rock_type: str) -> list[Formation]:
        """Return formations that list *rock_type* in their ``rock_types``."""
        q = rock_type.lower()
        return [
            f
            for f in self._store.values()
            if any(q in rt.lower() for rt in f.rock_types)
        ]

    def all_scenarios(self) -> dict[str, Formation]:
        """Return a copy of the internal store (name → Formation)."""
        return dict(self._store)

    # ------------------------------------------------------------------
    # Compatibility
    # ------------------------------------------------------------------

    def to_prior(self, name: str) -> dict[str, Any]:
        """Return a ``GEOLOGY_PRIORS``-compatible dict for *name*.

        This is the adapter used by
        :meth:`~pycsamt.forward.synthetic.LayeredModel.from_geology`.
        """
        return self.get(name).to_prior()

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------

    def to_dataframe(self, *, api: bool | None = None) -> Any:
        """Return a :class:`pandas.DataFrame` with one row per formation."""
        try:
            import pandas as pd  # noqa: PLC0415
        except ImportError as exc:
            raise ImportError(
                "pandas is required for to_dataframe()"
            ) from exc

        rows = []
        for f in self._store.values():
            rows.append(
                {
                    "name": f.name,
                    "rho_min": f.resistivity_range[0],
                    "rho_max": f.resistivity_range[1],
                    "rho_mid": round(f.rho_mid, 2),
                    "depth_min_m": f.depth_range[0],
                    "depth_max_m": f.depth_range[1],
                    "n_layers_min": f.n_layers_range[0],
                    "n_layers_max": f.n_layers_range[1],
                    "description": f.description,
                    "rock_types": ", ".join(f.rock_types),
                }
            )
        df = pd.DataFrame(rows).sort_values("name").reset_index(drop=True)

        return maybe_wrap_frame(
            df,
            api=api,
            name="geology_catalog",
            kind="metadata.geology",
            source=self.__class__.__name__,
            description=("Resistivity and lithology ranges by formation."),
        )

    def __repr__(self) -> str:
        return (
            f"GeologyCatalog({len(self._store)} formations: "
            + ", ".join(sorted(self._store)[:5])
            + ("…" if len(self._store) > 5 else "")
            + ")"
        )

    # ------------------------------------------------------------------
    # Private
    # ------------------------------------------------------------------

    def _load_builtins(self) -> None:
        for d in _BUILTIN_FORMATIONS:
            f = Formation.from_dict(d)
            self._store[f.name] = f


# ---------------------------------------------------------------------------
# Module-level singleton + compatibility function
# ---------------------------------------------------------------------------

CATALOG = GeologyCatalog()


def geology_prior(name: str) -> dict[str, Any]:
    """Return a ``GEOLOGY_PRIORS``-compatible dict for *name*.

    This is a drop-in replacement for direct access to the old
    ``GEOLOGY_PRIORS`` dict in ``forward/synthetic.py``.

    Parameters
    ----------
    name : str
        Formation name (case-insensitive).

    Returns
    -------
    dict
        Keys: ``n_layers``, ``log_rho_range``, ``depth_max_range``,
        ``description``.
    """
    return CATALOG.to_prior(name)
