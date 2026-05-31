# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Synthetic 1-D earth model generators.

The :class:`LayeredModel` dataclass is the central data object shared by
all forward solvers and the ML training pipeline.  It stores the
subsurface resistivity distribution as a vector of layer resistivities
and thicknesses, and provides several class-method constructors for
generating geologically realistic random models.

Typical usage
-------------
>>> from pycsamt.forward.synthetic import LayeredModel
>>> m = LayeredModel.random(n_layers=5, seed=42)
>>> m.resistivity          # array of Ω·m
>>> m.thickness            # array of metres
>>> x = m.to_vector()      # flat parameter vector for ML
>>> m2 = LayeredModel.from_vector(x, n_layers=5)

Geological priors
-----------------
``LayeredModel.from_geology(name)`` provides pre-set parameter ranges
for common geological scenarios:

``'sedimentary'``
    Alternating conductive (clay/shale) and resistive (sand/carbonate)
    layers, depth 0–3 km.

``'crystalline'``
    Resistive upper crust (1000–10 000 Ω·m) grading to conductive
    lower crust, depth 0–30 km.

``'geothermal'``
    Resistive cap rock over a conductive heat source/fluid zone.

``'marine'``
    Seawater (0.3 Ω·m) over resistive hydrocarbon reservoir, typical
    for marine CSEM scenarios.

``'permafrost'``
    Very resistive frozen layer (> 1000 Ω·m) over conductive unfrozen
    sediments.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional, Sequence, Tuple, Union

import numpy as np

__all__ = [
    "LayeredModel",
    "GEOLOGY_PRIORS",
]

# ─────────────────────────────────────────────────────────────────────────────
# Geological prior definitions
# ─────────────────────────────────────────────────────────────────────────────

#: Ready-to-use parameter bounds for ``LayeredModel.from_geology``.
#: Each entry: (n_layers_range, log10_rho_range, depth_max_range, description)
GEOLOGY_PRIORS: Dict[str, dict] = {
    "sedimentary": dict(
        n_layers=(3, 7),
        log_rho_range=(0.5, 3.5),   # 3 – 3162 Ω·m
        depth_max_range=(500, 3000),
        description="Alternating clay/shale and sand/carbonate layers",
    ),
    "crystalline": dict(
        n_layers=(3, 6),
        log_rho_range=(2.0, 4.5),   # 100 – 31 623 Ω·m
        depth_max_range=(5000, 30000),
        description="Resistive upper to conductive lower crust",
    ),
    "geothermal": dict(
        n_layers=(3, 5),
        log_rho_range=(0.3, 4.0),
        depth_max_range=(500, 5000),
        description="Resistive cap over conductive geothermal reservoir",
    ),
    "marine": dict(
        n_layers=(3, 6),
        log_rho_range=(-0.5, 3.0),  # seawater 0.3 Ω·m to resistive sand
        depth_max_range=(100, 2000),
        description="Seawater over possible HC reservoir (CSEM)",
    ),
    "permafrost": dict(
        n_layers=(3, 5),
        log_rho_range=(1.0, 4.5),
        depth_max_range=(50, 500),
        description="Frozen resistive layer over conductive unfrozen sediments",
    ),
}


# ─────────────────────────────────────────────────────────────────────────────
# LayeredModel
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class LayeredModel:
    """
    1-D layered earth model.

    The earth has ``n_layers`` layers.  The last layer is the halfspace
    (infinite thickness).  Thicknesses apply to layers 0 … n−2.

    Parameters
    ----------
    resistivity : array-like, shape (n_layers,)
        Layer resistivities in Ω·m, ordered top → bottom.
    thickness : array-like, shape (n_layers-1,)
        Layer thicknesses in metres.  The halfspace has no thickness.
    depth : ndarray or None
        Top-of-layer depths [m]; computed from thickness if not given.
    name : str
        Optional label (e.g. geological scenario name).

    Examples
    --------
    >>> from pycsamt.forward.synthetic import LayeredModel
    >>> m = LayeredModel(
    ...     resistivity=[100, 10, 500],
    ...     thickness=[300, 800],
    ... )
    >>> m.n_layers
    3
    >>> m.depth
    array([   0.,  300., 1100.])
    """

    resistivity: np.ndarray
    thickness: np.ndarray
    depth: np.ndarray = field(default=None, repr=False)
    name: str = ""

    def __post_init__(self):
        self.resistivity = np.asarray(self.resistivity, dtype=float)
        self.thickness = np.asarray(self.thickness, dtype=float)
        n = len(self.resistivity)
        if len(self.thickness) != n - 1:
            raise ValueError(
                f"len(thickness) must be n_layers - 1 = {n - 1}, "
                f"got {len(self.thickness)}"
            )
        if np.any(self.resistivity <= 0.0):
            raise ValueError("All resistivities must be strictly positive.")
        if np.any(self.thickness <= 0.0):
            raise ValueError("All layer thicknesses must be strictly positive.")
        if self.depth is None:
            self.depth = np.concatenate(
                [[0.0], np.cumsum(self.thickness)]
            )

    # ─── read-only helpers ────────────────────────────────────────────────

    @property
    def n_layers(self) -> int:
        """Number of layers including the halfspace."""
        return len(self.resistivity)

    @property
    def conductivity(self) -> np.ndarray:
        """Layer conductivities σ = 1/ρ  [S/m]."""
        return 1.0 / self.resistivity

    # ─── ML serialisation ─────────────────────────────────────────────────

    def to_vector(self, *, log_rho: bool = True) -> np.ndarray:
        """
        Flatten to a 1-D parameter vector suitable for ML targets.

        The vector layout is  ``[ρ₀, ρ₁, …, ρ_{n-1}, h₀, h₁, …, h_{n-2}]``
        where ρ values are log₁₀(Ω·m) when *log_rho* is ``True``.

        Parameters
        ----------
        log_rho : bool
            If ``True`` (default), resistivities are stored as
            log₁₀(ρ).  This normalises the dynamic range and improves
            neural network convergence.

        Returns
        -------
        v : ndarray, shape (2 * n_layers - 1,)
        """
        rho = np.log10(self.resistivity) if log_rho else self.resistivity.copy()
        return np.concatenate([rho, self.thickness])

    @classmethod
    def from_vector(
        cls,
        v: np.ndarray,
        n_layers: int,
        *,
        log_rho: bool = True,
        name: str = "",
    ) -> "LayeredModel":
        """
        Reconstruct a :class:`LayeredModel` from a flat parameter vector.

        Parameters
        ----------
        v : ndarray, shape (2 * n_layers - 1,)
            Parameter vector as produced by :meth:`to_vector`.
        n_layers : int
            Number of layers.
        log_rho : bool
            Whether resistivities in *v* are in log₁₀ space.
        """
        rho = v[:n_layers]
        thick = v[n_layers:]
        if log_rho:
            rho = 10.0 ** rho
        return cls(resistivity=rho, thickness=thick, name=name)

    # ─── constructors ─────────────────────────────────────────────────────

    @classmethod
    def random(
        cls,
        n_layers: int = 5,
        *,
        rho_range: Tuple[float, float] = (1.0, 10_000.0),
        depth_max: float = 2000.0,
        seed: Optional[Union[int, np.random.Generator]] = None,
        name: str = "random",
    ) -> "LayeredModel":
        """
        Generate a random layered model with log-uniform resistivities.

        Layer thicknesses are drawn from a Dirichlet-like distribution
        that sums to *depth_max* so that the model spans the full
        target depth range.

        Parameters
        ----------
        n_layers : int
            Number of layers (including halfspace).
        rho_range : (low, high)
            Resistivity bounds in Ω·m.
        depth_max : float
            Total depth spanned by the first (n-1) layers [m].
        seed : int, Generator, or None
            Random seed for reproducibility.
        """
        rng = _ensure_rng(seed)

        log_lo, log_hi = np.log10(rho_range[0]), np.log10(rho_range[1])
        log_rho = rng.uniform(log_lo, log_hi, n_layers)
        rho = 10.0 ** log_rho

        # Dirichlet-like split: n_layers-1 thicknesses summing to depth_max
        # n_layers-2 breaks → n_layers-1 fractions
        breaks = np.sort(rng.uniform(0.0, 1.0, max(1, n_layers - 2)))
        fracs = np.diff(np.concatenate([[0.0], breaks, [1.0]]))
        fracs = np.maximum(fracs, 1e-3)
        fracs /= fracs.sum()
        thick = fracs * depth_max

        return cls(resistivity=rho, thickness=thick, name=name)

    @classmethod
    def blocky(
        cls,
        n_layers: int = 4,
        *,
        rho_background: float = 100.0,
        rho_anomaly: float = 5.0,
        anomaly_layer: int = 1,
        depth_max: float = 1000.0,
        equal_thickness: bool = True,
        seed: Optional[Union[int, np.random.Generator]] = None,
        name: str = "blocky",
    ) -> "LayeredModel":
        """
        Build a model with a single conductive (or resistive) anomaly
        embedded in a resistive (or conductive) background.

        Parameters
        ----------
        rho_background : float
            Background layer resistivity [Ω·m].
        rho_anomaly : float
            Anomaly layer resistivity [Ω·m].
        anomaly_layer : int
            Zero-based index of the anomalous layer.
        """
        rng = _ensure_rng(seed)
        rho = np.full(n_layers, rho_background)
        rho[anomaly_layer % n_layers] = rho_anomaly

        if equal_thickness:
            thick = np.full(n_layers - 1, depth_max / (n_layers - 1))
        else:
            thick = cls.random(n_layers, depth_max=depth_max, seed=rng).thickness

        return cls(resistivity=rho, thickness=thick, name=name)

    @classmethod
    def smooth(
        cls,
        n_layers: int = 10,
        *,
        rho_surface: float = 100.0,
        rho_deep: float = 10.0,
        depth_max: float = 5000.0,
        perturbation: float = 0.2,
        seed: Optional[Union[int, np.random.Generator]] = None,
        name: str = "smooth",
    ) -> "LayeredModel":
        """
        Build a model with a smooth gradient from surface to depth,
        plus small random perturbations.

        Parameters
        ----------
        rho_surface : float
            Resistivity at the surface [Ω·m].
        rho_deep : float
            Resistivity at *depth_max* [Ω·m].
        perturbation : float
            Fractional standard deviation of log₁₀(ρ) noise added to
            the gradient.  0 = pure gradient.
        """
        rng = _ensure_rng(seed)
        log_rho = np.linspace(
            np.log10(rho_surface), np.log10(rho_deep), n_layers
        )
        log_rho += rng.normal(0.0, perturbation, n_layers)
        rho = 10.0 ** log_rho
        rho = np.maximum(rho, 0.01)

        thick = np.full(n_layers - 1, depth_max / (n_layers - 1))

        return cls(resistivity=rho, thickness=thick, name=name)

    @classmethod
    def from_geology(
        cls,
        name: str,
        *,
        seed: Optional[Union[int, np.random.Generator]] = None,
    ) -> "LayeredModel":
        """
        Generate a random model using a predefined geological prior.

        Parameters
        ----------
        name : str
            One of: ``'sedimentary'``, ``'crystalline'``,
            ``'geothermal'``, ``'marine'``, ``'permafrost'``.

        Raises
        ------
        KeyError
            If *name* is not in :data:`GEOLOGY_PRIORS`.
        """
        if name not in GEOLOGY_PRIORS:
            raise KeyError(
                f"Unknown geology {name!r}. "
                f"Available: {list(GEOLOGY_PRIORS)}"
            )
        rng = _ensure_rng(seed)
        p = GEOLOGY_PRIORS[name]

        n_lo, n_hi = p["n_layers"]
        n = int(rng.integers(n_lo, n_hi + 1))
        d_lo, d_hi = p["depth_max_range"]
        depth_max = rng.uniform(d_lo, d_hi)
        lr_lo, lr_hi = p["log_rho_range"]
        log_rho = rng.uniform(lr_lo, lr_hi, n)
        rho = 10.0 ** log_rho

        breaks = np.sort(rng.uniform(0.0, 1.0, max(1, n - 2)))
        fracs = np.diff(np.concatenate([[0.0], breaks, [1.0]]))
        fracs = np.maximum(fracs, 1e-3)
        fracs /= fracs.sum()
        thick = fracs * depth_max

        return cls(resistivity=rho, thickness=thick, name=name)

    # ─── visualisation ────────────────────────────────────────────────────

    def plot(
        self,
        ax=None,
        *,
        log_scale: bool = True,
        depth_max: Optional[float] = None,
        label: Optional[str] = None,
        **kwargs,
    ):
        """
        Plot the 1-D resistivity–depth profile.

        Parameters
        ----------
        ax : Axes or None
            Target axes.  Created if not given.
        log_scale : bool
            Use log₁₀ scale on the resistivity axis.
        depth_max : float or None
            Maximum depth shown.  Defaults to 1.2× the deepest interface.
        label : str or None
            Line label for legend.

        Returns
        -------
        ax : Axes
        """
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots(figsize=(3.5, 5))

        rho = self.resistivity
        # Build a staircase profile
        depths_top = self.depth
        depths_bot = np.concatenate(
            [depths_top[1:], [depths_top[-1] + self.thickness[-1]
                               if len(self.thickness) else depths_top[-1] * 2]]
        )
        rho_steps = np.repeat(rho, 2)
        depth_steps = np.zeros(2 * len(rho))
        depth_steps[0::2] = depths_top
        depth_steps[1::2] = depths_bot

        kw = dict(color="steelblue", linewidth=1.5)
        kw.update(kwargs)
        if label is not None:
            kw["label"] = label

        ax.plot(rho_steps, depth_steps, **kw)
        ax.invert_yaxis()
        if log_scale:
            ax.set_xscale("log")

        dmax = depth_max or depths_bot[-1] * 1.05
        ax.set_ylim(dmax, 0.0)
        ax.set_xlabel(r"Resistivity (Ω·m)")
        ax.set_ylabel("Depth (m)")
        ax.set_title(self.name or "1-D earth model")
        ax.grid(True, which="both", alpha=0.3)
        return ax

    # ─── repr ─────────────────────────────────────────────────────────────

    def __repr__(self) -> str:
        rho_str = np.array2string(self.resistivity, precision=1, max_line_width=60)
        thick_str = np.array2string(self.thickness, precision=1, max_line_width=60)
        return (
            f"LayeredModel(n_layers={self.n_layers}, "
            f"rho={rho_str}, thick={thick_str})"
        )


# ─────────────────────────────────────────────────────────────────────────────
# Helper
# ─────────────────────────────────────────────────────────────────────────────

def _ensure_rng(seed) -> np.random.Generator:
    """Return a ``numpy.random.Generator`` from any seed-like input."""
    if isinstance(seed, np.random.Generator):
        return seed
    return np.random.default_rng(seed)
