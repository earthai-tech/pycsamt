# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Parallelised synthetic EM dataset generator.

:func:`generate_dataset` is the single entry point.  It creates a
large collection of (data, model) pairs suitable for training 1-D
neural network inverters.

Quick start
-----------
>>> import numpy as np
>>> from pycsamt.forward.batch import generate_dataset
>>> ds = generate_dataset(
...     solver="mt1d",
...     n_samples=5_000,
...     freqs=np.logspace(-3, 4, 30),
...     n_layers=(3, 7),
...     noise_level=0.05,
...     seed=42,
...     n_jobs=4,
...     output="mt1d_train.npz",
... )
>>> ds.X.shape          # (5000, 60) — log10(rho_a) + phase at 30 freqs
>>> ds.y.shape          # (5000, 9) — log10(rho) and thickness for up to 7 layers

Dataset layout
--------------
``X``
    Feature matrix, shape ``(n_samples, n_features)``.
    For MT: ``[log10(rho_a_0), …, log10(rho_a_{nf-1}), phi_0, …, phi_{nf-1}]``.
    For TEM: ``[log10(|dBz/dt|_0), …, log10(|dBz/dt|_{nt-1})]``.

``y``
    Target matrix, shape ``(n_samples, n_params)``.
    ``[log10(rho_0), …, log10(rho_{nl-1}), thick_0, …, thick_{nl-2}]``.
    Padded with ``NaN`` for samples with fewer than ``max_layers`` layers.

``meta``
    Structured array with per-sample metadata:
    ``n_layers``, ``depth_max``, ``noise_level``, ``seed``.
"""

from __future__ import annotations

import warnings
from concurrent.futures import (
    ProcessPoolExecutor,
    as_completed,
)
from dataclasses import dataclass

import numpy as np

__all__ = [
    "generate_dataset",
    "ForwardDataset",
    "generate_dataset_3d",
    "SurveyDataset3D",
]


# ─────────────────────────────────────────────────────────────────────────────
# Dataset container
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class ForwardDataset:
    """
    Numpy container for a batch of (features, targets) EM samples.

    Parameters
    ----------
    X : ndarray, shape (n_samples, n_features)
        Model responses (log-scaled data).
    y : ndarray, shape (n_samples, n_params)
        Model parameters (log10 resistivity + thickness).
    freqs : ndarray or None
        Frequencies used (MT/CSAMT).
    times : ndarray or None
        Times used (TEM).
    meta : structured ndarray
        Per-sample metadata (n_layers, noise_level, …).
    solver : str
        Solver used to generate this dataset.
    """

    X: np.ndarray
    y: np.ndarray
    freqs: np.ndarray | None = None
    times: np.ndarray | None = None
    meta: np.ndarray | None = None
    solver: str = "mt1d"

    def save(self, path: str) -> None:
        """Save to a compressed ``.npz`` file."""
        arrays = dict(X=self.X, y=self.y, solver=np.array(self.solver))
        if self.freqs is not None:
            arrays["freqs"] = self.freqs
        if self.times is not None:
            arrays["times"] = self.times
        if self.meta is not None:
            # structured arrays need special handling
            arrays["meta_n_layers"] = self.meta["n_layers"]
            arrays["meta_noise"] = self.meta["noise_level"]
        np.savez_compressed(path, **arrays)

    @classmethod
    def load(cls, path: str) -> ForwardDataset:
        """Load from a ``.npz`` file produced by :meth:`save`."""
        d = np.load(path, allow_pickle=False)
        freqs = d["freqs"] if "freqs" in d else None
        times = d["times"] if "times" in d else None
        meta = None
        if "meta_n_layers" in d:
            n = len(d["meta_n_layers"])
            meta = np.zeros(
                n, dtype=[("n_layers", "i4"), ("noise_level", "f4")]
            )
            meta["n_layers"] = d["meta_n_layers"]
            meta["noise_level"] = d["meta_noise"]
        return cls(
            X=d["X"],
            y=d["y"],
            freqs=freqs,
            times=times,
            meta=meta,
            solver=str(d["solver"]),
        )

    def split(
        self,
        val_frac: float = 0.1,
        test_frac: float = 0.1,
        seed: int | None = None,
    ) -> tuple[ForwardDataset, ForwardDataset, ForwardDataset]:
        """
        Split into train / validation / test sets.

        Returns
        -------
        (train, val, test) : tuple of ForwardDataset
        """
        rng = np.random.default_rng(seed)
        n = len(self.X)
        idx = rng.permutation(n)
        n_test = int(n * test_frac)
        n_val = int(n * val_frac)
        test_idx = idx[:n_test]
        val_idx = idx[n_test : n_test + n_val]
        train_idx = idx[n_test + n_val :]

        def _subset(indices):
            m = self.meta[indices] if self.meta is not None else None
            return ForwardDataset(
                X=self.X[indices],
                y=self.y[indices],
                freqs=self.freqs,
                times=self.times,
                meta=m,
                solver=self.solver,
            )

        return _subset(train_idx), _subset(val_idx), _subset(test_idx)

    def __len__(self) -> int:
        return len(self.X)

    def __repr__(self) -> str:
        return (
            f"ForwardDataset(n={len(self.X)}, "
            f"n_features={self.X.shape[1]}, "
            f"n_params={self.y.shape[1]}, "
            f"solver='{self.solver}')"
        )


# ─────────────────────────────────────────────────────────────────────────────
# Worker (must be module-level for pickle in ProcessPoolExecutor)
# ─────────────────────────────────────────────────────────────────────────────


def _worker(args):
    """
    Generate one (X, y) sample.

    Parameters are passed as a single tuple to satisfy ProcessPoolExecutor.
    """
    (
        solver_name,
        n_layers,
        rho_range,
        depth_max,
        freqs,
        times,
        loop_radius,
        noise_level,
        noise_type,
        geology,
        include_phase,
        seed,
    ) = args

    from pycsamt.forward.noise import add_noise
    from pycsamt.forward.synthetic import LayeredModel

    rng = np.random.default_rng(seed)
    n_lay = (
        int(rng.integers(n_layers[0], n_layers[1] + 1))
        if isinstance(n_layers, tuple)
        else n_layers
    )

    if geology is not None:
        model = LayeredModel.from_geology(geology, seed=rng)
    else:
        model = LayeredModel.random(
            n_layers=n_lay,
            rho_range=rho_range,
            depth_max=depth_max,
            seed=rng,
        )

    # Run forward solver
    if solver_name == "mt1d":
        from pycsamt.forward.em1d import MT1DForward

        resp = MT1DForward(freqs).run(model)
    elif solver_name == "tem1d":
        from pycsamt.forward.em1d import TEM1DForward

        resp = TEM1DForward(times, loop_radius=loop_radius).run(model)
    elif solver_name == "csamt1d":
        from pycsamt.forward.em1d import CSAMT1DForward

        resp = CSAMT1DForward(freqs).run(model)
    else:
        raise ValueError(f"Unknown solver: {solver_name!r}")

    # Apply noise
    if noise_level > 0.0:
        resp = add_noise(
            resp, noise_type, level=noise_level, seed=rng.integers(2**31)
        )

    # Feature vector
    x_vec = resp.to_array(log_rho=True, include_phase=include_phase)

    # Target vector (fixed size, NaN-padded)
    y_vec = model.to_vector(log_rho=True)

    return x_vec, y_vec, model.n_layers, noise_level


# ─────────────────────────────────────────────────────────────────────────────
# Main entry point
# ─────────────────────────────────────────────────────────────────────────────


def generate_dataset(
    solver: str = "mt1d",
    n_samples: int = 10_000,
    *,
    freqs: np.ndarray | None = None,
    times: np.ndarray | None = None,
    n_layers: int | tuple[int, int] = (3, 7),
    rho_range: tuple[float, float] = (1.0, 10_000.0),
    depth_max: float = 2000.0,
    loop_radius: float = 50.0,
    noise_level: float = 0.05,
    noise_type: str = "gaussian",
    geology: str | None = None,
    include_phase: bool = True,
    seed: int | None = None,
    n_jobs: int = 1,
    output: str | None = None,
    verbose: bool = True,
) -> ForwardDataset:
    """
    Generate a batch of synthetic (data, model) pairs for ML training.

    Parameters
    ----------
    solver : {'mt1d', 'tem1d', 'csamt1d'}
        Forward solver to use.
    n_samples : int
        Total number of samples.
    freqs : ndarray or None
        Frequencies [Hz] for MT/CSAMT.  If None, uses
        ``np.logspace(-3, 4, 30)`` for MT.
    times : ndarray or None
        Times [s] for TEM.  If None, uses ``np.logspace(-6, -2, 25)``.
    n_layers : int or (lo, hi)
        Fixed number of layers, or a range from which the count is
        drawn uniformly at random per sample.
    rho_range : (low, high)
        Resistivity bounds in Ω·m.
    depth_max : float
        Maximum depth of the model [m].
    loop_radius : float
        TEM transmitter loop radius [m].
    noise_level : float
        Relative noise standard deviation.  0 = noise-free.
    noise_type : str
        Noise model: ``'gaussian'``, ``'multiplicative'``, ``'field'``.
    geology : str or None
        If given, models are drawn from :func:`LayeredModel.from_geology`
        using this geological scenario name.
    include_phase : bool
        Include impedance phase in the MT feature vector.
    seed : int or None
        Base random seed.  Worker seeds are derived deterministically.
    n_jobs : int
        Number of parallel worker processes.  ``-1`` uses all CPU cores.
    output : str or None
        If given, save the dataset to a ``.npz`` file at this path.
    verbose : bool
        Print progress.

    Returns
    -------
    ForwardDataset

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.forward.batch import generate_dataset
    >>> ds = generate_dataset(n_samples=100, seed=0)
    >>> ds.X.shape[0]
    100
    """
    import os

    solver = solver.lower().strip()
    if solver not in ("mt1d", "tem1d", "csamt1d"):
        raise ValueError(
            f"Unknown solver {solver!r}. Use 'mt1d', 'tem1d', 'csamt1d'."
        )

    # Default grids
    if freqs is None and solver in ("mt1d", "csamt1d"):
        freqs = np.logspace(-3, 4, 30)
    if times is None and solver == "tem1d":
        times = np.logspace(-6, -2, 25)
        warnings.warn(
            "TEM1D forward uses scipy numerical integration — "
            "generation may be slow for large n_samples.  "
            "Consider n_jobs > 1 or using a smaller n_samples for Phase 1.",
            UserWarning,
            stacklevel=2,
        )

    n_jobs_eff = os.cpu_count() if n_jobs == -1 else max(1, n_jobs)

    # Build argument list (deterministic seeds from base seed)
    rng_base = np.random.default_rng(seed)
    child_seeds = rng_base.integers(0, 2**31, n_samples)

    args_list = [
        (
            solver,
            n_layers,
            rho_range,
            depth_max,
            freqs,
            times,
            loop_radius,
            noise_level,
            noise_type,
            geology,
            include_phase,
            int(child_seeds[i]),
        )
        for i in range(n_samples)
    ]

    # Run workers
    results = []
    if n_jobs_eff == 1:
        for i, args in enumerate(args_list):
            results.append(_worker(args))
            if verbose and (i + 1) % max(1, n_samples // 10) == 0:
                print(f"  {i + 1}/{n_samples} samples generated")
    else:
        with ProcessPoolExecutor(max_workers=n_jobs_eff) as pool:
            futures = {
                pool.submit(_worker, a): i for i, a in enumerate(args_list)
            }
            done = 0
            for fut in as_completed(futures):
                results.append(fut.result())
                done += 1
                if verbose and done % max(1, n_samples // 10) == 0:
                    print(f"  {done}/{n_samples} samples generated")

    # Assemble into arrays
    x_list = [r[0] for r in results]
    y_list = [r[1] for r in results]
    nl_list = [r[2] for r in results]
    nv_list = [r[3] for r in results]

    # Pad y to equal length (max param vector size)
    max_y = max(len(y) for y in y_list)
    y_pad = np.full((n_samples, max_y), np.nan)
    for i, y in enumerate(y_list):
        y_pad[i, : len(y)] = y

    X = np.vstack(x_list).astype(np.float32)
    y_arr = y_pad.astype(np.float32)

    meta_dt = np.dtype([("n_layers", "i4"), ("noise_level", "f4")])
    meta = np.zeros(n_samples, dtype=meta_dt)
    meta["n_layers"] = nl_list
    meta["noise_level"] = nv_list

    ds = ForwardDataset(
        X=X,
        y=y_arr,
        freqs=freqs,
        times=times,
        meta=meta,
        solver=solver,
    )

    if output is not None:
        ds.save(output)
        if verbose:
            print(f"Dataset saved to {output}")

    return ds


# ─────────────────────────────────────────────────────────────────────────────
# Pseudo-3D survey dataset container
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class SurveyDataset3D:
    """
    Numpy container for a pseudo-3D multi-station survey dataset.

    Designed for training :class:`~pycsamt.ai.inversion.inv3d.GCNInverter3D`.
    All surveys share the same fixed station grid so that a single
    pre-computed adjacency matrix covers the whole dataset.

    Parameters
    ----------
    X : ndarray, shape ``(n_surveys, n_stations, n_features)``
        Per-station MT/CSAMT feature vectors (log-scaled data).
    y : ndarray, shape ``(n_surveys, n_stations, n_params)``
        Per-station model parameters ``[log10(ρ), thickness]``.
    coords : ndarray, shape ``(n_stations, 2)``
        Station (easting, northing) coordinates in metres.
    freqs : ndarray or None
        Frequencies [Hz] used by the forward solver.
    meta : structured ndarray or None
        Per-survey metadata (``corr_length``, ``noise_level``).
    solver : str
        Solver used to generate this dataset.
    """

    X: np.ndarray
    y: np.ndarray
    coords: np.ndarray
    freqs: np.ndarray | None = None
    meta: np.ndarray | None = None
    solver: str = "mt1d"

    @property
    def n_surveys(self) -> int:
        """Number of synthetic survey realisations."""
        return int(self.X.shape[0])

    @property
    def n_stations(self) -> int:
        """Number of stations per survey."""
        return int(self.X.shape[1])

    @property
    def n_features(self) -> int:
        """Feature dimensionality per station."""
        return int(self.X.shape[2])

    @property
    def n_params(self) -> int:
        """Parameter dimensionality per station (``2 × n_layers − 1``)."""
        return int(self.y.shape[2])

    def save(self, path: str) -> None:
        """Save to a compressed ``.npz`` file."""
        arrays = dict(
            X=self.X,
            y=self.y,
            coords=self.coords,
            solver=np.array(self.solver),
        )
        if self.freqs is not None:
            arrays["freqs"] = self.freqs
        if self.meta is not None:
            arrays["meta_corr_length"] = self.meta["corr_length"]
            arrays["meta_noise"] = self.meta["noise_level"]
        np.savez_compressed(path, **arrays)

    @classmethod
    def load(cls, path: str) -> SurveyDataset3D:
        """Load from a ``.npz`` file produced by :meth:`save`."""
        d = np.load(path, allow_pickle=False)
        freqs = d["freqs"] if "freqs" in d else None
        meta = None
        if "meta_corr_length" in d:
            n = len(d["meta_corr_length"])
            meta = np.zeros(
                n, dtype=[("corr_length", "f4"), ("noise_level", "f4")]
            )
            meta["corr_length"] = d["meta_corr_length"]
            meta["noise_level"] = d["meta_noise"]
        return cls(
            X=d["X"],
            y=d["y"],
            coords=d["coords"],
            freqs=freqs,
            meta=meta,
            solver=str(d["solver"]),
        )

    def split(
        self,
        val_frac: float = 0.1,
        test_frac: float = 0.1,
        seed: int | None = None,
    ) -> tuple[SurveyDataset3D, SurveyDataset3D, SurveyDataset3D]:
        """
        Split into train / validation / test sets along the survey axis.

        Returns
        -------
        (train, val, test) : tuple of SurveyDataset3D
        """
        rng = np.random.default_rng(seed)
        n = self.n_surveys
        idx = rng.permutation(n)
        n_test = int(n * test_frac)
        n_val = int(n * val_frac)
        test_idx = idx[:n_test]
        val_idx = idx[n_test : n_test + n_val]
        train_idx = idx[n_test + n_val :]

        def _subset(indices):
            m = self.meta[indices] if self.meta is not None else None
            return SurveyDataset3D(
                X=self.X[indices],
                y=self.y[indices],
                coords=self.coords,
                freqs=self.freqs,
                meta=m,
                solver=self.solver,
            )

        return _subset(train_idx), _subset(val_idx), _subset(test_idx)

    def __len__(self) -> int:
        return self.n_surveys

    def __repr__(self) -> str:
        return (
            f"SurveyDataset3D(n_surveys={self.n_surveys}, "
            f"n_stations={self.n_stations}, "
            f"n_features={self.n_features}, "
            f"n_params={self.n_params}, "
            f"solver='{self.solver}')"
        )


# ─────────────────────────────────────────────────────────────────────────────
# Gaussian Random Field helper
# ─────────────────────────────────────────────────────────────────────────────


def _grf_cholesky(
    coords: np.ndarray, corr_length: float, nugget: float = 1e-6
) -> np.ndarray:
    """
    Cholesky factor of a squared-exponential covariance matrix.

    The covariance between stations *i* and *j* is

    .. math::

        C_{ij} = \\exp\\!\\left(-\\frac{\\|\\mathbf{x}_i - \\mathbf{x}_j\\|^2}
                                      {2\\,\\ell^2}\\right)

    where :math:`\\ell` is the correlation length.  A small nugget is
    added to the diagonal for numerical stability.

    Parameters
    ----------
    coords : ndarray, shape (n_stations, 2)
        Station (x, y) positions in metres.
    corr_length : float
        Spatial correlation length in metres.
    nugget : float
        Diagonal regularisation term (default 1e-6).

    Returns
    -------
    L : ndarray, shape (n_stations, n_stations)
        Lower-triangular Cholesky factor.
    """
    d2 = np.sum(
        (coords[:, None, :] - coords[None, :, :]) ** 2,
        axis=-1,
    )
    C = np.exp(-d2 / (2.0 * corr_length**2))
    C += nugget * np.eye(len(coords))
    return np.linalg.cholesky(C)


# ─────────────────────────────────────────────────────────────────────────────
# Per-survey worker (module-level for ProcessPoolExecutor pickling)
# ─────────────────────────────────────────────────────────────────────────────


def _worker_3d(args):
    """
    Generate one pseudo-3D survey (all stations, fixed station grid).

    Parameters are packed into a single tuple so that
    ``ProcessPoolExecutor.submit`` can pickle the call.
    """
    (
        solver_name,
        n_layers,
        n_stations,
        log_rho_mean,
        log_rho_std,
        thickness_range,
        L_chol,
        freqs,
        noise_level,
        noise_type,
        include_phase,
        seed,
    ) = args

    from pycsamt.forward.em1d import MT1DForward
    from pycsamt.forward.noise import add_noise
    from pycsamt.forward.synthetic import LayeredModel

    rng = np.random.default_rng(seed)

    # ── spatially correlated log-resistivity fields ───────────────────────────
    # For each layer k draw a GRF: z_k ~ N(0, I), then L @ z_k gives a sample
    # from N(0, C).  Scale by log_rho_std and shift by log_rho_mean.
    log_rho_fields = np.empty((n_layers, n_stations))
    for k in range(n_layers):
        z = rng.standard_normal(n_stations)
        log_rho_fields[k] = log_rho_mean + log_rho_std * (L_chol @ z)

    # ── thicknesses: log-uniform per layer, spatially constant per survey ─────
    log_t_lo = np.log10(float(thickness_range[0]))
    log_t_hi = np.log10(float(thickness_range[1]))
    thick = 10.0 ** rng.uniform(log_t_lo, log_t_hi, n_layers - 1)

    # ── run forward solver at every station ───────────────────────────────────
    fwd = MT1DForward(freqs)
    x_list = []
    y_list = []
    for s in range(n_stations):
        rho_s = 10.0 ** log_rho_fields[:, s]
        # Clip to physically valid range before constructing the model
        rho_s = np.maximum(rho_s, 1e-3)
        model = LayeredModel(resistivity=rho_s, thickness=thick.copy())
        resp = fwd.run(model)
        if noise_level > 0.0:
            resp = add_noise(
                resp,
                noise_type,
                level=noise_level,
                seed=int(rng.integers(2**31)),
            )
        x_list.append(
            resp.to_array(log_rho=True, include_phase=include_phase)
        )
        y_list.append(model.to_vector(log_rho=True))

    X_survey = np.array(x_list, dtype=np.float32)  # (n_sta, n_feat)
    y_survey = np.array(y_list, dtype=np.float32)  # (n_sta, n_params)
    return X_survey, y_survey


# ─────────────────────────────────────────────────────────────────────────────
# Main 3-D entry point
# ─────────────────────────────────────────────────────────────────────────────


def generate_dataset_3d(
    solver: str = "mt1d",
    n_surveys: int = 1000,
    *,
    n_stations: int = 25,
    n_layers: int = 4,
    freqs: np.ndarray | None = None,
    extent: float = 10_000.0,
    corr_length: float = 2_000.0,
    log_rho_mean: float = 2.0,
    log_rho_std: float = 0.5,
    thickness_range: tuple[float, float] = (100.0, 2_000.0),
    station_layout: str = "grid",
    noise_level: float = 0.03,
    noise_type: str = "gaussian",
    include_phase: bool = True,
    seed: int | None = None,
    n_jobs: int = 1,
    output: str | None = None,
    verbose: bool = True,
) -> SurveyDataset3D:
    """
    Generate a pseudo-3D synthetic dataset for :class:`GCNInverter3D` training.

    Each survey realisation consists of *n_stations* MT soundings whose
    per-layer log-resistivities are drawn from a spatially correlated
    Gaussian random field (GRF) with a squared-exponential covariance:

    .. math::

        C(\\mathbf{x}_i,\\mathbf{x}_j)
            = \\exp\\!\\left(-\\frac{\\|\\mathbf{x}_i-\\mathbf{x}_j\\|^2}
                                   {2\\,\\ell^2}\\right)

    where :math:`\\ell` is *corr_length*.  Layer thicknesses are drawn
    log-uniformly per survey but kept spatially constant (flat-layer
    assumption within one survey).  All surveys share the same fixed
    station grid so that a single adjacency matrix covers the whole
    training set.

    Parameters
    ----------
    solver : str
        Forward solver.  Currently only ``'mt1d'`` is supported.
    n_surveys : int
        Number of independent survey realisations.
    n_stations : int
        Number of MT stations per survey.
    n_layers : int
        Fixed number of earth layers (including halfspace).
    freqs : ndarray or None
        Frequencies [Hz].  Defaults to ``np.logspace(-3, 4, 30)``.
    extent : float
        Side length of the square survey area in metres.
    corr_length : float
        Spatial correlation length of the GRF in metres.  Set smaller
        than station spacing for uncorrelated models; larger for smooth
        lateral variation.
    log_rho_mean : float
        Mean of log₁₀(ρ) for all layers (default 2 → 100 Ω·m).
    log_rho_std : float
        Marginal standard deviation of the GRF in log₁₀(ρ) units.
    thickness_range : (lo, hi)
        Bounds [m] for log-uniform layer thickness draws.
    station_layout : ``'grid'`` or ``'random'``
        ``'grid'`` places stations on a regular grid of
        ``⌈√n_stations⌉²`` points (default).  ``'random'`` draws
        uniform positions within ``[0, extent]²``.
    noise_level : float
        Relative noise standard deviation applied to each response.
    noise_type : str
        Noise model: ``'gaussian'``, ``'multiplicative'``, ``'field'``.
    include_phase : bool
        Include impedance phase in the feature vector.
    seed : int or None
        Base random seed.  Worker seeds are derived deterministically.
    n_jobs : int
        Number of parallel worker processes.  ``-1`` uses all CPU cores.
    output : str or None
        If given, save the dataset to a compressed ``.npz`` file.
    verbose : bool
        Print progress.

    Returns
    -------
    SurveyDataset3D
        Container with arrays ``X`` (n_surveys, n_stations, n_features),
        ``y`` (n_surveys, n_stations, n_params), and
        ``coords`` (n_stations, 2).

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.forward.batch import generate_dataset_3d
    >>> from pycsamt.ai.nets.gcn import build_adjacency
    >>> ds = generate_dataset_3d(n_surveys=500, n_stations=16,
    ...                          n_layers=4, corr_length=2000., seed=0)
    >>> ds.X.shape
    (500, 16, 60)
    >>> A = build_adjacency(ds.coords, radius=3000.)
    >>> from pycsamt.ai.inversion.inv3d import GCNInverter3D
    >>> inv = GCNInverter3D(n_features=ds.n_features, n_layers=4)
    >>> inv.fit(ds.X, ds.y, adjacency=A, epochs=5, verbose=False)
    """
    import os

    solver = solver.lower().strip()
    if solver != "mt1d":
        raise ValueError(
            f"generate_dataset_3d currently supports solver='mt1d'; "
            f"got {solver!r}."
        )

    if freqs is None:
        freqs = np.logspace(-3, 4, 30)
    freqs = np.asarray(freqs, dtype=float)

    # ── station coordinates (fixed for the whole dataset) ────────────────────
    rng_base = np.random.default_rng(seed)
    if station_layout == "grid":
        side = int(np.ceil(np.sqrt(n_stations)))
        gx = np.linspace(0.0, float(extent), side)
        gy = np.linspace(0.0, float(extent), side)
        gX, gY = np.meshgrid(gx, gy)
        raw = np.column_stack([gX.ravel(), gY.ravel()])
        coords = raw[:n_stations].astype(np.float32)
    elif station_layout == "random":
        coords = rng_base.uniform(0.0, float(extent), (n_stations, 2)).astype(
            np.float32
        )
    else:
        raise ValueError(
            f"station_layout must be 'grid' or 'random'; got {station_layout!r}."
        )

    # ── pre-compute Cholesky factor (shared across all surveys) ───────────────
    L_chol = _grf_cholesky(coords.astype(float), float(corr_length))

    # ── per-survey worker arguments ───────────────────────────────────────────
    child_seeds = rng_base.integers(0, 2**31, n_surveys)
    args_list = [
        (
            solver,
            int(n_layers),
            int(n_stations),
            float(log_rho_mean),
            float(log_rho_std),
            thickness_range,
            L_chol,
            freqs,
            float(noise_level),
            noise_type,
            bool(include_phase),
            int(child_seeds[i]),
        )
        for i in range(n_surveys)
    ]

    n_jobs_eff = os.cpu_count() if n_jobs == -1 else max(1, n_jobs)

    results = []
    if n_jobs_eff == 1:
        for i, a in enumerate(args_list):
            results.append(_worker_3d(a))
            if verbose and (i + 1) % max(1, n_surveys // 10) == 0:
                print(f"  {i + 1}/{n_surveys} surveys generated")
    else:
        with ProcessPoolExecutor(max_workers=n_jobs_eff) as pool:
            futures = {
                pool.submit(_worker_3d, a): i for i, a in enumerate(args_list)
            }
            done = 0
            for fut in as_completed(futures):
                results.append(fut.result())
                done += 1
                if verbose and done % max(1, n_surveys // 10) == 0:
                    print(f"  {done}/{n_surveys} surveys generated")

    # ── assemble into 3-D arrays ──────────────────────────────────────────────
    X = np.stack([r[0] for r in results], axis=0).astype(np.float32)
    y = np.stack([r[1] for r in results], axis=0).astype(np.float32)

    meta_dt = np.dtype([("corr_length", "f4"), ("noise_level", "f4")])
    meta = np.zeros(n_surveys, dtype=meta_dt)
    meta["corr_length"] = float(corr_length)
    meta["noise_level"] = float(noise_level)

    ds = SurveyDataset3D(
        X=X,
        y=y,
        coords=coords,
        freqs=freqs,
        meta=meta,
        solver=solver,
    )

    if output is not None:
        ds.save(output)
        if verbose:
            print(f"3D dataset saved to {output}")

    return ds
