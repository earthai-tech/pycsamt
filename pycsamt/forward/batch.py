# -*- coding: utf-8 -*-
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
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Optional, Sequence, Tuple, Union

import numpy as np

__all__ = [
    "generate_dataset",
    "ForwardDataset",
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
    freqs: Optional[np.ndarray] = None
    times: Optional[np.ndarray] = None
    meta: Optional[np.ndarray] = None
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
    def load(cls, path: str) -> "ForwardDataset":
        """Load from a ``.npz`` file produced by :meth:`save`."""
        d = np.load(path, allow_pickle=False)
        freqs = d["freqs"] if "freqs" in d else None
        times = d["times"] if "times" in d else None
        meta = None
        if "meta_n_layers" in d:
            n = len(d["meta_n_layers"])
            meta = np.zeros(n, dtype=[("n_layers", "i4"), ("noise_level", "f4")])
            meta["n_layers"] = d["meta_n_layers"]
            meta["noise_level"] = d["meta_noise"]
        return cls(
            X=d["X"], y=d["y"],
            freqs=freqs, times=times,
            meta=meta,
            solver=str(d["solver"]),
        )

    def split(
        self,
        val_frac: float = 0.1,
        test_frac: float = 0.1,
        seed: Optional[int] = None,
    ) -> Tuple["ForwardDataset", "ForwardDataset", "ForwardDataset"]:
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
        val_idx = idx[n_test: n_test + n_val]
        train_idx = idx[n_test + n_val:]

        def _subset(indices):
            m = self.meta[indices] if self.meta is not None else None
            return ForwardDataset(
                X=self.X[indices], y=self.y[indices],
                freqs=self.freqs, times=self.times,
                meta=m, solver=self.solver,
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
        solver_name, n_layers, rho_range, depth_max,
        freqs, times, loop_radius, noise_level,
        noise_type, geology, include_phase, seed,
    ) = args

    from pycsamt.forward.synthetic import LayeredModel
    from pycsamt.forward.noise import add_noise

    rng = np.random.default_rng(seed)
    n_lay = int(rng.integers(n_layers[0], n_layers[1] + 1)) if isinstance(n_layers, tuple) else n_layers

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
        resp = add_noise(resp, noise_type, level=noise_level, seed=rng.integers(2**31))

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
    freqs: Optional[np.ndarray] = None,
    times: Optional[np.ndarray] = None,
    n_layers: Union[int, Tuple[int, int]] = (3, 7),
    rho_range: Tuple[float, float] = (1.0, 10_000.0),
    depth_max: float = 2000.0,
    loop_radius: float = 50.0,
    noise_level: float = 0.05,
    noise_type: str = "gaussian",
    geology: Optional[str] = None,
    include_phase: bool = True,
    seed: Optional[int] = None,
    n_jobs: int = 1,
    output: Optional[str] = None,
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
        raise ValueError(f"Unknown solver {solver!r}. Use 'mt1d', 'tem1d', 'csamt1d'.")

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
            solver, n_layers, rho_range, depth_max,
            freqs, times, loop_radius, noise_level,
            noise_type, geology, include_phase,
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
            futures = {pool.submit(_worker, a): i for i, a in enumerate(args_list)}
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
        X=X, y=y_arr,
        freqs=freqs, times=times,
        meta=meta, solver=solver,
    )

    if output is not None:
        ds.save(output)
        if verbose:
            print(f"Dataset saved to {output}")

    return ds
