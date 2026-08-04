# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
TensorFlow implementations of the PINN optimisation
kernels.

All ``import tensorflow`` calls are deferred to
function bodies — importing this module never fails
on systems without TensorFlow.

Complex arithmetic notes
------------------------
The Wait (1954) recursion requires complex128
arithmetic.  TensorFlow supports ``tf.complex128``
and tracks gradients through it via Wirtinger
derivatives.  Since the loss is always real-valued,
the gradient w.r.t. real parameters (log_rho,
log_thick) is real-valued and computed correctly.

Gradient clamping is done with
``tf.clip_by_norm`` after ``tape.gradient``.
Parameter clamping uses ``var.assign``.
"""

from __future__ import annotations

import math

import numpy as np

__all__ = [
    "_mt1d_tf",
    "_mt1d_tf_batch",
    "_fit_station_tf",
    "_fit_2d_joint_tf",
    "_fit_3d_joint_tf",
]

_MU0 = 4e-7 * math.pi
_PI = math.pi
_LOG10E = math.log10(math.e)  # 1 / ln(10)


# ── single-station Wait recursion ──────────────────


def _mt1d_tf(freqs_t, log_rho, log_thick):
    r"""
    Wait (1954) 1-D MT recursion in TensorFlow.

    Parameters
    ----------
    freqs_t : Tensor (F,) float64
        Frequencies in Hz.
    log_rho : Variable (L,) float64
        Base-10 log-resistivity per layer.
    log_thick : Variable (L-1,) float64
        Base-10 log-thickness per layer in m.

    Returns
    -------
    log10_rho_a : Tensor (F,) float64
    phase : Tensor (F,) float64
    """
    import tensorflow as tf

    rho = tf.pow(10.0, tf.cast(log_rho, tf.float64))
    thick = tf.pow(10.0, tf.cast(log_thick, tf.float64))
    omega = 2.0 * _PI * freqs_t  # (F,)

    j_mu_w = tf.complex(
        tf.zeros_like(omega),
        tf.cast(_MU0 * omega, tf.float64),
    )  # (F,)

    rho_c = tf.cast(rho, tf.complex128)  # (L,)
    k_n = tf.sqrt(j_mu_w / rho_c[-1])
    Z = j_mu_w / k_n  # (F,)

    L = log_rho.shape[0]
    j_unit = tf.constant(1j, dtype=tf.complex128)

    for idx in range(L - 2, -1, -1):
        k_j = tf.sqrt(j_mu_w / rho_c[idx])
        Z0_j = j_mu_w / k_j
        h = tf.cast(thick[idx], tf.complex128)
        kh = k_j * h
        tanh_v = tf.math.tanh(j_unit * kh)
        denom = Z0_j + Z * tanh_v
        Z = Z0_j * (Z + Z0_j * tanh_v) / denom

    mu0_omega = tf.cast(_MU0 * omega, tf.float64)
    abs_Z2 = tf.cast(tf.abs(Z) ** 2, tf.float64)
    rho_a = abs_Z2 / mu0_omega  # (F,)
    phase = tf.cast(tf.math.angle(Z), tf.float64) * (180.0 / _PI)

    log10_rho_a = tf.math.log(tf.maximum(rho_a, 1e-30)) * _LOG10E
    return log10_rho_a, phase


# ── batched Wait recursion ─────────────────────────


def _mt1d_tf_batch(freqs_t, log_rho, log_thick):
    r"""
    Batched Wait (1954) MT-1D forward in TensorFlow.

    Parameters
    ----------
    freqs_t : Tensor (F,) float64
    log_rho : Variable (S, L) float64
    log_thick : Variable (S, L-1) float64

    Returns
    -------
    log10_rho_a : Tensor (S, F) float64
    phase : Tensor (S, F) float64
    """
    import tensorflow as tf

    _S, L = (
        log_rho.shape[0],
        log_rho.shape[1],
    )
    rho = tf.pow(10.0, tf.cast(log_rho, tf.float64))  # (S, L)
    thick = tf.pow(10.0, tf.cast(log_thick, tf.float64))  # (S, L-1)
    omega = 2.0 * _PI * freqs_t  # (F,)

    # j_mu_w shape (1, F) — broadcasts over stations
    j_mu_w = tf.complex(
        tf.zeros_like(omega),
        tf.cast(_MU0 * omega, tf.float64),
    )[tf.newaxis, :]  # (1, F)

    rho_c = tf.cast(rho, tf.complex128)  # (S, L)

    k_n = tf.sqrt(j_mu_w / rho_c[:, -1:])  # (S, F)
    Z = j_mu_w / k_n  # (S, F)

    j_unit = tf.constant(1j, dtype=tf.complex128)

    for idx in range(L - 2, -1, -1):
        k_j = tf.sqrt(j_mu_w / rho_c[:, idx : idx + 1])
        Z0_j = j_mu_w / k_j
        h = tf.cast(thick[:, idx : idx + 1], tf.complex128)
        kh = k_j * h
        tanh_v = tf.math.tanh(j_unit * kh)
        denom = Z0_j + Z * tanh_v
        Z = Z0_j * (Z + Z0_j * tanh_v) / denom

    mu0_omega = tf.cast(_MU0 * omega[tf.newaxis, :], tf.float64)
    abs_Z2 = tf.cast(tf.abs(Z) ** 2, tf.float64)  # (S, F)
    rho_a = abs_Z2 / mu0_omega  # (S, F)
    phase = tf.cast(tf.math.angle(Z), tf.float64) * (180.0 / _PI)

    log10_rho_a = tf.math.log(tf.maximum(rho_a, 1e-30)) * _LOG10E
    return log10_rho_a, phase


# ── helpers ────────────────────────────────────────


def _clip_grads_tf(grads, max_norm: float = 5.0):
    """Clip a list of TF gradients by global norm."""
    import tensorflow as tf

    return [
        tf.clip_by_norm(g, max_norm) if g is not None else g for g in grads
    ]


# ── 1-D per-station optimisation ───────────────────


def _fit_station_tf(
    obs,
    n_layers: int,
    depth_max: float,
    lam: float,
    lr: float,
    epochs: int,
    device: str,
    log_every: int,
    init_log_rho: np.ndarray | None = None,
    init_log_thick: np.ndarray | None = None,
) -> dict:
    """
    TensorFlow 1-D station optimisation.

    Parameters match :func:`_fit_station_torch`.

    Returns
    -------
    dict : log_rho, log_thick, history
    """
    import tensorflow as tf

    freq = obs.freq
    rho_obs = obs.rho_obs
    ph_obs = obs.phase_obs

    freqs_t = tf.constant(freq, dtype=tf.float64)
    lr_obs_t = tf.constant(
        np.log10(np.maximum(rho_obs, 1e-12)),
        dtype=tf.float64,
    )
    ph_obs_t = tf.constant(ph_obs, dtype=tf.float64)
    valid = tf.math.is_finite(lr_obs_t) & tf.math.is_finite(ph_obs_t)

    if init_log_rho is None:
        init_lr = float(np.mean(np.log10(np.maximum(rho_obs, 1e-3))))
        init_lr_v = np.full(n_layers, init_lr, dtype=np.float64)
    else:
        init_lr_v = np.asarray(init_log_rho, dtype=np.float64)

    if init_log_thick is None:
        dz = depth_max / n_layers
        init_lt_v = np.full(
            n_layers - 1,
            np.log10(max(dz, 1.0)),
            dtype=np.float64,
        )
    else:
        init_lt_v = np.asarray(init_log_thick, dtype=np.float64)

    log_rho_var = tf.Variable(init_lr_v, dtype=tf.float64)
    log_thick_var = tf.Variable(init_lt_v, dtype=tf.float64)
    opt = tf.keras.optimizers.Adam(learning_rate=lr)
    history: list[float] = []

    for ep in range(1, epochs + 1):
        with tf.GradientTape() as tape:
            lr_pred, ph_pred = _mt1d_tf(freqs_t, log_rho_var, log_thick_var)
            lr_safe = tf.where(
                valid,
                lr_obs_t,
                tf.stop_gradient(lr_pred),
            )
            ph_safe = tf.where(
                valid,
                ph_obs_t,
                tf.stop_gradient(ph_pred),
            )
            n_v = tf.maximum(
                tf.reduce_sum(tf.cast(valid, tf.float64)),
                1.0,
            )
            data_loss = (
                tf.reduce_sum(
                    (lr_pred - lr_safe) ** 2
                    + ((ph_pred - ph_safe) / 90.0) ** 2
                )
                / n_v
            )
            smooth = tf.reduce_mean((log_rho_var[1:] - log_rho_var[:-1]) ** 2)
            loss = data_loss + lam * smooth

        grads = tape.gradient(loss, [log_rho_var, log_thick_var])
        clipped = _clip_grads_tf(grads)
        opt.apply_gradients(zip(clipped, [log_rho_var, log_thick_var]))
        log_thick_var.assign(tf.clip_by_value(log_thick_var, 0.0, 5.0))
        val = float(loss.numpy())
        history.append(val)

        if log_every > 0 and (ep == 1 or ep % log_every == 0):
            print(f"  {obs.name}  ep {ep:>5d}/{epochs}  loss={val:.5f}")

    return {
        "log_rho": log_rho_var.numpy(),
        "log_thick": log_thick_var.numpy(),
        "history": history,
    }


# ── 2-D joint optimisation ─────────────────────────


def _fit_2d_joint_tf(
    log_rho_obs,
    ph_obs,
    freqs_t,
    n_layers,
    depth_max,
    lam_z,
    lam_x,
    lr,
    epochs,
    device,
    log_every,
    init_log_rho=None,
    init_log_thick=None,
    verbose=None,
):
    r"""
    TensorFlow 2-D joint optimisation via Adam.

    Parameters match :func:`_fit_2d_joint_torch`.

    Returns
    -------
    dict : log_rho (S,L), log_thick (S,L-1), history
    """
    import tensorflow as tf

    from pycsamt.api.view.progress import get_progress_bar

    S, F = log_rho_obs.shape
    L = n_layers

    lr_obs_t = tf.constant(log_rho_obs, dtype=tf.float64)
    ph_obs_t = tf.constant(ph_obs, dtype=tf.float64)
    valid = tf.math.is_finite(lr_obs_t) & tf.math.is_finite(ph_obs_t)

    if init_log_rho is None:
        mean_lr = float(np.nanmean(log_rho_obs))
        if not np.isfinite(mean_lr):
            mean_lr = 2.0
        init_lr_v = np.full((S, L), mean_lr)
    else:
        init_lr_v = np.asarray(init_log_rho, dtype=float)

    if init_log_thick is None:
        dz = max(depth_max / L, 1.0)
        init_lt_v = np.full((S, L - 1), np.log10(dz))
    else:
        init_lt_v = np.asarray(init_log_thick, dtype=float)

    freqs_c = tf.constant(freqs_t, dtype=tf.float64)
    log_rho_var = tf.Variable(init_lr_v, dtype=tf.float64)
    log_thick_var = tf.Variable(init_lt_v, dtype=tf.float64)
    opt = tf.keras.optimizers.Adam(learning_rate=lr)
    history: list[float] = []

    _verbose = verbose if verbose is not None else (2 if log_every > 0 else 0)
    with get_progress_bar(
        total=epochs,
        desc="Stage 2 (2-D joint)",
        unit="epoch",
        verbose=_verbose,
        log_every=log_every or None,
    ) as bar:
        for ep in range(1, epochs + 1):
            with tf.GradientTape() as tape:
                lr_pred, ph_pred = _mt1d_tf_batch(
                    freqs_c, log_rho_var, log_thick_var
                )
                lr_safe = tf.where(
                    valid,
                    lr_obs_t,
                    tf.stop_gradient(lr_pred),
                )
                ph_safe = tf.where(
                    valid,
                    ph_obs_t,
                    tf.stop_gradient(ph_pred),
                )
                n_v = tf.maximum(
                    tf.reduce_sum(tf.cast(valid, tf.float64)),
                    1.0,
                )
                data_loss = (
                    tf.reduce_sum(
                        (lr_pred - lr_safe) ** 2
                        + ((ph_pred - ph_safe) / 90.0) ** 2
                    )
                    / n_v
                )

                vert = tf.reduce_mean(
                    (log_rho_var[:, 1:] - log_rho_var[:, :-1]) ** 2
                )

                if S > 1:
                    lat = tf.reduce_mean(
                        (log_rho_var[1:, :] - log_rho_var[:-1, :]) ** 2
                    )
                else:
                    lat = tf.constant(0.0, dtype=tf.float64)

                loss = data_loss + lam_z * vert + lam_x * lat

            grads = tape.gradient(loss, [log_rho_var, log_thick_var])
            clipped = _clip_grads_tf(grads)
            opt.apply_gradients(zip(clipped, [log_rho_var, log_thick_var]))
            log_thick_var.assign(tf.clip_by_value(log_thick_var, 0.0, 5.0))
            val = float(loss.numpy())
            history.append(val)
            bar.update(1, metrics={"loss": val})

    return {
        "log_rho": log_rho_var.numpy(),
        "log_thick": log_thick_var.numpy(),
        "history": history,
    }


# ── 3-D joint optimisation ─────────────────────────


def _fit_3d_joint_tf(
    log_rho_obs,
    ph_obs,
    freqs_t,
    adjacency,
    n_layers,
    depth_max,
    lam_z,
    lam_g,
    lr,
    epochs,
    device,
    log_every,
    init_log_rho=None,
    init_log_thick=None,
    verbose=None,
):
    r"""
    TensorFlow 3-D joint optimisation via Adam.

    Parameters match :func:`_fit_3d_joint_torch`.

    Returns
    -------
    dict : log_rho (S,L), log_thick (S,L-1), history
    """
    import tensorflow as tf

    from pycsamt.api.view.progress import get_progress_bar

    S, F = log_rho_obs.shape
    L = n_layers

    lr_obs_t = tf.constant(log_rho_obs, dtype=tf.float64)
    ph_obs_t = tf.constant(ph_obs, dtype=tf.float64)
    valid = tf.math.is_finite(lr_obs_t) & tf.math.is_finite(ph_obs_t)
    A_t = tf.constant(adjacency, dtype=tf.float64)
    d_vec = tf.reduce_sum(A_t, axis=1)  # (S,)

    if init_log_rho is None:
        mean_lr = float(np.nanmean(log_rho_obs))
        if not np.isfinite(mean_lr):
            mean_lr = 2.0
        init_lr_v = np.full((S, L), mean_lr)
    else:
        init_lr_v = np.asarray(init_log_rho, dtype=float)

    if init_log_thick is None:
        dz = max(depth_max / L, 1.0)
        init_lt_v = np.full((S, L - 1), np.log10(dz))
    else:
        init_lt_v = np.asarray(init_log_thick, dtype=float)

    freqs_c = tf.constant(freqs_t, dtype=tf.float64)
    log_rho_var = tf.Variable(init_lr_v, dtype=tf.float64)
    log_thick_var = tf.Variable(init_lt_v, dtype=tf.float64)
    opt = tf.keras.optimizers.Adam(learning_rate=lr)
    history: list[float] = []

    _verbose = verbose if verbose is not None else (2 if log_every > 0 else 0)
    with get_progress_bar(
        total=epochs,
        desc="Stage 2 (3-D joint)",
        unit="epoch",
        verbose=_verbose,
        log_every=log_every or None,
    ) as bar:
        for ep in range(1, epochs + 1):
            with tf.GradientTape() as tape:
                lr_pred, ph_pred = _mt1d_tf_batch(
                    freqs_c, log_rho_var, log_thick_var
                )
                lr_safe = tf.where(
                    valid,
                    lr_obs_t,
                    tf.stop_gradient(lr_pred),
                )
                ph_safe = tf.where(
                    valid,
                    ph_obs_t,
                    tf.stop_gradient(ph_pred),
                )
                n_v = tf.maximum(
                    tf.reduce_sum(tf.cast(valid, tf.float64)),
                    1.0,
                )
                data_loss = (
                    tf.reduce_sum(
                        (lr_pred - lr_safe) ** 2
                        + ((ph_pred - ph_safe) / 90.0) ** 2
                    )
                    / n_v
                )

                vert = tf.reduce_mean(
                    (log_rho_var[:, 1:] - log_rho_var[:, :-1]) ** 2
                )

                # Graph Laplacian: L = D - A
                # d_vec: (S,); log_rho_var: (S,L)
                d_col = tf.expand_dims(d_vec, 1)  # (S, 1)
                Lm = d_col * log_rho_var - tf.linalg.matmul(A_t, log_rho_var)
                spatial = tf.reduce_sum(Lm * log_rho_var) / tf.cast(
                    tf.size(log_rho_var),
                    tf.float64,
                )

                loss = data_loss + lam_z * vert + lam_g * spatial

            grads = tape.gradient(loss, [log_rho_var, log_thick_var])
            clipped = _clip_grads_tf(grads)
            opt.apply_gradients(zip(clipped, [log_rho_var, log_thick_var]))
            log_thick_var.assign(tf.clip_by_value(log_thick_var, 0.0, 5.0))
            val = float(loss.numpy())
            history.append(val)
            bar.update(1, metrics={"loss": val})

    return {
        "log_rho": log_rho_var.numpy(),
        "log_thick": log_thick_var.numpy(),
        "history": history,
    }
