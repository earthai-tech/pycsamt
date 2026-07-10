# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Training convergence visualisation.

:func:`plot_convergence` renders training and validation loss curves
with a shaded 1-σ band (when multiple runs are provided), a vertical
marker at the early-stopping epoch, and an LR schedule indicator.

Usage
-----
>>> from pycsamt.ai.plot.convergence import plot_convergence
>>> fig = plot_convergence(trainer.history)
>>> fig.savefig("convergence.png", dpi=300)
"""
from __future__ import annotations

import numpy as np

__all__ = ["plot_convergence", "plot_lr_schedule"]


def plot_convergence(
    history: dict[str, list] | list[dict[str, list]],
    *,
    ax=None,
    log_scale: bool = True,
    best_epoch: int | None = None,
    smoothing: float = 0.0,
    show_lr: bool = True,
    title: str = "Training convergence",
    style: bool = True,
) -> Figure:  # noqa: F821
    """
    Plot train / validation loss curves from a trainer history dict.

    Parameters
    ----------
    history : dict or list of dict
        A dict with keys ``'train_loss'`` and ``'val_loss'`` (and
        optionally ``'lr'``), as returned by
        :attr:`~pycsamt.ai.training.trainer.EMTrainer.history`.
        A list of such dicts (from multiple runs) activates the
        mean ± 1-σ band mode.
    ax : Axes or None
        Target axes.  If ``None``, a new figure/axes is created.
    log_scale : bool
        Log₁₀ y-axis for the loss.
    best_epoch : int or None
        If given, draw a vertical dashed line at this epoch.
    smoothing : float in [0, 1)
        Exponential moving average smoothing coefficient.
        0 = no smoothing.
    show_lr : bool
        Overlay learning rate on a twin y-axis (if ``'lr'`` in history).
    title : str
        Axes title.
    style : bool
        Apply :class:`~pycsamt.ai.plot._style.EMStyle`.

    Returns
    -------
    fig : Figure
    """
    import matplotlib.pyplot as plt

    from ._style import EM_COLORS, EM_FIGSIZE, EMStyle

    ctx = EMStyle() if style else _NullCtx()
    with ctx:
        if ax is None:
            fig, ax = plt.subplots(figsize=EM_FIGSIZE["double"])
        else:
            fig = ax.get_figure()

        # Normalise to list of histories
        histories = history if isinstance(history, list) else [history]
        n_runs = len(histories)

        epochs = np.arange(1, len(histories[0]["train_loss"]) + 1)

        train_mat = np.array([h["train_loss"] for h in histories])
        val_mat = np.array([h["val_loss"] for h in histories])

        if smoothing > 0:
            train_mat = np.apply_along_axis(_ema, 1, train_mat, smoothing)
            val_mat = np.apply_along_axis(_ema, 1, val_mat, smoothing)

        train_mean = train_mat.mean(axis=0)
        val_mean = val_mat.mean(axis=0)

        # ── Loss curves ──────────────────────────────────────────────────
        ax.plot(epochs, train_mean,
                color=EM_COLORS["primary"], lw=1.8, label="Train loss")
        ax.plot(epochs, val_mean,
                color=EM_COLORS["secondary"], lw=1.8, label="Val loss")

        if n_runs > 1:
            train_std = train_mat.std(axis=0)
            val_std = val_mat.std(axis=0)
            ax.fill_between(epochs,
                            train_mean - train_std, train_mean + train_std,
                            color=EM_COLORS["primary"], alpha=0.15)
            ax.fill_between(epochs,
                            val_mean - val_std, val_mean + val_std,
                            color=EM_COLORS["secondary"], alpha=0.15)

        # Best epoch marker
        if best_epoch is None:
            best_epoch = int(np.argmin(val_mean)) + 1
        ax.axvline(best_epoch, color="#555555", lw=1.0,
                   linestyle="--", label=f"Best epoch ({best_epoch})")

        ax.set_xlabel("Epoch", fontsize=11)
        ax.set_ylabel("MSE loss", fontsize=11)
        if log_scale:
            ax.set_yscale("log")
        ax.set_title(title, fontsize=13)
        ax.legend(fontsize=9, loc="upper right")

        # ── LR overlay ───────────────────────────────────────────────────
        if show_lr and "lr" in histories[0]:
            lr_vals = np.array(histories[0]["lr"])
            ax2 = ax.twinx()
            ax2.plot(epochs, lr_vals,
                     color="#999999", lw=1.0, linestyle=":",
                     label="LR")
            ax2.set_ylabel("Learning rate", fontsize=9, color="#777777")
            ax2.tick_params(axis="y", labelsize=8, colors="#777777")
            ax2.set_yscale("log")
            ax2.spines["right"].set_visible(True)

        fig.tight_layout()
    return fig


def plot_lr_schedule(
    lr_history: list,
    *,
    ax=None,
    title: str = "Learning rate schedule",
    style: bool = True,
):
    """
    Standalone learning-rate schedule plot.

    Parameters
    ----------
    lr_history : list of float
        Per-epoch LR values (from ``trainer.history['lr']``).
    ax : Axes or None
    title : str
    style : bool

    Returns
    -------
    ax : Axes
    """
    import matplotlib.pyplot as plt

    from ._style import EM_FIGSIZE, EMStyle

    ctx = EMStyle() if style else _NullCtx()
    with ctx:
        if ax is None:
            _, ax = plt.subplots(figsize=EM_FIGSIZE["single"])

        epochs = np.arange(1, len(lr_history) + 1)
        ax.semilogy(epochs, lr_history, color="#555555", lw=1.4)
        ax.set_xlabel("Epoch")
        ax.set_ylabel("Learning rate")
        ax.set_title(title, fontsize=13)
    return ax


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _ema(arr: np.ndarray, alpha: float) -> np.ndarray:
    """Exponential moving average, forward-pass."""
    out = np.empty_like(arr)
    out[0] = arr[0]
    for i in range(1, len(arr)):
        out[i] = alpha * out[i - 1] + (1.0 - alpha) * arr[i]
    return out


class _NullCtx:
    def __enter__(self): return self
    def __exit__(self, *_): pass
