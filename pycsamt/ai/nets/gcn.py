# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Graph Convolutional Network for spatial EM inversion.

:class:`GCNNet` implements Kipf & Welling (2017) spectral graph
convolutions adapted for geophysical data:

    :math:`H^{(l+1)} = \\sigma\\!\\left(\\tilde{D}^{-1/2}
    \\tilde{A}\\tilde{D}^{-1/2} H^{(l)} W^{(l)}\\right)`

where :math:`\\tilde{A} = A + I` (self-loops),
:math:`\\tilde{D}_{ii} = \\sum_j \\tilde{A}_{ij}`.

No external graph library (PyTorch Geometric, DGL) is required; the
normalised adjacency matrix is pre-computed once from station
coordinates and passed as a dense tensor to every forward call.
"""
from __future__ import annotations

from typing import Sequence, Tuple

import numpy as np

__all__ = ["GCNNet", "build_adjacency"]


# ─────────────────────────────────────────────────────────────────────────────
# Adjacency utilities (framework-free)
# ─────────────────────────────────────────────────────────────────────────────

def build_adjacency(
    coords: np.ndarray,
    radius: float,
    *,
    self_loops: bool = True,
    normalise: bool = True,
) -> np.ndarray:
    """
    Build a symmetric adjacency matrix from 2-D station coordinates.

    Parameters
    ----------
    coords : ndarray, shape (n_stations, 2)
        Station (x, y) positions in any consistent unit (metres, degrees).
    radius : float
        Maximum inter-station distance for an edge to exist.  Uses the
        same unit as *coords*.
    self_loops : bool
        If ``True`` (default), add :math:`\\tilde{A} = A + I`.
    normalise : bool
        If ``True`` (default), apply symmetric normalisation
        :math:`\\tilde{D}^{-1/2}\\tilde{A}\\tilde{D}^{-1/2}`.

    Returns
    -------
    A : ndarray, shape (n_stations, n_stations), float32
        Adjacency matrix, optionally normalised.
    """
    coords = np.asarray(coords, dtype=np.float64)
    n = len(coords)
    diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]   # (n, n, 2)
    dist = np.sqrt((diff ** 2).sum(axis=-1))                       # (n, n)

    A = (dist <= radius).astype(np.float32)
    if not self_loops:
        np.fill_diagonal(A, 0.0)
    else:
        np.fill_diagonal(A, 1.0)

    if normalise:
        deg = A.sum(axis=1)                          # (n,)
        deg_inv_sqrt = np.where(deg > 0, deg ** -0.5, 0.0)
        D = np.diag(deg_inv_sqrt)
        A = D @ A @ D

    return A.astype(np.float32)


# ─────────────────────────────────────────────────────────────────────────────
# PyTorch GCN
# ─────────────────────────────────────────────────────────────────────────────

class GCNNet:
    """
    Factory that builds a PyTorch or TensorFlow GCN.

    Call :meth:`build` after import to obtain the framework module.

    Parameters
    ----------
    n_features : int
        Dimensionality of per-node input features (e.g. 2 × n_freqs).
    n_out : int
        Per-node output size (2n-1 for an n-layer model:
        n resistivities + n-1 thicknesses).
    hidden : sequence of int
        Hidden-layer widths for each GCN message-passing step.
    dropout : float
        Dropout probability applied after each hidden GCN layer.
    """

    def __init__(
        self,
        n_features: int,
        n_out: int,
        hidden: Sequence[int] = (256, 128, 64),
        dropout: float = 0.1,
    ) -> None:
        self.n_features = int(n_features)
        self.n_out      = int(n_out)
        self.hidden     = tuple(int(h) for h in hidden)
        self.dropout    = float(dropout)

    # ── PyTorch ──────────────────────────────────────────────────────────────

    def build(self) -> "torch.nn.Module":
        """Return a ``torch.nn.Module`` for this architecture."""
        import torch
        import torch.nn as nn
        import torch.nn.functional as F

        hidden  = self.hidden
        n_in    = self.n_features
        n_out   = self.n_out
        dropout = self.dropout

        class _GCNLayer(nn.Module):
            def __init__(self, in_f: int, out_f: int) -> None:
                super().__init__()
                self.W = nn.Linear(in_f, out_f, bias=True)
                self.bn = nn.BatchNorm1d(out_f)

            def forward(
                self, H: "torch.Tensor", A: "torch.Tensor"
            ) -> "torch.Tensor":
                # H : (n_nodes, in_f)   A : (n_nodes, n_nodes)
                agg = A @ H                        # (n_nodes, in_f)
                out = self.W(agg)                  # (n_nodes, out_f)
                out = self.bn(out)
                return F.relu(out)

        class _GCNModule(nn.Module):
            def __init__(self) -> None:
                super().__init__()
                dims = [n_in] + list(hidden)
                self.layers = nn.ModuleList(
                    [_GCNLayer(dims[i], dims[i + 1])
                     for i in range(len(dims) - 1)]
                )
                self.drop = nn.Dropout(dropout)
                self.head = nn.Linear(dims[-1], n_out)

            def forward(
                self,
                H: "torch.Tensor",
                A: "torch.Tensor",
            ) -> "torch.Tensor":
                """
                Parameters
                ----------
                H : Tensor (n_nodes, n_features)  — node feature matrix
                A : Tensor (n_nodes, n_nodes)      — normalised adjacency

                Returns
                -------
                out : Tensor (n_nodes, n_out)      — per-node predictions
                """
                for layer in self.layers:
                    H = self.drop(layer(H, A))
                return self.head(H)

        return _GCNModule()

    # ── TensorFlow / Keras ────────────────────────────────────────────────────

    def build_tf(self) -> "tf.keras.Model":
        """Return a ``tf.keras.Model`` for this architecture."""
        import tensorflow as tf

        n_in    = self.n_features
        n_out   = self.n_out
        hidden  = self.hidden
        dropout = self.dropout

        class GCNLayer(tf.keras.layers.Layer):
            def __init__(self, units: int, **kw) -> None:
                super().__init__(**kw)
                self.dense = tf.keras.layers.Dense(units, activation=None)
                self.bn    = tf.keras.layers.BatchNormalization()

            def call(
                self,
                inputs: Tuple["tf.Tensor", "tf.Tensor"],
                training: bool = False,
            ) -> "tf.Tensor":
                H, A = inputs          # H: (n, f)  A: (n, n)
                agg = tf.linalg.matmul(A, H)
                out = self.dense(agg)
                out = self.bn(out, training=training)
                return tf.nn.relu(out)

        H_in = tf.keras.Input(shape=(n_in,),    name="node_features")
        A_in = tf.keras.Input(shape=(None,),    name="adj_row")

        # We build the Keras model as a pure function (adjacency is data input)
        # H_in : (n_nodes, n_in)  A_in : (n_nodes, n_nodes) — passed as batches
        H = H_in
        for units in hidden:
            layer = GCNLayer(units)
            H = layer([H, A_in])
            H = tf.keras.layers.Dropout(dropout)(H)
        out = tf.keras.layers.Dense(n_out, name="output")(H)

        return tf.keras.Model(inputs=[H_in, A_in], outputs=out, name="GCNNet")
