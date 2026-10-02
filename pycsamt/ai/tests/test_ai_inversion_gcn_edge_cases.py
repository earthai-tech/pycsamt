# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Edge-case coverage for pycsamt.ai.nets.gcn.build_adjacency and
pycsamt.ai.inversion.inv3d.GCNInverter3D._prepare_adjacency.

test_inversion.py's TestGCNInverter3DInterface/WithBackend cover the
default (self-loops + normalised) adjacency path and the coords-based
fit path thoroughly. This file covers the branches those miss: the
``self_loops=False`` / ``normalise=False`` toggles on
:func:`build_adjacency`, and the validation errors in
``_prepare_adjacency`` for malformed ``adjacency``/``coords`` inputs.

Note: :class:`~pycsamt.ai.nets.gcn.GCNNet.build_tf` and the
TensorFlow-backend fit/predict branches of ``GCNInverter3D`` are not
exercised anywhere in this environment -- TensorFlow's native runtime
fails to import here (broken DLL, confirmed via a subprocess probe),
matching the project-wide convention of never importing tensorflow
in-process to test for availability.
"""

from __future__ import annotations

import numpy as np
import pytest


def test_build_adjacency_no_self_loops():
    from pycsamt.ai.nets.gcn import build_adjacency

    coords = np.array([[0.0, 0.0], [1.0, 0.0], [10.0, 10.0]])
    A = build_adjacency(coords, radius=2.0, self_loops=False, normalise=False)
    assert np.all(np.diag(A) == 0.0)
    # stations 0 and 1 are within radius -> edge; station 2 is isolated
    assert A[0, 1] == 1.0
    assert A[2, 0] == 0.0


def test_build_adjacency_unnormalised_is_binary():
    from pycsamt.ai.nets.gcn import build_adjacency

    coords = np.array([[0.0, 0.0], [1.0, 0.0], [10.0, 10.0]])
    A = build_adjacency(coords, radius=2.0, self_loops=True, normalise=False)
    assert set(np.unique(A)).issubset({0.0, 1.0})
    assert np.all(np.diag(A) == 1.0)


def test_prepare_adjacency_rejects_non_square_adjacency():
    from pycsamt.ai.inversion.inv3d import GCNInverter3D

    bad_A = np.zeros((3, 4), dtype=np.float32)
    with pytest.raises(ValueError, match="must be square"):
        GCNInverter3D._prepare_adjacency(bad_A, None, 5_000.0, 3)


def test_prepare_adjacency_rejects_bad_coords_shape():
    from pycsamt.ai.inversion.inv3d import GCNInverter3D

    bad_coords = np.zeros((3, 3), dtype=np.float64)  # must be (n, 2)
    with pytest.raises(ValueError, match=r"coords must be shape"):
        GCNInverter3D._prepare_adjacency(None, bad_coords, 5_000.0, 3)


def test_prepare_adjacency_accepts_precomputed_adjacency():
    from pycsamt.ai.inversion.inv3d import GCNInverter3D

    A = np.eye(4, dtype=np.float32)
    out = GCNInverter3D._prepare_adjacency(A, None, 5_000.0, 4)
    np.testing.assert_array_equal(out, A)
