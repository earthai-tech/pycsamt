# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.ai.nets — CNN1DNet, ResNet1DNet, FCN1DNet,
UNet2DNet, DRCNNNet, GCNNet, and build_adjacency.

Structure
---------
* Framework-free interface tests (always run): construction, repr,
  attribute types, build_adjacency (pure-NumPy).
* PyTorch-dependent tests (skip when torch absent): build(), forward
  pass shapes and dtypes, gradient flow, parameter counts.
"""
from __future__ import annotations

import importlib.util
import unittest

import numpy as np

_HAS_TORCH = importlib.util.find_spec("torch") is not None


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _torch_tensor(arr):
    import torch
    return torch.tensor(arr, dtype=torch.float32)


def _param_count(module) -> int:
    return sum(p.numel() for p in module.parameters())


# ─────────────────────────────────────────────────────────────────────────────
# build_adjacency (framework-free)
# ─────────────────────────────────────────────────────────────────────────────

class TestBuildAdjacency(unittest.TestCase):

    def _grid_coords(self, n: int = 4, spacing: float = 1000.0) -> np.ndarray:
        """Return a regular n×n grid of station (x, y) coordinates."""
        xs, ys = np.meshgrid(
            np.arange(n) * spacing,
            np.arange(n) * spacing,
        )
        return np.column_stack([xs.ravel(), ys.ravel()]).astype(np.float32)

    def test_output_shape(self):
        from pycsamt.ai.nets.gcn import build_adjacency
        coords = self._grid_coords(4)
        A = build_adjacency(coords, radius=1500.0)
        self.assertEqual(A.shape, (16, 16))

    def test_dtype_float32(self):
        from pycsamt.ai.nets.gcn import build_adjacency
        coords = self._grid_coords(3)
        A = build_adjacency(coords, radius=2000.0)
        self.assertEqual(A.dtype, np.float32)

    def test_symmetric(self):
        from pycsamt.ai.nets.gcn import build_adjacency
        coords = self._grid_coords(5)
        A = build_adjacency(coords, radius=1200.0)
        np.testing.assert_allclose(A, A.T, atol=1e-6)

    def test_self_loops_default(self):
        """Diagonal should be non-zero when self_loops=True (default)."""
        from pycsamt.ai.nets.gcn import build_adjacency
        coords = self._grid_coords(4)
        A = build_adjacency(coords, radius=1500.0, self_loops=True)
        self.assertTrue(np.all(np.diag(A) > 0))

    def test_no_self_loops(self):
        from pycsamt.ai.nets.gcn import build_adjacency
        coords = self._grid_coords(4)
        A = build_adjacency(coords, radius=1500.0,
                            self_loops=False, normalise=False)
        np.testing.assert_array_equal(np.diag(A), np.zeros(16))

    def test_normalised_row_sums_leq_one(self):
        """D^{-1/2} A D^{-1/2} normalisation — each element is in [0, 1]."""
        from pycsamt.ai.nets.gcn import build_adjacency
        coords = self._grid_coords(4)
        A = build_adjacency(coords, radius=1500.0, normalise=True)
        self.assertTrue(np.all(A >= 0), "Negative entry after normalisation")
        self.assertTrue(np.all(A <= 1.0 + 1e-6), "Entry > 1 after normalisation")

    def test_unnormalised_binary(self):
        from pycsamt.ai.nets.gcn import build_adjacency
        coords = self._grid_coords(4)
        A = build_adjacency(coords, radius=1500.0,
                            self_loops=False, normalise=False)
        unique_vals = np.unique(A)
        np.testing.assert_array_equal(unique_vals, [0.0, 1.0])

    def test_zero_radius_gives_identity(self):
        """Radius 0 with self_loops: only self-edges remain."""
        from pycsamt.ai.nets.gcn import build_adjacency
        coords = self._grid_coords(3)
        A = build_adjacency(coords, radius=0.0,
                            self_loops=True, normalise=False)
        np.testing.assert_array_equal(A, np.eye(9, dtype=np.float32))

    def test_large_radius_fully_connected(self):
        from pycsamt.ai.nets.gcn import build_adjacency
        coords = self._grid_coords(4)
        A = build_adjacency(coords, radius=1e9,
                            self_loops=False, normalise=False)
        # Every off-diagonal entry should be 1 (fully connected)
        self.assertTrue(np.all(A[~np.eye(16, dtype=bool)] == 1.0))

    def test_single_station(self):
        from pycsamt.ai.nets.gcn import build_adjacency
        coords = np.array([[0.0, 0.0]])
        A = build_adjacency(coords, radius=100.0)
        self.assertEqual(A.shape, (1, 1))

    def test_values_non_negative(self):
        from pycsamt.ai.nets.gcn import build_adjacency
        coords = self._grid_coords(5)
        A = build_adjacency(coords, radius=1500.0)
        self.assertTrue(np.all(A >= 0))


# ─────────────────────────────────────────────────────────────────────────────
# Interface tests (no backend required)
# ─────────────────────────────────────────────────────────────────────────────

class TestCNN1DNetInterface(unittest.TestCase):

    def test_importable(self):
        from pycsamt.ai.nets import CNN1DNet
        self.assertTrue(callable(CNN1DNet))

    def test_construction_defaults(self):
        from pycsamt.ai.nets import CNN1DNet
        net = CNN1DNet(n_features=40, n_out=9)
        self.assertEqual(net.n_features, 40)
        self.assertEqual(net.n_out, 9)
        self.assertEqual(net.channels, (32, 64, 128))
        self.assertEqual(net.kernel_size, 5)

    def test_construction_custom(self):
        from pycsamt.ai.nets import CNN1DNet
        net = CNN1DNet(n_features=20, n_out=5, channels=(16, 32), kernel_size=3, dropout=0.1)
        self.assertEqual(net.channels, (16, 32))
        self.assertEqual(net.kernel_size, 3)
        self.assertAlmostEqual(net.dropout, 0.1)

    def test_repr_contains_class(self):
        from pycsamt.ai.nets import CNN1DNet
        r = repr(CNN1DNet(20, 5))
        self.assertIn("CNN1DNet", r)
        self.assertIn("20", r)

    def test_top_level_export(self):
        from pycsamt.ai import CNN1DNet
        self.assertTrue(callable(CNN1DNet))


class TestResNet1DNetInterface(unittest.TestCase):

    def test_importable(self):
        from pycsamt.ai.nets import ResNet1DNet
        self.assertTrue(callable(ResNet1DNet))

    def test_construction_defaults(self):
        from pycsamt.ai.nets import ResNet1DNet
        net = ResNet1DNet(n_features=40, n_out=9)
        self.assertEqual(net.n_features, 40)
        self.assertEqual(net.n_out, 9)

    def test_repr_not_empty(self):
        from pycsamt.ai.nets import ResNet1DNet
        self.assertIn("ResNet1DNet", repr(ResNet1DNet(20, 5)))

    def test_top_level_export(self):
        from pycsamt.ai import ResNet1DNet
        self.assertTrue(callable(ResNet1DNet))


class TestFCN1DNetInterface(unittest.TestCase):

    def test_importable(self):
        from pycsamt.ai.nets import FCN1DNet
        self.assertTrue(callable(FCN1DNet))

    def test_construction(self):
        from pycsamt.ai.nets import FCN1DNet
        net = FCN1DNet(n_features=30, n_out=7)
        self.assertEqual(net.n_features, 30)
        self.assertEqual(net.n_out, 7)

    def test_repr(self):
        from pycsamt.ai.nets import FCN1DNet
        self.assertIn("FCN1DNet", repr(FCN1DNet(30, 7)))

    def test_top_level_export(self):
        from pycsamt.ai import FCN1DNet
        self.assertTrue(callable(FCN1DNet))


class TestUNet2DNetInterface(unittest.TestCase):

    def test_importable(self):
        from pycsamt.ai.nets import UNet2DNet
        self.assertTrue(callable(UNet2DNet))

    def test_construction_defaults(self):
        from pycsamt.ai.nets import UNet2DNet
        net = UNet2DNet(n_in=4, n_out=1)
        self.assertEqual(net.n_in, 4)
        self.assertEqual(net.n_out, 1)

    def test_repr(self):
        from pycsamt.ai.nets import UNet2DNet
        self.assertIn("UNet2DNet", repr(UNet2DNet(4, 1)))

    def test_top_level_export(self):
        from pycsamt.ai import UNet2DNet
        self.assertTrue(callable(UNet2DNet))


class TestDRCNNNetInterface(unittest.TestCase):

    def test_importable(self):
        from pycsamt.ai.nets import DRCNNNet
        self.assertTrue(callable(DRCNNNet))

    def test_construction(self):
        from pycsamt.ai.nets import DRCNNNet
        net = DRCNNNet(n_features_list=(120, 48), n_out=9)
        self.assertEqual(net.n_out, 9)
        self.assertIn(120, net.n_features_list)

    def test_repr(self):
        from pycsamt.ai.nets import DRCNNNet
        self.assertIn("DRCNNNet", repr(DRCNNNet((120, 48), 9)))

    def test_top_level_export(self):
        from pycsamt.ai import DRCNNNet
        self.assertTrue(callable(DRCNNNet))


class TestGCNNetInterface(unittest.TestCase):

    def test_importable(self):
        from pycsamt.ai.nets import GCNNet
        self.assertTrue(callable(GCNNet))

    def test_construction_defaults(self):
        from pycsamt.ai.nets import GCNNet
        net = GCNNet(n_features=20, n_out=5)
        self.assertEqual(net.n_features, 20)
        self.assertEqual(net.n_out, 5)
        self.assertEqual(net.hidden, (256, 128, 64))

    def test_custom_hidden(self):
        from pycsamt.ai.nets import GCNNet
        net = GCNNet(n_features=10, n_out=3, hidden=(64, 32))
        self.assertEqual(net.hidden, (64, 32))

    def test_repr(self):
        from pycsamt.ai.nets import GCNNet
        self.assertIn("GCNNet", repr(GCNNet(10, 3)))

    def test_top_level_export(self):
        from pycsamt.ai import GCNNet
        self.assertTrue(callable(GCNNet))

    def test_build_adjacency_exported(self):
        from pycsamt.ai.nets import build_adjacency
        self.assertTrue(callable(build_adjacency))

    def test_build_adjacency_top_level(self):
        from pycsamt.ai import build_adjacency
        self.assertTrue(callable(build_adjacency))


class TestNetsModuleExports(unittest.TestCase):
    """All public names in __all__ should be importable."""

    def test_all_exported(self):
        import pycsamt.ai.nets as nets_mod
        for name in nets_mod.__all__:
            self.assertTrue(hasattr(nets_mod, name), f"Missing: {name}")

    def test_top_level_ai_exports(self):
        import pycsamt.ai as ai
        for name in ("CNN1DNet", "ResNet1DNet", "FCN1DNet",
                     "UNet2DNet", "DRCNNNet", "GCNNet", "build_adjacency"):
            self.assertTrue(hasattr(ai, name), f"Missing top-level: {name}")


# ─────────────────────────────────────────────────────────────────────────────
# PyTorch-dependent forward-pass tests
# ─────────────────────────────────────────────────────────────────────────────

@unittest.skipUnless(_HAS_TORCH, "PyTorch not installed")
class TestCNN1DNetBuild(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        from pycsamt.ai.nets import CNN1DNet
        cls.net = CNN1DNet(n_features=40, n_out=9).build()

    def test_is_module(self):
        import torch.nn as nn
        self.assertIsInstance(self.net, nn.Module)

    def test_forward_shape(self):
        import torch
        x = torch.randn(8, 40)
        y = self.net(x)
        self.assertEqual(y.shape, (8, 9))

    def test_output_dtype_float32(self):
        import torch
        x = torch.randn(4, 40)
        y = self.net(x)
        self.assertEqual(y.dtype, torch.float32)

    def test_batch_size_one(self):
        import torch
        x = torch.randn(1, 40)
        y = self.net(x)
        self.assertEqual(y.shape, (1, 9))

    def test_gradients_flow(self):
        import torch
        x = torch.randn(4, 40, requires_grad=False)
        y = self.net(x).sum()
        y.backward()
        for p in self.net.parameters():
            if p.requires_grad:
                self.assertIsNotNone(p.grad)
                break

    def test_param_count_positive(self):
        self.assertGreater(_param_count(self.net), 0)

    def test_custom_channels(self):
        import torch

        from pycsamt.ai.nets import CNN1DNet
        net = CNN1DNet(n_features=24, n_out=5, channels=(16, 32)).build()
        y = net(torch.randn(3, 24))
        self.assertEqual(y.shape, (3, 5))

    def test_single_feature(self):
        import torch

        from pycsamt.ai.nets import CNN1DNet
        net = CNN1DNet(n_features=1, n_out=1, channels=(4,)).build()
        y = net(torch.randn(2, 1))
        self.assertEqual(y.shape, (2, 1))


@unittest.skipUnless(_HAS_TORCH, "PyTorch not installed")
class TestResNet1DNetBuild(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        from pycsamt.ai.nets import ResNet1DNet
        cls.net = ResNet1DNet(n_features=40, n_out=9).build()

    def test_is_module(self):
        import torch.nn as nn
        self.assertIsInstance(self.net, nn.Module)

    def test_forward_shape(self):
        import torch
        y = self.net(torch.randn(8, 40))
        self.assertEqual(y.shape, (8, 9))

    def test_output_dtype(self):
        import torch
        self.assertEqual(self.net(torch.randn(2, 40)).dtype, torch.float32)

    def test_batch_size_one(self):
        import torch
        self.assertEqual(self.net(torch.randn(1, 40)).shape, (1, 9))

    def test_gradients_flow(self):
        import torch
        self.net(torch.randn(4, 40)).sum().backward()
        grads = [p.grad for p in self.net.parameters() if p.requires_grad]
        self.assertTrue(any(g is not None for g in grads))

    def test_param_count_positive(self):
        self.assertGreater(_param_count(self.net), 0)

    def test_different_input_size(self):
        import torch

        from pycsamt.ai.nets import ResNet1DNet
        net = ResNet1DNet(n_features=20, n_out=5).build()
        self.assertEqual(net(torch.randn(3, 20)).shape, (3, 5))


@unittest.skipUnless(_HAS_TORCH, "PyTorch not installed")
class TestFCN1DNetBuild(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        from pycsamt.ai.nets import FCN1DNet
        cls.net = FCN1DNet(n_features=30, n_out=7).build()

    def test_forward_shape(self):
        import torch
        self.assertEqual(self.net(torch.randn(6, 30)).shape, (6, 7))

    def test_output_dtype(self):
        import torch
        self.assertEqual(self.net(torch.randn(2, 30)).dtype, torch.float32)

    def test_batch_size_one(self):
        import torch
        self.assertEqual(self.net(torch.randn(1, 30)).shape, (1, 7))

    def test_gradients_flow(self):
        import torch
        self.net(torch.randn(4, 30)).sum().backward()
        grads = [p.grad for p in self.net.parameters() if p.requires_grad]
        self.assertTrue(any(g is not None for g in grads))

    def test_variable_input_size(self):
        """FCN is length-agnostic — should work for different n_features."""
        import torch

        from pycsamt.ai.nets import FCN1DNet
        net = FCN1DNet(n_features=50, n_out=9).build()
        self.assertEqual(net(torch.randn(5, 50)).shape, (5, 9))


@unittest.skipUnless(_HAS_TORCH, "PyTorch not installed")
class TestUNet2DNetBuild(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        from pycsamt.ai.nets import UNet2DNet
        # n_in=4 channels (e.g. 4 EM components), n_out=1 (resistivity)
        cls.net = UNet2DNet(n_in=4, n_out=1).build()

    def test_is_module(self):
        import torch.nn as nn
        self.assertIsInstance(self.net, nn.Module)

    def test_forward_shape_spatial_preserved(self):
        """U-Net output H×W must equal input H×W."""
        import torch
        x = torch.randn(2, 4, 16, 16)
        y = self.net(x)
        self.assertEqual(y.shape[:2], (2, 1))
        self.assertEqual(y.shape[2:], x.shape[2:])

    def test_output_dtype(self):
        import torch
        # batch of 1: run in eval mode, BatchNorm cannot normalise a
        # single sample per channel in training mode
        x = torch.randn(1, 4, 16, 16)
        self.net.eval()
        try:
            with torch.no_grad():
                self.assertEqual(self.net(x).dtype, torch.float32)
        finally:
            self.net.train()

    def test_gradients_flow(self):
        import torch
        # Default UNet has 4 stages → min spatial dim = 2^4 = 16
        self.net(torch.randn(2, 4, 16, 16)).sum().backward()
        grads = [p.grad for p in self.net.parameters() if p.requires_grad]
        self.assertTrue(any(g is not None for g in grads))

    def test_param_count_positive(self):
        self.assertGreater(_param_count(self.net), 0)

    def test_multi_output_channels(self):
        import torch

        from pycsamt.ai.nets import UNet2DNet
        net = UNet2DNet(n_in=2, n_out=3).build()
        net.eval()  # single-sample batch: see test_output_dtype
        # default UNet has 4 stages -> minimum spatial dim is 16
        x = torch.randn(1, 2, 16, 16)
        with torch.no_grad():
            y = net(x)
        self.assertEqual(y.shape[1], 3)


@unittest.skipUnless(_HAS_TORCH, "PyTorch not installed")
class TestDRCNNNetBuild(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        from pycsamt.ai.nets import DRCNNNet
        # Two modalities: 40-feature MT + 20-feature TEM
        cls.net = DRCNNNet(n_features_list=(40, 20), n_out=9).build()

    def test_is_module(self):
        import torch.nn as nn
        self.assertIsInstance(self.net, nn.Module)

    def test_forward_shape_two_modalities(self):
        import torch
        x1 = torch.randn(4, 40)
        x2 = torch.randn(4, 20)
        y = self.net(x1, x2)
        self.assertEqual(y.shape, (4, 9))

    def test_output_dtype(self):
        import torch
        y = self.net(torch.randn(2, 40), torch.randn(2, 20))
        self.assertEqual(y.dtype, torch.float32)

    def test_wrong_modality_count_raises(self):
        import torch
        with self.assertRaises(ValueError):
            self.net(torch.randn(2, 40))   # expects 2 inputs, got 1

    def test_three_modalities(self):
        import torch

        from pycsamt.ai.nets import DRCNNNet
        net = DRCNNNet(n_features_list=(30, 20, 10), n_out=5).build()
        y = net(torch.randn(3, 30), torch.randn(3, 20), torch.randn(3, 10))
        self.assertEqual(y.shape, (3, 5))

    def test_gradients_flow(self):
        import torch
        x1, x2 = torch.randn(4, 40), torch.randn(4, 20)
        self.net(x1, x2).sum().backward()
        grads = [p.grad for p in self.net.parameters() if p.requires_grad]
        self.assertTrue(any(g is not None for g in grads))

    def test_param_count_positive(self):
        self.assertGreater(_param_count(self.net), 0)


@unittest.skipUnless(_HAS_TORCH, "PyTorch not installed")
class TestGCNNetBuild(unittest.TestCase):
    """GCN forward pass: H (n_nodes, n_features), A (n_nodes, n_nodes)."""

    _N_NODES = 9
    _N_FEAT  = 20
    _N_OUT   = 5

    @classmethod
    def setUpClass(cls):
        import torch

        from pycsamt.ai.nets import GCNNet, build_adjacency
        cls.net = GCNNet(n_features=cls._N_FEAT, n_out=cls._N_OUT,
                         hidden=(64, 32)).build()
        # Precompute a 3×3 grid adjacency
        xs, ys = np.meshgrid(np.arange(3) * 1000.0, np.arange(3) * 1000.0)
        coords = np.column_stack([xs.ravel(), ys.ravel()])
        A_np = build_adjacency(coords, radius=1500.0)
        cls.A = torch.tensor(A_np, dtype=torch.float32)

    def test_is_module(self):
        import torch.nn as nn
        self.assertIsInstance(self.net, nn.Module)

    def test_forward_shape(self):
        import torch
        H = torch.randn(self._N_NODES, self._N_FEAT)
        y = self.net(H, self.A)
        self.assertEqual(y.shape, (self._N_NODES, self._N_OUT))

    def test_output_dtype(self):
        import torch
        H = torch.randn(self._N_NODES, self._N_FEAT)
        y = self.net(H, self.A)
        self.assertEqual(y.dtype, torch.float32)

    def test_gradients_flow(self):
        import torch
        H = torch.randn(self._N_NODES, self._N_FEAT)
        self.net(H, self.A).sum().backward()
        grads = [p.grad for p in self.net.parameters() if p.requires_grad]
        self.assertTrue(any(g is not None for g in grads))

    def test_param_count_positive(self):
        self.assertGreater(_param_count(self.net), 0)

    def test_identity_adjacency(self):
        """With identity adjacency (no message passing) output is still valid."""
        import torch
        H = torch.randn(self._N_NODES, self._N_FEAT)
        A_eye = torch.eye(self._N_NODES)
        y = self.net(H, A_eye)
        self.assertEqual(y.shape, (self._N_NODES, self._N_OUT))
        self.assertTrue(torch.isfinite(y).all())

    def test_custom_hidden(self):
        import torch

        from pycsamt.ai.nets import GCNNet
        net = GCNNet(n_features=10, n_out=3, hidden=(32, 16)).build()
        H = torch.randn(5, 10)
        A = torch.eye(5)
        y = net(H, A)
        self.assertEqual(y.shape, (5, 3))

    def test_batch_via_loop(self):
        """Processing a batch survey-by-survey must give consistent shapes."""
        import torch
        batch_size = 4
        results = []
        for _ in range(batch_size):
            H = torch.randn(self._N_NODES, self._N_FEAT)
            y = self.net(H, self.A)
            results.append(y.shape)
        self.assertTrue(all(s == (self._N_NODES, self._N_OUT) for s in results))


# ─────────────────────────────────────────────────────────────────────────────
# Integration: generate_dataset_3d → build_adjacency → GCNNet forward
# ─────────────────────────────────────────────────────────────────────────────

@unittest.skipUnless(_HAS_TORCH, "PyTorch not installed")
class TestGCNPipelineIntegration(unittest.TestCase):
    """Verify the full 3-D dataset → adjacency → GCN forward pipeline."""

    @classmethod
    def setUpClass(cls):
        import torch

        from pycsamt.ai.nets import GCNNet, build_adjacency
        from pycsamt.forward.batch import generate_dataset_3d
        cls.ds = generate_dataset_3d(
            n_surveys=6,
            n_stations=9,
            n_layers=3,
            freqs=np.logspace(0, 3, 10),
            corr_length=2000.0,
            seed=0,
            verbose=False,
        )
        A_np = build_adjacency(cls.ds.coords, radius=5000.0)
        cls.A = torch.tensor(A_np, dtype=torch.float32)
        cls.net = GCNNet(
            n_features=cls.ds.n_features,
            n_out=cls.ds.n_params,
            hidden=(64, 32),
        ).build()

    def test_adjacency_shape(self):
        self.assertEqual(self.A.shape, (9, 9))

    def test_forward_per_survey(self):
        import torch
        ds = self.ds
        for i in range(len(ds)):
            H = torch.tensor(ds.X[i], dtype=torch.float32)
            y = self.net(H, self.A)
            self.assertEqual(y.shape, (9, ds.n_params))

    def test_output_finite(self):
        import torch
        H = torch.tensor(self.ds.X[0], dtype=torch.float32)
        y = self.net(H, self.A)
        self.assertTrue(torch.isfinite(y).all())


if __name__ == "__main__":
    unittest.main()
