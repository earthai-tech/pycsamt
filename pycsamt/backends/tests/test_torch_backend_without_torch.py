"""Coverage for :mod:`pycsamt.backends._torch` without PyTorch installed.

These tests use small protocol fakes rather than importing the optional
framework.  They keep backend dispatch and wrapper logic covered in minimal
CI environments; numerical integration remains in ``test_torch_backend``.
"""

from __future__ import annotations

import sys
import types
from contextlib import nullcontext
from unittest import mock

import numpy as np
import pytest

from pycsamt.backends._torch import TorchBackend


class _Builder:
    calls = []

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs
        self.__class__.calls.append((args, kwargs))

    def build(self):
        return self.args, self.kwargs


@pytest.mark.parametrize(
    "arch,method",
    [
        ("cnn", "_build_cnn1d"),
        ("resnet", "_build_resnet1d"),
        ("fcn", "_build_fcn1d"),
        ("unet", "_build_unet2d"),
        ("drcnn", "_build_drcnn"),
    ],
)
def test_build_dispatch_does_not_need_torch(arch, method):
    backend = TorchBackend()
    with (
        mock.patch.object(backend, "require"),
        mock.patch.object(backend, method, return_value=arch) as builder,
    ):
        assert backend.build({"arch": arch}) == arch
        builder.assert_called_once_with({"arch": arch})


def test_build_rejects_unknown_arch_without_importing_torch():
    with mock.patch.object(TorchBackend, "require"):
        with pytest.raises(ValueError, match="Unknown arch"):
            TorchBackend().build({"arch": "unknown"})


@pytest.mark.parametrize(
    "module_name,class_name,method,spec",
    [
        ("cnn1d", "CNN1DNet", "_build_cnn1d", {"n_features": 4, "n_out": 2}),
        ("resnet", "ResNet1DNet", "_build_resnet1d", {"n_features": 4, "n_out": 2}),
        ("fcn", "FCN1DNet", "_build_fcn1d", {"n_features": 4, "n_out": 2}),
        ("unet", "UNet2DNet", "_build_unet2d", {"n_in": 1}),
        (
            "drcnn",
            "DRCNNNet",
            "_build_drcnn",
            {"n_features_list": [2, 3], "n_out": 1},
        ),
    ],
)
def test_architecture_helpers_use_lazy_builder_modules(
    monkeypatch, module_name, class_name, method, spec
):
    module = types.ModuleType(f"pycsamt.ai.nets.{module_name}")
    setattr(module, class_name, _Builder)
    monkeypatch.setitem(sys.modules, module.__name__, module)
    result = getattr(TorchBackend(), method)(spec)
    assert isinstance(result, tuple)


def test_availability_and_device_resolution_with_fake_torch(monkeypatch):
    fake = types.ModuleType("torch")
    fake.cuda = types.SimpleNamespace(is_available=lambda: False)
    fake.backends = types.SimpleNamespace(
        mps=types.SimpleNamespace(is_available=lambda: True)
    )
    monkeypatch.setitem(sys.modules, "torch", fake)
    with mock.patch("pycsamt.backends._torch.probe_backend", return_value=True):
        assert TorchBackend.is_available()
    assert TorchBackend().resolve_device() == "mps"
    assert TorchBackend().resolve_device("cuda:7") == "cuda:7"


class _Tensor:
    def __init__(self, value):
        self.value = np.asarray(value)
        self.device = "cpu"

    def to(self, device):
        self.device = device
        return self

    def cpu(self):
        return self

    def numpy(self):
        return self.value

    def __len__(self):
        return len(self.value)


class _Loss:
    def __init__(self, value):
        self.value = value

    def backward(self):
        return None

    def item(self):
        return self.value


class _Model:
    def __init__(self):
        self.weight = _Tensor([1.0])
        self.loaded = None

    def parameters(self):
        return iter([self.weight])

    def to(self, device):
        self.weight.to(device)
        return self

    def train(self):
        return None

    def eval(self):
        return None

    def __call__(self, tensor):
        return _Tensor(tensor.value * 2)

    def state_dict(self):
        return {"weight": self.weight}

    def load_state_dict(self, state):
        self.loaded = state


def _fake_torch_modules(monkeypatch):
    torch = types.ModuleType("torch")
    nn = types.ModuleType("torch.nn")
    utils = types.ModuleType("torch.utils")
    data = types.ModuleType("torch.utils.data")

    class TensorDataset:
        def __init__(self, *tensors):
            self.tensors = tensors

    class DataLoader:
        def __init__(self, dataset, **kwargs):
            self.dataset = dataset

        def __iter__(self):
            yield self.dataset.tensors

    class Adam:
        def __init__(self, parameters, lr):
            list(parameters)

        def zero_grad(self):
            return None

        def step(self):
            return None

    class Scheduler:
        def __init__(self, optimizer, **kwargs):
            pass

        def step(self, value):
            return None

    data.TensorDataset = TensorDataset
    data.DataLoader = DataLoader
    nn.utils = types.SimpleNamespace(clip_grad_norm_=lambda *args: None)
    nn.MSELoss = type("MSELoss", (), {})
    nn.CrossEntropyLoss = type("CrossEntropyLoss", (), {})
    torch.nn = nn
    torch.from_numpy = _Tensor
    torch.no_grad = nullcontext
    torch.optim = types.SimpleNamespace(
        Adam=Adam,
        lr_scheduler=types.SimpleNamespace(ReduceLROnPlateau=Scheduler),
    )
    for name, module in {
        "torch": torch,
        "torch.nn": nn,
        "torch.utils": utils,
        "torch.utils.data": data,
    }.items():
        monkeypatch.setitem(sys.modules, name, module)
    return torch


def test_train_loop_with_fake_torch_covers_wrapper_logic(monkeypatch, capsys):
    _fake_torch_modules(monkeypatch)
    backend = TorchBackend()
    values = iter([0.5, 0.4, 0.6, 0.7])

    def loss_fn(prediction, target):
        return _Loss(next(values))

    model = _Model()
    data = np.ones((2, 2), dtype=np.float32)
    with (
        mock.patch.object(backend, "require"),
        mock.patch.object(backend, "_get_loss_fn", return_value=loss_fn),
    ):
        history = backend.train(
            model,
            data,
            data,
            data,
            data,
            epochs=2,
            batch_size=2,
            patience=1,
            grad_clip=1.0,
            device="cpu",
            verbose=True,
        )
    assert history == {"train_loss": [0.5, 0.6], "val_loss": [0.4, 0.7]}
    assert model.loaded is not None
    assert "Early stop at epoch 2" in capsys.readouterr().out


def test_predict_and_weight_roundtrip_with_fake_torch(monkeypatch):
    _fake_torch_modules(monkeypatch)
    backend = TorchBackend()
    model = _Model()
    with mock.patch.object(backend, "require"):
        result = backend.predict(model, np.arange(6).reshape(3, 2), batch_size=2)
    np.testing.assert_array_equal(result, np.arange(6).reshape(3, 2) * 2)

    weights = backend.get_weights(model)
    np.testing.assert_array_equal(weights["weight"], [1.0])
    backend.set_weights(model, {"weight": np.array([3.0])})
    np.testing.assert_array_equal(model.loaded["weight"].numpy(), [3.0])


def test_loss_dispatch_with_fake_nn(monkeypatch):
    _fake_torch_modules(monkeypatch)
    assert TorchBackend._get_loss_fn("mse").__class__.__name__ == "MSELoss"
    assert (
        TorchBackend._get_loss_fn("cross_entropy").__class__.__name__
        == "CrossEntropyLoss"
    )
    with pytest.raises(ValueError, match="Unknown loss"):
        TorchBackend._get_loss_fn("bogus")
