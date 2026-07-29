"""API-contract tests for :mod:`pycsamt.ai.inversion`.

These tests intentionally avoid fitting neural networks.  They protect the
public AI inversion configuration and estimator interfaces on machines where a
deep-learning backend may not be installed.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from pycsamt.ai.inversion import (
    EMInverter1D,
    EMInverter2D,
    EnsembleInverter,
    GCNInverter3D,
    InversionConfig,
    JointInverter,
    RunConfig,
)
from pycsamt.forward.config import ForwardConfig


def test_ai_inversion_package_exports_public_api():
    """The package-level imports should expose the documented user API."""
    assert EMInverter1D.__name__ == "EMInverter1D"
    assert EMInverter2D.__name__ == "EMInverter2D"
    assert GCNInverter3D.__name__ == "GCNInverter3D"
    assert JointInverter.__name__ == "JointInverter"
    assert EnsembleInverter.__name__ == "EnsembleInverter"
    assert InversionConfig.__name__ == "InversionConfig"
    assert RunConfig.__name__ == "RunConfig"


def test_inversion_config_builds_inverter_and_fit_kwargs(tmp_path):
    cfg = InversionConfig(
        arch="cnn1d",
        n_layers=4,
        solver="csamt1d",
        device="cpu",
        include_phase=False,
        log_thickness=False,
        augment_noise=0.0,
        epochs=12,
        batch_size=8,
        lr=2e-3,
        patience=3,
        val_frac=0.25,
        grad_clip=None,
        seed=42,
        checkpoint_dir=str(tmp_path),
        checkpoint_name="csamt_cnn",
        verbose=False,
    )

    cfg.validate()
    inv = cfg.to_inverter()

    assert isinstance(inv, EMInverter1D)
    assert inv.arch == "cnn1d"
    assert inv.n_layers == 4
    assert inv.solver == "csamt1d"
    assert inv.device == "cpu"
    assert inv.include_phase is False
    assert inv.log_thickness is False
    assert inv.augment_noise == 0.0

    assert cfg.to_fit_kwargs() == {
        "epochs": 12,
        "batch_size": 8,
        "lr": 2e-3,
        "patience": 3,
        "val_frac": 0.25,
        "grad_clip": None,
        "seed": 42,
        "verbose": False,
    }
    assert cfg.checkpoint_path() == tmp_path / "csamt_cnn.npz"


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("arch", "transformer", "arch"),
        ("solver", "modem", "solver"),
        ("n_layers", 1, "n_layers"),
        ("device", "gpu", "device"),
        ("epochs", 0, "epochs"),
        ("batch_size", 0, "batch_size"),
        ("lr", 0.0, "lr"),
        ("weight_decay", -1.0, "weight_decay"),
        ("patience", 0, "patience"),
        ("min_delta", -1e-3, "min_delta"),
        ("val_frac", 1.0, "val_frac"),
        ("grad_clip", 0.0, "grad_clip"),
        ("augment_noise", -0.01, "augment_noise"),
    ],
)
def test_inversion_config_validate_rejects_bad_values(field, value, message):
    cfg = InversionConfig()
    setattr(cfg, field, value)

    with pytest.raises(ValueError, match=message):
        cfg.validate()


def test_inversion_config_json_roundtrip_and_strict_unknown_keys(tmp_path):
    cfg = InversionConfig(
        arch="fcn",
        n_layers=6,
        solver="tem1d",
        include_phase=False,
        checkpoint_dir=str(tmp_path / "ckpt"),
        checkpoint_name="tem_fcn",
        seed=7,
    )

    path = cfg.to_template(tmp_path / "inv_config.json")
    loaded = InversionConfig.from_file(path)

    assert loaded.arch == "fcn"
    assert loaded.n_layers == 6
    assert loaded.solver == "tem1d"
    assert loaded.include_phase is False
    assert loaded.checkpoint_path() == tmp_path / "ckpt" / "tem_fcn.npz"
    assert "tem1d" in loaded.summary()

    bad_path = tmp_path / "bad_inv_config.json"
    bad_path.write_text(
        json.dumps({"config": {"arch": "resnet", "unknown_ai_key": 1}}),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="unknown_ai_key"):
        InversionConfig.from_file(bad_path)

    relaxed = InversionConfig.from_file(bad_path, strict=False)
    assert relaxed.arch == "resnet"


def test_inversion_config_from_inverter_snapshots_architecture():
    inv = EMInverter1D(
        arch="cnn1d",
        n_layers=3,
        solver="mt1d",
        device="cpu",
        include_phase=False,
        log_thickness=False,
        augment_noise=0.05,
    )

    cfg = InversionConfig.from_inverter(inv)

    assert cfg.arch == "cnn1d"
    assert cfg.n_layers == 3
    assert cfg.solver == "mt1d"
    assert cfg.device == "cpu"
    assert cfg.include_phase is False
    assert cfg.log_thickness is False
    assert cfg.augment_noise == 0.05


def test_run_config_delegates_and_roundtrips_json(tmp_path):
    forward = ForwardConfig(
        solver="mt1d",
        n_layers_min=5,
        n_layers_max=5,
        n_samples=16,
        n_freqs=8,
        include_phase=True,
        output_dir=None,
        verbose=False,
    )
    inversion = InversionConfig(
        solver="mt1d",
        n_layers=5,
        include_phase=True,
        epochs=5,
        batch_size=4,
        checkpoint_dir=str(tmp_path / "models"),
        checkpoint_name="mt1d_resnet",
        verbose=False,
    )
    run = RunConfig(
        forward=forward,
        inversion=inversion,
        name="mt1d_contract",
        description="small API contract test",
    )

    run.validate()
    dataset_kwargs = run.to_dataset_kwargs()
    fit_kwargs = run.to_fit_kwargs()

    assert dataset_kwargs["solver"] == "mt1d"
    assert dataset_kwargs["n_layers"] == 5
    assert dataset_kwargs["include_phase"] is True
    assert dataset_kwargs["output"] is None
    assert fit_kwargs["epochs"] == 5
    assert fit_kwargs["batch_size"] == 4
    assert isinstance(run.to_inverter(), EMInverter1D)
    assert run.checkpoint_path() == tmp_path / "models" / "mt1d_resnet.npz"
    assert "mt1d_contract" in run.summary()

    path = run.to_template(tmp_path / "run_config.json")
    loaded = RunConfig.from_file(path)
    loaded.validate()

    assert loaded.forward.solver == "mt1d"
    assert loaded.forward.n_layers_min == 5
    assert loaded.forward.n_layers_max == 5
    assert loaded.inversion.n_layers == 5
    assert loaded.inversion.epochs == 5


@pytest.mark.parametrize(
    ("forward", "inversion", "message"),
    [
        (
            ForwardConfig(solver="mt1d"),
            InversionConfig(solver="tem1d"),
            "forward.solver",
        ),
        (
            ForwardConfig(solver="mt1d", include_phase=False),
            InversionConfig(solver="mt1d", include_phase=True),
            "include_phase",
        ),
        (
            ForwardConfig(solver="mt1d", n_layers_min=4, n_layers_max=4),
            InversionConfig(solver="mt1d", n_layers=5),
            "fixed layer count",
        ),
    ],
)
def test_run_config_validate_rejects_incompatible_subconfigs(
    forward,
    inversion,
    message,
):
    with pytest.raises(ValueError, match=message):
        RunConfig(forward=forward, inversion=inversion).validate()


def test_run_config_json_strict_unknown_keys(tmp_path):
    path = tmp_path / "bad_run_config.json"
    payload = {
        "forward": {"solver": "mt1d", "unexpected_forward_key": True},
        "inversion": {"solver": "mt1d"},
    }
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="unexpected_forward_key"):
        RunConfig.from_file(path)

    loaded = RunConfig.from_file(path, strict=False)
    assert loaded.forward.solver == "mt1d"
    assert loaded.inversion.solver == "mt1d"


@pytest.mark.parametrize("suffix", [".py", ".json", ".yml"])
def test_run_config_templates_are_utf8_roundtrips(tmp_path, suffix):
    path = tmp_path / f"unicode_run{suffix}"
    written = RunConfig.write_template(
        path,
        name="résistivité",
        description="Frequency–depth contract in Ω·m",
    )

    text = written.read_text(encoding="utf-8")
    assert "résistivité" in text
    loaded = RunConfig.from_file(written)
    loaded.validate()


def test_em_inverter_2d_interface_without_backend():
    inv = EMInverter2D(
        n_components=3,
        n_depth=12,
        n_stations=7,
        n_freqs=9,
        device="cpu",
        log_rho_out=False,
        channels=(8, 16),
    )

    assert inv.arch == "unet"
    assert inv.solver == "mt2d"
    assert inv.n_layers == 12
    assert inv.n_components == 3
    assert inv.n_depth == 12
    assert inv.n_stations == 7
    assert inv.n_freqs == 9
    assert inv.log_rho_out is False
    assert "unfitted" in repr(inv)

    params = inv._get_params()
    assert params["n_components"] == 3
    assert params["n_depth"] == 12
    assert params["n_stations"] == 7
    assert params["n_freqs"] == 9
    # channels may be normalised to a list (JSON-safe params)
    assert tuple(params["channels"]) == (8, 16)

    with pytest.raises(RuntimeError, match=r"fit\(\) before predict"):
        inv.predict(np.ones((2, 3, 9, 7), dtype=np.float32))


def test_joint_inverter_interface_without_backend():
    inv = JointInverter(
        n_features_list=(12, 5),
        n_layers=4,
        growth_rate=16,
        n_dense_layers=3,
        hidden_dim=64,
        dropout=0.1,
        device="cpu",
        log_thickness=False,
    )

    assert inv.arch == "drcnn"
    assert inv.solver == "joint"
    assert inv.n_features_list == (12, 5)
    assert inv.n_layers == 4
    assert inv._n_out == 7
    assert inv.growth_rate == 16
    assert inv.n_dense_layers == 3
    assert inv.hidden_dim == 64
    assert inv.dropout == 0.1
    assert inv.log_thickness is False
    assert "modalities=2" in repr(inv)
    assert "unfitted" in repr(inv)

    params = inv._get_params()
    assert params["n_features_list"] == [12, 5]
    assert params["n_layers"] == 4
    assert params["hidden_dim"] == 64

    with pytest.raises(RuntimeError, match=r"fit\(\) before predict"):
        inv.predict(
            [
                np.ones((3, 12), dtype=np.float32),
                np.ones((3, 5), dtype=np.float32),
            ]
        )
