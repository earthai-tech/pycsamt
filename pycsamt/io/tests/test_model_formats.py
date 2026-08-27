# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.io.model_formats — inversion-model format registry."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.io import (
    ModelFormatError,
    detect_model_format,
    get_model_format,
    get_model_format_for_target,
    list_model_formats,
    register_model_format,
)


def test_pcsf_is_registered():
    formats = list_model_formats()
    assert "pcsf" in formats
    assert formats["pcsf"]["readable"]
    assert formats["pcsf"]["writable"]
    assert ".pcsf" in formats["pcsf"]["extensions"]


def test_get_model_format_for_target_resolves_pcsf(tmp_path):
    spec = get_model_format_for_target(tmp_path / "model.pcsf")
    assert spec.name == "pcsf"


def test_get_model_format_for_target_rejects_unknown_extension(tmp_path):
    with pytest.raises(ModelFormatError):
        get_model_format_for_target(tmp_path / "model.unknownext")


def test_detect_model_format_reads_a_real_pcsf_file(tmp_path):
    from pycsamt.format import Grid2DGeometry, PCSFModel, write_pcsf

    model = PCSFModel(
        geometry=Grid2DGeometry(x=np.array([0.0, 1.0]), z=np.array([0.0, 1.0])),
        resistivity=np.ones((2, 2)) * 100.0,
    )
    path = write_pcsf(model, tmp_path / "detect.pcsf")
    assert detect_model_format(path) == "pcsf"


def test_detect_model_format_rejects_non_hdf5_file(tmp_path):
    path = tmp_path / "not_hdf5.pcsf"
    path.write_text("this is plain text, not HDF5")
    with pytest.raises(ModelFormatError):
        detect_model_format(path)


def test_pcsf_reader_writer_round_trip_through_registry(tmp_path):
    from pycsamt.format import Grid2DGeometry, PCSFModel

    spec = get_model_format("pcsf")
    model = PCSFModel(
        geometry=Grid2DGeometry(x=np.array([0.0, 1.0]), z=np.array([0.0, 1.0])),
        resistivity=np.ones((2, 2)) * 50.0,
    )
    path = tmp_path / "registry_roundtrip.pcsf"
    spec.writer(model, path)
    restored = spec.reader(path)
    np.testing.assert_allclose(restored.resistivity, model.resistivity)


def test_register_model_format_rejects_non_callable_reader():
    with pytest.raises(TypeError):
        register_model_format("bogus", reader="not-callable")


def _pcsm_model():
    from pycsamt.format import Grid2DGeometry, PCSFModel

    return PCSFModel(
        geometry=Grid2DGeometry(x=np.array([0.0, 1.0]), z=np.array([0.0, 1.0])),
        resistivity=np.ones((2, 2)) * 100.0,
    )


def test_pcsm_is_registered():
    formats = list_model_formats()
    assert "pcsm" in formats
    assert formats["pcsm"]["readable"]
    assert formats["pcsm"]["writable"]
    assert ".pcsm" in formats["pcsm"]["extensions"]


def test_get_model_format_for_target_resolves_pcsm(tmp_path):
    spec = get_model_format_for_target(tmp_path / "model.pcsm")
    assert spec.name == "pcsm"


def test_get_model_format_for_target_resolves_gzipped_pcsm(tmp_path):
    # Path("model.pcsm.gz").suffix is ".gz" alone (only the last
    # suffix) -- the registry's ".gz" extension hint covers this.
    spec = get_model_format_for_target(tmp_path / "model.pcsm.gz")
    assert spec.name == "pcsm"


def test_detect_model_format_reads_a_real_pcsm_file(tmp_path):
    from pycsamt.format import write_pcsm

    path = write_pcsm(_pcsm_model(), tmp_path / "detect.pcsm")
    assert detect_model_format(path) == "pcsm"


def test_detect_model_format_reads_a_real_gzipped_pcsm_file(tmp_path):
    from pycsamt.format import write_pcsm

    path = write_pcsm(_pcsm_model(), tmp_path / "detect.pcsm.gz")
    assert detect_model_format(path) == "pcsm"


def test_detect_model_format_does_not_confuse_pcsf_and_pcsm(tmp_path):
    from pycsamt.format import write_pcsf, write_pcsm

    model = _pcsm_model()
    pcsf_path = write_pcsf(model, tmp_path / "a.pcsf")
    pcsm_path = write_pcsm(model, tmp_path / "a.pcsm")
    assert detect_model_format(pcsf_path) == "pcsf"
    assert detect_model_format(pcsm_path) == "pcsm"


def test_pcsm_reader_writer_round_trip_through_registry(tmp_path):
    spec = get_model_format("pcsm")
    model = _pcsm_model()
    path = tmp_path / "registry_roundtrip.pcsm"
    spec.writer(model, path)
    restored = spec.reader(path)
    np.testing.assert_allclose(restored.resistivity, model.resistivity)
