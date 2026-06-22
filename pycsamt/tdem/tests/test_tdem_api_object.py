"""Tests for TDEM integration with common API object helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.api import PyCSAMTObject
from pycsamt.tdem import (
    TEMSounding,
    TEMtoEDI,
    read_temavg_survey,
    transform_temavg_survey,
)
from pycsamt.tdem.transform import FourierTransform, LateTimeTransform
from pycsamt.tdem.waveform import SquareWaveform

DATA_DIR = (
    Path(__file__).parents[3]
    / "data"
    / "TEMAVG"
    / "JIANGSU"
)
AVG_FILE = DATA_DIR / "TEM100.AVG"


def test_sounding_inherits_common_object_and_tracks_runtime_attrs():
    """TEMSounding should inherit the common API root."""
    sounding = TEMSounding(
        np.array([1e-3, 2e-3]),
        np.array([1.0, 0.5]),
        current=1.0,
        tx_area=100.0,
        verbose=2,
        logger=object(),
    )

    assert isinstance(sounding, PyCSAMTObject)
    assert sounding.verbose == 2
    assert sounding.logger is not None
    assert "verbose" not in repr(sounding)
    assert "logger" not in repr(sounding)


def test_transformers_and_waveforms_use_common_object_repr():
    """Runtime attributes should stay available but hidden."""
    tr = LateTimeTransform(verbose=3, logger=object())
    ft = FourierTransform(verbose=4, logger=object())
    conv = TEMtoEDI(verbose=5, logger=object())
    wave = SquareWaveform(verbose=6, logger=object())

    for obj in (tr, ft, conv, wave):
        assert isinstance(obj, PyCSAMTObject)
        assert obj.verbose > 0
        assert obj.logger is not None
        assert "verbose" not in repr(obj)
        assert "logger" not in repr(obj)


@pytest.mark.skipif(
    not AVG_FILE.exists(),
    reason="TEMAVG sample data not available",
)
def test_temavg_survey_propagates_verbose_and_logger():
    """Folder readers should pass runtime attrs to child objects."""
    logger = object()

    survey = read_temavg_survey(
        DATA_DIR,
        pattern="TEM100.AVG",
        verbose=7,
        logger=logger,
    )
    sounding = survey.to_soundings()[0]

    assert isinstance(survey, PyCSAMTObject)
    assert survey.verbose == 7
    assert survey.logger is logger
    assert survey.get("TEM100").verbose == 7
    assert sounding.verbose == 7
    assert sounding.logger is logger
    assert "verbose" not in repr(survey)
    assert "logger" not in repr(survey)


@pytest.mark.skipif(
    not AVG_FILE.exists(),
    reason="TEMAVG sample data not available",
)
def test_temavg_workflow_bundle_tracks_runtime_attrs():
    """Workflow output bundles should inherit common object behavior."""
    logger = object()

    bundle = transform_temavg_survey(
        DATA_DIR,
        stems=["TEM100"],
        verbose=8,
        logger=logger,
    )

    assert isinstance(bundle, PyCSAMTObject)
    assert bundle.verbose == 8
    assert bundle.logger is logger
    assert bundle.soundings[0].verbose == 8
    assert "verbose" not in repr(bundle)
    assert "logger" not in repr(bundle)

