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
from pycsamt.tdem.transform import (
    FourierTransform,
    LateTimeTransform,
)
from pycsamt.tdem.waveform import SquareWaveform

DATA_DIR = Path(__file__).parents[3] / "data" / "TEMAVG" / "JIANGSU"
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


def test_temtoedi_writes_a_valid_reloadable_edi_with_local_coordinates(
    tmp_path,
):
    """EDI files must reload even when x/y are local survey metres.

    Regression test: ``TEMSounding.x``/``.y`` are commonly local survey
    coordinates (tens to thousands of metres), not geographic decimal
    degrees. ``TEMtoEDI`` used to assign them straight to
    ``Head.lat``/``Head.long``; the resulting ``ValueError`` was caught
    by a bare ``except Exception: pass`` around the whole ``>HEAD``
    construction, so the file silently ended up with no ``>HEAD``
    section at all (and no ``>=DEFINEMEAS``/``>=MTSECT`` section either),
    making it unreadable by :func:`~pycsamt.emtools.ensure_sites`.
    """
    from pycsamt.emtools import ensure_sites

    sounding = TEMSounding(
        time_gates=np.logspace(-5, -2, 20),
        data=5e-5 * np.logspace(-5, -2, 20) ** (-5.0 / 2.0),
        current=8.0,
        tx_area=1.0e4,
        station_name="LOCAL01",
        x=250.0,  # local survey metres, not longitude
        y=1080.0,  # local survey metres, not latitude
        elevation=1090.0,
    )

    conv = TEMtoEDI(method="late_time", out_dir=str(tmp_path))
    written = conv.save(sounding)
    assert len(written) == 1

    text = Path(written[0]).read_text()
    assert text.startswith(">HEAD")
    assert ">=DEFINEMEAS" in text
    assert ">=MTSECT" in text or ">=EMAPSECT" in text

    sites = ensure_sites(tmp_path, recursive=False)
    assert len(sites) == 1
    assert sites[0].z.shape[0] == sounding.n_gates


def test_temtoedi_keeps_geographic_lat_lon_when_in_range():
    """Genuine geographic coordinates should still reach Head.lat/long."""
    sounding = TEMSounding(
        time_gates=np.logspace(-5, -2, 15),
        data=5e-5 * np.logspace(-5, -2, 15) ** (-5.0 / 2.0),
        current=8.0,
        tx_area=1.0e4,
        station_name="GEO01",
        x=113.487,  # longitude
        y=26.052,  # latitude
        elevation=574.5,
    )

    conv = TEMtoEDI(method="late_time")
    edi = conv._result_to_edifile(
        conv._transformer.transform(sounding)
    )
    head = edi.get_section("head")
    assert head.long == pytest.approx(113.487, abs=1e-9)
    assert head.lat == pytest.approx(26.052, abs=1e-9)
