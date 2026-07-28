# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Unit tests for :mod:`pycsamt.gis.config`.

Covers :class:`ElevationAPIConfig` lookup helpers and the module-level
``HAS_GDAL`` / ``EPSG_DICT`` exports built at import time.
"""

from __future__ import annotations

import importlib
from unittest.mock import patch

import pytest

from pycsamt.gis.config import (
    EPSG_DICT,
    HAS_GDAL,
    ElevationAPIConfig,
    GDALMissingError,
    GisError,
)


def test_has_gdal_is_bool():
    assert isinstance(HAS_GDAL, bool)


def test_epsg_dict_is_populated():
    # Either from pyproj's legacy 'epsg' file or the bundled epsg.npy
    # fallback -- either way it must be a non-empty int->str mapping.
    assert isinstance(EPSG_DICT, dict)
    assert len(EPSG_DICT) > 0
    a_key = next(iter(EPSG_DICT))
    assert isinstance(a_key, int)
    assert isinstance(EPSG_DICT[a_key], str)


def test_gis_error_and_gdal_missing_error_are_exceptions():
    assert issubclass(GisError, Exception)
    assert issubclass(GDALMissingError, RuntimeError)
    with __import__("pytest").raises(GisError):
        raise GisError("boom")


def test_elevation_api_config_default_api():
    cfg = ElevationAPIConfig.get_api_config()
    assert cfg == ElevationAPIConfig.APIS[ElevationAPIConfig.DEFAULT_API]
    assert cfg["url"] == "https://api.open-meteo.com/v1/elevation"


def test_elevation_api_config_named_apis():
    for name in ("open_meteo", "open_topo_data", "usgs_ned"):
        cfg = ElevationAPIConfig.get_api_config(name)
        assert cfg == ElevationAPIConfig.APIS[name]
        assert "url" in cfg
        assert "response_key" in cfg


def test_elevation_api_config_unknown_name_falls_back_to_default():
    cfg = ElevationAPIConfig.get_api_config("not-a-real-api")
    assert cfg == ElevationAPIConfig.APIS[ElevationAPIConfig.DEFAULT_API]


# ----------------- import-time branch coverage (reload-based) -----------------
#
# ``pycsamt.gis.config`` makes its GDAL/pyproj availability decisions at
# *import* time. This test environment has no GDAL, so the "GDAL is
# available" branch and the legacy pyproj ``epsg`` flat-file parsing branch
# are never exercised by a normal import. We force both by reloading the
# module under a patched environment and restore the pristine module
# afterwards so later tests see the real state.


@pytest.fixture()
def _reload_config_module():
    import pycsamt.gis.config as config_mod

    yield config_mod
    # Restore genuine environment state for subsequent tests/imports.
    importlib.reload(config_mod)


def test_config_reload_with_gdal_available_branch(_reload_config_module):
    config_mod = _reload_config_module
    with patch(
        "pycsamt.decorators.GdalDataCheck._locate_gdal_data",
        return_value=True,
    ):
        importlib.reload(config_mod)
    assert config_mod.HAS_GDAL is True


def test_config_reload_with_legacy_pyproj_epsg_file(tmp_path, _reload_config_module):
    config_mod = _reload_config_module
    epsg_file = tmp_path / "epsg"
    epsg_file.write_text(
        "# comment line, ignored\n"
        "no angle-bracket codes on this line, skipped\n"
        "<4326> +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs <>\n"
        "<32611> +proj=utm +zone=11 +datum=WGS84 +units=m +no_defs <>\n",
        encoding="utf-8",
    )
    with patch("pyproj.pyproj_datadir", str(tmp_path), create=True):
        importlib.reload(config_mod)
    assert config_mod.EPSG_DICT[4326].startswith("+proj=longlat")
    assert config_mod.EPSG_DICT[32611].startswith("+proj=utm")
