# -*- coding: utf-8 -*-
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pycsamt.site import edit as ed


def test_maybe_df_detects_common_names_from_dicts() -> None:
    table = [
        {"station": "S01", "lat": 10.0, "lon": 20.0, "elev": 100.0},
        {"station": "S02", "lat": 11.0, "lon": 21.0, "elev": 110.0},
    ]
    df, cols = ed._maybe_df(table)

    # resolved columns
    assert cols["station"] == "station"
    assert cols["lat"] == "lat"
    assert cols["lon"] == "lon"
    assert cols["elev"] == "elev"
    assert cols["easting"] is None
    assert cols["northing"] is None

    # dataframe content
    assert set(df.columns) >= {"station", "lat", "lon", "elev"}
    assert df.shape[0] == 2
    assert df["station"].tolist() == ["S01", "S02"]
    assert np.allclose(df["lat"].to_numpy(), [10.0, 11.0])
    assert np.allclose(df["lon"].to_numpy(), [20.0, 21.0])
    assert np.allclose(df["elev"].to_numpy(), [100.0, 110.0])


def test_maybe_df_detects_aliases_from_csv(tmp_path: Path) -> None:
    csvp = tmp_path / "coords.csv"
    # aliases: name, latitude, long, elevation
    csvp.write_text(
        "name,latitude,long,elevation\n"
        "A01,35.1,12.8,123.0\n"
        "A02,35.2,12.9,124.0\n",
        encoding="utf-8",
    )
    df, cols = ed._maybe_df(csvp)

    assert cols["station"] == "name"
    assert cols["lat"] == "latitude"
    assert cols["lon"] == "long"
    assert cols["elev"] == "elevation"
    assert df.shape[0] == 2
    assert df["name"].tolist() == ["A01", "A02"]


def test_frame_to_mapping_with_latlon() -> None:
    df = pd.DataFrame(
        {
            "station": ["P01", "P02"],
            "lat": [40.0, 40.1],
            "lon": [8.0, 8.1],
            "elev": [5.0, 10.0],
        }
    )
    cols = {
        "station": "station",
        "lat": "lat",
        "lon": "lon",
        "elev": "elev",
        "easting": None,
        "northing": None,
    }
    mp = ed._frame_to_mapping(df, cols)
    assert set(mp.keys()) == {"P01", "P02"}
    assert mp["P01"] == (40.0, 8.0, 5.0)
    assert mp["P02"] == (40.1, 8.1, 10.0)


def test_frame_to_mapping_requires_crs_for_en() -> None:
    df = pd.DataFrame(
        {
            "site": ["Q01"],
            "easting": [400000.0],
            "northing": [5750000.0],
        }
    )
    cols = {
        "station": "site",
        "lat": None,
        "lon": None,
        "elev": None,
        "easting": "easting",
        "northing": "northing",
    }
    with pytest.raises(ValueError):
        _ = ed._frame_to_mapping(df, cols, crs_from=None)


def test_project_en_to_lonlat_roundtrip_if_pyproj() -> None:
    pyproj = pytest.importorskip("pyproj")
    from pyproj import Transformer

    # choose a lon/lat inside UTM zone 31N (EPSG:32631)
    lon0, lat0 = 3.0, 45.0
    t_fwd = Transformer.from_crs("EPSG:4326", "EPSG:32631",
                                 always_xy=True)
    e0, n0 = t_fwd.transform(lon0, lat0)

    lon, lat = ed._project_en_to_lonlat(
        np.asarray([e0], float),
        np.asarray([n0], float),
        crs_from="EPSG:32631",
    )
    assert lon.shape == (1,)
    assert lat.shape == (1,)
    assert lon[0] == pytest.approx(lon0, rel=0, abs=1e-5)
    assert lat[0] == pytest.approx(lat0, rel=0, abs=1e-5)


def test_set_coords_from_table_invokes_set_coords_all(
    monkeypatch,
) -> None:
    df = pd.DataFrame(
        {
            "station": ["S1", "S2"],
            "lat": [1.0, 2.0],
            "lon": [10.0, 20.0],
            "elev": [100.0, 200.0],
        }
    )

    calls = {}

    def fake_set_coords_all(sites, mapping, *, inplace):
        calls["sites"] = sites
        calls["mapping"] = mapping
        calls["inplace"] = inplace
        return "OK"

    monkeypatch.setattr(ed, "set_coords_all", fake_set_coords_all)

    out = ed.set_coords_from_table("SITES_OBJ", df, inplace=True)
    assert out == "OK"
    assert calls["sites"] == "SITES_OBJ"
    assert calls["inplace"] is True
    assert set(calls["mapping"].keys()) == {"S1", "S2"}
    assert calls["mapping"]["S1"] == (1.0, 10.0, 100.0)
    assert calls["mapping"]["S2"] == (2.0, 20.0, 200.0)


def test_set_coords_from_en_invokes_set_coords(monkeypatch) -> None:
    # avoid pyproj dependency here by faking the projection step
    def fake_proj(e, n, crs_from, to_crs="EPSG:4326"):
        return np.asarray([20.0]), np.asarray([10.0])

    called = {}

    def fake_set_coords(site, *, lat, lon, elev, inplace):
        called["site"] = site
        called["lat"] = lat
        called["lon"] = lon
        called["elev"] = elev
        called["inplace"] = inplace
        return "DONE"

    monkeypatch.setattr(ed, "_project_en_to_lonlat", fake_proj)
    monkeypatch.setattr(ed, "set_coords", fake_set_coords)

    out = ed.set_coords_from_en(
        "EDI_OBJ",
        easting=400000.0,
        northing=5750000.0,
        crs_from="EPSG:32631",
        elev=250.0,
        inplace=False,
    )
    assert out == "DONE"
    # projected lon,lat -> 20,10 from fake
    assert called["site"] == "EDI_OBJ"
    assert called["lat"] == pytest.approx(10.0)
    assert called["lon"] == pytest.approx(20.0)
    assert called["elev"] == pytest.approx(250.0)
    assert called["inplace"] is False
