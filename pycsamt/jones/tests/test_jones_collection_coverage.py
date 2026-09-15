from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.jones.collection import JCollection
from pycsamt.jones.j import JFile


# ─────────────────────────────────────────────────────────────────────────
# JCollectionMixin: add_from / sort / paths / nf_stats
# ─────────────────────────────────────────────────────────────────────────


def test_add_from_populates_collection_and_logs_errors(
    j_single_file: Path, tmp_path: Path
):
    bad = tmp_path / "garbage.j"
    bad.write_text("not a real J file\n", encoding="utf-8")
    col = JCollection(verbose=0)
    col.add_from([j_single_file, bad])
    assert len(col) >= 1


def test_sort_by_station_default_key(jc_files):
    col = JCollection.from_sources(jc_files, verbose=0)
    out = col.sort()
    assert isinstance(out, JCollection)
    sites = [getattr(jf, "site", None) for jf in out]
    assert sites == sorted(sites)


def test_sort_by_callable_key(jc_files):
    col = JCollection.from_sources(jc_files, verbose=0)
    out = col.sort(key=lambda jf: -(getattr(jf, "lat", 0) or 0), reverse=False)
    assert isinstance(out, JCollection) and len(out) == len(col)


def test_sort_by_n_freq_path_lat_lon_keys(jc_files):
    col = JCollection.from_sources(jc_files, verbose=0)
    for key in ("n_freq", "path", "lat", "lon"):
        out = col.sort(key=key)
        assert len(out) == len(col)


def test_sort_by_n_freq_when_freq_has_no_size_attr():
    class _FakeJF:
        site = "S1"
        freq = None

    col = JCollection()
    col.add(_FakeJF())
    out = col.sort(key="n_freq")
    assert len(out) == 1


def test_sort_by_unknown_key_falls_back_to_getattr(j_single_file: Path):
    # A single-item collection avoids comparing two None keys (which
    # Python's sort cannot order) while still exercising the getattr
    # fallback branch.
    col = JCollection.from_sources([j_single_file], verbose=0)
    out = col.sort(key="nonexistent_attr")
    assert len(out) == len(col)


def test_paths_property(jc_files):
    col = JCollection.from_sources(jc_files, verbose=0)
    paths = col.paths
    assert isinstance(paths, list) and len(paths) == len(col)
    assert all(isinstance(p, str) for p in paths)


def test_nf_stats_empty_and_nonempty(jc_files):
    empty = JCollection()
    stats = empty.nf_stats()
    assert stats == {"min": 0, "max": 0, "mean": 0.0}

    col = JCollection.from_sources(jc_files, verbose=0)
    stats2 = col.nf_stats()
    assert stats2["min"] <= stats2["mean"] <= stats2["max"]


# ─────────────────────────────────────────────────────────────────────────
# JCollection.merge
# ─────────────────────────────────────────────────────────────────────────


def test_from_sources_logs_when_errors_present(tmp_path: Path):
    bad = tmp_path / "garbage.j"
    bad.write_text("not a real J file\n", encoding="utf-8")
    col = JCollection.from_sources([bad], verbose=0)
    assert isinstance(col, JCollection)
    assert len(col) == 0


def test_merge_default_replace(j_single_file: Path, jc_files):
    a = JCollection.from_sources([j_single_file], verbose=0)
    b = JCollection.from_sources(jc_files, verbose=0)
    merged = a.merge(b)
    assert len(merged) == len(a) + len(b)


def test_merge_keep_skips_duplicates(j_single_file: Path):
    a = JCollection.from_sources([j_single_file], verbose=0)
    b = JCollection.from_sources([j_single_file], verbose=0)
    merged = a.merge(b, on_dup="keep")
    assert len(merged) == len(a)


def test_merge_invalid_on_dup_raises(j_single_file: Path):
    a = JCollection.from_sources([j_single_file], verbose=0)
    b = JCollection.from_sources([j_single_file], verbose=0)
    with pytest.raises(ValueError):
        a.merge(b, on_dup="bogus")


# ─────────────────────────────────────────────────────────────────────────
# JCollection._resolve: index / iterate / path fallbacks
# ─────────────────────────────────────────────────────────────────────────


def test_resolve_by_index_case_insensitive(j_single_file: Path):
    col = JCollection.from_sources([j_single_file], verbose=0)
    site = list(col)[0].site
    jf = col._resolve(site.lower())
    assert jf.site == site


def test_resolve_falls_back_to_iteration_when_index_stale(j_single_file: Path):
    col = JCollection.from_sources([j_single_file], verbose=0)
    site = list(col)[0].site
    col._index.clear()  # force fallback to full iteration
    jf = col._resolve(site)
    assert jf.site == site


def test_resolve_falls_back_to_path_matching(tmp_path: Path):
    class _FakeJF:
        site = None
        path = tmp_path / "weird_name.j"

    col = JCollection()
    fake = _FakeJF()
    col.add(fake)
    col._index.clear()  # force fallback past the index and site loops
    resolved = col._resolve("weird_name")
    assert resolved is fake


def test_resolve_raises_keyerror_when_not_found(j_single_file: Path):
    col = JCollection.from_sources([j_single_file], verbose=0)
    with pytest.raises(KeyError):
        col._resolve("NO_SUCH_SITE")


def test_resolve_path_fallback_skips_items_with_no_path():
    class _NoPathJF:
        site = None
        path = None

    col = JCollection()
    col.add(_NoPathJF())
    col._index.clear()
    with pytest.raises(KeyError):
        col._resolve("anything")


# ─────────────────────────────────────────────────────────────────────────
# JCollection.get: every 'what' branch
# ─────────────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def loaded_collection(jc_files):
    return JCollection.from_sources(jc_files, verbose=0)


@pytest.fixture
def a_site(loaded_collection):
    return list(loaded_collection)[0].site


def test_get_unknown_site_returns_default(loaded_collection):
    assert loaded_collection.get("NOPE", "freq", default="DEF") == "DEF"


def test_get_freq(loaded_collection, a_site):
    out = loaded_collection.get(a_site, "freq")
    assert out is not None


def test_get_z_and_components(loaded_collection, a_site):
    z = loaded_collection.get(a_site, "z")
    assert z is not None
    for comp in ("zxx", "zxy", "zyx", "zyy"):
        arr = loaded_collection.get(a_site, comp)
        assert arr is not None and len(arr) == len(z)


def test_get_z_component_default_when_no_z(loaded_collection, a_site):
    jf = loaded_collection._resolve(a_site)
    real_z = jf.Z
    jf.Z = None
    try:
        out = loaded_collection.get(a_site, "zxy", default="NOZ")
        assert out == "NOZ"
    finally:
        jf.Z = real_z


def test_get_tipper_variants():
    class _FakeTip:
        tipper = np.zeros((3, 1, 2), complex)

    class _FakeJF:
        site = "TIPSITE"
        Tip = _FakeTip()

    col = JCollection()
    col.add(_FakeJF())
    tip = col.get("TIPSITE", "tip")
    tx = col.get("TIPSITE", "tx")
    ty = col.get("TIPSITE", "ty")
    assert tip.shape == (3, 2)
    assert tx.shape == (3,) and ty.shape == (3,)


def test_get_tipper_default_when_tipper_attr_none():
    class _FakeTip:
        tipper = None

    class _FakeJF:
        site = "TIPSITE2"
        Tip = _FakeTip()

    col = JCollection()
    col.add(_FakeJF())
    assert col.get("TIPSITE2", "tip", default="NOARR") == "NOARR"


def test_get_tipper_default_when_missing(loaded_collection, a_site):
    jf = loaded_collection._resolve(a_site)
    real_tip = jf.Tip
    jf.Tip = None
    try:
        assert loaded_collection.get(a_site, "tip", default="NOTIP") == "NOTIP"
    finally:
        jf.Tip = real_tip


def test_get_resphase_variants(loaded_collection, a_site):
    jf = loaded_collection._resolve(a_site)
    if jf.Res is None:
        pytest.skip("no resphase data on this fixture site")
    for w in ("rxy", "phixy"):
        loaded_collection.get(a_site, w)  # just exercise the branch


def test_get_resphase_default_when_missing(loaded_collection, a_site):
    jf = loaded_collection._resolve(a_site)
    real_res = jf.Res
    jf.Res = None
    try:
        assert loaded_collection.get(a_site, "rxy", default="NORES") == "NORES"
    finally:
        jf.Res = real_res


def test_get_site_metadata_fields(loaded_collection, a_site):
    assert loaded_collection.get(a_site, "station") == a_site
    assert loaded_collection.get(a_site, "site") == a_site
    loaded_collection.get(a_site, "name")
    loaded_collection.get(a_site, "lat")
    loaded_collection.get(a_site, "lon")
    loaded_collection.get(a_site, "elev")
    loaded_collection.get(a_site, "az")


def test_get_path_and_filename(loaded_collection, a_site):
    p = loaded_collection.get(a_site, "path")
    fn = loaded_collection.get(a_site, "filename")
    assert isinstance(p, str) and isinstance(fn, str)


def test_get_path_default_when_no_path(loaded_collection, a_site):
    jf = loaded_collection._resolve(a_site)
    real_path = jf.path
    jf.path = None
    try:
        assert loaded_collection.get(a_site, "path", default="NOPATH") == "NOPATH"
    finally:
        jf.path = real_path


def test_get_unknown_what_returns_default(loaded_collection, a_site):
    assert loaded_collection.get(a_site, "bogus_field", default="X") == "X"


# ─────────────────────────────────────────────────────────────────────────
# JCollection.set
# ─────────────────────────────────────────────────────────────────────────


def test_set_replaces_with_jfile(j_single_file: Path):
    col = JCollection.from_sources([j_single_file], verbose=0)
    site = list(col)[0].site
    other = JFile.from_file(j_single_file)
    out = col.set(site, jfile=other)
    assert out is other


def test_set_update_station_lat_lon_elev_az(j_single_file: Path):
    col = JCollection.from_sources([j_single_file], verbose=0)
    site = list(col)[0].site
    jf = col.set(
        site,
        update={
            "lat": 12.5,
            "lon": -33.1,
            "elev": 100.0,
            "az": 15.0,
        },
    )
    assert jf.heads.info.items["LATITUDE"] == "12.5"
    assert jf.heads.info.items["LONGITUDE"] == "-33.1"
    assert jf.heads.info.items["ELEVATION"] == "100.0"
    assert jf.heads.info.items["AZIMUTH"] == "15.0"


def test_set_update_station_rename(j_single_file: Path):
    col = JCollection.from_sources([j_single_file], verbose=0)
    site = list(col)[0].site
    jf = col.set(site, update={"station": "RENAMED"})
    assert jf.station == "RENAMED" or jf.heads.head.station == "RENAMED"


def test_set_jfile_replace_falls_back_to_index_when_add_raises(
    j_single_file: Path, monkeypatch
):
    col = JCollection.from_sources([j_single_file], verbose=0)
    site = list(col)[0].site
    other = JFile.from_file(j_single_file)

    def _boom(_jf):
        raise RuntimeError("boom")

    monkeypatch.setattr(col, "add", _boom)
    out = col.set(site, jfile=other)
    assert out is other
    assert col._index[site] is other


def test_set_station_rename_falls_back_to_head_exception_swallowed(
    j_single_file: Path,
):
    col = JCollection.from_sources([j_single_file], verbose=0)
    site = list(col)[0].site
    jf = col._resolve(site)
    jf.heads.head = None  # both jf.station and heads.head.station now fail
    out = col.set(site, update={"station": "X"})
    assert out is jf  # no crash; exception swallowed


def test_set_info_updates_are_noop_when_info_missing(j_single_file: Path):
    col = JCollection.from_sources([j_single_file], verbose=0)
    site = list(col)[0].site
    jf = col._resolve(site)
    jf.heads.info = None
    out = col.set(site, update={"lat": 1.0})
    assert out is jf


def test_set_info_item_assignment_exception_swallowed(j_single_file: Path):
    col = JCollection.from_sources([j_single_file], verbose=0)
    site = list(col)[0].site
    jf = col._resolve(site)
    jf.heads.info.items = None  # subscript assignment now raises
    out = col.set(site, update={"lat": 1.0})
    assert out is jf


def test_set_resphase_assignment_exception_swallowed(j_single_file: Path):
    class _ReadOnlyRes:
        @property
        def resistivity(self):
            return None

    col = JCollection.from_sources([j_single_file], verbose=0)
    site = list(col)[0].site
    jf = col._resolve(site)
    jf.Res = _ReadOnlyRes()  # no setter -> assignment raises, swallowed
    out = col.set(site, update={"resphase": [1.0]})
    assert out is jf


def test_set_update_freq_z_tip_resphase(j_single_file: Path):
    col = JCollection.from_sources([j_single_file], verbose=0)
    site = list(col)[0].site
    jf = col._resolve(site)
    if jf.Z is None:
        pytest.skip("fixture site has no Z object")
    n = jf.Z.z.shape[0]
    new_freq = np.linspace(1.0, 2.0, n)
    new_z = np.zeros_like(jf.Z.z)
    updated = col.set(site, update={"freq": new_freq, "z": new_z})
    assert np.allclose(updated.Z._freq, new_freq)
    assert np.allclose(updated.Z._z, new_z)

    if jf.Tip is not None:
        new_tip = np.zeros_like(np.asarray(jf.Tip.tipper))
        updated2 = col.set(site, update={"tip": new_tip})
        assert np.allclose(updated2.Tip._tipper, new_tip)

    if jf.Res is not None:
        new_rp = np.ones_like(np.asarray(jf.Res.resistivity))
        col.set(site, update={"resphase": new_rp})


# ─────────────────────────────────────────────────────────────────────────
# JCollection.adjust
# ─────────────────────────────────────────────────────────────────────────


def test_adjust_set_lat_lon_elev_and_rename(j_single_file: Path):
    col = JCollection.from_sources([j_single_file], verbose=0)
    site = list(col)[0].site
    jf = col.adjust(site, lat=10.0, lon=20.0, elev=5.0, rename="NEWNAME")
    assert jf.heads.info.items["LATITUDE"] == "10.0"
    assert jf.heads.info.items["LONGITUDE"] == "20.0"
    assert jf.heads.info.items["ELEVATION"] == "5.0"


def test_adjust_dlat_dlon_shifts_existing_value(j_single_file: Path):
    col = JCollection.from_sources([j_single_file], verbose=0)
    site = list(col)[0].site
    jf = col._resolve(site)
    jf.heads.info.items["LATITUDE"] = "10.0"
    jf.heads.info.items["LONGITUDE"] = "20.0"
    jf.heads.info._site_cache = None
    out = col.adjust(site, dlat=1.0, dlon=-2.0)
    assert out.heads.info.items["LATITUDE"] == "11.0"
    assert out.heads.info.items["LONGITUDE"] == "18.0"


def test_adjust_rename_falls_back_to_head_exception_swallowed(
    j_single_file: Path,
):
    col = JCollection.from_sources([j_single_file], verbose=0)
    site = list(col)[0].site
    jf = col._resolve(site)
    jf.heads.head = None  # both jf.station and heads.head.station now fail
    out = col.adjust(site, rename="X")
    assert out is jf  # no crash; exception swallowed


def test_adjust_getf_returns_default_on_bad_existing_value(j_single_file: Path):
    col = JCollection.from_sources([j_single_file], verbose=0)
    site = list(col)[0].site
    jf = col._resolve(site)
    jf.heads.info.items["LATITUDE"] = "not-a-number"
    jf.heads.info._site_cache = None
    out = col.adjust(site, dlat=1.0)
    # _getf falls back to 0.0 default, so result is "0.0 + 1.0"
    assert out.heads.info.items["LATITUDE"] == "1.0"


def test_adjust_raises_when_no_info(j_single_file: Path):
    col = JCollection.from_sources([j_single_file], verbose=0)
    site = list(col)[0].site
    jf = col._resolve(site)
    jf.heads.info = None
    with pytest.raises(ValueError):
        col.adjust(site, lat=1.0)


# ─────────────────────────────────────────────────────────────────────────
# _site_of / _lat_of / _lon_of / sites / latitude / longitude
# ─────────────────────────────────────────────────────────────────────────


def test_site_of_prefers_direct_attrs():
    class _F:
        site = "S1"

    assert JCollection._site_of(_F()) == "S1"


def test_site_of_falls_back_to_nested_head():
    class _Head:
        station = "S2"

    class _Heads:
        head = _Head()

    class _F:
        heads = _Heads()

    assert JCollection._site_of(_F()) == "S2"


def test_site_of_returns_none_when_nothing_found():
    class _F:
        pass

    assert JCollection._site_of(_F()) is None


def test_lat_of_and_lon_of_fallback_chain():
    class _Heads:
        latitude = 5.0
        longitude = 6.0

    class _F:
        heads = _Heads()

    assert JCollection._lat_of(_F()) == 5.0
    assert JCollection._lon_of(_F()) == 6.0


def test_lat_of_and_lon_of_none_when_not_numeric():
    class _F:
        lat = "not-a-number"
        lon = "not-a-number"

    assert JCollection._lat_of(_F()) is None
    assert JCollection._lon_of(_F()) is None


def test_sites_latitude_longitude_properties(jc_files):
    col = JCollection.from_sources(jc_files, verbose=0)
    assert len(col.sites) == len(col)
    assert len(col.latitude) == len(col)
    assert len(col.longitude) == len(col)


# ─────────────────────────────────────────────────────────────────────────
# export()
# ─────────────────────────────────────────────────────────────────────────


def test_export_writes_files_and_summary(tmp_path: Path, j_single_file: Path):
    col = JCollection.from_sources([j_single_file], verbose=0)
    out_dir = tmp_path / "export_out"
    result = col.export(
        out_dir,
        export_summary=True,
        datatype="R",
        overwrite=True,
    )
    assert "successful" in result and "failed" in result
    assert len(result["successful"]) >= 1
    assert (out_dir / "summary.csv").exists()


def test_export_falls_back_when_tqdm_unavailable(
    tmp_path: Path, j_single_file: Path, monkeypatch
):
    import sys

    monkeypatch.setitem(sys.modules, "tqdm", None)
    col = JCollection.from_sources([j_single_file], verbose=0)
    result = col.export(tmp_path / "no_tqdm", datatype="R", overwrite=True)
    assert len(result["successful"]) >= 1


def test_export_summary_failure_recorded(tmp_path: Path, j_single_file: Path):
    col = JCollection.from_sources([j_single_file], verbose=0)
    col.summary = lambda: (_ for _ in ()).throw(RuntimeError("boom"))
    result = col.export(
        tmp_path / "summary_fail", export_summary=True, datatype="R"
    )
    assert any(name == "summary.csv" for name, _ in result["failed"])


def test_export_records_failures_gracefully(tmp_path: Path, j_single_file: Path):
    col = JCollection.from_sources([j_single_file], verbose=0)
    jf = list(col)[0]
    jf.write = lambda **kwargs: (_ for _ in ()).throw(RuntimeError("boom"))
    out_dir = tmp_path / "export_fail"
    result = col.export(out_dir)
    assert result["failed"]
    assert result["failed"][0][0] == jf.site


# ─────────────────────────────────────────────────────────────────────────
# fetch()
# ─────────────────────────────────────────────────────────────────────────


def test_fetch_by_site_case_insensitive(loaded_collection, a_site):
    out = loaded_collection.fetch(site=a_site.lower(), first=True)
    assert out is not None and out.site == a_site


def test_fetch_by_site_no_match_returns_none_first(loaded_collection):
    assert loaded_collection.fetch(site="NOPE", first=True) is None


def test_fetch_by_site_no_match_returns_empty_list(loaded_collection):
    assert loaded_collection.fetch(site="NOPE", first=False) == []


def test_fetch_by_lat_lon_tolerance(loaded_collection, a_site):
    jf = loaded_collection._resolve(a_site)
    if jf.lat is None or jf.lon is None:
        pytest.skip("fixture site has no lat/lon")
    out = loaded_collection.fetch(lat=jf.lat, lon=jf.lon, tol=1e-6, first=True)
    assert out is not None and out.site == a_site


def test_fetch_by_lat_outside_tolerance_excludes(loaded_collection, a_site):
    jf = loaded_collection._resolve(a_site)
    if jf.lat is None:
        pytest.skip("fixture site has no lat")
    out = loaded_collection.fetch(lat=jf.lat + 100.0, tol=0.001)
    assert all(o.site != a_site for o in out)


def test_fetch_by_lon_outside_tolerance_excludes(loaded_collection, a_site):
    jf = loaded_collection._resolve(a_site)
    if jf.lon is None:
        pytest.skip("fixture site has no lon")
    out = loaded_collection.fetch(lon=jf.lon + 100.0, tol=0.001)
    assert all(o.site != a_site for o in out)


def test_fetch_by_arbitrary_kwarg_on_jfile_or_nested(loaded_collection, a_site):
    jf = loaded_collection._resolve(a_site)
    out = loaded_collection.fetch(site=a_site, azimuth=jf.azimuth)
    assert any(o.site == a_site for o in out) or jf.azimuth is None


def test_fetch_kwarg_mismatch_excludes(loaded_collection, a_site):
    out = loaded_collection.fetch(site=a_site, azimuth=-99999.0)
    assert all(o.site != a_site for o in out)


def test_fetch_kwarg_falls_back_to_nested_info_attribute(loaded_collection, a_site):
    # "latitude" is not a direct JFile or Head attribute (only "lat" is),
    # so this exercises the fallback to jf.heads.info.latitude.
    jf = loaded_collection._resolve(a_site)
    if jf.lat is None:
        pytest.skip("fixture site has no lat")
    out = loaded_collection.fetch(site=a_site, latitude=jf.lat)
    assert any(o.site == a_site for o in out)


def test_fetch_kwarg_string_comparison_case_insensitive(loaded_collection, a_site):
    out = loaded_collection.fetch(site=a_site, station=a_site.lower())
    assert any(o.site == a_site for o in out)


def test_fetch_kwarg_string_comparison_mismatch_excludes(loaded_collection, a_site):
    out = loaded_collection.fetch(site=a_site, station="DEFINITELY_WRONG")
    assert all(o.site != a_site for o in out)


# ─────────────────────────────────────────────────────────────────────────
# _summary_stats
# ─────────────────────────────────────────────────────────────────────────


def test_summary_stats_empty_collection_message():
    col = JCollection()
    assert "No statistics" in col._summary_stats([])


def test_summary_stats_nonempty_collection(jc_files):
    col = JCollection.from_sources(jc_files, verbose=0)
    text = col._summary_stats(col.summary())
    assert "Statistical Summary" in text
    assert "Total Sites" in text
