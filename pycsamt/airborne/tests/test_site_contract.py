from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.airborne import (
    AirborneEMRecord,
    AirborneSite,
    AirborneSites,
    NavigationTrack,
    ensure_asites,
)
from pycsamt.airborne.afmag import (
    build_airmt_line,
    build_original_afmag_line,
)
from pycsamt.airborne.mobilemt import (
    build_mobilemt_dataset,
    build_mobilemt_line,
)
from pycsamt.airborne.ztem import (
    ZTEMReferenceStation,
    build_ztem_dataset,
    build_ztem_line,
)
from pycsamt.metadata import LocationMeta, SiteMeta

_DATA_ROOT = Path(__file__).resolve().parents[3] / "data"
_GOLD_SPRINGS = _DATA_ROOT / "ZTEM" / "gold_springs_nv"
_FLAMMEFJELD = _DATA_ROOT / "mobileMT" / "flammefjeld_greenland"


def _tipper(nf: int = 4) -> np.ndarray:
    base = np.arange(nf * 2, dtype=float).reshape(nf, 2)
    return base + 1j * (base + 0.25)


def _admittance(n_samples: int, nf: int = 3) -> np.ndarray:
    data = np.zeros((n_samples, nf, 3, 2), dtype=complex)
    for i in range(n_samples):
        data[i, :, 0, 1] = 1.0 + 0.1j * (i + 1)
        data[i, :, 1, 0] = -1.0 - 0.1j * (i + 1)
    return data


def _ztem_line(line_id: str = "L001", n: int = 3):
    nav = NavigationTrack(
        sample_ids=tuple(f"P{i:03d}" for i in range(n)),
        latitude=[5.0 + 0.001 * i for i in range(n)],
        longitude=[-3.0 - 0.001 * i for i in range(n)],
    )
    data = np.stack([_tipper(4) + index for index in range(n)])
    ref = ZTEMReferenceStation(
        station_id="BASE01",
        site=SiteMeta(
            site_id="BASE01",
            location=LocationMeta(latitude=5.2, longitude=-3.1),
        ),
    )
    return build_ztem_line(
        line_id,
        nav,
        data,
        frequency=[30.0, 45.0, 90.0, 180.0],
        reference_station=ref,
    )


def _mobilemt_line(line_id: str = "M001", n: int = 3):
    nav = NavigationTrack(
        sample_ids=tuple(f"Q{i:03d}" for i in range(n)),
        latitude=[10.0 + 0.001 * i for i in range(n)],
        longitude=[20.0 + 0.001 * i for i in range(n)],
    )
    return build_mobilemt_line(
        line_id,
        nav,
        _admittance(n),
        frequency=[25.0, 100.0, 1000.0],
    )


def _interstation_tensor(n_samples: int, nf: int = 3) -> np.ndarray:
    data = np.zeros((n_samples, nf, 3, 2), dtype=complex)
    for i in range(n_samples):
        data[i, :, 0, 0] = 1.0
        data[i, :, 1, 1] = 1.0
        data[i, :, 2, 0] = 0.1 + 0.02j * (i + 1)
        data[i, :, 2, 1] = -0.05 - 0.01j * (i + 1)
    return data


def _airmt_line(line_id: str = "AM001", n: int = 3):
    nav = NavigationTrack(
        sample_ids=tuple(f"R{i:03d}" for i in range(n)),
        latitude=[15.0 + 0.001 * i for i in range(n)],
        longitude=[25.0 + 0.001 * i for i in range(n)],
    )
    return build_airmt_line(
        line_id,
        nav,
        _interstation_tensor(n),
        frequency=[20.0, 100.0, 500.0],
    )


def _original_afmag_line(line_id: str = "OA001", n: int = 3):
    nav = NavigationTrack(
        sample_ids=tuple(f"W{i:03d}" for i in range(n)),
        latitude=[20.0 + 0.001 * i for i in range(n)],
        longitude=[30.0 + 0.001 * i for i in range(n)],
    )
    tilt = np.array([[5.0 + i, -3.0 - i] for i in range(n)])
    return build_original_afmag_line(
        line_id, nav, tilt, frequency=[150.0, 510.0],
    )


# ─────────────────────────────────────────────────────────────────────────
# AirborneSite
# ─────────────────────────────────────────────────────────────────────────


def test_airborne_site_reads_tipper_directly_no_edi():
    line = _ztem_line()
    record = line.get_record("P000")
    site = AirborneSite(record, line_id="L001", technology="ztem")
    assert site.tipper.shape == (4, 1, 2)
    assert np.allclose(site.freq, [30.0, 45.0, 90.0, 180.0])
    assert site.z is None
    assert site.admittance is None
    assert site.technology == "ztem"
    assert site.line_id == "L001"
    assert site.sample_id == "P000"


def test_airborne_site_exposes_admittance_for_mobilemt():
    line = _mobilemt_line()
    record = line.get_record("Q000")
    site = AirborneSite(record, technology="mobilemt")
    assert site.tipper is None
    assert site.admittance.shape == (3, 3, 2)
    assert site.has_component("admittance")
    assert not site.has_component("tipper")


def test_airborne_site_exposes_interstation_tensor_for_airmt():
    line = _airmt_line()
    record = line.get_record("R000")
    site = AirborneSite(record, technology="afmag_airmt")
    assert site.tipper is None
    assert site.admittance is None
    assert site.afmag_tilt_deg is None
    assert site.interstation_tensor.shape == (3, 3, 2)
    assert site.has_component("interstation_tensor")
    assert not site.has_component("tipper")
    # build_airmt_line defaults include_amplification_parameter=True
    assert site.afmag_amplification_parameter is not None
    assert site.has_component("ap")


def test_airborne_site_exposes_scalar_tilt_for_original_afmag():
    line = _original_afmag_line()
    record = line.get_record("W000")
    site = AirborneSite(record, technology="afmag_original")
    assert site.tipper is None
    assert site.interstation_tensor is None
    assert site.afmag_amplification_parameter is None
    tilt = site.afmag_tilt_deg
    assert tilt.shape == (2,)
    assert np.allclose(tilt, [5.0, -3.0])
    assert site.has_component("tilt")
    assert not site.has_component("interstation_tensor")


def test_airborne_site_rejects_non_record():
    with pytest.raises(TypeError):
        AirborneSite("not-a-record")


def test_airborne_site_coords_prefer_explicit_override():
    line = _ztem_line()
    record = line.get_record("P000")
    site = AirborneSite(record, coords=(1.0, 2.0, 3.0))
    assert site.coords == (1.0, 2.0, 3.0)


def test_airborne_site_coords_nan_when_unknown():
    record = AirborneEMRecord(sample_id="S1")
    site = AirborneSite(record)
    lat, lon, elev = site.coords
    assert np.isnan(lat) and np.isnan(lon) and np.isnan(elev)


def test_airborne_site_name_falls_back_to_sample_id():
    record = AirborneEMRecord(sample_id="S1")
    site = AirborneSite(record)
    assert site.name == "S1"
    assert site.station == "S1"


def test_airborne_site_summary_and_repr():
    line = _ztem_line()
    site = AirborneSite(
        line.get_record("P000"), line_id="L001", technology="ztem",
    )
    s = site.summary()
    assert s["nfreq"] == 4
    assert s["tipper"] is True
    assert s["admittance"] is False
    assert "AirborneSite(" in repr(site)


def test_airborne_site_to_dataframe_tipper_and_admittance():
    ztem_site = AirborneSite(_ztem_line().get_record("P000"))
    df = ztem_site.to_dataframe("tipper")
    assert list(df.columns) == ["Tx", "Ty"]
    assert len(df) == 4

    mmt_site = AirborneSite(_mobilemt_line().get_record("Q000"))
    df2 = mmt_site.to_dataframe("admittance")
    assert list(df2.columns) == ["Yxx", "Yxy", "Yyx", "Yyy", "Yhzx", "Yhzy"]

    with pytest.raises(ValueError):
        ztem_site.to_dataframe("bogus")


def test_airborne_site_to_xml_roundtrip(tmp_path):
    site = AirborneSite(_ztem_line().get_record("P000"), technology="ztem")
    xml_str = site.to_xml()
    assert "<T " in xml_str

    target = tmp_path / "P000.xml"
    site.to_xml(target)
    reloaded = AirborneSite.from_xml(target)
    assert np.allclose(reloaded.tipper, site.tipper)


def test_airborne_site_to_xml_without_emtf_raises():
    record = AirborneEMRecord(sample_id="S1")
    site = AirborneSite(record)
    with pytest.raises(ValueError):
        site.to_xml()


# ─────────────────────────────────────────────────────────────────────────
# AirborneSites
# ─────────────────────────────────────────────────────────────────────────


def test_airborne_sites_from_line_preserves_navigation_order():
    line = _ztem_line(n=3)
    asites = AirborneSites.from_line(line)
    assert len(asites) == 3
    assert [s.sample_id for s in asites] == ["P000", "P001", "P002"]
    assert all(s.line_id == "L001" for s in asites)
    assert all(s.technology == "ztem" for s in asites)
    lat0, lon0, _ = asites[0].coords
    assert np.isclose(lat0, 5.0) and np.isclose(lon0, -3.0)


def test_airborne_sites_from_dataset_flattens_all_lines():
    line1 = _ztem_line("L001", n=2)
    line2 = _ztem_line("L002", n=2)
    # rebuild line2 with distinct sample ids to avoid collisions
    nav2 = NavigationTrack(
        sample_ids=("R000", "R001"),
        latitude=[6.0, 6.001],
        longitude=[-4.0, -4.001],
    )
    from pycsamt.airborne.ztem import build_ztem_line as _build

    line2 = _build(
        "L002", nav2, np.stack([_tipper(4), _tipper(4) + 1]),
        frequency=[30.0, 45.0, 90.0, 180.0],
    )
    dataset = build_ztem_dataset("SURVEY", [line1, line2])
    asites = AirborneSites.from_dataset(dataset)
    assert len(asites) == 4
    assert set(asites.line_ids) == {"L001", "L002"}
    assert asites.technologies == ("ztem",)


def test_airborne_sites_container_protocol():
    asites = AirborneSites.from_line(_ztem_line(n=3))
    assert len(asites) == 3
    assert asites[0].sample_id == "P000"
    assert asites["p001"].sample_id == "P001"
    assert asites.by_index(2).sample_id == "P002"
    assert asites.get("missing") is None
    with pytest.raises(KeyError):
        asites["missing"]


def test_airborne_sites_select_and_map():
    asites = AirborneSites.from_line(_ztem_line(n=3))
    subset = asites.select(names=["P001"])
    assert len(subset) == 1
    subset2 = asites.select(predicate=lambda s: s.has_component("tipper"))
    assert len(subset2) == 3
    names = asites.map(lambda s: s.name)
    assert names == ["P000", "P001", "P002"]


def test_airborne_sites_closest():
    asites = AirborneSites.from_line(_ztem_line(n=3))
    near = asites.closest(5.0, -3.0)
    assert near.sample_id == "P000"
    far = asites.closest(-80.0, 170.0, tol=1000.0)
    assert far is None


def test_airborne_sites_write_xml_roundtrip(tmp_path):
    asites = AirborneSites.from_line(_ztem_line(n=3))
    paths = asites.write_xml(tmp_path)
    assert len(paths) == 3
    assert all(p.exists() for p in paths)
    reloaded = AirborneSites.from_xml_dir(tmp_path)
    assert len(reloaded) == 3


def test_airborne_sites_init_accepts_mixed_items(tmp_path):
    line = _ztem_line(n=2)
    record0 = line.get_record("P000")
    site1 = AirborneSite(line.get_record("P001"), line_id="L001")
    asites = AirborneSites([record0, site1])
    assert len(asites) == 2
    assert asites[0].sample_id == "P000"


def test_airborne_sites_init_rejects_bad_item():
    with pytest.raises(TypeError):
        AirborneSites([123])


# ─────────────────────────────────────────────────────────────────────────
# ensure_asites
# ─────────────────────────────────────────────────────────────────────────


def test_ensure_asites_passthrough():
    asites = AirborneSites.from_line(_ztem_line())
    assert ensure_asites(asites) is asites


def test_ensure_asites_from_dataset_and_line():
    line = _ztem_line(n=2)
    dataset = build_ztem_dataset("SURVEY", [line])
    a1 = ensure_asites(dataset)
    a2 = ensure_asites(line)
    assert len(a1) == 2
    assert len(a2) == 2


def test_ensure_asites_from_directory(tmp_path):
    AirborneSites.from_line(_ztem_line(n=3)).write_xml(tmp_path)
    asites = ensure_asites(tmp_path)
    assert len(asites) == 3


def test_ensure_asites_none_raises():
    with pytest.raises(ValueError):
        ensure_asites(None)


def test_ensure_asites_strict_empty_raises(tmp_path):
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    with pytest.raises(ValueError):
        ensure_asites(empty_dir, strict=True)


def test_ensure_asites_non_strict_empty_returns_empty(tmp_path):
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    asites = ensure_asites(empty_dir, strict=False)
    assert len(asites) == 0


def test_ensure_asites_on_dup_policies():
    line = _ztem_line(n=2)
    site_a = AirborneSite(line.get_record("P000"))
    site_b = AirborneSite(line.get_record("P000"))  # same name
    with pytest.raises(ValueError):
        ensure_asites([site_a, site_b], on_dup="raise")
    kept_first = ensure_asites([site_a, site_b], on_dup="keep_first")
    assert len(kept_first) == 1
    replaced = ensure_asites([site_a, site_b], on_dup="replace")
    assert len(replaced) == 1


def test_ensure_asites_invalid_on_dup_raises():
    line = _ztem_line(n=1)
    with pytest.raises(ValueError):
        ensure_asites(line, on_dup="bogus")


def test_ensure_asites_mixed_list_of_lines_and_paths(tmp_path):
    AirborneSites.from_line(_ztem_line(n=2)).write_xml(tmp_path)
    line2 = _mobilemt_line(n=2)
    combined = ensure_asites([tmp_path, line2])
    assert len(combined) == 4
    assert set(combined.technologies) == {"ztem", "mobilemt"}


# ─────────────────────────────────────────────────────────────────────────
# Real synthetic data (skipped if not present locally)
# ─────────────────────────────────────────────────────────────────────────


@pytest.mark.skipif(
    not _GOLD_SPRINGS.exists(), reason="synthetic ZTEM sample data not found"
)
def test_real_synthetic_ztem_xml_loads():
    # gold_springs_nv is a 7-line, 105-station block survey.
    asites = ensure_asites(_GOLD_SPRINGS)
    assert len(asites) == 105
    assert asites.technologies == ("ztem",)
    s0 = asites[0]
    assert s0.tipper is not None
    assert s0.freq.size == 6
    lat, lon, elev = s0.coords
    assert np.isfinite(lat) and np.isfinite(lon) and np.isfinite(elev)


@pytest.mark.skipif(
    not _FLAMMEFJELD.exists(),
    reason="synthetic MobileMT sample data not found",
)
def test_real_synthetic_mobilemt_xml_loads():
    asites = ensure_asites(_FLAMMEFJELD)
    assert len(asites) == 12
    assert asites.technologies == ("mobilemt",)
    s0 = asites[0]
    assert s0.admittance is not None
    assert s0.admittance.shape[1:] == (3, 2)


# ─────────────────────────────────────────────────────────────────────────
# _navigation_coords elevation fallback branches
# ─────────────────────────────────────────────────────────────────────────


def test_navigation_coords_uses_platform_elevation_when_present():
    from pycsamt.airborne.site import _navigation_coords

    nav = NavigationTrack(
        sample_ids=("A",), latitude=(1.0,), longitude=(2.0,),
        platform_elevation=(150.0,),
    )
    _lat, _lon, elev = _navigation_coords(nav, 0)
    assert elev == 150.0


def test_navigation_coords_falls_back_to_terrain_elevation():
    from pycsamt.airborne.site import _navigation_coords

    nav = NavigationTrack(
        sample_ids=("A",), latitude=(1.0,), longitude=(2.0,),
        terrain_elevation=(90.0,),
    )
    _lat, _lon, elev = _navigation_coords(nav, 0)
    assert elev == 90.0


# ─────────────────────────────────────────────────────────────────────────
# AirborneSite.from_xml identifier resolution + _stem_from_source
# ─────────────────────────────────────────────────────────────────────────


def test_from_xml_uses_station_when_no_site_id():
    doc = _ztem_line().get_record("P000").emtf
    doc.station = "STATION42"
    site = AirborneSite.from_xml(doc)
    assert site.sample_id == "STATION42"


def test_from_xml_falls_back_to_literal_site_when_source_not_path():
    doc = _ztem_line().get_record("P000").emtf
    doc.station = None
    site = AirborneSite.from_xml(doc)
    assert site.sample_id == "site"


# ─────────────────────────────────────────────────────────────────────────
# Simple accessor properties
# ─────────────────────────────────────────────────────────────────────────


def test_airborne_site_simple_accessors():
    record = _ztem_line().get_record("P000")
    site = AirborneSite(record)
    assert site.record is record
    assert site.emtf is record.emtf
    assert site.tf is record.emtf
    assert isinstance(site.quality, dict)
    assert isinstance(site.fields, dict)
    assert site.site_meta is record.emtf.site
    assert site.site_layout is record.emtf.site_layout
    assert site.provenance is record.emtf.provenance
    assert site.processing is record.emtf.processing
    assert site.copyright is record.emtf.copyright
    assert site.quality_meta is record.emtf.quality


def test_airborne_site_technology_none_when_unresolvable():
    record = AirborneEMRecord(sample_id="S1")
    site = AirborneSite(record)
    assert site.technology is None


def test_airborne_site_name_falls_back_to_doc_station():
    record = _ztem_line().get_record("P000")
    record.emtf.site = None
    record.emtf.station = "DOCSTATION"
    site = AirborneSite(record)
    assert site.name == "DOCSTATION"


def test_airborne_site_coords_nan_when_site_has_no_location():
    record = _ztem_line().get_record("P000")
    record.emtf.site = SiteMeta(site_id="X")
    site = AirborneSite(record)
    lat, lon, elev = site.coords
    assert np.isnan(lat) and np.isnan(lon) and np.isnan(elev)


def test_airborne_site_optional_accessors_none_without_emtf():
    record = AirborneEMRecord(sample_id="S1")
    site = AirborneSite(record)
    assert site.admittance is None
    assert site.interstation_tensor is None
    assert site.afmag_tilt_deg is None
    assert site.afmag_amplification_parameter is None


def test_airborne_site_has_component_default_falls_back_to_z():
    record = _ztem_line().get_record("P000")
    site = AirborneSite(record)
    assert site.has_component("z") is False
    assert site.has_component("something_unknown") is False


def test_airborne_site_to_dataframe_z_kind_empty_when_no_impedance():
    record = _ztem_line().get_record("P000")
    site = AirborneSite(record)
    df = site.to_dataframe("impedance")
    assert list(df.columns) == []


# ─────────────────────────────────────────────────────────────────────────
# AirborneSite.tipper setter
# ─────────────────────────────────────────────────────────────────────────


def test_tipper_setter_raises_without_emtf():
    record = AirborneEMRecord(sample_id="S1")
    site = AirborneSite(record)
    with pytest.raises(ValueError):
        site.tipper = np.zeros((2, 2), dtype=complex)


def test_tipper_setter_raises_without_tipper_tf():
    record = _mobilemt_line().get_record("Q000")
    site = AirborneSite(record)
    with pytest.raises(ValueError):
        site.tipper = np.zeros((2, 2), dtype=complex)


def test_tipper_setter_accepts_nf2_and_nf1_2_shapes():
    site = AirborneSite(_ztem_line().get_record("P000"))
    new_values = np.arange(8, dtype=complex).reshape(4, 2)
    site.tipper = new_values
    assert np.allclose(site.tipper[:, 0, :], new_values)

    site2 = AirborneSite(_ztem_line().get_record("P000"))
    canonical = np.arange(8, dtype=complex).reshape(4, 1, 2)
    site2.tipper = canonical
    assert np.allclose(site2.tipper, canonical)


# ─────────────────────────────────────────────────────────────────────────
# AirborneSites construction from a bare str/Path/EMTF, single-item wrapping
# ─────────────────────────────────────────────────────────────────────────


def test_airborne_sites_accepts_single_path_item(tmp_path):
    site = AirborneSite(_ztem_line().get_record("P000"))
    target = tmp_path / "P000.xml"
    site.to_xml(target)
    asites = AirborneSites(target)
    assert len(asites) == 1
    assert asites[0].sample_id == "P000"


def test_airborne_sites_accepts_single_emtf_item():
    doc = _ztem_line().get_record("P000").emtf
    asites = AirborneSites(doc)
    assert len(asites) == 1


def test_airborne_sites_accepts_single_non_iterable_unsupported_item():
    with pytest.raises(TypeError):
        AirborneSites(123)


# ─────────────────────────────────────────────────────────────────────────
# AirborneSites.from_xml_dir error handling
# ─────────────────────────────────────────────────────────────────────────


def test_from_xml_dir_skips_malformed_file_by_default(tmp_path):
    AirborneSites.from_line(_ztem_line(n=2)).write_xml(tmp_path)
    bad = tmp_path / "broken.xml"
    bad.write_text("not valid xml <<<", encoding="utf-8")
    asites = AirborneSites.from_xml_dir(tmp_path)
    assert len(asites) == 2


def test_from_xml_dir_strict_raises_on_malformed_file(tmp_path):
    AirborneSites.from_line(_ztem_line(n=1)).write_xml(tmp_path)
    bad = tmp_path / "broken.xml"
    bad.write_text("not valid xml <<<", encoding="utf-8")
    with pytest.raises(Exception):
        AirborneSites.from_xml_dir(tmp_path, strict=True)


# ─────────────────────────────────────────────────────────────────────────
# AirborneSites misc: to_emtf_list / select() / closest()
# ─────────────────────────────────────────────────────────────────────────


def test_airborne_sites_to_emtf_list():
    asites = AirborneSites.from_line(_ztem_line(n=2))
    docs = asites.to_emtf_list()
    assert len(docs) == 2
    assert all(doc is not None for doc in docs)


def test_airborne_sites_select_with_no_args_returns_copy():
    asites = AirborneSites.from_line(_ztem_line(n=2))
    copy = asites.select()
    assert len(copy) == 2
    assert copy is not asites


def test_airborne_sites_closest_skips_nan_coords_and_returns_none_when_empty():
    nav_site = AirborneSites.from_line(_ztem_line(n=1))[0]
    nan_site = AirborneSite(AirborneEMRecord(sample_id="NAN1"))
    combined = AirborneSites([nan_site, nav_site])
    nearest = combined.closest(5.0, -3.0)
    assert nearest.sample_id == "P000"

    only_nan = AirborneSites([AirborneEMRecord(sample_id="NAN2")])
    assert only_nan.closest(5.0, -3.0) is None


# ─────────────────────────────────────────────────────────────────────────
# ensure_asites: dataset-in-iterable, TypeError skip/raise, verbose warning
# ─────────────────────────────────────────────────────────────────────────


def test_ensure_asites_iterable_containing_a_dataset():
    line = _ztem_line(n=2)
    dataset = build_ztem_dataset("SURVEY", [line])
    combined = ensure_asites([dataset])
    assert len(combined) == 2


def test_ensure_asites_iterable_skips_unsupported_items_when_not_strict():
    line = _ztem_line(n=1)
    combined = ensure_asites([line, 123])
    assert len(combined) == 1


def test_ensure_asites_iterable_raises_on_unsupported_item_when_strict():
    with pytest.raises(TypeError):
        ensure_asites([123], strict=True)


def test_ensure_asites_verbose_warns_on_empty_result():
    with pytest.warns(RuntimeWarning):
        result = ensure_asites([123], verbose=1)
    assert len(result) == 0


def test_ensure_asites_empty_iterable_strict_raises():
    with pytest.raises(ValueError):
        ensure_asites([], strict=True)


def test_ensure_asites_passthrough_invalid_on_dup_raises():
    asites = AirborneSites.from_line(_ztem_line(n=1))
    with pytest.raises(ValueError):
        ensure_asites(asites, on_dup="bogus")


def test_ensure_asites_passthrough_empty_strict_raises():
    empty = AirborneSites([])
    with pytest.raises(ValueError):
        ensure_asites(empty, strict=True)


# ─────────────────────────────────────────────────────────────────────────
# from_line / from_dataset type checks and AirborneSites.__repr__
# ─────────────────────────────────────────────────────────────────────────


def test_from_line_rejects_non_line():
    with pytest.raises(TypeError):
        AirborneSites.from_line("not-a-line")


def test_from_dataset_rejects_non_dataset():
    with pytest.raises(TypeError):
        AirborneSites.from_dataset("not-a-dataset")


def test_airborne_sites_repr():
    asites = AirborneSites.from_line(_ztem_line(n=1))
    text = repr(asites)
    assert "AirborneSites(n=1" in text
    assert "ztem" in text

    empty = AirborneSites([])
    assert "?" in repr(empty)
