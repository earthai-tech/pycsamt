from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pycsamt.seg.edi import EDIFile
from pycsamt.site.base import Site, Sites
from pycsamt.site.metadata import (
    SiteMetadataEditor,
    rename_sites,
    update_metadata,
    update_metadata_all,
)

_PROJECT_ROOT = Path(__file__).parents[3]
_WILLY_ROOT = _PROJECT_ROOT / "data" / "AMT" / "WILLY_DATA"


def _real_line(name: str) -> list[Path]:
    paths = sorted((_WILLY_ROOT / name).glob("*.edi"))
    if not paths:
        pytest.skip(f"{name} integration data is not available")
    return paths


def _copy_edi(tmp_path: Path, source: Path, stem: str) -> EDIFile:
    path = tmp_path / f"{stem}.edi"
    path.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")
    edi = EDIFile(path)
    edi.name = stem
    edi.station = stem
    return edi


def _sites(tmp_path: Path, source: Path) -> Sites:
    return Sites(
        [
            _copy_edi(tmp_path, source, "A01"),
            _copy_edi(tmp_path, source, "A02"),
        ]
    )


def test_update_one_site_copy_and_metadata(simulated_edi: Path) -> None:
    original = Site(EDIFile(simulated_edi))
    old_name = original.name

    updated = update_metadata(
        original,
        {
            "name": "L01_001",
            "lat": 5.25,
            "lon": -3.75,
            "elev": 120.0,
            "head": {"project": "LINE_01"},
            "info": {"project": "LINE_01", "reviewed_by": "operator"},
        },
    )

    assert isinstance(updated, Site)
    assert updated is not original
    assert original.name == old_name
    assert updated.name == "L01_001"
    assert tuple(map(float, updated.coords)) == (5.25, -3.75, 120.0)
    assert updated.edi.get_section("head").project == "LINE_01"
    info = updated.edi.get_section("info")
    assert info.Source.project == "LINE_01"
    assert info.reviewed_by == "operator"
    assert "REVIEWED_BY=operator" in info.info_text


def test_dictionary_rename_preserves_order_and_audits(tmp_path: Path, simulated_edi: Path) -> None:
    sites = _sites(tmp_path, simulated_edi)
    editor = SiteMetadataEditor({"A01": {"name": "L22_010"}, "A02": {"name": "L22_011"}})

    renamed = editor.apply(sites)

    assert isinstance(renamed, Sites)
    assert [site.name for site in renamed] == ["L22_010", "L22_011"]
    assert [site.name for site in sites] == ["A01", "A02"]
    audit = editor.audit()
    assert isinstance(audit, pd.DataFrame)
    assert audit["old_name"].tolist() == ["A01", "A02"]
    assert audit["new_name"].tolist() == ["L22_010", "L22_011"]
    assert audit["status"].tolist() == ["updated", "updated"]


def test_sequence_and_callable_rename(tmp_path: Path, simulated_edi: Path) -> None:
    sites = _sites(tmp_path, simulated_edi)
    sequenced = rename_sites(sites, ["P01", "P02"])
    generated = rename_sites(sites, lambda _edi, index: f"G{index + 1:02d}")

    assert [site.name for site in sequenced] == ["P01", "P02"]
    assert [site.name for site in generated] == ["G01", "G02"]


def test_missing_and_duplicate_names_are_validated_before_mutation(
    tmp_path: Path, simulated_edi: Path
) -> None:
    sites = _sites(tmp_path, simulated_edi)

    with pytest.raises(KeyError, match="UNKNOWN"):
        rename_sites(sites, {"UNKNOWN": "X01"})
    with pytest.raises(ValueError, match="duplicate"):
        rename_sites(sites, {"A01": "SAME", "A02": "same"})

    assert [site.name for site in sites] == ["A01", "A02"]


def test_sites_convenience_methods_and_inplace(tmp_path: Path, simulated_edi: Path) -> None:
    sites = _sites(tmp_path, simulated_edi)
    returned = sites.rename({"A01": "B01", "A02": "B02"}, inplace=True)
    updated = returned.update_metadata({"B01": {"elev": 101.0}, "B02": {"elev": 102.0}})

    assert returned is sites
    assert [site.name for site in sites] == ["B01", "B02"]
    assert [float(site.coords[2]) for site in updated] == [101.0, 102.0]


def test_apply_and_write_uses_updated_station_names(tmp_path: Path, simulated_edi: Path) -> None:
    sites = _sites(tmp_path, simulated_edi)
    editor = SiteMetadataEditor({"A01": "E01", "A02": "E02"})
    manifest = tmp_path / "manifest.csv"

    result = editor.apply_and_write(
        sites,
        tmp_path / "exported",
        manifest_csv=manifest,
    )

    assert isinstance(result, Sites)
    assert [path.name for path in editor.output_paths_] == ["E01.edi", "E02.edi"]
    assert all(path.exists() for path in editor.output_paths_)
    assert pd.read_csv(manifest)["station"].tolist() == ["E01", "E02"]


def test_plan_reports_requested_before_and_after_without_mutation(
    tmp_path: Path, simulated_edi: Path
) -> None:
    sites = _sites(tmp_path, simulated_edi)
    original_elevation = float(sites[0].coords[2])
    editor = SiteMetadataEditor(
        {"A01": {"name": "B01", "head": {"project": "P01"}}},
        missing="ignore",
    )

    audit = editor.plan(sites)

    assert sites[0].name == "A01"
    assert float(sites[0].coords[2]) == original_elevation
    row = audit.iloc[0]
    assert row["requested_fields"] == ["name", "head.project"]
    assert row["before"]["name"] == "A01"
    assert row["after"]["name"] == "B01"
    assert row["after"]["head.project"] == "P01"
    assert editor.output_paths_ == []


def test_dataframe_and_csv_sources_resolve_aliases_and_dotted_fields(
    tmp_path: Path, simulated_edi: Path
) -> None:
    sites = _sites(tmp_path, simulated_edi)
    table = pd.DataFrame(
        {
            "station": ["A01", "A02"],
            "new_name": ["T01", "T02"],
            "latitude": [5.1, 5.2],
            "longitude": [-3.1, -3.2],
            "elevation": [101.0, np.nan],
            "head.project": ["TABLE", "TABLE"],
            "info.processingtag": ["checked", "checked"],
        }
    )

    from_frame = update_metadata_all(sites, table)
    csv_path = tmp_path / "metadata.csv"
    table.to_csv(csv_path, index=False)
    from_csv = update_metadata_all(sites, csv_path)

    for result in (from_frame, from_csv):
        assert [site.name for site in result] == ["T01", "T02"]
        assert result[0].coords == (5.1, -3.1, 101.0)
        assert result[1].coords[:2] == (5.2, -3.2)
        assert result[0].edi.get_section("head").project == "TABLE"
        assert result[0].edi.get_section("info").Processing.processingtag == "checked"


def test_set_unset_transform_and_callable_values(tmp_path: Path, simulated_edi: Path) -> None:
    sites = _sites(tmp_path, simulated_edi)
    sites[0].edi.get_section("head").project = "old"
    sites[0].edi.get_section("head").county = "remove-me"
    old_elevation = float(sites[0].coords[2])
    editor = SiteMetadataEditor(
        {
            "A01": {
                "name": lambda current, _edi, index: f"{current}_R{index}",
                "set": {
                    "head.project": lambda current: current.upper(),
                    "info.processing.processingsoftware.name": "META-TEST",
                    "edi.review_state": "approved",
                },
                "transform": {"head.elev": lambda value: value + 2.5},
                "unset": ["head.county"],
            }
        },
        missing="ignore",
    )

    updated = editor.apply(sites)
    first = updated[0]

    assert first.name == "A01_R0"
    assert first.edi.get_section("head").project == "OLD"
    assert first.edi.get_section("head").county is None
    assert float(first.coords[2]) == old_elevation + 2.5
    assert first.edi.review_state == "approved"
    assert first.edi.get_section("info").Processing.ProcessingSoftware.name == "META-TEST"
    assert set(editor.records_[0].changed_fields) >= {
        "name",
        "head.project",
        "head.elev",
        "head.county",
        "info.processing.processingsoftware.name",
        "edi.review_state",
    }


def test_coords_mapping_partial_update_and_validation(tmp_path: Path, simulated_edi: Path) -> None:
    site = _sites(tmp_path, simulated_edi)[0]
    original = site.coords

    updated = update_metadata(
        site,
        {"coords": {"latitude": 6.25, "longitude": -4.5}},
    )

    assert updated.coords == (6.25, -4.5, original[2])
    with pytest.raises(ValueError, match="latitude"):
        update_metadata(site, {"lat": 91.0})
    with pytest.raises(ValueError, match="longitude"):
        update_metadata(site, {"lon": np.inf})
    with pytest.raises(ValueError, match="elevation"):
        update_metadata(site, {"elev": np.nan})


def test_arbitrary_sections_case_insensitive_paths_and_sectid_sync(
    tmp_path: Path, simulated_edi: Path
) -> None:
    site = _sites(tmp_path, simulated_edi)[0]
    updated = update_metadata(
        site,
        {
            "name": "SYNC01",
            "sections": {
                "definemeas": {"custom.review": "accepted"},
                "quality": {"score": 0.95},
            },
            "info": {"processing.processingsoftware.name": "REVIEWER"},
        },
    )

    assert updated.edi.get_section("definemeas").custom.review == "accepted"
    assert updated.edi.get_section("quality").score == 0.95
    assert updated.edi.get_section("info").Processing.ProcessingSoftware.name == "REVIEWER"
    assert updated.edi.get_section("mtsect").sectid == "SYNC01"


def test_custom_info_is_added_replaced_and_removed(
    simulated_edi: Path,
) -> None:
    original = Site(EDIFile(simulated_edi))
    added = update_metadata(original, {"info": {"reviewed_by": "alice"}})
    replaced = update_metadata(added, {"info": {"reviewed_by": "bob"}})
    removed = update_metadata(replaced, {"unset": ["info.reviewed_by"]})

    assert "REVIEWED_BY=alice" in added.edi.get_section("info").info_text
    assert "REVIEWED_BY=bob" in replaced.edi.get_section("info").info_text
    assert sum("REVIEWED_BY=" in line for line in replaced.edi.get_section("info").info_text) == 1
    assert not any("REVIEWED_BY=" in line for line in removed.edi.get_section("info").info_text)


def test_error_policies_and_atomic_inplace_commit(tmp_path: Path, simulated_edi: Path) -> None:
    sites = _sites(tmp_path, simulated_edi)
    first_edi = sites[0].edi
    editor = SiteMetadataEditor({"A01": {"name": "B01"}, "A02": {"lat": 100.0}})

    with pytest.raises(ValueError, match="latitude"):
        editor.apply(sites, inplace=True)
    assert [site.name for site in sites] == ["A01", "A02"]
    assert sites[0].edi is first_edi

    warning_editor = SiteMetadataEditor(
        {"A01": {"name": "B01"}, "A02": {"lat": 100.0}},
        on_error="warn",
    )
    with pytest.warns(UserWarning, match="latitude"):
        returned = warning_editor.apply(sites, inplace=True)
    assert returned is sites
    assert sites[0].edi is first_edi
    assert [site.name for site in sites] == ["B01", "A02"]
    assert [record.status for record in warning_editor.records_] == [
        "updated",
        "error",
    ]


def test_custom_validators_and_missing_policies(tmp_path: Path, simulated_edi: Path) -> None:
    sites = _sites(tmp_path, simulated_edi)

    def reject(edi: EDIFile) -> bool:
        return not str(edi.name).startswith("BAD")

    with pytest.raises(ValueError, match="validator rejected"):
        SiteMetadataEditor(
            {"A01": {"name": "BAD01"}},
            missing="ignore",
            validators=[reject],
        ).apply(sites)

    with pytest.warns(UserWarning, match="UNKNOWN"):
        unchanged = SiteMetadataEditor({"UNKNOWN": {"name": "X"}}, missing="warn").apply(sites)
    assert [site.name for site in unchanged] == ["A01", "A02"]


def test_invalid_sources_names_and_actions_fail_clearly(
    tmp_path: Path, simulated_edi: Path
) -> None:
    sites = _sites(tmp_path, simulated_edi)

    with pytest.raises(ValueError, match="duplicate case-insensitive key"):
        SiteMetadataEditor({"A01": {"name": "X"}, "a01": {"name": "Y"}}).apply(sites)
    with pytest.raises(ValueError, match="empty station"):
        rename_sites(sites, {"A01": ""}, missing="ignore")
    with pytest.raises(ValueError, match="cannot be unset"):
        update_metadata(sites[0], {"unset": ["head.dataid"]})
    with pytest.raises(TypeError, match="transform"):
        update_metadata(sites[0], {"transform": {"head.elev": 3.0}})
    with pytest.raises(TypeError, match="one-shot"):
        SiteMetadataEditor([{"name": "X01"}, {"name": "X02"}]).apply(
            (site.edi for site in sites), inplace=True
        )


def test_empty_audit_has_stable_schema() -> None:
    audit = SiteMetadataEditor({}).audit()

    assert list(audit.columns) == [
        "index",
        "old_name",
        "new_name",
        "changed_fields",
        "status",
        "error",
        "requested_fields",
        "before",
        "after",
    ]


@pytest.mark.integration
def test_real_l18_full_line_rename_preserves_data_and_source() -> None:
    paths = _real_line("L18PLT")
    sites = Sites([EDIFile(path) for path in paths])
    original_names = [site.name for site in sites]
    original_freq = [np.array(site.freq, copy=True) for site in sites]
    mapping = {name: f"L18_{index + 1:03d}" for index, name in enumerate(original_names)}

    editor = SiteMetadataEditor(
        {
            name: {
                "name": mapping[name],
                "head": {"project": "WILLY_L18"},
                "info": {"processingtag": "metadata-reviewed"},
            }
            for name in original_names
        }
    )
    updated = editor.apply(sites)

    assert len(updated) == len(paths) == 28
    assert [site.name for site in sites] == original_names
    assert [site.name for site in updated] == list(mapping.values())
    assert all(site.edi.get_section("head").project == "WILLY_L18" for site in updated)
    assert all(
        site.edi.get_section("info").Processing.processingtag == "metadata-reviewed"
        for site in updated
    )
    assert all(
        np.array_equal(site.freq, expected) for site, expected in zip(updated, original_freq)
    )
    assert all(site.edi.get_section("mtsect").sectid == site.name for site in updated)
    assert editor.audit()["status"].eq("updated").all()


@pytest.mark.integration
def test_real_l22_dataframe_plan_apply_export_and_reload(tmp_path: Path) -> None:
    paths = _real_line("L22PLT")
    # Three differently named stations are enough to exercise real headers,
    # INFO nesting, transfer-function sections, writer, and parser round-trip.
    selected_paths = [paths[0], paths[len(paths) // 2], paths[-1]]
    sites = Sites([EDIFile(path) for path in selected_paths])
    original_names = [site.name for site in sites]
    original_shapes = [np.asarray(site.z).shape for site in sites]
    table = pd.DataFrame(
        {
            "station": original_names,
            "new_name": [f"L22_META_{index + 1:02d}" for index in range(3)],
            "head.project": ["WILLY_L22"] * 3,
            "info.processingtag": ["roundtrip"] * 3,
            "section.mtsect.reviewcode": ["QC-PASS"] * 3,
        }
    )
    editor = SiteMetadataEditor(table)

    preview = editor.plan(sites)
    assert [site.name for site in sites] == original_names
    assert preview["status"].eq("updated").all()

    manifest = tmp_path / "l22_metadata_manifest.csv"
    updated = editor.apply_and_write(
        sites,
        tmp_path / "l22_metadata",
        manifest_csv=manifest,
    )
    reloaded = [EDIFile(path) for path in editor.output_paths_]

    assert [site.name for site in updated] == [
        "L22_META_01",
        "L22_META_02",
        "L22_META_03",
    ]
    assert [path.name for path in editor.output_paths_] == [
        "L22_META_01.edi",
        "L22_META_02.edi",
        "L22_META_03.edi",
    ]
    assert pd.read_csv(manifest)["station"].tolist() == [
        "L22_META_01",
        "L22_META_02",
        "L22_META_03",
    ]
    assert [edi.get_section("head").dataid for edi in reloaded] == [
        "L22_META_01",
        "L22_META_02",
        "L22_META_03",
    ]
    assert [np.asarray(edi.Z.z).shape for edi in reloaded] == original_shapes
