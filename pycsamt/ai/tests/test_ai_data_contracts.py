"""Contract tests for :mod:`pycsamt.ai.data`."""

from __future__ import annotations

import json

import numpy as np
import pytest

from pycsamt.ai.data import (
    ArtifactRecord,
    ComplexZScore,
    DatasetManifest,
    ImpedanceConvention,
    NormalizedSurvey,
    RealizationSplit,
    SurveyCoverage,
    SurveyData,
    canonical_hash,
    merge_surveys,
    realization_folds,
    sha256_file,
    split_realizations,
)


def _survey() -> SurveyData:
    z = np.arange(24, dtype=float).reshape(3, 4, 2)
    z = z + 1j * (z + 0.5)
    z[1, 2, 0] = np.nan + 1j * np.nan
    errors = np.full(z.shape, 0.05)
    return SurveyData(
        impedance=z,
        frequencies_hz=[1000.0, 100.0, 10.0, 1.0],
        station_names=("S01", "S02", "S03"),
        components=("xy", "yx"),
        coordinates_m=[[0, 0], [100, 0], [200, 10]],
        impedance_error=errors,
        tipper=np.ones((3, 4, 2), dtype=complex),
        crs="EPSG:32630",
        metadata={"survey": "test", "rotation_deg": 0.0},
    )


def test_survey_contract_shape_mask_coordinates_and_readonly():
    survey = _survey()
    assert survey.shape == (3, 4, 2)
    assert survey.coordinates_m.shape == (3, 3)
    assert np.all(np.isnan(survey.coordinates_m[:, 2]))
    assert survey.n_valid == survey.impedance.size - 1
    assert not survey.valid[1, 2, 0]
    assert not survey.impedance.flags.writeable
    with pytest.raises(ValueError):
        survey.impedance[0, 0, 0] = 1j


@pytest.mark.parametrize(
    "change, match",
    [
        ({"frequencies_hz": [1, 10, 5, 100]}, "strictly monotonic"),
        ({"station_names": ("same", "same", "third")}, "unique"),
        ({"coordinates_m": [[0, 0], [1, np.nan], [2, 0]]}, "x/y"),
        ({"valid": np.ones((2, 2), dtype=bool)}, "same shape"),
    ],
)
def test_survey_contract_rejects_ambiguous_axes(change, match):
    kwargs = dict(
        impedance=np.ones((3, 4, 2), dtype=complex),
        frequencies_hz=[1000, 100, 10, 1],
        station_names=("a", "b", "c"),
        components=("xy", "yx"),
        coordinates_m=np.zeros((3, 2)),
    )
    kwargs.update(change)
    with pytest.raises(ValueError, match=match):
        SurveyData(**kwargs)


def test_survey_select_and_npz_roundtrip(tmp_path):
    selected = _survey().select(stations=[2, 0], frequencies=[1, 2], components=[1])
    assert selected.shape == (2, 2, 1)
    assert selected.station_names == ("S03", "S01")
    assert selected.components == ("yx",)
    path = selected.to_npz(tmp_path / "survey.npz")
    restored = SurveyData.from_npz(path)
    np.testing.assert_array_equal(restored.impedance, selected.impedance)
    np.testing.assert_array_equal(restored.valid, selected.valid)
    np.testing.assert_array_equal(restored.tipper, selected.tipper)
    assert restored.station_names == selected.station_names
    assert dict(restored.metadata) == dict(selected.metadata)
    assert restored.convention == selected.convention


def test_impedance_convention_validates_normalizes_and_roundtrips():
    convention = ImpedanceConvention(time_dependence="exp(-iwt)", rotation_deg=375.0)
    assert convention.rotation_deg == pytest.approx(15.0)
    assert ImpedanceConvention.from_dict(convention.to_dict()) == convention
    with pytest.raises(ValueError, match="time_dependence"):
        ImpedanceConvention(time_dependence="unknown")
    with pytest.raises(ValueError, match="units"):
        ImpedanceConvention(units="mV/km/nT")


def test_coverage_named_access_selection_summary_and_metadata_copy():
    survey = _survey()
    coverage = survey.coverage()
    assert isinstance(coverage, SurveyCoverage)
    assert coverage.overall == pytest.approx(survey.n_valid / survey.impedance.size)
    assert not coverage.complete
    assert coverage.tipper_overall == 1.0
    assert survey.station_index("S02") == 1
    assert survey.component_index("yx") == 1
    values, errors, valid = survey.component_data("xy")
    assert values.shape == errors.shape == valid.shape == (3, 4)
    selected = survey.select_names(
        stations=["S03", "S01"],
        components=["yx"],
        frequency_min_hz=10.0,
        frequency_max_hz=100.0,
    )
    assert selected.shape == (2, 2, 1)
    assert selected.station_names == ("S03", "S01")
    tagged = survey.with_metadata({"line": "L18"})
    assert tagged.metadata["survey"] == "test"
    assert tagged.metadata["line"] == "L18"
    assert "line" not in survey.metadata
    assert survey.summary()["frequency_order"] == "descending"


def test_compatibility_and_merge_are_strict_and_lossless():
    left = _survey().select(stations=[0, 1])
    right = _survey().select(stations=[2])
    left.assert_compatible(right)
    merged = merge_surveys([left, right])
    assert merged.shape == _survey().shape
    assert merged.station_names == _survey().station_names
    np.testing.assert_array_equal(merged.impedance, _survey().impedance)
    assert merged.metadata["source_survey_count"] == 2

    incompatible = SurveyData(
        np.ones((1, 4, 2), complex),
        [1000, 100, 10, 1],
        ["other"],
        ["xy", "yx"],
        [[0, 0]],
        crs="EPSG:4326",
    )
    with pytest.raises(ValueError, match="reference systems"):
        left.assert_compatible(incompatible)


def test_named_selection_and_lookup_fail_loudly():
    survey = _survey()
    with pytest.raises(KeyError, match="unknown station"):
        survey.station_index("missing")
    with pytest.raises(KeyError, match="unknown component"):
        survey.component_index("xx")
    with pytest.raises(ValueError, match="select no observations"):
        survey.select_names(frequency_min_hz=1e9)
    with pytest.raises(ValueError, match="cannot exceed"):
        survey.select_names(frequency_min_hz=100, frequency_max_hz=10)


def test_complex_zscore_uses_valid_training_values_and_roundtrips():
    survey = _survey()
    state = ComplexZScore.fit(survey)
    features, mask = state.transform(survey, fill_value=0.0)
    assert features.shape == (3, 4, 2, 2)
    assert mask.shape == features.shape
    assert np.all(features[~mask] == 0.0)
    restored = state.inverse_transform(features)
    np.testing.assert_allclose(restored[survey.valid], survey.impedance[survey.valid])
    state2 = ComplexZScore.from_dict(json.loads(json.dumps(state.to_dict())))
    np.testing.assert_allclose(state2.mean, state.mean)
    assert state2.state_hash == state.state_hash
    assert state.training_survey_count == 1
    assert state.training_station_count == survey.n_stations
    assert state.count[2, 0, 0] == 2
    with pytest.raises(ValueError, match="frequency grid"):
        state.transform(survey.select(frequencies=[1, 2]))


def test_normalized_survey_propagates_errors_flattens_and_is_readonly():
    survey = _survey()
    state = ComplexZScore.fit(survey)
    normalized = state.transform_survey(survey)
    assert isinstance(normalized, NormalizedSurvey)
    assert normalized.shape == (3, 4, 2, 2)
    assert normalized.errors.shape == normalized.shape
    assert normalized.n_valid_observations == survey.n_valid
    values, valid = normalized.flatten()
    assert values.shape == valid.shape == (3, 16)
    assert not normalized.values.flags.writeable
    restored = state.inverse_transform(normalized.values, valid=normalized.valid)
    np.testing.assert_allclose(restored[survey.valid], survey.impedance[survey.valid])
    assert np.all(np.isnan(restored[~survey.valid]))


def test_complex_zscore_inverse_variance_weighting_and_ddof():
    z = np.array([[[0 + 0j]], [[10 + 10j]]])
    errors = np.array([[[1.0]], [[10.0]]])
    survey = SurveyData(
        z,
        [1],
        ["precise", "noisy"],
        ["xy"],
        [[0, 0], [1, 0]],
        impedance_error=errors,
    )
    uniform = ComplexZScore.fit(survey)
    weighted = ComplexZScore.fit(survey, weighting="inverse_variance")
    assert weighted.mean[0, 0, 0] < uniform.mean[0, 0, 0]
    assert weighted.weighting == "inverse_variance"
    with pytest.raises(ValueError, match="ddof must be zero"):
        ComplexZScore.fit(survey, weighting="inverse_variance", ddof=1)
    sample = ComplexZScore.fit(survey, ddof=1)
    assert sample.scale[0, 0, 0] == pytest.approx(np.sqrt(50.0))


def test_complex_zscore_clipping_convention_and_error_contracts():
    survey = _survey()
    state = ComplexZScore.fit(survey)
    features, _ = state.transform(survey, clip=0.5)
    assert np.max(np.abs(features)) <= 0.5
    errors, mask = state.transform_errors(survey)
    assert errors.shape == mask.shape == features.shape
    incompatible = SurveyData(
        survey.impedance,
        survey.frequencies_hz,
        survey.station_names,
        survey.components,
        survey.coordinates_m,
        impedance_error=survey.impedance_error,
        convention=ImpedanceConvention(time_dependence="exp(-iwt)"),
    )
    with pytest.raises(ValueError, match="convention differs"):
        state.transform(incompatible)
    no_errors = SurveyData(
        survey.impedance,
        survey.frequencies_hz,
        survey.station_names,
        survey.components,
        survey.coordinates_m,
    )
    with pytest.raises(ValueError, match="does not contain"):
        state.transform_errors(no_errors)


def test_schema_one_normalization_state_remains_readable():
    legacy = {
        "schema_version": 1,
        "kind": "complex_cartesian_zscore",
        "mean": [[[1.0, 2.0]]],
        "scale": [[[3.0, 4.0]]],
        "frequencies_hz": [1.0],
        "components": ["xy"],
        "eps": 1e-8,
    }
    restored = ComplexZScore.from_dict(legacy)
    assert restored.convention is None
    assert restored.count.tolist() == [[[1, 1]]]
    assert ComplexZScore.from_dict(restored.to_dict()).state_hash == restored.state_hash


def test_realization_split_is_deterministic_disjoint_and_serializable():
    ids = [f"model-{i}" for i in range(20)]
    first = split_realizations(ids, validation_fraction=0.2, test_fraction=0.2, seed=7)
    second = split_realizations(ids, validation_fraction=0.2, test_fraction=0.2, seed=7)
    assert first == second
    assert set(first.all_ids) == set(ids)
    assert first.partition(first.train[0]) == "train"
    assert RealizationSplit.from_dict(first.to_dict()) == first
    assert len(first.split_hash) == 64
    assert first.sizes == {"train": 12, "validation": 4, "test": 4}
    assert first.fractions["train"] == pytest.approx(0.6)
    np.testing.assert_array_equal(
        first.mask([first.train[0], first.test[0]], "train"),
        [True, False],
    )
    first.assert_complete(list(reversed(ids)))
    with pytest.raises(ValueError, match="disjoint"):
        RealizationSplit(("a",), ("a",), ())


def test_split_is_input_order_independent_and_lineage_safe():
    ids = [f"r{i}" for i in range(12)]
    forward = split_realizations(
        ids, validation_fraction=0.25, test_fraction=0.25, seed=9
    )
    reverse = split_realizations(
        list(reversed(ids)),
        validation_fraction=0.25,
        test_fraction=0.25,
        seed=9,
    )
    assert forward == reverse

    lineage = {
        "a-clean": "a",
        "a-noisy": "a",
        "b-clean": "b",
        "b-noisy": "b",
        "c": "c",
        "d": "d",
    }
    grouped = split_realizations(
        list(lineage),
        validation_fraction=0.2,
        test_fraction=0.2,
        seed=3,
        lineage=lineage,
    )
    assert grouped.partition("a-clean") == grouped.partition("a-noisy")
    assert grouped.partition("b-clean") == grouped.partition("b-noisy")
    grouped.assert_no_lineage_leakage()
    assert grouped.strategy == "lineage_random"
    grouped_reversed = split_realizations(
        list(reversed(lineage)),
        validation_fraction=0.2,
        test_fraction=0.2,
        seed=3,
        lineage=lineage,
    )
    assert grouped_reversed == grouped


def test_lineage_leakage_partial_reassignment_and_unknown_ids_fail():
    lineage = {"a1": "a", "a2": "a", "b": "b", "c": "c"}
    split = RealizationSplit(("a1", "a2", "b"), (), ("c",), lineage=lineage)
    with pytest.raises(ValueError, match="affected lineage"):
        split.reassign(["a1"], "validation")
    moved = split.reassign(["a1", "a2"], "validation")
    assert moved.validation == ("a1", "a2")
    assert moved.strategy.endswith("+manual")
    with pytest.raises(KeyError, match="unknown realization"):
        split.partition("missing")
    with pytest.raises(KeyError, match="unknown realization IDs"):
        split.mask(["missing"], "train")
    np.testing.assert_array_equal(
        split.mask(["missing"], "train", unknown="false"), [False]
    )


def test_realization_folds_cover_every_group_once_without_leakage():
    lineage = {
        "a1": "a",
        "a2": "a",
        "b": "b",
        "c": "c",
        "d1": "d",
        "d2": "d",
    }
    folds = realization_folds(list(lineage), n_splits=3, seed=5, lineage=lineage)
    assert len(folds) == 3
    assert sorted(item for fold in folds for item in fold.test) == sorted(lineage)
    for fold in folds:
        assert not fold.validation
        fold.assert_no_lineage_leakage()
        assert set(fold.train).isdisjoint(fold.test)


def test_split_schema_one_and_validation_edges():
    legacy = {
        "schema_version": 1,
        "train": ["a"],
        "validation": [],
        "test": ["b"],
        "seed": 0,
    }
    restored = RealizationSplit.from_dict(legacy)
    assert restored.strategy == "random"
    assert not restored.lineage
    with pytest.raises(ValueError, match="exactly match"):
        split_realizations(
            ["a", "b", "c"],
            lineage={"a": "parent", "b": "parent"},
        )
    with pytest.raises(ValueError, match="too few independent"):
        split_realizations(
            ["a1", "a2"],
            validation_fraction=0.2,
            test_fraction=0.2,
            lineage={"a1": "a", "a2": "a"},
        )
    with pytest.raises(ValueError, match="n_splits"):
        realization_folds(["a", "b"], n_splits=3)


def test_manifest_hash_and_json_roundtrip_detect_tampering(tmp_path):
    split = RealizationSplit(("r1", "r2"), ("r3",), ("r4",), seed=1)
    manifest = DatasetManifest(
        dataset_id="willy-synthetic-v1",
        generator="correlated2d",
        generator_version="0.1.0",
        configuration={"seed": 1, "correlation_m": [1000.0, 100.0]},
        split=split,
        sample_count=4,
        created_utc="2026-07-26T00:00:00Z",
        artifacts={"data.npz": "a" * 64},
    )
    assert len(canonical_hash(manifest.configuration)) == 64
    with pytest.raises(TypeError):
        manifest.configuration["seed"] = 2
    restored = DatasetManifest.read_json(
        manifest.write_json(tmp_path / "manifest.json")
    )
    assert restored.to_dict() == manifest.to_dict()
    bad = manifest.to_dict()
    bad["configuration"]["seed"] = 99
    with pytest.raises(ValueError, match="does not match"):
        DatasetManifest.from_dict(bad)


def test_manifest_recursively_freezes_configuration_and_normalizes_time():
    manifest = DatasetManifest(
        dataset_id="dataset.v2",
        generator="package.generator",
        generator_version="git:abc123",
        configuration={"nested": {"values": [1, 2]}},
        split=RealizationSplit(("r1",), (), ()),
        sample_count=1,
        created_utc="2026-07-26T02:00:00+02:00",
    )
    assert manifest.created_utc == "2026-07-26T00:00:00Z"
    assert manifest.configuration["nested"]["values"] == (1, 2)
    with pytest.raises(TypeError):
        manifest.configuration["nested"]["new"] = 3
    assert len(manifest.manifest_hash) == 64
    assert manifest.realization_count == 1


def test_artifact_record_hash_verification_and_immutable_update(tmp_path):
    artifact = tmp_path / "models.npz"
    artifact.write_bytes(b"synthetic-models")
    record = ArtifactRecord.from_file(
        artifact, media_type="application/x-npz", role="models"
    )
    assert record.sha256 == sha256_file(artifact)
    assert record.size_bytes == artifact.stat().st_size
    assert ArtifactRecord.from_dict(record.to_dict()) == record

    manifest = DatasetManifest(
        "dataset",
        "generator",
        "1.0",
        {},
        RealizationSplit(("r1",), (), ()),
        1,
    ).with_artifact("models.npz", record)
    assert manifest.verify_artifacts(tmp_path) == {"models.npz": True}
    assert not DatasetManifest(
        "dataset",
        "generator",
        "1.0",
        {},
        manifest.split,
        1,
    ).artifacts

    artifact.write_bytes(b"tampered")
    assert manifest.verify_artifacts(tmp_path) == {"models.npz": False}
    with pytest.raises(ValueError, match="integrity verification failed"):
        manifest.verify_artifacts(tmp_path, raise_on_error=True)


def test_manifest_rejects_unsafe_paths_invalid_hashes_and_naive_time():
    split = RealizationSplit(("r1",), (), ())
    with pytest.raises(ValueError, match="64 hexadecimal"):
        ArtifactRecord("not-a-hash")
    with pytest.raises(ValueError, match="relative paths"):
        DatasetManifest(
            "dataset",
            "generator",
            "1",
            {},
            split,
            1,
            artifacts={"../x": "a" * 64},
        )
    with pytest.raises(ValueError, match="timezone"):
        DatasetManifest(
            "dataset", "generator", "1", {}, split, 1, created_utc="2026-07-26"
        )
    with pytest.raises(ValueError, match="portable"):
        DatasetManifest("bad id", "generator", "1", {}, split, 1)
    with pytest.raises(ValueError, match="colliding JSON key"):
        canonical_hash({1: "integer", "1": "string"})


def test_schema_one_manifest_is_upgraded_on_read():
    split = RealizationSplit(("r1",), (), ())
    configuration = {"seed": 4}
    legacy = {
        "schema_version": 1,
        "dataset_id": "legacy",
        "generator": "generator",
        "generator_version": "1",
        "configuration": configuration,
        "configuration_hash": canonical_hash(configuration),
        "split": split.to_dict(),
        "sample_count": 1,
        "created_utc": None,
        "artifacts": {"data.npz": "b" * 64},
    }
    restored = DatasetManifest.from_dict(legacy)
    assert restored.schema_version == 2
    assert restored.artifacts["data.npz"].sha256 == "b" * 64
