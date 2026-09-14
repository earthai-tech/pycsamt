from __future__ import annotations

import numpy as np
import pytest

from pycsamt.airborne import validation as v
from pycsamt.airborne.mobilemt import MobileMTReferenceStation
from pycsamt.metadata import ProcessingMeta, RemoteReferenceMeta, SiteMeta


class _Boom(ValueError):
    pass


# ─────────────────────────────────────────────────────────────────────────
# clean_identifier / normalize_optional_identifier
# ─────────────────────────────────────────────────────────────────────────


def test_clean_identifier_strips_and_validates():
    assert v.clean_identifier(" abc ", name="x") == "abc"
    with pytest.raises(ValueError):
        v.clean_identifier("   ", name="x")
    with pytest.raises(_Boom):
        v.clean_identifier("", name="x", error_cls=_Boom)


def test_normalize_optional_identifier():
    assert v.normalize_optional_identifier(None) is None
    assert v.normalize_optional_identifier("   ") is None
    assert v.normalize_optional_identifier(" abc ") == "abc"


# ─────────────────────────────────────────────────────────────────────────
# normalize_numeric_vector / normalize_object_vector
# ─────────────────────────────────────────────────────────────────────────


def test_normalize_numeric_vector_none_and_scalar():
    assert v.normalize_numeric_vector(None, name="x", size=3) is None
    arr = v.normalize_numeric_vector(5.0, name="x", size=1)
    assert arr.shape == (1,)


def test_normalize_numeric_vector_allows_nan_but_not_inf():
    arr = v.normalize_numeric_vector([1.0, np.nan, 3.0], name="x", size=3)
    assert np.isnan(arr[1])
    with pytest.raises(ValueError):
        v.normalize_numeric_vector([1.0, np.inf], name="x", size=2)


def test_normalize_numeric_vector_bad_ndim_and_size():
    with pytest.raises(ValueError):
        v.normalize_numeric_vector([[1.0, 2.0]], name="x", size=2)
    with pytest.raises(_Boom):
        v.normalize_numeric_vector([1.0, 2.0], name="x", size=3, error_cls=_Boom)


def test_normalize_object_vector():
    assert v.normalize_object_vector(None, name="ts", size=2) is None
    out = v.normalize_object_vector(["a", "b"], name="ts", size=2)
    assert out == ("a", "b")
    with pytest.raises(ValueError):
        v.normalize_object_vector(["a"], name="ts", size=2)


# ─────────────────────────────────────────────────────────────────────────
# normalize_frequency / resolve_frequency_or_periods
# ─────────────────────────────────────────────────────────────────────────


def test_normalize_frequency_scalar_promotion():
    arr = v.normalize_frequency(90.0)
    assert arr.shape == (1,)


def test_normalize_frequency_rejects_bad_shapes_and_values():
    with pytest.raises(ValueError):
        v.normalize_frequency([[1.0, 2.0]])
    with pytest.raises(ValueError):
        v.normalize_frequency([])
    with pytest.raises(ValueError):
        v.normalize_frequency([1.0, np.nan])
    with pytest.raises(ValueError):
        v.normalize_frequency([1.0, -1.0])
    with pytest.raises(_Boom):
        v.normalize_frequency([0.0], error_cls=_Boom)


def test_resolve_frequency_or_periods_requires_exactly_one():
    with pytest.raises(ValueError):
        v.resolve_frequency_or_periods(frequency=None, periods=None)
    with pytest.raises(ValueError):
        v.resolve_frequency_or_periods(frequency=[1.0], periods=[1.0])


def test_resolve_frequency_or_periods_from_frequency():
    freq, periods = v.resolve_frequency_or_periods(frequency=[10.0], periods=None)
    assert np.allclose(freq, [10.0])
    assert np.allclose(periods, [0.1])


def test_resolve_frequency_or_periods_from_periods():
    freq, periods = v.resolve_frequency_or_periods(frequency=None, periods=[0.1])
    assert np.allclose(periods, [0.1])
    assert np.allclose(freq, [10.0])


# ─────────────────────────────────────────────────────────────────────────
# resolve_line_frequency_grid
# ─────────────────────────────────────────────────────────────────────────


def test_resolve_line_frequency_grid_shared_vector():
    common, rows = v.resolve_line_frequency_grid(
        [10.0, 20.0], n_samples=5, n_frequency=2,
    )
    assert rows is None
    assert np.allclose(common, [10.0, 20.0])


def test_resolve_line_frequency_grid_shared_vector_length_mismatch():
    with pytest.raises(ValueError):
        v.resolve_line_frequency_grid([10.0], n_samples=5, n_frequency=2)


def test_resolve_line_frequency_grid_per_sample_matrix():
    grid = np.ones((3, 2))
    common, rows = v.resolve_line_frequency_grid(
        grid, n_samples=3, n_frequency=2,
    )
    assert common is None
    assert rows.shape == (3, 2)


def test_resolve_line_frequency_grid_invalid_shape():
    with pytest.raises(ValueError):
        v.resolve_line_frequency_grid(
            np.ones((3, 3)), n_samples=5, n_frequency=2,
        )


# ─────────────────────────────────────────────────────────────────────────
# normalize_positive_float / normalize_count_range / normalize_frequency_range
# ─────────────────────────────────────────────────────────────────────────


def test_normalize_positive_float():
    assert v.normalize_positive_float(2.5, name="rate") == 2.5
    with pytest.raises(ValueError):
        v.normalize_positive_float(0.0, name="rate")
    with pytest.raises(ValueError):
        v.normalize_positive_float(float("nan"), name="rate")


def test_normalize_count_range():
    assert v.normalize_count_range((1, 5), name="count") == (1, 5)
    assert v.normalize_count_range((3, 3), name="count") == (3, 3)
    with pytest.raises(ValueError):
        v.normalize_count_range((1, 2, 3), name="count")
    with pytest.raises(ValueError):
        v.normalize_count_range((0, 5), name="count")
    with pytest.raises(ValueError):
        v.normalize_count_range((5, 1), name="count")


def test_normalize_frequency_range():
    assert v.normalize_frequency_range((1.0, 10.0), name="band") == (1.0, 10.0)
    with pytest.raises(ValueError):
        v.normalize_frequency_range((1.0,), name="band")
    with pytest.raises(ValueError):
        v.normalize_frequency_range((float("nan"), 10.0), name="band")
    with pytest.raises(ValueError):
        v.normalize_frequency_range((0.0, 10.0), name="band")
    with pytest.raises(ValueError):
        v.normalize_frequency_range((10.0, 1.0), name="band")


# ─────────────────────────────────────────────────────────────────────────
# normalize_fixed_channels / normalize_estimate_array / normalize_sample_axis_array
# ─────────────────────────────────────────────────────────────────────────


def test_normalize_fixed_channels():
    out = v.normalize_fixed_channels(
        (" Ex ", "Ey"), expected=("Ex", "Ey"), name="channels",
    )
    assert out == ("Ex", "Ey")
    with pytest.raises(ValueError):
        v.normalize_fixed_channels(
            ("Hx", "Hy"), expected=("Ex", "Ey"), name="channels",
        )


def test_normalize_estimate_array_promotes_2d():
    arr = v.normalize_estimate_array(
        np.ones((3, 2)), n_frequency=1, tail=(3, 2), name="VAR",
    )
    assert arr.shape == (1, 3, 2)


def test_normalize_estimate_array_rejects_bad_shape_and_dtype():
    with pytest.raises(ValueError):
        v.normalize_estimate_array(
            np.ones((2, 3, 2)), n_frequency=1, tail=(3, 2), name="VAR",
        )
    with pytest.raises(ValueError):
        v.normalize_estimate_array(
            np.array([["a", "b"]]), n_frequency=1, tail=(1, 2), name="VAR",
        )


def test_normalize_sample_axis_array_promotes_single_sample():
    arr = v.normalize_sample_axis_array(
        np.ones((2, 2)), name="admittance", n_samples=1, expected=(1, 2, 2),
    )
    assert arr.shape == (1, 2, 2)


def test_normalize_sample_axis_array_rejects_bad_shape():
    with pytest.raises(ValueError):
        v.normalize_sample_axis_array(
            np.ones((2, 2)), name="admittance", n_samples=3, expected=(3, 2, 2),
        )


# ─────────────────────────────────────────────────────────────────────────
# normalize_record_mask
# ─────────────────────────────────────────────────────────────────────────


def test_normalize_record_mask_default_all_true():
    mask = v.normalize_record_mask(None, n_samples=3)
    assert mask.tolist() == [True, True, True]


def test_normalize_record_mask_validates_shape():
    mask = v.normalize_record_mask([True, False, True], n_samples=3)
    assert mask.tolist() == [True, False, True]
    with pytest.raises(ValueError):
        v.normalize_record_mask([True, False], n_samples=3)


# ─────────────────────────────────────────────────────────────────────────
# reference_station_mapping
# ─────────────────────────────────────────────────────────────────────────


def test_reference_station_mapping_none():
    assert v.reference_station_mapping(None, channel_fields=()) is None


def test_reference_station_mapping_minimal():
    ref = MobileMTReferenceStation(station_id="BASE01")
    out = v.reference_station_mapping(
        ref, channel_fields=("electric_channels",),
    )
    assert out["station_id"] == "BASE01"
    assert out["electric_channels"] == ["Ex", "Ey"]
    assert "site" not in out
    assert "attrs" not in out


def test_reference_station_mapping_with_site_and_attrs():
    ref = MobileMTReferenceStation(
        station_id="BASE01",
        site=SiteMeta(site_id="BASE01"),
        attrs={"note": "fixed"},
    )
    out = v.reference_station_mapping(
        ref, channel_fields=("electric_channels",),
    )
    assert out["site"]["site_id"] == "BASE01"
    assert out["attrs"] == {"note": "fixed"}


# ─────────────────────────────────────────────────────────────────────────
# merge_remote_reference_processing
# ─────────────────────────────────────────────────────────────────────────


def test_merge_remote_reference_processing_none_reference_returns_processing():
    processing = ProcessingMeta(processed_by="Alice")
    assert v.merge_remote_reference_processing(
        None, processing, reference_type="fixed_ground_magnetic",
        technology="ztem",
    ) is processing
    assert v.merge_remote_reference_processing(
        None, None, reference_type="fixed_ground_magnetic", technology="ztem",
    ) is None


def test_merge_remote_reference_processing_rejects_bad_processing_type():
    ref = MobileMTReferenceStation(station_id="BASE01")
    with pytest.raises(TypeError):
        v.merge_remote_reference_processing(
            ref, "not-a-processing-meta",
            reference_type="fixed_ground_magnetic", technology="ztem",
        )


def test_merge_remote_reference_processing_synthesizes_when_processing_none():
    ref = MobileMTReferenceStation(station_id="BASE01")
    processing = v.merge_remote_reference_processing(
        ref, None, reference_type="fixed_ground_magnetic",
        technology="ztem", extra={"channels": ["Hx", "Hy"]},
    )
    assert processing.remote_reference.site == "BASE01"
    assert processing.remote_reference.extra["technology"] == "ztem"
    assert processing.remote_reference.extra["channels"] == ["Hx", "Hy"]


def test_merge_remote_reference_processing_merges_when_no_existing_remote():
    ref = MobileMTReferenceStation(station_id="BASE01")
    processing = ProcessingMeta(
        processed_by="Alice", run_list=["run1"],
    )
    merged = v.merge_remote_reference_processing(
        ref, processing, reference_type="fixed_ground_magnetic",
        technology="ztem",
    )
    assert merged.processed_by == "Alice"
    assert merged.run_list == ["run1"]
    assert merged.remote_reference.site == "BASE01"


def test_merge_remote_reference_processing_matching_existing_site_passthrough():
    ref = MobileMTReferenceStation(station_id="BASE01")
    processing = ProcessingMeta(
        remote_reference=RemoteReferenceMeta(
            reference_type="fixed_ground_magnetic", site="BASE01",
        ),
    )
    merged = v.merge_remote_reference_processing(
        ref, processing, reference_type="fixed_ground_magnetic",
        technology="ztem",
    )
    assert merged is processing


def test_merge_remote_reference_processing_existing_remote_without_site_passthrough():
    ref = MobileMTReferenceStation(station_id="BASE01")
    processing = ProcessingMeta(
        remote_reference=RemoteReferenceMeta(
            reference_type="fixed_ground_magnetic", site=None,
        ),
    )
    merged = v.merge_remote_reference_processing(
        ref, processing, reference_type="fixed_ground_magnetic",
        technology="ztem",
    )
    assert merged is processing


def test_merge_remote_reference_processing_conflicting_site_raises():
    ref = MobileMTReferenceStation(station_id="BASE02")
    processing = ProcessingMeta(
        remote_reference=RemoteReferenceMeta(
            reference_type="fixed_ground_magnetic", site="BASE01",
        ),
    )
    with pytest.raises(ValueError):
        v.merge_remote_reference_processing(
            ref, processing, reference_type="fixed_ground_magnetic",
            technology="ztem",
        )
    with pytest.raises(_Boom):
        v.merge_remote_reference_processing(
            ref, processing, reference_type="fixed_ground_magnetic",
            technology="ztem", error_cls=_Boom,
        )


def test_emtf_class_returns_emtf_type():
    from pycsamt.emtf.document import EMTF

    assert v.emtf_class() is EMTF
