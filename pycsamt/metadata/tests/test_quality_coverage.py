from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from pycsamt.metadata.quality import (
    ComponentQuality,
    DataQuality,
    QualityComment,
    QualityFlag,
    TransferFunctionQuality,
    _bar,
    _extract_component,
    _safe_array,
    assess_collection,
    quality_dataframe,
)


@pytest.mark.parametrize(
    ("coverage", "flag"),
    [
        (0.0, QualityFlag.MISSING),
        (0.2, QualityFlag.POOR),
        (0.5, QualityFlag.PARTIAL),
        (0.9, QualityFlag.GOOD),
    ],
)
def test_quality_flag_thresholds_and_rank(coverage, flag):
    assert QualityFlag.from_coverage(coverage) is flag
    assert QualityFlag.GOOD.rank == 3
    assert QualityFlag.worst([]) is QualityFlag.MISSING
    assert QualityFlag.best([]) is QualityFlag.MISSING
    flags = [QualityFlag.GOOD, QualityFlag.POOR, QualityFlag.PARTIAL]
    assert QualityFlag.worst(flags) is QualityFlag.POOR
    assert QualityFlag.best(flags) is QualityFlag.GOOD


def test_component_quality_handles_missing_real_complex_and_snr():
    for value in (None, []):
        missing = ComponentQuality.from_array("Zxx", value)
        assert missing.flag is QualityFlag.MISSING and missing.n_total == 0

    real = ComponentQuality.from_array(
        "Zxy", [1.0, np.nan, 3.0], snr=[10.0, np.nan, 20.0]
    )
    assert real.coverage == pytest.approx(2 / 3)
    assert real.snr_mean == 15.0 and real.snr_std == 5.0
    assert real.pct_str == "67%"
    assert real.to_dict()["snr_mean"] == 15.0
    assert "SNR=15.0" in repr(real)

    complex_quality = ComponentQuality.from_array(
        "Zyx", [1 + 2j, complex(np.nan, 1)]
    )
    assert complex_quality.n_valid == 1
    no_finite_snr = ComponentQuality.from_array("Zyy", [1], snr=[np.nan])
    assert no_finite_snr.snr_mean is None


def _site(name="S1"):
    z = np.arange(12, dtype=float).reshape(3, 2, 2).astype(complex)
    z[1, 0, 0] = np.nan + 1j * np.nan
    return SimpleNamespace(
        name=name,
        freq=np.array([1.0, 10.0, 100.0]),
        z=z,
        tipper=np.ones((3, 1, 2), dtype=complex),
    )


def test_data_quality_from_site_accessors_output_and_collection():
    quality = DataQuality.from_site(_site())
    assert quality.station == "S1" and quality.n_freq == 3
    assert quality.freq_min == 1.0 and quality.freq_max == 100.0
    assert quality.get("zXX").n_valid == 2
    assert quality.get("missing") is None
    assert len(quality.z_components) == 4
    assert quality.has_tipper
    assert quality.mean_coverage > 0
    assert "Frequencies" in quality.summary()
    assert quality.to_dict()["station"] == "S1"
    assert "overall=" in repr(quality)
    assert assess_collection([_site("A"), _site("B")])[1].station == "B"

    frame = quality_dataframe([_site()], api=False)
    assert frame.loc[0, "station"] == "S1"
    assert "cov_Zxx" in frame.columns


def test_data_quality_missing_site_content():
    quality = DataQuality.from_site(SimpleNamespace(name="empty", freq=[]))
    assert quality.n_freq == 0
    assert quality.overall is QualityFlag.MISSING
    assert quality.mean_coverage == 0.0
    assert not quality.has_tipper
    assert "MISSING" in quality.summary()


def test_quality_comment_and_transfer_quality_validation():
    with pytest.raises(ValueError, match="non-empty"):
        QualityComment(" ")
    comment = QualityComment(" ok ", author=" ")
    assert comment.text == "ok" and comment.author is None

    for kwargs, message in [
        ({"rating": 6}, "rating"),
        ({"good_from_period": 0}, "finite and positive"),
        ({"good_to_period": np.inf}, "finite and positive"),
        ({"warning_flag": 2}, "warning_flag"),
        ({"comments": ["bad"]}, "QualityComment"),
        ({"warnings": ["bad"]}, "QualityComment"),
    ]:
        with pytest.raises((ValueError, TypeError), match=message):
            TransferFunctionQuality(**kwargs)

    assert TransferFunctionQuality().rating_label is None
    labels = {
        rating: TransferFunctionQuality(rating=rating).rating_label
        for rating in range(6)
    }
    assert labels[5] == "great" and labels[0] == "not_assessed"
    assert not TransferFunctionQuality(warning_flag=0).has_warning
    assert TransferFunctionQuality(warnings=[QualityComment("warn")]).has_warning


def test_quality_array_helpers_and_progress_bar():
    assert _safe_array(None) is None and _safe_array([]) is None
    np.testing.assert_array_equal(_safe_array([1]), [1])
    cube = np.arange(8).reshape(2, 2, 2)
    np.testing.assert_array_equal(_extract_component(cube, 3), cube[:, 1, 1])
    flat4 = np.arange(8).reshape(2, 4)
    np.testing.assert_array_equal(_extract_component(flat4, 2), flat4[:, 2])
    np.testing.assert_array_equal(_extract_component(np.eye(2), 0), np.eye(2).ravel())
    assert len(_bar(-1, width=4)) == 4
    assert len(_bar(2, width=4)) == 4
