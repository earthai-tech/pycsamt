# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.models.mare2dem.merge`.

Every :class:`EMDataFile` here is built in-memory with
:class:`MTConfig`/:class:`CSEMConfig` — nothing depends on the
gitignored ``data/mare2dem/`` example data, so the coverage gain is
real in CI (the existing ``TestMergeDataFiles.test_merge_real_csem_files``
in ``test_mare2dem_data.py`` is the only CSEM-merge test and it is
gated behind a bundled-data fixture that skips there).
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.models.mare2dem.iotools.emdata import (
    CSEMConfig,
    DCConfig,
    EMDataFile,
    MTConfig,
    UTMOrigin,
)
from pycsamt.models.mare2dem.merge import (
    _check_utm,
    _merge_arrays,
    _merge_csem,
    _merge_lists,
    _merge_mt,
    merge_data_files,
    merge_emdata,
)


# ---------------------------------------------------------------------------
# _merge_arrays / _merge_lists
# ---------------------------------------------------------------------------


def test_merge_arrays_empty_a_returns_b_copy():
    a = np.empty((0, 2))
    b = np.array([[1.0, 2.0], [3.0, 4.0]])
    combined, ib = _merge_arrays(a, b)
    np.testing.assert_array_equal(combined, b)
    np.testing.assert_array_equal(ib, [1, 2])
    combined[0, 0] = 999.0
    assert b[0, 0] == 1.0  # copy, not a view


def test_merge_arrays_dedup_existing_and_new_rows():
    a = np.array([[1.0, 2.0], [3.0, 4.0]])
    b = np.array([[1.0, 2.0], [5.0, 6.0]])
    combined, ib = _merge_arrays(a, b)
    assert combined.shape == (3, 2)
    assert ib[0] == 1  # duplicate of a's row 0
    assert ib[1] == 3  # newly appended


def test_merge_lists_dedup_existing_and_new():
    a = ["S1", "S2"]
    b = ["S1", "S3"]
    combined, ib = _merge_lists(a, b)
    assert combined == ["S1", "S2", "S3"]
    assert list(ib) == [1, 3]


# ---------------------------------------------------------------------------
# _check_utm
# ---------------------------------------------------------------------------


def test_check_utm_passes_when_equal():
    em1 = EMDataFile(utm=UTMOrigin(grid=12, north0=1.0, east0=2.0, theta=0.0))
    em2 = EMDataFile(utm=UTMOrigin(grid=12, north0=1.0, east0=2.0, theta=0.0))
    _check_utm(em1, em2)  # no raise


@pytest.mark.parametrize(
    "field", ["north0", "east0", "theta"]
)
def test_check_utm_raises_when_any_field_differs(field):
    em1 = EMDataFile(utm=UTMOrigin(grid=12, north0=1.0, east0=2.0, theta=0.0))
    kwargs = {"grid": 12, "north0": 1.0, "east0": 2.0, "theta": 0.0}
    kwargs[field] = kwargs[field] + 100.0
    em2 = EMDataFile(utm=UTMOrigin(**kwargs))
    with pytest.raises(ValueError, match="UTM"):
        _check_utm(em1, em2)


# ---------------------------------------------------------------------------
# Helpers to build synthetic MT / CSEM EMDataFile objects
# ---------------------------------------------------------------------------


def _mt_em(freqs, rx_y, comment="", name_prefix="S"):
    em = EMDataFile()
    n_rx = len(rx_y)
    em.mt = MTConfig(
        frequencies=np.asarray(freqs, dtype=float),
        receivers=np.zeros((n_rx, 8)),
        receiver_name=[f"{name_prefix}{i:03d}" for i in range(n_rx)],
    )
    em.mt.receivers[:, 1] = rx_y
    rows = []
    for ifreq in range(1, len(freqs) + 1):
        for irx in range(1, n_rx + 1):
            rows.append([123, ifreq, irx, irx, 1.5, 0.1])
            rows.append([104, ifreq, irx, irx, 45.0, 2.0])
    em.data = np.array(rows, dtype=float)
    em.comment = comment
    return em


def _csem_em(freqs, tx_x, rx_x, phase_convention="lag"):
    em = EMDataFile()
    n_tx = len(tx_x)
    n_rx = len(rx_x)
    em.csem = CSEMConfig(
        phase_convention=phase_convention,
        frequencies=np.asarray(freqs, dtype=float),
        transmitters=np.zeros((n_tx, 7)),
        transmitter_type=["edipole"] * n_tx,
        transmitter_name=[f"Tx{i}" for i in range(n_tx)],
        receivers=np.zeros((n_rx, 8)),
        receiver_name=[f"Rx{i}" for i in range(n_rx)],
    )
    em.csem.transmitters[:, 0] = tx_x
    em.csem.receivers[:, 0] = rx_x
    rows = []
    for ifreq in range(1, len(freqs) + 1):
        for itx in range(1, n_tx + 1):
            for irx in range(1, n_rx + 1):
                rows.append([21, ifreq, itx, irx, 1.0, 0.1])
                rows.append([22, ifreq, itx, irx, 10.0, 1.0])
    em.data = np.array(rows, dtype=float)
    return em


# ---------------------------------------------------------------------------
# _merge_mt
# ---------------------------------------------------------------------------


def test_merge_mt_src_has_no_mt_section_returns_unchanged():
    out = _mt_em([1.0], [0.0])
    src = EMDataFile()  # src.mt is None
    result = _merge_mt(out, src, np.empty((0, 6)))
    assert result is out


def test_merge_mt_out_mt_is_none_takes_src_directly():
    out = EMDataFile()  # out.mt is None
    out.data = np.empty((0, 6))
    src = _mt_em([1.0, 2.0], [0.0, 100.0])
    result = _merge_mt(out, src, src.data)
    assert result.mt is src.mt
    assert len(result.data) == len(src.data)


def test_merge_mt_full_merge_with_frequency_and_receiver_dedup():
    out = _mt_em([1.0, 10.0], [0.0, 1000.0], name_prefix="A")
    src = _mt_em([10.0, 100.0], [1000.0, 2000.0], name_prefix="B")
    result = _merge_mt(out, src, src.data)
    assert len(result.mt.frequencies) == 3
    assert len(result.mt.receivers) == 3
    # receiver_name must stay positionally aligned with receivers: the
    # duplicate row (y=1000) keeps out's own name "A001", not src's
    # "B000", and the genuinely-new row (y=2000) appends "B001".
    assert result.mt.receiver_name == ["A000", "A001", "B001"]


# ---------------------------------------------------------------------------
# _merge_csem
# ---------------------------------------------------------------------------


def test_merge_csem_src_has_no_csem_section_returns_unchanged():
    out = _csem_em([1.0], [0.0], [100.0])
    src = EMDataFile()
    result = _merge_csem(out, src, np.empty((0, 6)), keep_duplicate_rx=False)
    assert result is out


def test_merge_csem_out_csem_is_none_takes_src_directly():
    out = EMDataFile()
    out.data = np.empty((0, 6))
    src = _csem_em([1.0], [0.0], [100.0, 200.0])
    result = _merge_csem(out, src, src.data, keep_duplicate_rx=False)
    assert result.csem is src.csem


def test_merge_csem_phase_convention_mismatch_raises():
    out = _csem_em([1.0], [0.0], [100.0], phase_convention="lag")
    src = _csem_em([1.0], [0.0], [100.0], phase_convention="lead")
    with pytest.raises(ValueError, match="phase convention"):
        _merge_csem(out, src, src.data, keep_duplicate_rx=False)


def test_merge_csem_dedup_transmitters_and_receivers():
    out = _csem_em([1.0], [0.0, 500.0], [100.0, 200.0])
    src = _csem_em([1.0], [500.0, 1000.0], [200.0, 300.0])  # shares tx@500, rx@200
    result = _merge_csem(out, src, src.data, keep_duplicate_rx=False)
    assert len(result.csem.transmitters) == 3
    assert len(result.csem.receivers) == 3
    assert set(result.csem.transmitter_type) == {"edipole"}
    # names must stay positionally aligned with transmitters/receivers
    # (same reasoning as the MT case above): out's own 2 names are kept
    # unchanged (index 1 is the row that coincides with src's tx@500),
    # and exactly 1 new name is appended for the genuinely-new tx@1000.
    assert result.csem.transmitter_name[:2] == ["Tx0", "Tx1"]
    assert len(result.csem.transmitter_name) == 3
    assert len(result.csem.receiver_name) == 3


def test_merge_csem_keep_duplicate_rx_true_does_not_dedup():
    out = _csem_em([1.0], [0.0], [100.0, 200.0])
    src = _csem_em([1.0], [0.0], [100.0, 200.0])  # identical receivers
    result = _merge_csem(out, src, src.data, keep_duplicate_rx=True)
    assert len(result.csem.receivers) == 4


def test_merge_csem_bdipole_type_roundtrips():
    out = _csem_em([1.0], [0.0], [100.0])
    out.csem.transmitter_type = ["bdipole"]
    src = _csem_em([1.0], [500.0], [200.0])
    src.csem.transmitter_type = ["bdipole"]
    result = _merge_csem(out, src, src.data, keep_duplicate_rx=False)
    assert set(result.csem.transmitter_type) == {"bdipole"}


# ---------------------------------------------------------------------------
# merge_emdata / merge_data_files
# ---------------------------------------------------------------------------


def test_merge_emdata_default_comment_generated():
    em1 = _mt_em([1.0], [0.0])
    em2 = _mt_em([2.0], [100.0])
    merged = merge_emdata([em1, em2])
    assert "Merged 2 data files with pycsamt" in merged.comment


def test_merge_emdata_custom_comment_used():
    em1 = _mt_em([1.0], [0.0])
    em2 = _mt_em([2.0], [100.0])
    merged = merge_emdata([em1, em2], comment="my custom merge")
    assert merged.comment == "my custom merge"


def test_merge_emdata_skips_src_with_no_data():
    em1 = _mt_em([1.0], [0.0])
    em2 = EMDataFile()
    em2.data = np.empty((0, 6))
    merged = merge_emdata([em1, em2])
    assert merged.n_data == em1.n_data


def test_merge_emdata_raises_for_dc_data():
    em1 = _mt_em([1.0], [0.0])
    em2 = _mt_em([2.0], [100.0])
    em2.dc = DCConfig()
    with pytest.raises(NotImplementedError, match="DC resistivity"):
        merge_emdata([em1, em2])


def test_merge_emdata_mt_and_csem_together():
    em1 = _mt_em([1.0], [0.0])
    em2 = _csem_em([1.0], [0.0], [100.0])
    merged = merge_emdata([em1, em2])
    assert merged.mt is not None
    assert merged.csem is not None
    assert merged.n_data == em1.n_data + em2.n_data


def test_merge_emdata_requires_two_files():
    em = _mt_em([1.0], [0.0])
    with pytest.raises(ValueError, match="least two"):
        merge_emdata([em])


def test_merge_data_files_requires_two_paths(tmp_path):
    with pytest.raises(ValueError, match="least two"):
        merge_data_files([tmp_path / "only_one.emdata"], tmp_path / "out.emdata")


def test_merge_data_files_roundtrip(tmp_path):
    from pycsamt.models.mare2dem.iotools.emdata import write_emdata

    em1 = _mt_em([1.0], [0.0, 100.0])
    em2 = _mt_em([2.0], [100.0, 200.0])
    write_emdata(em1, tmp_path / "a.emdata")
    write_emdata(em2, tmp_path / "b.emdata")
    out = tmp_path / "merged.emdata"
    merged = merge_data_files([tmp_path / "a.emdata", tmp_path / "b.emdata"], out)
    assert out.exists()
    assert len(merged.mt.frequencies) == 2
