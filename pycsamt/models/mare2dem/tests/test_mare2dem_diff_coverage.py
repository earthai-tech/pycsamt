# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage tests for models/mare2dem/diff.py.

Builds small synthetic ``.resistivity`` files with ``write_resistivity``
so the test does not depend on the (gitignored, not present in CI)
bundled ``data/mare2dem`` example directory.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.models.mare2dem.diff import diff_resistivity
from pycsamt.models.mare2dem.iotools.resistivity import (
    ResistivityFile,
    read_resistivity,
    write_resistivity,
)


def _make_resistivity_file(path, resistivity, *, poly_file="demo.poly"):
    rf = ResistivityFile(
        resistivity_file=str(path),
        poly_file=poly_file,
        data_file="demo.emdata",
        settings_file="demo.settings",
        anisotropy="isotropic",
    )
    n = len(resistivity)
    rf.resistivity = np.asarray(resistivity, dtype=float).reshape(n, 1)
    rf.free_parameter = np.ones((n, 1), dtype=float)
    rf.bounds = np.tile([-2.0, 5.0], (n, 1))
    rf.prejudice = np.zeros((n, 2), dtype=float)
    write_resistivity(rf, path)
    return rf


@pytest.fixture
def two_models(tmp_path):
    f1 = tmp_path / "iter00.resistivity"
    f2 = tmp_path / "iter20.resistivity"
    _make_resistivity_file(f1, [10.0, 100.0, 1000.0])
    _make_resistivity_file(f2, [20.0, 50.0, 500.0])
    return f1, f2


def test_diff_resistivity_default_log10(two_models, tmp_path):
    f1, f2 = two_models
    out = tmp_path / "diff.resistivity"
    result = diff_resistivity(f1, f2, out)

    assert isinstance(result, ResistivityFile)
    assert out.exists()
    assert result.num_regions == 3

    rf1 = read_resistivity(f1)
    rf2 = read_resistivity(f2)
    expected = np.log10(rf1.resistivity) - np.log10(rf2.resistivity)
    np.testing.assert_allclose(result.resistivity, expected, rtol=1e-5)

    # header/config fields copied from file1, not recomputed
    assert result.poly_file == rf1.poly_file
    assert result.data_file == rf1.data_file
    assert result.settings_file == rf1.settings_file
    assert result.version == rf1.version
    assert result.anisotropy == rf1.anisotropy
    assert result.target_misfit == rf1.target_misfit
    assert result.max_iterations == rf1.max_iterations
    np.testing.assert_allclose(result.global_bounds, rf1.global_bounds)
    np.testing.assert_allclose(result.free_parameter, rf1.free_parameter)

    reread = read_resistivity(out)
    np.testing.assert_allclose(reread.resistivity, expected, rtol=1e-4)


def test_diff_resistivity_custom_diff_fn(two_models, tmp_path):
    f1, f2 = two_models
    out = tmp_path / "pct_change.resistivity"
    result = diff_resistivity(
        f1,
        f2,
        out,
        diff_fn=lambda A, B: np.abs((A - B) / A * 100),
    )
    rf1 = read_resistivity(f1)
    rf2 = read_resistivity(f2)
    expected = np.abs((rf1.resistivity - rf2.resistivity) / rf1.resistivity * 100)
    np.testing.assert_allclose(result.resistivity, expected, rtol=1e-5)


def test_diff_resistivity_linear_difference(two_models, tmp_path):
    f1, f2 = two_models
    out = tmp_path / "linear_diff.resistivity"
    result = diff_resistivity(f1, f2, out, diff_fn=lambda A, B: A - B)
    rf1 = read_resistivity(f1)
    rf2 = read_resistivity(f2)
    np.testing.assert_allclose(
        result.resistivity, rf1.resistivity - rf2.resistivity, rtol=1e-5
    )


def test_diff_resistivity_mismatched_regions_raises(tmp_path):
    f1 = tmp_path / "a.resistivity"
    f2 = tmp_path / "b.resistivity"
    _make_resistivity_file(f1, [10.0, 100.0])
    _make_resistivity_file(f2, [10.0, 100.0, 1000.0])
    out = tmp_path / "diff.resistivity"
    with pytest.raises(ValueError, match="Cannot diff"):
        diff_resistivity(f1, f2, out)


def test_diff_resistivity_none_bounds_and_prejudice_default_empty(
    tmp_path, monkeypatch
):
    f1 = tmp_path / "iter00.resistivity"
    f2 = tmp_path / "iter20.resistivity"
    _make_resistivity_file(f1, [10.0, 100.0])
    _make_resistivity_file(f2, [20.0, 50.0])
    # simulate a header-only read leaving bounds/prejudice None
    read1 = read_resistivity(f1)
    read1.bounds = None
    read1.prejudice = None

    import pycsamt.models.mare2dem.diff as diff_mod

    original_read = diff_mod.read_resistivity

    def _patched_read(path):
        if str(path) == str(f1):
            return read1
        return original_read(path)

    monkeypatch.setattr(diff_mod, "read_resistivity", _patched_read)
    out = tmp_path / "diff.resistivity"
    result = diff_resistivity(f1, f2, out)

    assert result.bounds.size == 0
    assert result.prejudice.size == 0
