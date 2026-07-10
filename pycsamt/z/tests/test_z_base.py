
import logging

import numpy as np
import pytest

from pycsamt.exceptions import ZError
from pycsamt.z.base import (
    BaseEM,
)


class Dummy(BaseEM):
    def __init__(self, name=None, verbose=0):
        super().__init__(name=name, verbose=verbose)


def test_freq_property_ok():
    d = Dummy()
    d.freq = [10.0, 1.0, 0.1]
    assert d.has_freq
    assert d.n_freq == 3
    assert isinstance(d.freq, np.ndarray)


def test_freq_property_invalid_shape():
    d = Dummy()
    with pytest.raises(ZError):
        d.freq = [[1.0, 2.0], [3.0, 4.0]]


def test_freq_property_nonfinite_or_nonpos():
    d = Dummy()
    with pytest.raises(ZError):
        d.freq = [1.0, np.nan]
    with pytest.raises(ZError):
        d.freq = [10.0, -1.0]


def test_n_freq_infers_from_arrays_when_no_freq():
    d = Dummy()
    d._z = np.zeros((3, 2, 2), dtype=complex)
    assert not d.has_freq
    assert d.n_freq == 3


def test_has_errors_flag():
    d = Dummy()
    d._z = np.zeros((2, 2, 2), dtype=complex)
    assert not d.has_errors
    d._z_err = np.zeros_like(d._z, dtype=float)
    assert d.has_errors
    d._z_err = None
    assert not d.has_errors


def test_copy_and_deepcopy_semantics():
    d = Dummy()
    d.freq = [10.0, 1.0]
    arr = np.zeros((2, 2, 2), dtype=complex)
    d._z = arr

    c = d.copy()
    assert c is not d
    assert c._z is arr

    dc = d.deepcopy()
    assert dc is not d
    assert dc._z is not arr
    assert np.array_equal(dc._z, arr)


def test_subset_slices_all_freq_aligned_arrays():
    d = Dummy()
    f = np.array([10.0, 1.0, 0.1], float)
    z = np.arange(3 * 2 * 2, dtype=float).reshape(3, 2, 2)
    tip = np.arange(3 * 1 * 2, dtype=float).reshape(3, 1, 2)

    d.freq = f
    d._z = z.copy()
    d._tipper = tip.copy()
    d._noslice = np.zeros((2,), float)

    s = d.subset([1, 2])
    assert s.n_freq == 2
    assert np.all(s.freq == f[[1, 2]])
    assert s._z.shape == (2, 2, 2)
    assert np.all(s._z == z[[1, 2]])
    assert s._tipper.shape == (2, 1, 2)
    assert np.all(s._tipper == tip[[1, 2]])
    # untouched non-freq-first arrays
    assert s._noslice.shape == (2,)


def test_subset_boolean_mask():
    d = Dummy()
    f = np.array([10.0, 1.0, 0.1], float)
    z = np.arange(3 * 2 * 2, dtype=float).reshape(3, 2, 2)
    d.freq = f
    d._z = z.copy()

    m = np.array([True, False, True])
    s = d.subset(m)
    assert s.n_freq == 2
    assert np.all(s.freq == f[m])
    assert np.all(s._z == z[m])


def test_subset_when_n_is_zero_returns_deepcopy():
    d = Dummy()
    d._aux = np.array([1, 2, 3], int)
    s = d.subset([0, 1])  # no freq → deepcopy
    assert s is not d
    assert np.array_equal(s._aux, d._aux)
    # mutate original and ensure copy unaffected
    d._aux[0] = 99
    assert s._aux[0] != d._aux[0]


def test_select_alias():
    d = Dummy()
    d.freq = [10.0, 1.0, 0.1]
    d._z = np.zeros((3, 2, 2), complex)
    a = d.subset([0, 2])
    b = d.select([0, 2])
    # shallow structural equality checks
    assert a.n_freq == b.n_freq
    assert np.all(a.freq == b.freq)
    assert np.all(a._z == b._z)


def test_validate_shapes_ok_and_error():
    d = Dummy()
    d.freq = [10.0, 1.0, 0.1]
    d._z = np.zeros((3, 2, 2), complex)
    d._phase = np.zeros((3, 2, 2), float)
    # should not raise
    d.validate_shapes()

    d._resistivity = np.zeros((2, 2, 2), float)
    with pytest.raises(ZError):
        d.validate_shapes()


def test_summary_and_repr_include_key_info():
    d = Dummy(name="site01")
    d.freq = [10.0, 1.0]
    d._z = np.zeros((2, 2, 2), complex)
    s = d.summary()
    r = repr(d)
    assert "n_freq: 2" in s
    assert "arrays:" in s
    assert "site01" in s
    assert "n_freq=2" in r


def test_verbose_levels_and_logger_present():
    d = Dummy(verbose=0)
    assert hasattr(d, "log")
    d.set_verbose(2)
    assert d.verbose == 2
    # compare numeric level against DEBUG
    assert d.log.level == logging.DEBUG
    d.set_verbose(1)
    assert d.log.level == logging.INFO
    d.set_verbose(0)
    assert d.log.level == logging.WARNING
