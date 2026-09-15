# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations

import copy

import pytest

from pycsamt.zonge.config import Zonge


class _Sub(Zonge):
    def __init__(self, verbose: bool = False, alpha: int = 1, **kws):
        super().__init__(verbose=verbose, **kws)
        self.alpha = alpha


class _Host(Zonge):
    def __init__(self, verbose: bool = False, child: Zonge | None = None, **kws):
        super().__init__(verbose=verbose, **kws)
        self.child = child if child is not None else _Sub()


def test_init_sets_verbose_and_logger():
    z = Zonge(verbose=True)
    assert z.verbose is True
    assert z._logger is not None


def test_get_param_names_excludes_self_and_var_keyword():
    names = Zonge._get_param_names()
    assert names == ["verbose"]

    names2 = _Sub._get_param_names()
    assert names2 == ["alpha", "verbose"]


def test_get_params_shallow():
    s = _Sub(verbose=True, alpha=5)
    params = s.get_params(deep=False)
    assert params == {"verbose": True, "alpha": 5}


def test_get_params_deep_includes_nested_params():
    host = _Host(verbose=False, child=_Sub(verbose=True, alpha=7))
    params = host.get_params(deep=True)
    assert params["child__verbose"] is True
    assert params["child__alpha"] == 7
    assert params["child"] is host.child


def test_set_params_updates_attributes_and_returns_self():
    s = _Sub(alpha=1)
    out = s.set_params(alpha=99, verbose=True)
    assert out is s
    assert s.alpha == 99
    assert s.verbose is True


def test_set_params_noop_when_no_kwargs():
    s = _Sub(alpha=1)
    out = s.set_params()
    assert out is s
    assert s.alpha == 1


def test_set_params_raises_for_unknown_param():
    s = _Sub()
    with pytest.raises(ValueError):
        s.set_params(bogus_param=1)


def test_shallow_copy_shares_nested_objects():
    child = _Sub(alpha=3)
    host = _Host(child=child)
    dup = copy.copy(host)
    assert dup is not host
    assert dup.child is host.child  # shallow: same nested object


def test_deep_copy_duplicates_nested_objects():
    child = _Sub(alpha=3)
    host = _Host(child=child)
    dup = copy.deepcopy(host)
    assert dup is not host
    assert dup.child is not host.child
    assert dup.child.alpha == host.child.alpha


def test_str_hides_private_attrs():
    s = _Sub(verbose=True, alpha=42)
    text = str(s)
    assert "alpha=42" in text
    assert "verbose=True" in text
    assert "_logger" not in text


def test_repr_shows_class_name_and_verbose():
    z = Zonge(verbose=True)
    assert repr(z) == "Zonge(verbose=True)"
