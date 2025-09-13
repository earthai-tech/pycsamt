# -*- coding: utf-8 -*-
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.core import transformers as tr
from pycsamt.seg.collection import EDICollection
from pycsamt.seg.edi import EDIFile
from pycsamt.jones.j import JFile
from pycsamt.jones.collection import JCollection


@pytest.mark.usefixtures("legacy_data_file", "j_single_file")
def test_public_api_all():
    assert set(tr.__all__) == {"AVGtoEDI", "JtoEDI"}


def test_avg_to_edi_true_file(legacy_data_file: Path):
    out = tr.AVGtoEDI().transform(legacy_data_file)
    assert isinstance(out, EDICollection)
    assert len(out) >= 1
    first = next(iter(out))
    assert isinstance(first, EDIFile)
    assert getattr(first.Z, "_freq", None) is not None
    has_z = getattr(first.Z, "_z", None) is not None
    if not has_z:
        assert hasattr(first.Z, "_freq")
    # default ordering is descending
    f = np.asarray(first.Z._freq)
    assert np.all(np.diff(f) <= 0)


def test_avg_to_edi_accepts_str_and_path(legacy_data_file: Path):
    p = legacy_data_file
    # as Path
    coll1 = tr.AVGtoEDI().transform(p)
    # as str
    coll2 = tr.AVGtoEDI().transform(str(p))
    assert isinstance(coll1, EDICollection)
    assert isinstance(coll2, EDICollection)
    assert len(coll1) >= 1 and len(coll2) >= 1


def test_j_to_edi_true_file(j_single_file: Path):
    out = tr.JtoEDI().transform(j_single_file)
    assert isinstance(out, EDIFile)
    assert getattr(out.Z, "_freq", None) is not None
    f = np.asarray(out.Z._freq)
    assert np.all(np.diff(f) <= 0)


def test_j_to_edi_accepts_str_and_path(j_single_file: Path):
    p = j_single_file
    ed1 = tr.JtoEDI().transform(p)
    ed2 = tr.JtoEDI().transform(str(p))
    assert isinstance(ed1, EDIFile)
    assert isinstance(ed2, EDIFile)


def test_j_collection_roundtrip_if_constructible(
    jc_files: list[Path],
):
    paths = jc_files[:2]
    jfs = [JFile.from_file(p) for p in paths]
    coll = None
    try:
        coll = JCollection(items=jfs)
    except TypeError:
        try:
            coll = JCollection(jfs)
        except Exception:
            pytest.skip("Cannot construct JCollection programmatically")
    out = tr.JtoEDI().transform(coll)
    assert isinstance(out, EDICollection)
    assert len(out) == len(jfs)


def test_avg_transform_is_consistent_length(legacy_data_file: Path):
    c1 = tr.AVGtoEDI().transform(legacy_data_file)
    c2 = tr.AVGtoEDI().transform(legacy_data_file)
    n1 = [len(ed.Z._freq) for ed in c1]
    n2 = [len(ed.Z._freq) for ed in c2]
    assert n1 == n2 and all(n > 0 for n in n1)


def test_j_transform_freq_and_z_lengths_match(j_single_file: Path):
    ed = tr.JtoEDI().transform(j_single_file)
    f = np.asarray(ed.Z._freq)
    z = getattr(ed.Z, "_z", None)
    if z is not None:
        assert z.shape[0] == f.size
