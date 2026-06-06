"""Tests for common pyCSAMT API object helpers."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pycsamt.api import MetadataMixin, PyCSAMTObject
from pycsamt.core.base import CoreObject


def test_pycsamt_object_repr_excludes_logger_and_verbose():
    """Automatic repr should skip noisy runtime attributes."""

    class Obj(PyCSAMTObject):
        def __init__(self):
            self.name = "S001"
            self.verbose = 3
            self.logger = object()
            self.values = [1, 2, 3]

    text = repr(Obj())

    assert text.startswith("Obj(")
    assert "name='S001'" in text
    assert "values=list(len=3" in text
    assert "verbose" not in text
    assert "logger" not in text


def test_pycsamt_object_summarizes_arrays_and_mappings():
    """Large numerical and mapping values should be compact."""

    class Obj(PyCSAMTObject):
        __repr_fields__ = ("data", "attrs")

        def __init__(self):
            self.data = np.zeros((2, 3))
            self.attrs = {"a": 1, "b": 2, "c": 3, "d": 4}

    text = repr(Obj())

    assert "ndarray(shape=(2, 3), dtype=float64)" in text
    assert "dict(len=4" in text


def test_pycsamt_object_supports_dataclasses():
    """Dataclass fields should be discovered automatically."""

    @dataclass(repr=False)
    class Obj(PyCSAMTObject):
        name: str
        count: int
        verbose: int = 0

    text = repr(Obj("A", 2, verbose=5))

    assert text == "Obj(name='A', count=2)"


def test_pycsamt_object_supports_slots_and_to_dict():
    """Slotted classes should display and serialize public slots."""

    class Obj(PyCSAMTObject):
        __slots__ = ("name", "value", "verbose")

        def __init__(self):
            self.name = "slotty"
            self.value = 4
            self.verbose = 1

    obj = Obj()

    assert repr(obj) == "Obj(name='slotty', value=4)"
    assert obj.to_dict() == {"name": "slotty", "value": 4}


def test_metadata_mixin_attaches_free_form_metadata():
    """MetadataMixin should manage metadata without base coupling."""

    class Obj(MetadataMixin, PyCSAMTObject):
        pass

    obj = Obj()
    obj.update_metadata(project="demo", station="S001")

    assert obj.metadata_dict() == {
        "project": "demo",
        "station": "S001",
    }


def test_core_object_inherits_pycsamt_object_and_exclusions():
    """Existing CoreObject subclasses should inherit the new root."""

    class Obj(CoreObject):
        def __init__(self):
            self.name = "core"
            self.verbose = 2

    obj = Obj()

    assert isinstance(obj, PyCSAMTObject)
    assert repr(obj) == "Obj(name='core')"
