from __future__ import annotations

from pycsamt.api import PyCSAMTObject
from pycsamt.loc import (
    BaseLoc,
    Bounds,
    GeoPath,
    Location,
    UTMZone,
)
from pycsamt.property import (
    FieldAliases,
    FileRecognizer,
    TermDefinitions,
)


def test_location_helpers_inherit_pycsamt_object():
    assert issubclass(BaseLoc, PyCSAMTObject)
    assert isinstance(Location(), PyCSAMTObject)
    assert isinstance(UTMZone("31N"), PyCSAMTObject)
    assert isinstance(Bounds(0.0, 0.0, 1.0, 1.0), PyCSAMTObject)
    assert isinstance(GeoPath([(0.0, 0.0), (0.0, 1.0)]), PyCSAMTObject)


def test_property_helpers_inherit_pycsamt_object():
    assert isinstance(TermDefinitions(), PyCSAMTObject)
    assert isinstance(FieldAliases(), PyCSAMTObject)
    assert isinstance(FileRecognizer(), PyCSAMTObject)


def test_session_exports_are_lazy_at_package_root():
    import pycsamt

    assert pycsamt.Session.__name__ == "Session"
    assert pycsamt.work_session.__name__ == "work_session"
    assert pycsamt.Normalize.__name__ == "Normalize"
    assert pycsamt.normalize_session.__name__ == "normalize_session"
