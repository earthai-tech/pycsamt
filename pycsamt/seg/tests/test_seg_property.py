# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

import pytest

from pycsamt.seg.property import (
    PlainMeta,
    References,
    Person,
    Software,
    Source,
    Processing,
    Copyright,
    PropertiesMixin,
)

__all__ = [
    "PlainMeta",
    "References",
    "Person",
    "Software",
    "Source",
    "Processing",
    "Copyright",
    "PropertiesMixin",
]


# ------------------------
# Basic PlainMeta helpers
# ------------------------
def test_plainmeta_update_and_clone_uses_validate():
    p = Person(name="A")
    # update with harmless field (email None so ok)
    p.update(name="B")
    assert p.name == "B"

    # clone with override and independence
    q = p.clone(name="C")
    assert q.name == "C"
    assert p.name == "B"
    # dict-like export
    d = q.to_dict()
    assert isinstance(d, dict) and d["name"] == "C"


def test_plainmeta_validate_called_on_update_error():
    p = Person()
    with pytest.raises(ValueError):
        # invalid email triggers validate
        p.update(email="no-at-symbol")


# ------------------------
# Value objects validation
# ------------------------
def test_references_year_validation():
    r = References(author="X", year=2020)
    r.validate()  # ok
    with pytest.raises(ValueError):
        r.update(year=-1)


def test_person_email_validation():
    ok = Person(email="u@e.com")
    ok.validate()  # ok
    bad = Person()
    with pytest.raises(ValueError):
        bad.update(email="ue.com")


def test_software_validation_fields_and_url():
    s = Software(name="pkg", version="1.0")
    s.validate()  # ok

    with pytest.raises(ValueError):
        Software(name="pkg", version="1.0",
                 url="ftp://x").validate()

    with pytest.raises(ValueError):
        Software(name="", version="1.0").validate()

    with pytest.raises(ValueError):
        Software(name="pkg", version="").validate()


def test_source_defaults_and_validation():
    s = Source()
    # has ISO-like Z suffix
    assert isinstance(s.creationdate, str)
    assert s.creationdate.endswith("Z")
    assert s.creatingsoftware == "pyCSAMT"
    s.validate()  # ok

    with pytest.raises(ValueError):
        s.update(creatingsoftware="").validate()


def test_processing_signconv_and_runlist_validation():
    p = Processing()

    # normalize +i w t
    p.update(signconvention="exp(+i w t)")
    assert p.signconvention == "exp(+i ω t)"

    # normalize -i w t
    p.update(signconvention="exp(-i w t)")
    assert p.signconvention == "exp(-i ω t)"

    # runlist must be list[str]
    with pytest.raises(ValueError):
        p.update(runlist="not-a-list")

    with pytest.raises(ValueError):
        p.update(runlist=[1, 2])  # not str

    p.update(runlist=["run1", "run2"])
    assert p.runlist == ["run1", "run2"]


def test_copyright_release_status_validation():
    c = Copyright()
    for st in ("open", "public", "proprietary",
               "OPEN", "Public"):
        c.update(release_status=st)
        c.validate()

    with pytest.raises(ValueError):
        c.update(release_status="unknown").validate()


# ------------------------
# PropertiesMixin behavior
# ------------------------
class Host(PropertiesMixin):
    """Minimal host that mixes properties in."""
    pass


def test_properties_mixin_ensures_holders_on_init():
    h = Host()
    assert hasattr(h, "Source")
    assert hasattr(h, "Processing")
    assert hasattr(h, "Copyright")
    assert isinstance(h.Source, Source)
    assert isinstance(h.Processing, Processing)
    assert isinstance(h.Copyright, Copyright)


def test_properties_mixin_attach_and_update_properties():
    h = Host()

    # attach explicit holders
    s = Source(project="P")
    pr = Processing(processedby="Ann")
    cp = Copyright(owner="Org")
    h.attach_properties(source=s, processing=pr, copyright=cp)
    assert h.Source.project == "P"
    assert h.Processing.processedby == "Ann"
    assert h.Copyright.owner == "Org"

    # update via dicts and flat keys
    h.update_properties(
        source={"survey": "S1", "sitename": "SiteA"},
        processedby="Bob",  # flat routed to Processing
        processingsoftware={"name": "Tool",
                            "version": "2.0"},
        release_status="open",
    )
    assert h.Source.survey == "S1"
    assert h.Source.sitename == "SiteA"
    assert h.Processing.processedby == "Bob"
    assert h.Processing.ProcessingSoftware.name == "Tool"
    assert h.Processing.ProcessingSoftware.version == "2.0"
    assert h.Copyright.release_status.lower() == "open"


def test_properties_mixin_processingsoftware_name_string():
    h = Host()
    h.update_properties(processingsoftware="PackX")
    assert h.Processing.ProcessingSoftware.name == "PackX"


def test_properties_mixin_properties_as_dict_structure():
    h = Host()
    h.update_properties(
        source={"project": "PX"},
        processedby="Jim",
        release_status="public",
    )
    d = h.properties_as_dict()
    assert set(d.keys()) == {"source", "processing", "copyright"}
    assert d["source"]["project"] == "PX"
    assert d["processing"]["processedby"] == "Jim"
    assert d["copyright"]["release_status"] == "public"


def test_properties_mixin_copy_properties_returns_clones():
    h = Host()
    h.update_properties(
        source={"project": "A"},
        processedby="M",
        owner="Co",
    )
    copies = h.copy_properties()
    copies["source"].update(project="B")
    copies["processing"].update(processedby="N")
    copies["copyright"].update(owner="Co2")

    # original unchanged
    assert h.Source.project == "A"
    assert h.Processing.processedby == "M"
    assert h.Copyright.owner == "Co"


def test_properties_mixin_summary_contains_key_bits():
    h = Host()
    h.update_properties(
        source={"project": "P", "sitename": "S"},
        processedby="X",
        release_status="public",
        owner="Org",
    )
    s = h.property_summary()
    # relaxed checks; presence of markers
    assert "project='P'" in s
    assert "by='X'" in s or "Processing(" in s
    assert "status='public'" in s or "owner='Org'" in s
