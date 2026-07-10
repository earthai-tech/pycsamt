# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations

from pycsamt.seg.components import ComponentsMixin


class _Host(ComponentsMixin):
    """Minimal host that exposes a sections dict and the
    trio add/get/has required by ComponentsMixin."""
    def __init__(self) -> None:
        self.sections: dict[str, object] = {}

    def add_section(self, key: str, obj: object) -> None:
        self.sections[str(key).lower()] = obj

    def get_section(self, key: str) -> object | None:
        return self.sections.get(str(key).lower())

    def has_section(self, key: str) -> bool:
        return str(key).lower() in self.sections


class _WList:
    def __init__(self, tag: str) -> None:
        self.tag = tag

    def write(self) -> list[str]:
        return [f">{self.tag}\n", "  X=1\n"]


class _WStr:
    def __init__(self, tag: str) -> None:
        self.tag = tag

    def write(self) -> str:
        return f">{self.tag}\n  Y=2\n"


class _WBad:
    def write(self) -> list[str]:
        raise RuntimeError("boom")


class _NoWrite:
    pass


def test_basic_set_get_drop_and_presence() -> None:
    h = _Host()
    assert h.cget("head") is None
    assert not h.chas("head")

    h.cset("head", {"dataid": "S1"})
    assert h.chas("head")
    assert h.cget("HEAD") == {"dataid": "S1"}  # key norm

    h.cdrop("HEAD")
    assert not h.chas("head")
    assert h.cget("head", default=123) == 123


def test_typed_accessors_core_and_alias() -> None:
    h = _Host()

    h.set_head({"dataid": "ST01"})
    h.set_info({"project": "P"})
    h.set_definemeasurement({"operator": "OP"})

    # spectra/time-series typed setters
    h.set_spectra_sect({"sectid": "SS1"})
    h.set_spectra_io({"n": 3})
    h.set_spectra({"blocks": [1, 2]})
    h.set_timeseries_sect({"sectid": "TS1"})
    h.set_timeseries_io({"m": 2})
    h.set_timeseries({"ids": ["HX"]})

    assert h.get_head()["dataid"] == "ST01"
    assert h.get_info()["project"] == "P"
    assert h.get_definemeasurement()["operator"] == "OP"

    # spectra getters should see what we stored
    assert h.get_spectra_sect()["sectid"] == "SS1"
    assert h.get_spectra_io()["n"] == 3
    assert h.get_spectra()["blocks"] == [1, 2]

    # time-series getters
    assert h.get_timeseries_sect()["sectid"] == "TS1"
    assert h.get_timeseries_io()["m"] == 2
    assert h.get_timeseries()["ids"] == ["HX"]


def test_definemeasurement_fallback_definemeas() -> None:
    h = _Host()
    # store under the short key, fetch via canonical getter
    h.cset("definemeas", {"ok": True})
    dm = h.get_definemeasurement()
    assert isinstance(dm, dict) and dm["ok"] is True


def test_snapshot_selected_and_all() -> None:
    h = _Host()
    h.cset("head", 1)
    h.cset("info", 2)
    h.cset("other", 3)

    sub = h.snapshot(keys=["HEAD", "other", "missing"])
    assert sub == {"head": 1, "other": 3, "missing": None}

    allsnap = h.snapshot()
    # should include all currently registered keys
    assert all(k in allsnap for k in ("head", "info", "other"))


def test_compose_headers_from_order_and_tolerance() -> None:
    h = _Host()
    h.cset("head", _WList("HEAD"))
    h.cset("info", _WStr("INFO"))
    h.cset("bad", _WBad())       # should be tolerated
    h.cset("noval", _NoWrite())  # ignored (no write())

    out = h.compose_headers_from(
        keys_order=["info", "head", "noval", "bad"]
    )
    # INFO should come before HEAD due to keys_order
    pos_info = out.find(">INFO")
    pos_head = out.find(">HEAD")
    assert pos_info != -1 and pos_head != -1
    assert pos_info < pos_head

    # bad writer is skipped; noval ignored
    assert "boom" not in out
    assert ">NOVAL" not in out

    # content from writers is present
    assert "Y=2" in out  # from _WStr
    assert "X=1" in out  # from _WList


def test_alias_getters_when_only_section_key_exists() -> None:
    h = _Host()
    # Simulate how iter_sections may store spectra under "spectra"
    h.cset("spectra", {"sectid": "S-ALIAS"})
    # get_spectra_sect falls back to "spectra" if "spectra_sect" missing
    ssect = h.get_spectra_sect()
    assert isinstance(ssect, dict) and ssect["sectid"] == "S-ALIAS"

    h.cset("timeseries", {"sectid": "T-ALIAS"})
    tsect = h.get_timeseries_sect()
    assert isinstance(tsect, dict) and tsect["sectid"] == "T-ALIAS"
