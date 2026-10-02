# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Desktop → Agent Master survey hand-off (``?handoff=<token>``)."""

from __future__ import annotations

import glob
import json
from pathlib import Path

import pytest

pytest.importorskip("dash", reason="dash required")

from dash.exceptions import PreventUpdate

from pycsamt.app.agent_master import _handoff as ho
from pycsamt.app.agent_master._ids import IDs

_ROOT = Path(__file__).resolve().parents[4]
_L22 = sorted(glob.glob(str(_ROOT / "data" / "AMT" / "WILLY_DATA" / "L22PLT"
                            / "*.edi")))[:3]
_L18 = sorted(glob.glob(str(_ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
                            / "*.edi")))[:2]
pytestmark = pytest.mark.skipif(len(_L22) < 3 or len(_L18) < 2,
                                reason="Baohuashan data missing")


@pytest.fixture
def root(tmp_path, monkeypatch):
    monkeypatch.setenv("PYCSAMT_AGENT_HANDOFF_DIR", str(tmp_path))
    return tmp_path


@pytest.fixture(scope="module")
def sites():
    from pycsamt.site.base import to_sites

    return to_sites(_L22 + _L18)


def _lines(sites):
    return {s.name: ("L" + s.name.split("-")[0]) for s in sites}


def test_write_then_read_round_trip(root, sites):
    token = ho.write_handoff(sites, _lines(sites), edited=True,
                             label="WILLY")
    session = ho.read_handoff(token)
    assert session["n_edi"] == 5
    assert set(session["groups"]) == {"L18", "L22"}
    assert len(session["groups"]["L22"]) == 3
    assert all(Path(f).exists() for fs in session["groups"].values()
               for f in fs)
    assert session["edited"] is True and session["mode"] == "desktop"
    assert Path(session["path"]) == root / token


def test_written_edis_are_the_in_memory_survey(root, sites):
    """The hand-off writes what the desktop holds, not the source files."""
    from pycsamt.site.base import to_sites

    token = ho.write_handoff(sites, _lines(sites))
    f = ho.read_handoff(token)["groups"]["L22"][0]
    assert not f.startswith(str(_ROOT / "data"))
    back = to_sites([f])
    assert len(list(back)) == 1


def test_stations_without_line_go_to_survey(root, sites):
    token = ho.write_handoff(sites)
    assert list(ho.read_handoff(token)["groups"]) == ["Survey"]


def test_only_recent_handoffs_are_kept(root, sites):
    toks = [ho.write_handoff(sites, keep=2) for _ in range(4)]
    kept = sorted(p.name for p in root.iterdir())
    assert kept == sorted(toks[-2:])


def test_rejects_bad_tokens(root, sites):
    assert ho.read_handoff("") is None
    assert ho.read_handoff("../../etc") is None
    assert ho.read_handoff("0" * 20) is None  # well-formed but unknown
    assert ho.token_from_search("?handoff=../x") == ""
    assert ho.token_from_search("?other=1") == ""
    tok = ho.write_handoff(sites)
    assert ho.token_from_search(f"?handoff={tok}&x=1") == tok


def test_empty_survey_refused(root):
    with pytest.raises(ValueError):
        ho.write_handoff([])


def _handoff_fn(agent_app):
    entry = next(v for k, v in agent_app.callback_map.items()
                 if IDs.WELCOME_NOTE in k)
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def test_agent_master_starts_on_the_handoff(agent_app, root, sites):
    fn = _handoff_fn(agent_app)
    with pytest.raises(PreventUpdate):  # no token: the usual welcome
        fn("")
    token = ho.write_handoff(sites, _lines(sites), edited=True)
    store, badge_cls, badge, splash, note = fn(f"?handoff={token}")
    assert store["n_edi"] == 5 and set(store["groups"]) == {"L18", "L22"}
    assert "visible" in badge_cls and "desktop" in badge
    assert "wlc-gone" in splash  # the welcome splash is skipped
    text = json.dumps(note.to_plotly_json(), default=str)
    assert "5 stations" in text and "edited in the desktop" in text
