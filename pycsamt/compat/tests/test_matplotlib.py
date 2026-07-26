from __future__ import annotations

import matplotlib
from packaging.version import parse

from pycsamt.compat import matplotlib as mpl_compat


def test_installs_removed_cm_get_cmap(monkeypatch):
    """Modern Matplotlib still supports legacy optional GUI callers."""

    monkeypatch.setattr(mpl_compat, "_MPL_VERSION", parse("3.11"))
    monkeypatch.delattr(matplotlib.cm, "get_cmap", raising=False)

    mpl_compat._install_legacy_cm_get_cmap()

    cmap = matplotlib.cm.get_cmap("viridis", lut=4)
    assert cmap.name == "viridis"
    assert cmap.N == 4


def test_does_not_replace_existing_cm_get_cmap(monkeypatch):
    sentinel = object()
    monkeypatch.setattr(matplotlib.cm, "get_cmap", sentinel, raising=False)

    mpl_compat._install_legacy_cm_get_cmap()

    assert matplotlib.cm.get_cmap is sentinel
