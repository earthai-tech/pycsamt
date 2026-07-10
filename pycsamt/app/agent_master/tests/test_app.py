# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Smoke tests for agent_master."""

from __future__ import annotations

import pytest


def test_import():
    from pycsamt.app.agent_master import (
        create_app,
        launch,
    )
    assert callable(create_app)
    assert callable(launch)


def test_create_app():
    pytest.importorskip("dash")
    pytest.importorskip(
        "dash_bootstrap_components"
    )
    from pycsamt.app.agent_master import (
        create_app,
    )
    app = create_app(debug=False)
    assert app is not None
    assert app.layout is not None


def test_ids_unique():
    from pycsamt.app.agent_master._ids import IDs
    vals = [
        v
        for k, v in vars(IDs).items()
        if not k.startswith("_")
        and isinstance(v, str)
    ]
    assert len(vals) == len(set(vals)), (
        "Duplicate component IDs detected"
    )


def test_edi_detect_from_ids():
    from pycsamt.app.agent_master.callbacks.edi import (
        _detect_from_ids,
    )
    ids = [
        "22-001A",
        "22-002A",
        "26-001A",
        "L30-001",
    ]
    groups = _detect_from_ids(ids)
    assert "L22" in groups
    assert "L26" in groups
    assert len(groups["L22"]) == 2


def test_edi_folder_to_lines_empty(
    tmp_path,
):
    from pycsamt.app.agent_master.callbacks.edi import (
        _folder_to_lines,
    )
    # empty folder -> empty dict
    groups = _folder_to_lines(str(tmp_path))
    assert isinstance(groups, dict)


def test_settings_cfg_roundtrip(tmp_path):
    from unittest.mock import patch

    from pycsamt.app.agent_master.callbacks\
            .settings import _load_cfg, _save_cfg

    fake = tmp_path / "agent_master.json"
    cfg = {
        "provider": "claude",
        "export_fmt": "svg",
    }
    with patch(
        "pycsamt.app.agent_master"
        ".callbacks.settings._CFG_FILE",
        fake,
    ), patch(
        "pycsamt.app.agent_master"
        ".callbacks.settings._CFG_DIR",
        tmp_path,
    ):
        _save_cfg(cfg)
        loaded = _load_cfg()

    assert loaded["provider"] == "claude"
    assert loaded["export_fmt"] == "svg"
