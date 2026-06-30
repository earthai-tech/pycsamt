# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for optional map backend boundaries."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.map import overlays
from pycsamt.map.export import write_image


class _KaleidoMissingFigure:
    def write_image(self, *_args, **_kwargs) -> None:
        raise ValueError("Image export requires kaleido")


def test_write_image_explains_missing_kaleido() -> None:
    out = "map.png"
    with pytest.raises(ImportError, match="kaleido"):
        write_image(_KaleidoMissingFigure(), out)


def test_contour_overlay_falls_back_without_scipy(
    monkeypatch,
) -> None:
    monkeypatch.setattr(
        overlays,
        "scipy_griddata",
        lambda: None,
    )
    trace = overlays.build_contour_overlay(
        np.array([0.0, 1.0, 0.0]),
        np.array([0.0, 0.0, 1.0]),
        np.array([10.0, 20.0, 30.0]),
        grid_size=4,
    )
    assert trace.z is not None
