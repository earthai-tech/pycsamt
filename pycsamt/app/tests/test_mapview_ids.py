# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Unit tests for pycsamt.app.mapview._ids — component ID constants."""

from __future__ import annotations


class TestIDs:
    def test_all_values_are_non_empty_strings(self):
        from pycsamt.app.mapview._ids import IDs

        values = [
            v
            for k, v in vars(IDs).items()
            if not k.startswith("_")
        ]
        assert values, "IDs class should expose at least one constant"
        assert all(isinstance(v, str) and v for v in values)

    def test_all_ids_are_unique(self):
        from pycsamt.app.mapview._ids import IDs

        values = [
            v
            for k, v in vars(IDs).items()
            if not k.startswith("_")
        ]
        assert len(values) == len(set(values))

    def test_all_ids_share_mv_prefix(self):
        from pycsamt.app.mapview._ids import IDs

        values = [
            v
            for k, v in vars(IDs).items()
            if not k.startswith("_")
        ]
        assert all(v.startswith("mv-") for v in values)

    def test_expected_core_ids_present(self):
        from pycsamt.app.mapview._ids import IDs

        for attr in (
            "ROOT",
            "SESSION_ID",
            "STORE_DATA",
            "STORE_THEME",
            "STORE_VIEW",
            "CANVAS_GRAPH",
            "MODAL_LOAD",
            "BTN_LOAD",
        ):
            assert hasattr(IDs, attr), f"IDs.{attr} missing"
