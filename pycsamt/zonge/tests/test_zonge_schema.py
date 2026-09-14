# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations

from pycsamt.zonge.schema import (
    ALL_ALIASES,
    QC_ALIASES,
    UNIT_MAP,
    get_aliases,
    get_unit,
)


def test_get_unit_formatted_default():
    assert get_unit("rho") == UNIT_MAP["rho"][1]


def test_get_unit_simple():
    assert get_unit("rho", formatted=False) == UNIT_MAP["rho"][0]


def test_get_unit_unknown_returns_none():
    assert get_unit("not_a_real_name") is None


def test_get_aliases_qc_kind():
    out = get_aliases("pc_emag", kind="qc")
    assert set(out) == set(QC_ALIASES["pc_emag"])


def test_get_aliases_all_kind_default():
    out = get_aliases("rho")
    assert set(out) == set(ALL_ALIASES["rho"])


def test_get_aliases_unknown_canonical_returns_empty_tuple():
    assert get_aliases("no_such_canonical") == ()


def test_get_aliases_normalize_case_insensitive():
    out_lower = get_aliases("rho")
    out_upper = get_aliases("RHO")
    assert out_lower == out_upper


def test_get_aliases_normalize_false_is_case_sensitive():
    out = get_aliases("RHO", normalize=False)
    assert out == ()  # no exact "RHO" key in ALL_ALIASES


def test_get_aliases_custom_aliases_extends_all_map():
    out = get_aliases(
        "rho", custom_aliases={"rho": ("resistivity_ohm_m",)}
    )
    assert "resistivity_ohm_m" in out
    for a in ALL_ALIASES["rho"]:
        assert a in out


def test_get_aliases_custom_aliases_normalized_key():
    out = get_aliases(
        "rho", custom_aliases={"RHO": ("custom_alias",)}
    )
    assert "custom_alias" in out


def test_get_aliases_custom_aliases_extends_qc_map_when_kind_qc():
    out = get_aliases(
        "pc_emag",
        kind="qc",
        custom_aliases={"pc_emag": ("custom_pc_emag",)},
    )
    assert "custom_pc_emag" in out
    for a in QC_ALIASES["pc_emag"]:
        assert a in out


def test_get_aliases_custom_aliases_for_non_qc_canonical_ignored_in_qc_map():
    # "rho" is not a QC canonical name, so passing custom aliases for it
    # while kind="qc" should not create a new QC entry.
    out = get_aliases(
        "rho",
        kind="qc",
        custom_aliases={"rho": ("custom_alias",)},
    )
    assert out == ()


def test_get_aliases_does_not_mutate_module_level_maps():
    before = dict(ALL_ALIASES)
    get_aliases("rho", custom_aliases={"rho": ("temp_alias",)})
    assert ALL_ALIASES == before
    assert "temp_alias" not in ALL_ALIASES.get("rho", ())
