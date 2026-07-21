# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Unit tests for the vendored :mod:`pycsamt.utils.version` module (PEP 440
version parsing/comparison, lineage of ``packaging.version``).
"""

from __future__ import annotations

import pytest

from pycsamt.utils.version import (
    Infinity,
    InfinityType,
    InvalidVersion,
    LegacyVersion,
    NegativeInfinity,
    NegativeInfinityType,
    Version,
    _cmpkey,
    _legacy_cmpkey,
    _parse_letter_version,
    _parse_local_version,
    _parse_version_parts,
    parse,
)

# --------------------------------- parse --------------------------------


def test_parse_valid_pep440_returns_version():
    v = parse("1.2.3")
    assert isinstance(v, Version)
    assert not isinstance(v, LegacyVersion)


@pytest.mark.parametrize(
    "raw",
    ["not-a-version", "1.0-beta-that-is-weird", "", "abc", "1.0.0.0.0-x-y"],
)
def test_parse_invalid_pep440_returns_legacy_version(raw):
    with pytest.warns(DeprecationWarning):
        v = parse(raw)
    assert isinstance(v, LegacyVersion)
    assert str(v) == raw


# ------------------------- Version: component parsing -------------------


def test_version_epoch_default_and_explicit():
    assert Version("2.0").epoch == 0
    assert Version("1!2.0").epoch == 1


def test_version_release_tuple():
    assert Version("1").release == (1,)
    assert Version("1.0.0").release == (1, 0, 0)
    assert Version("1.2.3.4").release == (1, 2, 3, 4)


@pytest.mark.parametrize(
    "raw, expected_letter, expected_number",
    [
        ("1.0a1", "a", 1),
        ("1.0alpha1", "a", 1),
        ("1.0-alpha1", "a", 1),
        ("1.0b2", "b", 2),
        ("1.0beta2", "b", 2),
        ("1.0c1", "rc", 1),
        ("1.0rc1", "rc", 1),
        ("1.0pre1", "rc", 1),
        ("1.0preview1", "rc", 1),
        ("1.0a", "a", 0),  # implicit 0
    ],
)
def test_version_prerelease_normalizes_spelling(
    raw, expected_letter, expected_number
):
    v = Version(raw)
    assert v.pre == (expected_letter, expected_number)
    assert v.is_prerelease is True
    assert v.is_devrelease is False


@pytest.mark.parametrize(
    "raw, expected_post",
    [
        ("1.0.post1", 1),
        ("1.0-r1", 1),
        ("1.0-rev1", 1),
        ("1.0-post1", 1),
        ("1.0-1", 1),  # implicit post via bare "-N"
    ],
)
def test_version_postrelease_forms(raw, expected_post):
    v = Version(raw)
    assert v.post == expected_post
    assert v.is_postrelease is True


def test_version_no_postrelease_is_none():
    v = Version("1.0")
    assert v.post is None
    assert v.is_postrelease is False


@pytest.mark.parametrize(
    "raw, expected_dev",
    [
        ("1.0.dev1", 1),
        ("1.0-dev1", 1),
        ("1.0dev1", 1),
        ("1.0.dev", 0),  # implicit 0
    ],
)
def test_version_devrelease_forms(raw, expected_dev):
    v = Version(raw)
    assert v.dev == expected_dev
    assert v.is_devrelease is True
    assert v.is_prerelease is True  # dev counts as prerelease


def test_version_no_devrelease_is_none():
    assert Version("1.0").dev is None
    assert Version("1.0").is_devrelease is False


@pytest.mark.parametrize(
    "raw, expected_local",
    [
        ("1.0+local.123", "local.123"),
        ("1.0+abc_1", "abc.1"),
        ("1.0+abc-def.5", "abc.def.5"),
    ],
)
def test_version_local_segment(raw, expected_local):
    v = Version(raw)
    assert v.local == expected_local


def test_version_no_local_is_none():
    assert Version("1.0").local is None


# --------------------- Version: str/repr round-tripping ------------------


@pytest.mark.parametrize(
    "raw, normalized",
    [
        ("1.0", "1.0"),
        ("1!2.0", "1!2.0"),
        ("1.0-alpha1", "1.0a1"),
        ("1.0b2", "1.0b2"),
        ("1.0rc1", "1.0rc1"),
        ("1.0-r1", "1.0.post1"),
        ("1.0.post1", "1.0.post1"),
        ("1.0.dev1", "1.0.dev1"),
        ("1.0+local.123", "1.0+local.123"),
        ("1.0+abc_1", "1.0+abc.1"),
    ],
)
def test_version_str_normalizes(raw, normalized):
    assert str(Version(raw)) == normalized


def test_version_repr():
    assert repr(Version("1.0")) == "<Version('1.0')>"


# ------------------- Version: public / base_version / flags --------------


def test_version_public_strips_local():
    v = Version("1.0+local.123")
    assert v.public == "1.0"


def test_version_public_without_local_equals_str():
    v = Version("1.0.post1")
    assert v.public == str(v)


def test_version_base_version_strips_everything_but_epoch_release():
    v = Version("1!1.0.post1.dev1+local")
    assert v.base_version == "1!1.0"


def test_version_base_version_omits_epoch_when_zero():
    v = Version("1.0.post1.dev1+local")
    assert v.base_version == "1.0"


def test_version_is_prerelease_combinations():
    assert Version("1.0a1").is_prerelease is True
    assert Version("1.0.dev1").is_prerelease is True
    assert Version("1.0").is_prerelease is False
    assert Version("1.0.post1").is_prerelease is False


def test_version_major_minor_micro():
    v = Version("1.2.3")
    assert (v.major, v.minor, v.micro) == (1, 2, 3)
    assert Version("1").minor == 0
    assert Version("1").micro == 0
    assert Version("1.2").micro == 0


# --------------------------- Version: InvalidVersion ----------------------


@pytest.mark.parametrize(
    "raw",
    ["not-a-version", "", "abc", "1.0-beta-that-is-weird", "1.0.0-"],
)
def test_version_invalid_raises(raw):
    with pytest.raises(InvalidVersion):
        Version(raw)


# ------------------------ Version: comparison operators -------------------

# Canonical PEP 440 ordering (strictly increasing).
_ORDERED_VERSIONS = [
    "1.0.dev456",
    "1.0a1",
    "1.0a2.dev456",
    "1.0a12.dev456",
    "1.0a12",
    "1.0b1.dev456",
    "1.0b2",
    "1.0b2.post345.dev456",
    "1.0b2.post345",
    "1.0rc1.dev456",
    "1.0rc1",
    "1.0",
    "1.0+abc.5",
    "1.0+abc.7",
    "1.0+5",
    "1.0.post456.dev34",
    "1.0.post456",
    "1.1.dev1",
]


@pytest.mark.parametrize(
    "lower, higher",
    list(zip(_ORDERED_VERSIONS, _ORDERED_VERSIONS[1:])),
)
def test_version_ordering_adjacent_pairs(lower, higher):
    lo, hi = Version(lower), Version(higher)
    assert lo < hi
    assert lo <= hi
    assert hi > lo
    assert hi >= lo
    assert lo != hi
    assert not (lo == hi)


def test_version_full_chain_is_sorted():
    versions = [Version(v) for v in _ORDERED_VERSIONS]
    shuffled = list(reversed(versions))
    assert sorted(shuffled) == versions


def test_version_equality_ignores_trailing_release_zeros():
    assert Version("1.0.0") == Version("1.0")
    assert Version("1.0.0") == Version("1")
    assert hash(Version("1.0.0")) == hash(Version("1.0"))


def test_version_equality_and_inequality():
    assert Version("1.0") == Version("1.0")
    assert Version("1.0") != Version("2.0")
    assert not (Version("1.0") == "1.0")  # not a _BaseVersion -> NotImplemented -> False


def test_version_comparison_with_non_baseversion_is_notimplemented():
    v = Version("1.0")
    with pytest.raises(TypeError):
        v < "1.0"  # noqa: B015 - comparison intentionally triggers TypeError
    with pytest.raises(TypeError):
        v <= "1.0"  # noqa: B015
    with pytest.raises(TypeError):
        v > "1.0"  # noqa: B015
    with pytest.raises(TypeError):
        v >= "1.0"  # noqa: B015
    assert (v == "1.0") is False
    assert (v != "1.0") is True


def test_version_hashable_in_set_and_dict():
    versions = {Version("1.0"), Version("1.0.0"), Version("2.0")}
    assert len(versions) == 2  # "1.0" and "1.0.0" collide

    mapping = {Version("1.0"): "first"}
    mapping[Version("1.0.0")] = "second"
    assert mapping[Version("1.0")] == "second"


# ------------------------------ LegacyVersion -----------------------------


def test_legacy_version_construction_warns_and_stores_raw():
    with pytest.warns(DeprecationWarning):
        lv = LegacyVersion("not-a-version-at-all")
    assert str(lv) == "not-a-version-at-all"
    assert lv.public == "not-a-version-at-all"
    assert lv.base_version == "not-a-version-at-all"
    assert repr(lv) == "<LegacyVersion('not-a-version-at-all')>"


def test_legacy_version_attribute_defaults():
    with pytest.warns(DeprecationWarning):
        lv = LegacyVersion("weird!!version")
    assert lv.epoch == -1
    assert lv.release is None
    assert lv.pre is None
    assert lv.post is None
    assert lv.dev is None
    assert lv.local is None
    assert lv.is_prerelease is False
    assert lv.is_postrelease is False
    assert lv.is_devrelease is False


def test_legacy_version_equals_itself():
    with pytest.warns(DeprecationWarning):
        lv1 = LegacyVersion("foo-bar")
        lv2 = LegacyVersion("foo-bar")
    assert lv1 == lv2
    assert hash(lv1) == hash(lv2)


def test_legacy_version_sorts_before_any_pep440_version():
    with pytest.warns(DeprecationWarning):
        lv = LegacyVersion("not-a-version")
    v = Version("0.0.0.0.0.1")  # smallest plausible real version
    assert lv < v
    assert v > lv
    assert lv != v


def test_legacy_version_ordering_between_two_legacy():
    with pytest.warns(DeprecationWarning):
        lv_a = LegacyVersion("1.0")
        lv_b = LegacyVersion("1.1")
    assert lv_a < lv_b


# ------------------------- InfinityType / NegativeInfinityType -----------


def test_infinity_repr_and_hash():
    assert repr(Infinity) == "Infinity"
    assert repr(NegativeInfinity) == "-Infinity"
    assert hash(Infinity) == hash(repr(Infinity))


def test_infinity_comparisons_against_numbers():
    assert Infinity > 1e300
    assert Infinity >= 1e300
    assert not (Infinity < 1e300)
    assert not (Infinity <= 1e300)
    assert Infinity != 1e300
    assert not (Infinity == 1e300)


def test_negative_infinity_comparisons_against_numbers():
    assert NegativeInfinity < -1e300
    assert NegativeInfinity <= -1e300
    assert not (NegativeInfinity > -1e300)
    assert not (NegativeInfinity >= -1e300)
    assert NegativeInfinity != -1e300
    assert not (NegativeInfinity == -1e300)


def test_infinity_equality_with_own_type_only():
    other_infinity = InfinityType()
    assert Infinity == other_infinity
    assert not (Infinity != other_infinity)
    assert Infinity != NegativeInfinity
    assert not (Infinity == NegativeInfinity)


def test_negative_infinity_equality_with_own_type_only():
    other_neg = NegativeInfinityType()
    assert NegativeInfinity == other_neg
    assert not (NegativeInfinity != other_neg)


def test_infinity_gt_ge_and_negative_lt_le_always_true():
    assert Infinity > object()
    assert Infinity >= object()
    assert NegativeInfinity < object()
    assert NegativeInfinity <= object()


def test_infinity_negation():
    assert -Infinity is NegativeInfinity
    assert -NegativeInfinity is Infinity
    assert -(-Infinity) is Infinity


# ------------------------------ internal helpers ---------------------------


def test_parse_version_parts_pads_numeric_and_marks_alpha():
    parts = list(_parse_version_parts("1.0"))
    assert parts[0] == "00000001"
    assert parts[1] == "00000000"
    assert parts[-1] == "*final"


def test_parse_version_parts_replacement_map():
    # "pre"/"preview" -> "c", "rc" -> "c", "dev" -> "@", "-" -> "final-"
    parts = list(_parse_version_parts("1.0-dev"))
    assert "@" in [p.lstrip("*") for p in parts] or any(
        p == "*@" for p in parts
    )


def test_legacy_cmpkey_epoch_is_always_negative_one():
    epoch, parts = _legacy_cmpkey("1.0")
    assert epoch == -1
    assert isinstance(parts, tuple)


def test_legacy_cmpkey_strips_trailing_zero_parts():
    # trailing ".0" components are dropped from the numeric tail
    _, parts_with_zero = _legacy_cmpkey("1.0.0")
    _, parts_without = _legacy_cmpkey("1")
    assert parts_with_zero == parts_without


def test_legacy_cmpkey_removes_dash_before_prerelease_tag():
    # "1.0-alpha" and "1.0alpha" should normalize to the same key because the
    # trailing "final-" marker inserted for the literal "-" is popped again
    # once a prerelease tag ("*alpha" < "*final") is encountered.
    assert _legacy_cmpkey("1.0-alpha") == _legacy_cmpkey("1.0alpha")


@pytest.mark.parametrize(
    "letter, number, expected",
    [
        ("Alpha", "1", ("a", 1)),
        ("BETA", "2", ("b", 2)),
        ("c", None, ("rc", 0)),
        ("pre", "3", ("rc", 3)),
        ("preview", "4", ("rc", 4)),
        ("rev", "5", ("post", 5)),
        ("r", "6", ("post", 6)),
        ("rc", "7", ("rc", 7)),
        (None, "8", ("post", 8)),  # implicit post
        (None, None, None),
        ("", "9", ("post", 9)),  # empty letter is falsy -> implicit post path
        ("", None, None),  # empty letter, no number -> falls through to None
    ],
)
def test_parse_letter_version_direct(letter, number, expected):
    assert _parse_letter_version(letter, number) == expected


def test_parse_local_version_direct():
    assert _parse_local_version(None) is None
    assert _parse_local_version("abc.1.twelve") == ("abc", 1, "twelve")
    assert _parse_local_version("a_b-c.d") == ("a", "b", "c", "d")


@pytest.mark.parametrize(
    "epoch, release, pre, post, dev, local",
    [
        (0, (1, 0, 0), None, None, None, None),
        (0, (1, 0), ("a", 1), None, None, None),
        (0, (1, 0), None, ("post", 1), None, None),
        (0, (1, 0), None, None, ("dev", 1), None),  # dev-only -> pre=-Inf
        (0, (1, 0), ("a", 1), ("post", 1), ("dev", 1), None),
        (0, (1, 0), None, None, None, (1, "abc")),  # local mixed types
    ],
)
def test_cmpkey_direct_shapes(epoch, release, pre, post, dev, local):
    key = _cmpkey(epoch, release, pre, post, dev, local)
    assert key[0] == epoch
    # release with no trailing zeros stays the same length or shrinks
    assert isinstance(key[1], tuple)


def test_cmpkey_dev_only_sorts_before_any_prerelease():
    # pre=None, post=None, dev=set -> _pre becomes NegativeInfinity
    dev_only = _cmpkey(0, (1, 0), None, None, ("dev", 1), None)
    with_pre = _cmpkey(0, (1, 0), ("a", 0), None, None, None)
    assert dev_only < with_pre


def test_cmpkey_local_none_sorts_before_local_present():
    no_local = _cmpkey(0, (1, 0), None, None, None, None)
    with_local = _cmpkey(0, (1, 0), None, None, None, (1,))
    assert no_local < with_local


def test_cmpkey_local_alpha_sorts_before_numeric():
    alpha_local = _cmpkey(0, (1, 0), None, None, None, ("abc",))
    numeric_local = _cmpkey(0, (1, 0), None, None, None, (1,))
    assert alpha_local < numeric_local
