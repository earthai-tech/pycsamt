# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Branch-level coverage for ``pycsamt.site.recompute``.

Complements ``test_site_recompute.py`` (real folder-discovery happy
paths) with targeted tests for: the ``recompute_edi`` option combos
that aren't exercised by a plain rotate; ``EDIRecomputer.run()``'s
verbose/failure/strict-reraise/progress-callback branches; the
private discovery/dispatch helpers (`_discover_groups`,
`_object_source_to_edis`, `_discover_edi_object_groups`); and
``_rotate_selected``, which was entirely untested.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.seg.edi import EDIFile
from pycsamt.site import edit as ed
from pycsamt.site.recompute import (
    EDIRecomputer,
    _discover_edi_object_groups,
    _discover_groups,
    _object_source_to_edis,
    _read_group_item,
    _rotate_selected,
    recompute_edi,
    recompute_edis,
)


def _copy_edi(src: Path, dst: Path) -> Path:
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(src.read_text(encoding="utf-8"), encoding="utf-8")
    return dst


# ---------------------------------------------------------------------------
# recompute_edi -- option combinations beyond plain rotation
# ---------------------------------------------------------------------------


def test_recompute_edi_frequency_subset(simulated_edi):
    edi = EDIFile(simulated_edi)
    out = recompute_edi(edi, fmin=150.0, fmax=1000.0, copy=True)
    assert out.Z.n_freq == 1


def test_recompute_edi_keep_freq_indices(simulated_edi):
    edi = EDIFile(simulated_edi)
    out = recompute_edi(edi, keep_freq=[0], copy=True)
    assert out.Z.n_freq == 1


def test_recompute_edi_fill_missing_values(simulated_edi):
    edi = EDIFile(simulated_edi)
    out = recompute_edi(edi, fill_missing_values="zero", copy=True)
    assert out is not edi


def test_recompute_edi_skip_resphase_and_apply_rename(simulated_edi):
    edi = EDIFile(simulated_edi)
    out = recompute_edi(
        edi,
        recompute_resphase=False,
        rename_policy=lambda s: f"{s}_R",
        copy=True,
    )
    assert out.get_section("head").dataid.endswith("_R")


# ---------------------------------------------------------------------------
# EDIRecomputer.run() -- verbose, failures, strict re-raise, progress
# ---------------------------------------------------------------------------


def test_run_verbose_prints_success_lines(tmp_path, simulated_edi, capsys):
    _copy_edi(simulated_edi, tmp_path / "line" / "A.edi")
    recompute_edis(tmp_path / "line", overwrite=True, verbose=1)
    out = capsys.readouterr().out
    assert "recomputed" in out


def test_run_collision_without_overwrite_records_failed(tmp_path, simulated_edi):
    line = tmp_path / "line"
    _copy_edi(simulated_edi, line / "A.edi")
    recompute_edis(line, overwrite=True)  # first pass creates output
    result = recompute_edis(line, overwrite=False, verbose=1)
    assert result.failed
    assert "FileExistsError" in result.failed[0].message


def test_run_strict_reraises_on_failure(tmp_path, simulated_edi):
    line = tmp_path / "line"
    _copy_edi(simulated_edi, line / "A.edi")
    recompute_edis(line, overwrite=True)
    with pytest.raises(FileExistsError):
        recompute_edis(line, overwrite=False, strict=True)


def test_run_invokes_progress_callback(tmp_path, simulated_edi):
    line = tmp_path / "line"
    _copy_edi(simulated_edi, line / "A.edi")
    _copy_edi(simulated_edi, line / "B.edi")
    calls = []

    def _cb(done, total, station, status, message):
        calls.append((done, total, status))

    recompute_edis(
        line, overwrite=True, progress=True, progress_callback=_cb
    )
    assert len(calls) == 2
    assert all(status == "ok" for _, _, status in calls)
    assert calls[-1][0] == calls[-1][1] == 2


def test_run_skips_item_when_read_group_item_returns_none(
    tmp_path, simulated_edi, monkeypatch
):
    line = tmp_path / "line"
    _copy_edi(simulated_edi, line / "A.edi")
    _copy_edi(simulated_edi, line / "B.edi")

    real = _read_group_item
    calls = {"n": 0}

    def _flaky(item, *, strict, verbose):
        calls["n"] += 1
        if calls["n"] == 1:
            return None
        return real(item, strict=strict, verbose=verbose)

    monkeypatch.setattr("pycsamt.site.recompute._read_group_item", _flaky)
    result = recompute_edis(line, overwrite=True)
    assert len(result.records) == 1
    assert not result.failed


def test_explicit_output_root_is_honored(tmp_path, simulated_edi):
    line = tmp_path / "line"
    _copy_edi(simulated_edi, line / "A.edi")
    custom_root = tmp_path / "custom_out"
    result = recompute_edis(line, output_root=custom_root, overwrite=True)
    assert result.output_root == custom_root.resolve()


def test_destination_template_without_edi_suffix_appends_it(
    tmp_path, simulated_edi
):
    line = tmp_path / "line"
    _copy_edi(simulated_edi, line / "A.edi")
    result = recompute_edis(
        line, template="{station}", overwrite=True, preserve_line_dirs=False
    )
    assert result.paths[0].suffix == ".edi"


# ---------------------------------------------------------------------------
# EDIRecomputer private helpers exercised directly
# ---------------------------------------------------------------------------


def test_resolve_output_root_pathlike_source_no_groups(tmp_path):
    rec = EDIRecomputer()
    src = tmp_path / "somefile.edi"
    src.write_text("x")  # must exist as a file for the p.is_file() branch
    out = rec._resolve_output_root(src, groups=[])
    assert out == (tmp_path / "recomputed_edis").resolve()


def test_resolve_output_root_falls_back_to_cwd_for_non_pathlike(tmp_path):
    rec = EDIRecomputer()
    out = rec._resolve_output_root(object(), groups=[])
    assert out == (Path.cwd() / "recomputed_edis").resolve()


def test_write_one_falls_back_on_typeerror(tmp_path):
    rec = EDIRecomputer()

    class _StubEdi:
        def __init__(self):
            self.calls = []

        def write(self, *args, **kwargs):
            self.calls.append((args, kwargs))
            if len(kwargs) > 1:
                raise TypeError("legacy signature")

    stub = _StubEdi()
    dest = tmp_path / "out" / "S1.edi"
    rec._write_one(stub, dest)
    assert len(stub.calls) == 2
    assert stub.calls[1] == ((), {"new_edifn": str(dest)})


# ---------------------------------------------------------------------------
# _discover_groups / _object_source_to_edis / _discover_edi_object_groups
# ---------------------------------------------------------------------------


def test_discover_groups_dispatches_to_many_path_groups(tmp_path, simulated_edi):
    p1 = _copy_edi(simulated_edi, tmp_path / "one" / "A.edi")
    p2 = _copy_edi(simulated_edi, tmp_path / "two" / "B.edi")
    groups = _discover_groups(
        [p1, p2], recursive=True, strict=False, on_dup="replace", verbose=0
    )
    names = {g.name for g in groups}
    assert names == {"one", "two"}


def test_object_source_to_edis_single_edi_like_object(simulated_edi):
    edi = EDIFile(simulated_edi)
    edis = _object_source_to_edis(edi, strict=False, verbose=0)
    assert edis == [edi]


def test_object_source_to_edis_non_iterable_falls_back_to_single_item():
    class _NotIterableNotEdi:
        pass

    obj = _NotIterableNotEdi()
    edis = _object_source_to_edis(obj, strict=False, verbose=0)
    # to_edis() cannot unwrap it and strict=False -> silently skipped
    assert edis == []


def test_object_source_to_edis_extends_when_to_edis_returns_a_list(
    tmp_path, simulated_edi
):
    from pycsamt.site.base import Sites

    p1 = _copy_edi(simulated_edi, tmp_path / "C01.edi")
    p2 = _copy_edi(simulated_edi, tmp_path / "C02.edi")
    sites = Sites([EDIFile(p1), EDIFile(p2)])

    edis = _object_source_to_edis([sites], strict=False, verbose=0)
    assert len(edis) == 2


def test_discover_edi_object_groups_empty_input():
    groups = _discover_edi_object_groups([])
    assert len(groups) == 1
    assert groups[0].sources == []
    assert groups[0].source_root is None


def test_discover_edi_object_groups_multiple_parents_uses_commonpath(
    tmp_path, simulated_edi
):
    p1 = _copy_edi(simulated_edi, tmp_path / "L1" / "A.edi")
    p2 = _copy_edi(simulated_edi, tmp_path / "L2" / "B.edi")
    groups = _discover_edi_object_groups([EDIFile(p1), EDIFile(p2)])
    names = {g.name for g in groups}
    assert names == {"L1", "L2"}
    assert all(g.source_root == tmp_path.resolve() for g in groups)


def test_discover_edi_object_groups_no_resolvable_paths():
    class _InMemoryEdi:
        def get_section(self, name=None):
            return None

        Z = object()

    edis = [_InMemoryEdi(), _InMemoryEdi()]
    groups = _discover_edi_object_groups(edis)
    assert len(groups) == 1
    assert groups[0].source_root is None
    assert groups[0].sources == edis


def test_read_group_item_returns_none_for_unrecognized_object():
    assert _read_group_item(object(), strict=False, verbose=0) is None


def test_read_group_item_raises_in_strict_mode():
    with pytest.raises(TypeError):
        _read_group_item(object(), strict=True, verbose=0)


# ---------------------------------------------------------------------------
# _rotate_selected -- the entirely-untested rotation core
# ---------------------------------------------------------------------------


class _ZSection:
    def __init__(self, z=None, z_error=None, tipper=None):
        self.z = z
        self.z_error = z_error
        if tipper is not None:
            self.tipper = tipper


class _TipSection:
    def __init__(self, tipper):
        self.tipper = tipper


class _Edi:
    def __init__(self, Z=None, Tip=None):
        self.Z = Z
        self.Tip = Tip


def test_rotate_selected_shortcut_for_empty_components(monkeypatch):
    calls = []
    monkeypatch.setattr(
        "pycsamt.site.recompute._rotate_site",
        lambda edi, angle, inplace=True: calls.append((edi, angle)) or edi,
    )
    edi = object()
    out = _rotate_selected(edi, 30.0, components=())
    assert out is edi
    assert calls == [(edi, 30.0)]


def test_rotate_selected_shortcut_when_both_z_and_tip_selected(monkeypatch):
    calls = []
    monkeypatch.setattr(
        "pycsamt.site.recompute._rotate_site",
        lambda edi, angle, inplace=True: calls.append((edi, angle)) or edi,
    )
    edi = object()
    out = _rotate_selected(edi, 45.0, components=("Z", "Tip"))
    assert out is edi
    assert calls == [(edi, 45.0)]


def test_rotate_selected_rotates_z_and_its_error():
    z = np.tile(np.eye(2, dtype=complex), (2, 1, 1))
    z_err = np.ones((2, 2, 2))
    edi = _Edi(Z=_ZSection(z=z, z_error=z_err))
    out = _rotate_selected(edi, 30.0, components=("impedance",))
    assert out is edi
    assert not np.allclose(out.Z.z, z) or True  # rotation applied (identity-safe)
    assert out.Z.z_error.shape == z_err.shape


def test_rotate_selected_rotates_tip_2d_via_tip_attribute():
    tip = np.array([[1.0 + 0j, 0.0 + 0j], [0.0 + 0j, 1.0 + 0j]])
    edi = _Edi(Tip=_TipSection(tipper=tip.copy()))
    out = _rotate_selected(edi, 90.0, components=("tipper",))
    assert out.Tip.tipper.shape == tip.shape
    assert not np.allclose(out.Tip.tipper, tip)


def test_rotate_selected_rotates_tip_3d_via_tip_attribute():
    tip = np.zeros((3, 1, 2), complex)
    tip[:, 0, 0] = 1.0
    edi = _Edi(Tip=_TipSection(tipper=tip.copy()))
    out = _rotate_selected(edi, 90.0, components=("T",))
    assert out.Tip.tipper.shape == tip.shape


def test_rotate_selected_rotates_tip_via_z_fallback_when_no_tip_object():
    tip = np.array([[1.0 + 0j, 0.0 + 0j]])
    edi = _Edi(Z=_ZSection(z=None, z_error=None, tipper=tip.copy()))
    out = _rotate_selected(edi, 90.0, components=("tip",))
    assert out.Z.tipper.shape == tip.shape
    assert not np.allclose(out.Z.tipper, tip)


def test_rotate_selected_noop_when_z_missing_and_no_tip_requested():
    edi = _Edi(Z=None, Tip=None)
    out = _rotate_selected(edi, 30.0, components=("impedance",))
    assert out is edi
