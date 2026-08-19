# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Regression tests for the EDI/EMTF-XML dual-backed Site/Sites design."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.emtf.converters.edi import DataLossWarning
from pycsamt.emtf.document import EMTF
from pycsamt.seg.edi import EDIFile
from pycsamt.site.base import Site, Sites, to_sites


def _load_edi(p: Path) -> EDIFile:
    return EDIFile(p)  # type: ignore


def _dup_edi(tmp_path: Path, src: Path, stem: str) -> Path:
    dst = tmp_path / f"{stem}.edi"
    dst.write_text(src.read_text(encoding="utf-8"), encoding="utf-8")
    return dst


def _mk_xml(tmp_path: Path, edi: EDIFile, stem: str) -> Path:
    """Produce a real XML fixture via pycsamt's own EDI->EMTF pipeline."""
    tf = EMTF.from_edi(edi)
    dst = tmp_path / f"{stem}.xml"
    tf.write_xml(dst)
    return dst


# ---------------------------------------------------------------------------
# Site: construction and lazy dual backend
# ---------------------------------------------------------------------------


def test_site_from_edi_backend_and_lazy_tf(simulated_edi: Path) -> None:
    s = Site.from_edi(simulated_edi)
    assert s.backend == "edi"
    tf = s.tf
    assert isinstance(tf, EMTF)
    assert s.tf is tf  # cached, same object on repeated access


def test_site_from_xml_backend_and_lazy_edi(
    tmp_path: Path, simulated_edi: Path
) -> None:
    edi = _load_edi(simulated_edi)
    xml_path = _mk_xml(tmp_path, edi, "sim_xml")

    s = Site.from_xml(xml_path)
    assert s.backend == "xml"
    assert isinstance(s.tf, EMTF)

    e = s.edi
    assert hasattr(e, "get_section")
    assert s.edi is e  # cached, same object on repeated access


def test_site_from_xml_numeric_accessors_match_source(
    tmp_path: Path, simulated_edi: Path
) -> None:
    # Derive the XML fixture from the already-`Site`-normalized EDI (its
    # dataid is stamped from the file stem at `Site()` construction) so
    # both sides share the same station identity for a fair comparison.
    p = _dup_edi(tmp_path, simulated_edi, "NUMC")
    s_edi = Site(_load_edi(p))
    xml_path = _mk_xml(tmp_path, s_edi.edi, "NUMC_xml")

    s_xml = Site.from_xml(xml_path)

    np.testing.assert_allclose(np.sort(s_edi.freq), np.sort(s_xml.freq))
    assert s_xml.name == s_edi.name
    assert s_xml.coords == pytest.approx(s_edi.coords, nan_ok=True)


def test_site_from_xml_name_and_coords_do_not_force_edi_materialization(
    tmp_path: Path, simulated_edi: Path
) -> None:
    edi = _load_edi(simulated_edi)
    xml_path = _mk_xml(tmp_path, edi, "sim_lazy")

    s = Site.from_xml(xml_path)
    assert s.name  # touch name
    assert s.coords  # touch coords
    assert s.meta  # touch meta
    # None of the above should have materialized the EDI view.
    assert s._edi_obj is None


def test_site_from_xml_rich_metadata_properties(
    tmp_path: Path, simulated_edi: Path
) -> None:
    edi = _load_edi(simulated_edi)
    xml_path = _mk_xml(tmp_path, edi, "sim_meta")

    s = Site.from_xml(xml_path)
    assert s.site_meta is s.tf.site
    assert s.site_layout is s.tf.site_layout
    assert s.provenance is s.tf.provenance
    assert s.processing is s.tf.processing
    assert s.copyright is s.tf.copyright
    assert s.quality_meta is s.tf.quality


def test_site_to_xml_roundtrip(tmp_path: Path, simulated_edi: Path) -> None:
    edi = _load_edi(simulated_edi)
    s = Site(edi)
    out = tmp_path / "roundtrip.xml"
    s.to_xml(out)
    assert out.exists()

    s2 = Site.from_xml(out)
    np.testing.assert_allclose(np.sort(s2.freq), np.sort(s.freq))


def test_site_to_xml_without_target_returns_string(
    simulated_edi: Path,
) -> None:
    s = Site.from_edi(simulated_edi)
    text = s.to_xml()
    assert isinstance(text, str)
    assert "<EM_TF>" in text


# ---------------------------------------------------------------------------
# Sites: mixed EDI/XML construction
# ---------------------------------------------------------------------------


def test_sites_accepts_mixed_edi_and_xml_objects(
    tmp_path: Path, simulated_edi: Path
) -> None:
    p1 = _dup_edi(tmp_path, simulated_edi, "M01")
    p2 = _dup_edi(tmp_path, simulated_edi, "M02")
    e1, e2 = _load_edi(p1), _load_edi(p2)
    xml2 = _mk_xml(tmp_path, e2, "M02_xml")
    tf2 = EMTF.from_xml(xml2)

    sites = Sites([e1, tf2])
    assert len(sites) == 2
    backends = sorted(s.backend for s in sites)
    assert backends == ["edi", "xml"]


def test_sites_to_emtf_list_and_write_xml(
    tmp_path: Path, simulated_edi: Path
) -> None:
    p1 = _dup_edi(tmp_path, simulated_edi, "W01")
    p2 = _dup_edi(tmp_path, simulated_edi, "W02")
    sites = Sites([_load_edi(p1), _load_edi(p2)])

    docs = sites.to_emtf_list()
    assert len(docs) == 2
    assert all(isinstance(d, EMTF) for d in docs)

    outdir = tmp_path / "xml_out"
    paths = sites.write_xml(outdir)
    assert len(paths) == 2
    assert all(p.exists() for p in paths)


# ---------------------------------------------------------------------------
# to_sites / ensure_sites coercion: XML recognition
# ---------------------------------------------------------------------------


def test_to_sites_single_xml_path(tmp_path: Path, simulated_edi: Path) -> None:
    edi = _load_edi(simulated_edi)
    xml_path = _mk_xml(tmp_path, edi, "coerce_single")

    out = to_sites(xml_path)
    assert isinstance(out, Sites)
    assert len(out) == 1
    assert out[0].backend == "xml"


def test_to_sites_single_emtf_object(
    tmp_path: Path, simulated_edi: Path
) -> None:
    edi = _load_edi(simulated_edi)
    tf = EMTF.from_edi(edi)

    out = to_sites(tf)
    assert isinstance(out, Sites)
    assert len(out) == 1
    assert out[0].backend == "xml"


def test_to_sites_list_of_xml_paths(
    tmp_path: Path, simulated_edi: Path
) -> None:
    edi = _load_edi(simulated_edi)
    p1 = _mk_xml(tmp_path, edi, "L01")
    p2 = _mk_xml(tmp_path, edi, "L02")

    out = to_sites([p1, p2])
    assert isinstance(out, Sites)
    assert len(out) == 2
    assert all(s.backend == "xml" for s in out)


def test_to_sites_mixed_edi_and_xml_directory(
    tmp_path: Path, simulated_edi: Path
) -> None:
    survey = tmp_path / "survey"
    survey.mkdir()
    edi_path = survey / "D01.edi"
    edi_path.write_text(
        simulated_edi.read_text(encoding="utf-8"), encoding="utf-8"
    )
    edi = _load_edi(edi_path)
    tf = EMTF.from_edi(edi)
    tf.write_xml(survey / "D02.xml")

    out = to_sites(survey)
    assert isinstance(out, Sites)
    names = sorted(s.name for s in out)
    assert names == ["D01", "SIM01"]  # D02.xml keeps its own station id
    backends = {s.name: s.backend for s in out}
    assert backends["D01"] == "edi"


def test_to_sites_xml_backed_site_in_a_list(
    tmp_path: Path, simulated_edi: Path
) -> None:
    edi = _load_edi(simulated_edi)
    xml_path = _mk_xml(tmp_path, edi, "site_in_list")
    xml_site = Site.from_xml(xml_path)
    edi_site = Site(_load_edi(_dup_edi(tmp_path, simulated_edi, "EDIONE")))

    out = to_sites([edi_site, xml_site])
    assert isinstance(out, Sites)
    assert len(out) == 2
    assert {s.backend for s in out} == {"edi", "xml"}


# ---------------------------------------------------------------------------
# Regression: list-of-Site coercion (previously silently dropped)
# ---------------------------------------------------------------------------


def test_to_sites_list_of_edi_backed_sites_not_dropped(
    tmp_path: Path, simulated_edi: Path
) -> None:
    p1 = _dup_edi(tmp_path, simulated_edi, "K01")
    p2 = _dup_edi(tmp_path, simulated_edi, "K02")
    s1, s2 = Site(_load_edi(p1)), Site(_load_edi(p2))

    out = to_sites([s1, s2])
    assert isinstance(out, Sites)
    # Previously silently returned an empty Sites: a bare Site does
    # not look EDI-shaped itself (no get_section/Z), so it was
    # dropped by iter_edifiles's duck-type check.
    assert len(out) == 2
    assert sorted(s.name for s in out) == ["K01", "K02"]


def test_ensure_sites_list_of_sites_not_dropped(
    tmp_path: Path, simulated_edi: Path
) -> None:
    from pycsamt.emtools._core import ensure_sites

    p1 = _dup_edi(tmp_path, simulated_edi, "EN1")
    p2 = _dup_edi(tmp_path, simulated_edi, "EN2")
    s1, s2 = Site(_load_edi(p1)), Site(_load_edi(p2))

    out = ensure_sites([s1, s2])
    assert len(out) == 2


def test_ensure_sites_single_xml_path(
    tmp_path: Path, simulated_edi: Path
) -> None:
    from pycsamt.emtools._core import ensure_sites

    edi = _load_edi(simulated_edi)
    xml_path = _mk_xml(tmp_path, edi, "ensure_single")

    out = ensure_sites(xml_path)
    assert len(out) == 1
    assert out[0].backend == "xml"


def test_ensure_sites_list_of_xml_paths(
    tmp_path: Path, simulated_edi: Path
) -> None:
    from pycsamt.emtools._core import ensure_sites

    edi = _load_edi(simulated_edi)
    p1 = _mk_xml(tmp_path, edi, "EX1")
    p2 = _mk_xml(tmp_path, edi, "EX2")

    out = ensure_sites([p1, p2])
    assert len(out) == 2


def test_ensure_sites_sites_passthrough(
    tmp_path: Path, simulated_edi: Path
) -> None:
    from pycsamt.emtools._core import ensure_sites

    p1 = _dup_edi(tmp_path, simulated_edi, "PT1")
    sites = Sites([_load_edi(p1)])

    out = ensure_sites(sites)
    # `ensure_sites` always applies `.ordered(order_by)`, which returns a
    # fresh `Sites` container (unless mode is exactly "input") -- but it
    # must reuse the same underlying `Site` objects, not silently
    # re-derive them (that's the regression this test guards against).
    assert list(out) == list(sites)


# ---------------------------------------------------------------------------
# Data-loss policy on materialization
# ---------------------------------------------------------------------------


def test_edi_materialization_honors_on_loss_raise(
    tmp_path: Path, simulated_edi: Path
) -> None:
    from pycsamt.emtf.converters.edi import EMTFEDIConversionError

    edi = _load_edi(simulated_edi)
    tf = EMTF.from_edi(edi)
    # Attach metadata EDI cannot represent so materialization has
    # something to warn/raise about.
    from pycsamt.metadata import Person, ProvenanceMeta

    tf.provenance = ProvenanceMeta(
        creator=Person(name="Ada Lovelace", email="ada@example.org"),
        submitter=Person(name="Ada Lovelace", email="ada@example.org"),
    )

    s_ignore = Site(tf, on_loss="ignore")
    assert s_ignore.edi is not None  # no raise

    s_raise = Site(EMTF.from_edi(edi), on_loss="raise")
    s_raise.tf.provenance = ProvenanceMeta(
        creator=Person(name="Ada Lovelace"),
        submitter=Person(name="Ada Lovelace"),
    )
    with pytest.raises(EMTFEDIConversionError):
        _ = s_raise.edi


def test_edi_materialization_default_warns_not_silent(
    tmp_path: Path, simulated_edi: Path
) -> None:
    from pycsamt.metadata import Person, ProvenanceMeta

    edi = _load_edi(simulated_edi)
    tf = EMTF.from_edi(edi)
    tf.provenance = ProvenanceMeta(
        creator=Person(name="Ada Lovelace"),
        submitter=Person(name="Ada Lovelace"),
    )
    s = Site(tf)  # default on_loss="warn"
    with pytest.warns(DataLossWarning):
        _ = s.edi
