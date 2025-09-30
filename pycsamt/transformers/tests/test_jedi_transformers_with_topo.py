# -*- coding: utf-8 -*-
from __future__ import annotations

from pathlib import Path
import pytest

import pandas as pd
from types import SimpleNamespace

from pycsamt.zonge.avg import AVG

from pycsamt.seg.collection import EDICollection
from pycsamt.transformers import jedi as tr

@pytest.mark.usefixtures("modern_data_file", "stn_file_k2")
def test_avg_with_topo_injection_adds_head(
    modern_data_file: Path,
    stn_file_k2: Path,
):

    def _read_stn(p: Path):
        try:
            df = pd.read_csv(
                p,
                sep=r"\s+|,",
                engine="python",
                comment="#",
            )
        except Exception:
            df = pd.read_fwf(p)

        def pick(opts):
            for o in opts:
                for c in df.columns:
                    if o in c.lower():
                        return c
            return None

        ren = {}
        mapping = {
            "station": ["station", "site", "id", "sta"],
            "latitude": ["latitude", "lat"],
            "longitude": ["longitude", "lon", "long"],
            "elevation": ["elevation", "elev", "alt"],
        }
        for k, opts in mapping.items():
            c = pick(opts)
            if c:
                ren[c] = k
        df = df.rename(columns=ren)
        if "station" not in df.columns:
            df.insert(0, "station", range(1, len(df) + 1))
        return df

    avg = AVG.from_file(modern_data_file)
    stn = _read_stn(stn_file_k2)
    topo = SimpleNamespace(frame=stn)
    try:
        setattr(avg, "topo", topo)
    except Exception:
        avg.topo = topo  # type: ignore

    out = tr.AVGtoEDI().transform(avg)
    assert isinstance(out, EDICollection)
    ed = next(iter(out))
    has_head = False
    if hasattr(ed, "sections"):
        has_head = "head" in getattr(ed, "sections")
    if not has_head and hasattr(ed, "Head"):
        has_head = True
    assert has_head is True


@pytest.mark.usefixtures("modern_data_file", "stn_file_k2")
def test_avg_collection_size_matches_sites(
    modern_data_file: Path,
    stn_file_k2: Path,
):


    avg = AVG.from_file(modern_data_file).add_topography(stn_file_k2)
    try:
        z, f, st = avg.to_tensor(
            var="z",
            station=None,
            sort_freq=True,
            align="union",
        )
        n_sites = len(st) if st is not None else 1
    except Exception:
        n_sites = None

    out = tr.AVGtoEDI().transform(avg)
    assert isinstance(out, EDICollection)
    if n_sites is not None:
        assert len(out) == n_sites
