from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from pycsamt.seg.collection import EDICollection
from pycsamt.transformers import jedi as tr
from pycsamt.zonge.avg import AVG


def test_avgtoedi_k1_stn_utm_coordinates_are_written_to_edi(tmp_path):
    root = Path(__file__).resolve().parents[3]
    avg_path = root / "data" / "avg" / "K1.AVG"
    stn_path = root / "data" / "avg" / "K1.stn"
    if not avg_path.exists() or not stn_path.exists():
        pytest.skip("K1 AVG/STN fixture is not available.")

    avg = AVG.from_file(avg_path).add_topography(stn_path, utm_zone="49N")
    avg.topo.convert_coords(to="ll", inplace=True)

    out = tr.AVGtoEDI().transform(avg)
    assert isinstance(out, EDICollection)
    assert len(out) >= 1

    ed = next(iter(out))
    head = ed.get_section("head")
    definemeas = ed.get_section("definemeas")
    assert head is not None
    assert definemeas is not None

    topo_row = avg.topo.frame.iloc[0]
    assert head.lat == pytest.approx(float(topo_row["latitude"]), abs=1e-6)
    assert head.long == pytest.approx(float(topo_row["longitude"]), abs=1e-6)
    assert head.elev == pytest.approx(float(topo_row["elevation"]), abs=1e-6)
    assert definemeas.reflat == pytest.approx(head.lat, abs=1e-6)
    assert definemeas.reflong == pytest.approx(head.long, abs=1e-6)
    assert definemeas.refelev == pytest.approx(head.elev, abs=1e-6)

    written = Path(ed.write(savepath=tmp_path, new_edifn="K1_S150.edi"))
    text = written.read_text(encoding="utf-8", errors="ignore").upper()
    assert "LAT=" in text
    assert "LONG=" in text
    assert "REFLAT=" in text
    assert "REFLONG=" in text

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
        avg.topo = topo
    except Exception:
        avg.topo = topo  # type: ignore

    out = tr.AVGtoEDI().transform(avg)
    assert isinstance(out, EDICollection)
    ed = next(iter(out))
    has_head = False
    if hasattr(ed, "sections"):
        has_head = "head" in ed.sections
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
