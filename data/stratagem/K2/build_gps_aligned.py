"""Build k2-gps-aligned.csv from k2-gps-ref.csv + k2-field-report.docx.

Reconciles two real-world data-quality problems in the raw K2 delivery so
that pycsamt.stratagem.CoordinateInjector (which requires one coordinate
row per EDI/raw station) can be run against this dataset:

1. Station-count mismatch (87 EDI/raw files vs. 80 GPS rows). Cross-
   referencing the field report explains the gap: file 2HX.01 is a
   calibration/"parallel test" shot (chainage 1 m, off the survey line),
   and files 2HX.043/069/075/087 are "check points" — QC repeats shot at
   the same physical location as the immediately preceding station. The
   field-report chainage ("里程") also reveals a physical point (560 m)
   that was surveyed by GPS but never successfully recorded by the CSAMT
   instrument ("公路中间车很多测不了" — too much road traffic to measure);
   that orphan GPS row has no matching file and is simply unused here.

2. The `station` column in k2-gps-ref.csv is shifted one row relative to
   the `step` column (verify: station[i]'s numeric suffix == step[i+1]).
   `step` is treated as the trustworthy per-row chainage.

3. The last GPS row (Site=K79, step=1640.4 m) is a bad shot: its
   consecutive-point distance from the second-to-last row is ~1640 m
   instead of the ~20 m design spacing, and its coordinates sit right on
   top of the profile's *start* point rather than continuing the eastward/
   northward trend of the last ~10 stations. It is dropped, and the true
   endpoint (chainage 1640 m, stations 86 and 87) is instead extrapolated
   from a linear fit through the last 5 trusted points.

3 further chainage gaps (~1080, ~1260, ~1300 m) were never GPS-shot at all
(no ~1640 m outlier involved) and are filled by linear interpolation
between their chainage-adjacent neighbours.

Output columns (one row per raw/EDI station number, 1..87):
    station_num, edi_file, raw_x_file, chainage_design_m,
    northing, easting, elev, source, use_for_survey

`source` documents how each row's coordinate was obtained:
    gps_exact | interpolated | extrapolated | checkpoint_dup | test_shot_approx

`use_for_survey` is False only for station 1 (the calibration/test shot,
not a real profile position) — drop it before spatial pipeline steps such
as StaticShiftCorrector's AMA neighbour averaging.

Run from the repo root:
    python data/stratagem/K2/build_gps_aligned.py
"""

from __future__ import annotations

import re
from pathlib import Path

import docx
import numpy as np
import pandas as pd

HERE = Path(__file__).parent


def _load_field_report(path: Path) -> pd.DataFrame:
    d = docx.Document(str(path))
    rows = [[c.text.strip() for c in r.cells] for r in d.tables[0].rows[1:]]
    fr = pd.DataFrame(
        rows,
        columns=[
            "date",
            "weather",
            "fileno",
            "line",
            "point",
            "chainage",
            "terrain",
            "notes",
        ],
    )
    fr["chainage"] = fr["chainage"].astype(float)
    fr["is_checkpoint"] = fr["notes"].str.contains("检查点", na=False) | fr[
        "terrain"
    ].str.contains("检查点", na=False)

    def file_num(s):
        if pd.isna(s) or not str(s).strip():
            return None
        m = re.search(r"\.(\d+)$", str(s))
        return int(m.group(1)) if m else None

    fr["station_num"] = fr["fileno"].apply(file_num)
    filed = fr[fr["station_num"].notna()].copy()
    filed["station_num"] = filed["station_num"].astype(int)
    filed = filed.sort_values("station_num").reset_index(drop=True)
    assert filed["station_num"].tolist() == list(range(1, 88)), (
        "field report station numbering has changed — reconciliation "
        "logic needs re-checking"
    )
    return filed


def _load_gps(path: Path) -> pd.DataFrame:
    gps = pd.read_csv(path, encoding="utf-8-sig")
    gps.columns = [c.strip() for c in gps.columns]
    gps["step"] = gps["step"].astype(float)
    gps = gps.sort_values("step").reset_index(drop=True)

    seg_dist = np.sqrt(gps["lon"].diff() ** 2 + gps["lat"].diff() ** 2)
    bad = (seg_dist > 100).fillna(False)  # design spacing is ~20 m
    return gps[~bad].reset_index(drop=True)


def _lookup(
    gps_trusted: pd.DataFrame, chain: float
) -> tuple[float, float, float, str]:
    within = gps_trusted[(gps_trusted["step"] - chain).abs() <= 5.0]
    if not within.empty:
        r = within.iloc[(within["step"] - chain).abs().argmin()]
        return r["lon"], r["lat"], r["elev"], "gps_exact"

    below = gps_trusted[gps_trusted["step"] <= chain]
    above = gps_trusted[gps_trusted["step"] >= chain]
    if not below.empty and not above.empty:
        b, a = below.iloc[-1], above.iloc[0]
        f = (chain - b["step"]) / (a["step"] - b["step"])
        return (
            b["lon"] + f * (a["lon"] - b["lon"]),
            b["lat"] + f * (a["lat"] - b["lat"]),
            b["elev"] + f * (a["elev"] - b["elev"]),
            "interpolated",
        )

    tail = gps_trusted.tail(5)
    A = np.vstack([tail["step"], np.ones(len(tail))]).T
    lon_m, lon_b = np.linalg.lstsq(A, tail["lon"], rcond=None)[0]
    lat_m, lat_b = np.linalg.lstsq(A, tail["lat"], rcond=None)[0]
    elev_m, elev_b = np.linalg.lstsq(A, tail["elev"], rcond=None)[0]
    return (
        lon_m * chain + lon_b,
        lat_m * chain + lat_b,
        elev_m * chain + elev_b,
        "extrapolated",
    )


def build(out_path: Path | None = None) -> pd.DataFrame:
    filed = _load_field_report(HERE / "k2-field-report.docx")
    gps_trusted = _load_gps(HERE / "k2-gps-ref.csv")

    out_rows = []
    prev_by_chain: dict[float, tuple[float, float, float]] = {}
    for _, row in filed.iterrows():
        sn = int(row["station_num"])
        chain = row["chainage"]
        if sn == 1:
            northing, easting, elev, _ = _lookup(gps_trusted, 0.0)
            tag, use = "test_shot_approx", False
        elif row["is_checkpoint"] and chain in prev_by_chain:
            northing, easting, elev = prev_by_chain[chain]
            tag, use = "checkpoint_dup", True
        else:
            northing, easting, elev, tag = _lookup(gps_trusted, chain)
            use = True
            prev_by_chain[chain] = (northing, easting, elev)

        out_rows.append(
            {
                "station_num": sn,
                "edi_file": f"Z2HX{sn:03d}.edi",
                "raw_x_file": f"X2HX.{sn:03d}",
                "chainage_design_m": chain,
                "northing": northing,
                "easting": easting,
                "elev": elev,
                "source": tag,
                "use_for_survey": use,
            }
        )

    aligned = pd.DataFrame(out_rows)
    if out_path is not None:
        aligned.to_csv(out_path, index=False)
    return aligned


if __name__ == "__main__":
    df = build(HERE / "k2-gps-aligned.csv")
    print(df["source"].value_counts())
    print(f"wrote {HERE / 'k2-gps-aligned.csv'} ({len(df)} rows)")
