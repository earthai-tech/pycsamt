# L14 — CSAMT example survey (Hunan, China)

A real controlled-source audio-frequency magnetotelluric (CSAMT) field
line, used as the worked example for pyCSAMT's legacy Zonge AVG reader,
`.stn` topography reader, and AVG-to-EDI conversion pipeline.

## What's here

Only a minimal curated subset is tracked in git; everything else is
regenerable from it via the two scripts below and stays git-ignored to
keep the repo small.

| File                                       | Tracked          | Notes                                    |
|---------------------------------------------|:-----------------:|-------------------------------------------|
| `L14.avg`                                  | yes               | Raw legacy (kind-1) AMTAVG 7.40 file      |
| `L14.stn`                                  | yes               | Station topography, UTM zone 49N          |
| `avg_to_edi.py`                            | yes               | AVG + STN -> EDI, writes `edi/`           |
| `run_pipeline.py`                          | yes               | Runs the pipeline on `edi/`, writes `output/` |
| `edi/S1000.edi`, `S2120.edi`, `S3280.edi`  | yes               | 3 of 58 stations (first/middle/last)      |
| `edi/*.edi` (remaining 55)                 | no (local only)   | regenerate with `avg_to_edi.py`           |
| `output/pipeline.yaml`, `report.html`, `summary.txt` | yes     | small run artefacts                       |
| `output/processed/*.edi`, `output/plots/*.png` | no (local only) | regenerate with `run_pipeline.py`         |

## Survey characteristics

58 stations spaced 40 m apart (1000-3280 m along profile), single ExHy
dipole-dipole component (no cross-tensor or tipper channels), 40
logarithmically-spaced frequencies from 1.33 Hz to 9600 Hz. Station
coordinates are given in UTM zone 49N, placing the line in Hunan
Province, China.

## Usage

```console
# 1. Read the raw AVG + topography, convert to EDI
python data/avg/L14/avg_to_edi.py

# 2. Run the full_processing pipeline preset on the resulting EDIs
python data/avg/L14/run_pipeline.py
```

Or programmatically:

```python
from pycsamt.zonge import AVG

avg = AVG.from_file("data/avg/L14/L14.avg")
avg.add_topography("data/avg/L14/L14.stn", utm_zone="49N")
avg.topo.convert_coords(to="ll", inplace=True)
```
