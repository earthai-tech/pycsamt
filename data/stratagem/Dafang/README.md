# Dafang line 5 — Stratagem (Geometrics/EMI) CSAMT field line

A real controlled-source/audio-magnetotelluric field survey collected with
Stratagem hardware, line "DF5X" (site name **Dafang**, line 5), used as the
bundled example dataset for the field-to-inversion tutorial
(:doc:`../../../docs/source/tutorials/process_stratagem_dafang_to_inversion`).
22 stations at a nominal 20 m point spacing along a single ~400 m
centerline, near 25.76°N, 109.63°E (Guangxi).

## What's here

| Path | Tracked | Contents |
|------|:-------:|----------|
| `df5-edi/` | yes, 21 of 22 | WinGLink-exported EDI files (`ZDF5X001.edi` … `ZDF5X022.edi`, station 019 missing — see below). Coordinates are **not yet injected** — `LAT`/`LONG` are WinGLink placeholders (`0:00:00.00`). |
| `df5-hx/` | yes, all 3 components | Raw Stratagem hardware files: `XDF5X.001`–`.022`, `YDF5X.001`–`.022`, `ZDF5X.001`–`.022` (19-column ASCII, ~36 MB total — small enough to bundle every component, unlike the K2 dataset), plus calibration/sensor files (`SENSORS.TBL`, `@`, `AFEV3.*`, `BE26V2.*`, `S746_HX.*`, `S747_HY.*`). |
| `df5-coords.csv` | yes | The Dafang line-5 rows (`D5-00` … `D5-400`, 21 rows) extracted from a multi-line master coordinate table delivered separately from the raw hardware data. Columns: `stations, easting, northing, elev, sigma_step, step_cal, step, den, err, longitude, latitude` — already reconciled, one row per surviving EDI station, in acquisition order. |
| `df5-watex-edi/` | yes, 21 files | An independent prior correction of this same line by **WATEX** (an earlier package by the same author): coordinates injected, all 39 raw frequencies retained. Station names `S00`…`S20`. Used only as a cross-validation reference — not required for the pyCSAMT pipeline itself. |
| `df5-watex-edi-processed/` | yes, 21 files | The WATEX-frequency-filtered version of the same 21 stations (21 of the original 39 frequencies kept). Used as an independent "second opinion" to sanity-check pyCSAMT's own static-shift correction. |

## What's excluded

Field photographs, a `.docx` daily field-work report, and a duplicate
GBK-encoded raw survey-measurement CSV (superseded by `df5-coords.csv`,
which carries the same numbers already reconciled and UTF-8 encoded) are
not tracked here.

## A real gap in the WinGLink export

The raw hardware delivery (`df5-hx/`) has all 22 stations. The WinGLink
export (`df5-edi/`) is missing station **019** — `ZDF5X019.edi` simply does
not exist, even though `XDF5X.019`/`YDF5X.019`/`ZDF5X.019` are present in
the raw delivery. `df5-coords.csv` already reflects this: its 21 rows
(`D5-00` … `D5-400`) line up 1:1 with the 21 *surviving* EDI files in
acquisition order, not with all 22 raw stations. `EDIBatch("df5-edi").fit()`
therefore returns exactly 21 objects, and `df5-coords.csv` needs no
row-dropping before it is handed to `CoordinateInjector` — unlike the K2
dataset, where a calibration shot and four checkpoint repeats had to be
identified and excluded explicitly.

## Coordinate system and a column-naming quirk

The `easting`/`northing` column *values* in `df5-coords.csv` are swapped
relative to their column *names* — the same quirk documented for the K2
dataset. The column named `easting` (~2,850,000) is large enough to be a
Gauss-Kruger-style northing (it matches degrees latitude × 111 km/degree
almost exactly), and the column named `northing` (~362,000) is the real
easting. Feed them to `CoordinateInjector` with the names swapped:

```python
CoordinateInjector(epsg=32649, order="forward").fit(
    edi_objects, "df5-coords.csv",
    easting_col="northing", northing_col="easting",
    elev_col="elev", station_col="stations",
)
```

`epsg=32649` (UTM Zone 49N WGS84) was confirmed, not assumed: the raw
`(easting, northing)` pair for `D5-00` reprojects under EPSG 32649 to
`(25.761665, 109.630124)`, matching `df5-coords.csv`'s own precomputed
`longitude`/`latitude` columns for that row to six decimal places, and
matching the independently-injected coordinates already present in
`df5-watex-edi/ZDF5X001.edi` (`LAT=25:45:41.99`, `LONG=109:37:48.45`) once
converted from sexagesimal. `CoordinateInjector`'s own default,
`epsg=15921`, is not resolvable in every PROJ installation (it raised
`crs not found` in the environment this dataset was verified in) — use the
explicit `epsg=32649` shown above rather than relying on the default.

## Cross-validating against WATEX

`df5-watex-edi-processed/` is real, independently-processed output for
the same physical stations, from a different codebase's static-shift and
frequency-selection logic. It is not a ground truth (WATEX's own
correction has no more claim to "correct" than pyCSAMT's), but the two
should broadly agree in order of magnitude and period-dependence if
pyCSAMT's own correction is doing something reasonable. They will not
match exactly point-for-point: the two tools select different frequency
bands (WATEX keeps 21 of 39; pyCSAMT's `FrequencyFilter(fmin=10, fmax=10000)`
keeps a different count per station) and apply different static-shift and
noise-removal algorithms. See the tutorial's "Cross-Validating Against An
Independent Tool" section for a worked comparison.

## Real pyCSAMT bugs found while building the tutorial from this line

Building an actual, real, compiled Occam2D run from this line's corrected
output surfaced two confirmed bugs in `pycsamt.models.occam2d`, both now
fixed:

1. `OccamRunner.compile()` looked for a binary literally named `Occam2D`
   after a successful `make`, which never matches the `Occam2D.exe` that
   Windows/MinGW toolchains actually produce.
2. `OccamMesh.write()`/`OccamMesh.read()` treated the mesh control line's
   4th field as an air-layer count. In the real PW2D format that field is
   `nrfix` (the fixed-resistivity count, always `0` for pyCSAMT's own
   output) — writing the air-layer count there instead desynchronises
   every field the real Fortran reader expects afterward, and a real
   compiled Occam2D fails immediately with "Mesh file ended prematurely".
   K2's own tutorial never caught this because no session had a working
   Fortran toolchain available before; this dataset's smaller, faster
   build made it practical to compile the bundled solver and run it for
   real.

Neither bug affected users who never had `gfortran`/`make` available —
`InputBuilder.build()` itself, and every prior tutorial's dry-run
commands, were unaffected either way.

## Usage

```python
from pycsamt.stratagem import EDIBatch, StratagemRawReader, CoordinateInjector

batch = EDIBatch("data/stratagem/Dafang/df5-edi").fit()
rdr = StratagemRawReader("data/stratagem/Dafang/df5-hx", component="X").fit()
injector = CoordinateInjector(epsg=32649, order="forward").fit(
    batch.edi_objects_, "data/stratagem/Dafang/df5-coords.csv",
    easting_col="northing", northing_col="easting",
    elev_col="elev", station_col="stations",
)
```

See `docs/source/tutorials/process_stratagem_dafang_to_inversion.rst` for
the full field-to-inversion walkthrough, including static-shift
correction, frequency filtering, noise removal, the WATEX cross-check, and
a real Occam2D inversion run through convergence.
