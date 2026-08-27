# PCSF conversion demo — Occam2D / ModEM / MARE2DEM -> `.pcsf`

Converts real bundled inversion results from three backends into the
pyCSAMT Common Subsurface Format (PCSF) — the backend-neutral
inversion-result container described in
`PYCSAMT-PCSF-INVERSION-FORMAT-PLAN.md` (repository root). Each `.pcsf`
file is a real, self-describing HDF5 container: open one directly with
`h5py`, `h5dump`, or any HDF5 viewer.

Run from the repository root:

```bash
python examples/pcsf_conversion_demo/run_demo.py
```

## What gets converted

| Source | Real dataset | PCSF geometry kind |
|---|---|---|
| Occam2D | `data/occam2D` (47 stations, Tongkeng CSAMT) | `grid2d` |
| ModEM 3-D | `data/modem/willy_27freq_watex_line02_sample` (41x50x288 real inverted volume) | `grid3d` |
| MARE2DEM | `data/mare2dem/demo_mt_inversion` (6540 real regions) | `mesh_unstructured` |

**DUHI** (`pycsamt.ai.inversion`) has no separate adapter and so no
separate file here — its output only becomes a final resistivity model
once folded back into an Occam2D run, so a DUHI-prepared result
converts through the exact same `occam2d_to_pcsf()` path as the plain
Occam2D files below.

## What gets written

```
output/
  occam2d_no_topo.pcsf        grid2d, no topography group
  occam2d_with_topo.pcsf      grid2d, topo=.csv (2 illustrative stations, named match)
  occam2d_with_bln_topo.pcsf  grid2d, topo=.bln (47/47, positional match + UTM/EPSG)
  modem3d_no_topo.pcsf        grid3d, station z is ModEM's own flat 0.0
  modem3d_with_topo.pcsf      grid3d, topo=Sites (REAL lon/lat/elevation, 112/125 stations)
  mare2dem.pcsf                mesh_unstructured, real mesh + real per-region ρ
  demo-summary.json             machine-readable file-by-file summary
```

### Smart topo attachment (`topo=`)

Every adapter accepts a `topo=` argument, resolved by
`pycsamt.format.topo_source.resolve_topo` — see `SPEC.md` section 6.1
for the full design. This demo exercises all three of its source
kinds:

- **`occam2d_with_topo.pcsf`**: `topo/occam2d_topo.csv`, a **named**
  source (`.stn`/`.csv`/a plain dict/a Sites object all match this
  way) — matched by station id, with only 2 of the 47 real
  `data/occam2D` stations listed. Direct replacement for what used to
  be a hand-typed `station_elevations` dict.
- **`occam2d_with_bln_topo.pcsf`**: `topo/occam2d_topo.bln`, a
  **positional** source — a bare `.bln` carries no station identity,
  so `resolve_topo` first checks its point count against the 47
  expected stations (a count mismatch would raise), then attributes
  all 47 in survey order, converting the file's illustrative UTM
  zone 48N coordinates to lon/lat via `epsg=32648`.
- **`modem3d_with_topo.pcsf`**: `topo=` is a **Sites/MapData object**
  (`data/AMT/WILLY_DATA`, loaded via `pycsamt.map.load_lines`) passed
  directly, no intermediate file — real elevation *and* lon/lat this
  time (not elevation alone). ModEM's own station names carry a `23-`
  survey-year prefix `data/AMT/WILLY_DATA`'s own ids don't (same
  physical stations); a small renaming wrapper in `run_demo.py`
  bridges that (the one thing no generic tool can infer). `topo=`
  overrides the `.dat` file's own real `GG_Lat`/`GG_Lon` per station —
  the console output shows the `UserWarning` that precedence rule
  raises, deliberately let through rather than silenced.

All three `topo/*.csv`/`*.bln` files in this demo are explicitly
labelled illustrative in their own header comments — `data/occam2D`
has no independently-sourced real coordinate anywhere in the
repository, so these are fixtures demonstrating the mechanism, not a
claim of real Tongkeng CSAMT station positions.

MARE2DEM's `demo_mt_inversion` dataset has no per-station elevation
source, so `mare2dem.pcsf` carries no `topography/` group.

### The MARE2DEM mesh is real, not regridded

`InversionResult` never loads a mesh itself (only the per-region
`.resistivity` table), so this demo rebuilds the real triangular mesh
in-process from the run's own `demo.poly` PSLG using the `triangle`
Python package (already a hard pycsamt dependency — no external
Triangle binary needed). That reproduces the run's real region
partition **exactly**: all 6540 of 6540 unique region ids, no gaps —
confirming `demo.poly` really is the PSLG this run solved on, not an
approximation.

## Inspecting a file

The script itself prints every written file's full group/dataset tree
(shapes, dtypes, root attributes) right after writing it — scroll the
console output, or rerun and pipe to a file. Programmatically:

```python
from pycsamt.format import read_pcsf

model = read_pcsf("examples/pcsf_conversion_demo/output/modem3d_with_topo.pcsf")
model.kind                    # "grid3d"
model.resistivity.shape       # (41, 50, 288), linear ohm.m
model.geometry.origin         # real-world grid centre, metres
model.topography.station_id   # 112 real station ids with elevation
```

Or with `h5py` directly, no pycsamt import needed:

```python
import h5py
with h5py.File("examples/pcsf_conversion_demo/output/mare2dem.pcsf") as f:
    f.visititems(lambda name, obj: print(name, getattr(obj, "shape", "")))
```

## Round-trip checks

Every conversion is written then re-read, comparing resistivity (and,
where relevant, mesh connectivity / grid origin) for bit-exact equality
— printed as `round-trip check: ... bit-exact after write -> read` for
each backend.
