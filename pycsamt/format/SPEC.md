# PCSF — pyCSAMT Common Subsurface Format

**Version:** 0.1.0
**Status:** Stable reference implementation, pre-1.0 (see [Versioning policy](#versioning-policy))
**Container:** [HDF5](https://www.hdfgroup.org/solutions/hdf5/) (via [`h5py`](https://www.h5py.org/)); a lossless ASCII sibling encoding, **PCSM** (`.pcsm`, see [§9](#9-pcsm-the-ascii-sibling-encoding)), is also part of this specification
**File extension:** `.pcsf` (binary); `.pcsm` (ASCII sibling)
**Reference implementation:** [`pycsamt.format`](./__init__.py) (`schema.py`, `io.py`, `text.py`, `adapters/`)

## 1. Purpose

Electromagnetic inversion produces a resistivity model, but every solver
disagrees about how to represent one on disk:

| Backend | Geometry | Resistivity encoding | Axis order |
|---|---|---|---|
| Occam2D | rectilinear, single profile | `log10(rho)` | `(n_z, n_x)` |
| ModEM 3-D | rectilinear, native 3-D volume | `ln(rho)` (header-tagged; can be `LOG10`/`LINEAR`) | `(nz, ny, nx)` |
| MARE2DEM | unstructured triangular mesh, per region | linear `rho`, per region (not per cell) | none (mesh) |
| DUHI (AI-inversion) | rides Occam2D once folded back via `map_ai_grid_to_occam` | — | — |
| Any other AI/DL model (UNet, GCN, ResNet, third-party) | whatever the model predicts on — grid2d/grid3d/mesh_unstructured, via `adapters.generic` | caller-chosen (`linear`/`log10`/`ln`) | caller-chosen, matching the geometry kind |

No shared on-disk representation connects any of these to a viewer.
Each downstream consumer — `app/mapview`, `app/web`'s 3-D map view, a
plain analysis script — has had to special-case the producing backend
before it can draw a single pixel, and no artifact persisting a true
3-D volume existed anywhere in pyCSAMT prior to this format.

PCSF is a single, self-describing, backend-neutral container that all
four backends' results convert *to*. It is the scientific artifact an
inversion run produces, not a fifth backend-specific format: Occam2D,
ModEM, MARE2DEM, and DUHI remain the systems that solve the inverse
problem; PCSF is what they hand to everything downstream of that.

## 2. Design principles

These are binding constraints on the format, not merely
recommendations — a reader or writer that violates them is non-conformant.

1. **Canonical resistivity is always linear ohm-m.** Stored in
   `model/resistivity`, regardless of the source backend's native
   encoding. A reader must never assume a file's resistivity is in any
   other unit or scale.
2. **Native encoding is preserved as provenance, never silently
   assumed.** If a file carries `model/resistivity_native`, its
   `model` group also carries an `encoding` attribute
   (`"log10"` | `"ln"` | `"linear"`); a value in `resistivity_native`
   without a declared `encoding` is invalid.
3. **Geometry is real-world, not solver-local.** `grid2d.x` is real
   station chainage (metres), never Occam2D's internal mesh-local
   frame; `grid3d` carries an explicit `origin`/`rotation_deg` when the
   source has one, not an implicit assumption of axis-aligned,
   zero-origin coordinates.
4. **Geometry kind is declared, never inferred.** A reader must read
   `geometry` group's `kind` attribute before interpreting any array
   inside it; array shape alone is not a reliable discriminator (a
   `(n,)` array means something different under `grid2d` than under
   `mesh_unstructured`).
5. **No derived 3-D views are mandatory.** Fence panels, isosurfaces,
   and block slices are computed client-side from the base geometry +
   model, exactly as `pycsamt.map.volume` and
   `app/web/callbacks/map3d.py`'s Plotly builders already do. The one
   documented exception is `multiline/derived_volume`, which is
   optional, explicitly tagged `synthesized: true` plus its
   `derivation_method`, and must never be mistaken for a native 3-D
   inversion result by a conformant reader.
6. **Unstructured meshes are preserved, not regridded.** A
   `mesh_unstructured` file keeps MARE2DEM's true triangular mesh
   (nodes + connectivity + region ids) at its native resolution. A
   regridded `grid2d` convenience export
   (`pycsamt.format.regrid.mesh_to_grid2d`) is a valid *additional*
   file a writer may produce, never a replacement for the mesh-native
   one — the returned model is explicitly marked
   `metadata["synthesized"] = True` so it can never be mistaken for a
   native inversion output (§2, principle 5's same rule, applied here).
7. **Zero new mandatory dependencies beyond `h5py`.** `h5py` is
   already pyCSAMT's only hard array-container dependency. PCSF is
   deliberately not built on `xarray`, `netCDF4`, `zarr`, `vtk`,
   `pyvista`, `GDAL`, or `rasterio` — none of those are core
   dependencies of pyCSAMT today, and the format must not raise the
   install bar for every `app/mapview`/`app/web` user just to read a
   result file.

## 3. Container structure

A `.pcsf` file is a single HDF5 file with the following top-level
layout. Group and dataset names below are exact — they are the names
`pycsamt.format.io` reads and writes, not illustrative examples.

```text
/                                   (root attrs)
    pcsf_version        str    e.g. "0.1.0"
    source_backend       str    "occam2d" | "modem3d" | "mare2dem" | "duhi" | "generic"
    created_by            str    free text (tool + version)
    created_at             str    ISO-8601 timestamp
    resistivity_unit        str    always "ohm.m"
    crs                       str    optional, pyproj-compatible
    description                str    optional

geometry/                          (attrs: kind)
    ... one of the four shapes in §4

model/
    resistivity                 dataset  canonical, always linear ohm.m
    resistivity_native            dataset, optional  source-backend passthrough
        (attrs: encoding = "log10" | "ln" | "linear", required if present)
    resistivity_by_region          dataset, optional  (mesh_unstructured only)
    uncertainty                     dataset, optional  same shape as resistivity
    sensitivity                      dataset, optional  same shape as resistivity

stations/                          optional
    name                         dataset[str]  (n,)
    x, y, z                        dataset[float64]  (n,) each
    line_id                          dataset[str], optional  (n,)
    lon, lat                          dataset[float64], optional  (n,) each
        WGS84 decimal degrees -- the one explicit, unambiguous carrier
        of real-world position in this group; ``x``/``y``/``z`` stay
        geometry-local (along-profile chainage for ``grid2d``, the
        model grid's own frame for ``grid3d``, whatever frame the
        writer used for ``mesh_unstructured``) and must never be
        assumed to be geographic coordinates even when this file's
        root ``crs`` attribute is set. Set together or not at all.

topography/                        optional  (attrs: kind = "per_station" | "raster")
    per_station/                       (kind == "per_station")
        station_id                 dataset[str]
        elevation                    dataset[float64]
    raster/                             (kind == "raster")
        x, y                        dataset[float64]  (n_x,), (n_y,)
        elevation                    dataset[float64]  (n_y, n_x)

survey_json                        dataset[str], optional  JSON blob
metadata_json                      dataset[str], optional  JSON blob

history/                           optional
    <key>                        dataset[float64]  e.g. "rms", "lambda", per iteration
```

`survey_json`/`metadata_json` are UTF-8 JSON strings (via
`json.dumps`/`json.loads`), not nested HDF5 groups — this keeps
free-form survey metadata (from `pycsamt.metadata`'s `SurveyMeta`,
`BBox`, `ProvenanceMeta`) and any adapter-specific extras
round-trippable without PCSF having to mirror every field of those
classes as its own typed schema.

### Compression

Arrays with `size >= 64` elements are written with gzip compression
(`compression_opts=4`); smaller arrays (a handful of stations, a short
history series) are stored uncompressed, since gzip's per-chunk
overhead outweighs the space saved below that threshold. This is an
implementation default, not a format requirement — a conformant reader
must not assume any particular compression state.

## 4. Geometry kinds

`geometry`'s `kind` attribute is one of four values. A reader must
branch on it before touching any dataset beneath `geometry/`.

### 4.1 `grid2d` — single-profile rectilinear geometry

Produced by Occam2D (and by DUHI, once folded back into an Occam2D
result — see §6).

| Path | Shape | Notes |
|---|---|---|
| `geometry/x` | `(n_x,)` | real station chainage, m |
| `geometry/z` | `(n_z,)` | depth cell centres, m, positive down |
| `geometry/x_nodes` | `(n_x+1,)`, optional | cell-edge coordinates |
| `geometry/z_nodes` | `(n_z+1,)`, optional | cell-edge coordinates |
| `geometry/origin` | `(2,)`, optional | real-world offset when `x` is locally referenced |
| attr `azimuth_deg` | scalar, optional | profile bearing |

`model/resistivity` shape: `(n_z, n_x)`.

### 4.2 `grid3d` — native 3-D tensor volume

Produced by ModEM 3-D. The only PCSF geometry kind representing a true
native 3-D inversion volume (as opposed to `multiline/derived_volume`,
which is a synthesized approximation — see §4.4).

| Path | Shape | Notes |
|---|---|---|
| `geometry/x`, `y`, `z` | `(n_x,)`, `(n_y,)`, `(n_z,)` | cell-centre coordinates, m |
| `geometry/x_nodes`, `y_nodes`, `z_nodes` | optional | cell-edge coordinates |
| `geometry/origin` | `(3,)`, optional | real-world grid origin |
| attr `rotation_deg` | scalar, default 0.0 | grid rotation about the vertical axis |
| attr `n_air` | int, default 0 | explicit air-layer count |

`model/resistivity` shape: **`(n_z, n_y, n_x)`** — deliberately kept
identical to ModEM's own native axis order rather than transposed to
some other convention, to avoid introducing a second axis-order bug on
top of the class of bugs a tensor-order mismatch can cause (cf. the
`em2d` TM-mode sign bug found elsewhere in this codebase — same
species of error, different subsystem).

### 4.3 `mesh_unstructured` — native triangular mesh

Produced by MARE2DEM. No forced regrid onto a tensor grid; the file
keeps the real element resolution the solver ran on.

| Path | Shape | Notes |
|---|---|---|
| `geometry/nodes` | `(n, 2)` or `(n, 3)` | node coordinates, m |
| `geometry/connectivity` | `(m, 3)`, int64 | triangle node indices |
| `geometry/region_ids` | `(m,)`, int32 | region id per triangle |
| attr `plane` | `"xz"` \| `"xy"` \| `"3d"`, default `"xz"` | physical plane; MARE2DEM profiles are conventionally `(y, z)` but stored generically as `"xz"` with `x` holding the profile's own along-line coordinate |

`model/resistivity` shape: **either** `(m,)` (per-triangle, expanded
from the per-region table via `region_ids`) **or** `(n_regions,)`
(per-region, collapsed) — both are valid; a reader must check which
shape it received rather than assuming per-triangle. The compact
per-region table, when the per-triangle array was expanded from one,
is additionally available at `model/resistivity_by_region`.

`model/resistivity_by_node`, shape `(n,)` (one value per mesh node),
is an optional companion for a graph-based source (e.g. a GCN) whose
native output is per-vertex rather than per-cell — the per-triangle
`model/resistivity` is then the arithmetic mean of each triangle's 3
vertex values (see `pycsamt.format.adapters.generic.mesh_to_pcsf`).
Only valid alongside a `mesh_unstructured` geometry.

### 4.4 `multiline` — formalized fence/stack of profiles

Formalizes what `app/web/callbacks/map3d.py` previously reconstructed
only at render time from a stack of independent 2-D sections.

```text
geometry/
    line_order              dataset[str]   ordered list of line_id values
    lines/<line_id>/
        geometry/            a nested grid2d geometry (§4.1), scoped to this line
        resistivity            dataset  (n_z, n_x) canonical linear ohm.m for this line
        attr offset_y            float     cross-line position, m
        attr offset_kind          "real" | "synthetic"
        attr azimuth_deg           float, optional
    derived_volume/            optional
        grid/                    a nested grid3d geometry (§4.2)
        resistivity                 dataset matching grid's (n_z, n_y, n_x)
        attr derivation_method       "linear_interp" | "kriging" | "idw"
        attr synthesized              always true
        derived_from                dataset[str]  line_id values the volume was built from
```

`model/resistivity` at the file root **must be absent** for
`multiline` files — each line carries its own resistivity under
`geometry/lines/<line_id>/resistivity` instead. A conformant writer
that emits both is non-conformant; a conformant reader must treat a
root-level `model/resistivity` alongside `kind="multiline"` as an
error, not silently prefer one.

### 4.5 Topography kinds

`topography` is independent of `geometry`'s kind — any geometry kind
may carry either topography kind, or none. A reader must branch on
`topography`'s own `kind` attribute the same way it branches on
`geometry/kind`.

- **`per_station`** — a scalar elevation per named station
  (`station_id` + `elevation`, both length `n`). Matches the
  pre-existing convention in `pycsamt.map.topo`.
- **`raster`** — a standalone gridded elevation surface, independent
  of any station table: `x` `(n_x,)`, `y` `(n_y,)`, and `elevation`
  `(n_y, n_x)` on the `(y, x)` meshgrid implied by `x`/`y` (the
  row-major convention `numpy.meshgrid(x, y)` produces by default).
  Construction is via `pycsamt.format.topography.topography_from_grid`,
  which takes plain arrays a caller already has — PCSF never parses a
  georeferenced raster *file format* (GeoTIFF, ASCII grid, ...)
  itself, so no GDAL/rasterio dependency is introduced by supporting
  this kind (see §2, principle 7): reading such a file into `x`/`y`/
  `elevation` arrays is the caller's own responsibility, by whatever
  means it likes, outside PCSF's dependency chain.

## 5. Versioning policy

`pcsf_version` is a semantic-versioning string (`MAJOR.MINOR.PATCH`)
written to every file's root attributes.

- **PATCH** (`0.1.0` -> `0.1.1`): documentation or implementation
  clarifications with no on-disk representation change. Existing files
  remain valid without modification.
- **MINOR** (`0.1.0` -> `0.2.0`): additive, backward-compatible
  changes only — a new optional dataset/attribute, a new geometry kind,
  a new `resistivity_native_encoding` value. A reader built against an
  older minor version must still be able to read the parts of the file
  it understands and ignore the rest. The `topography/kind="raster"`
  extension (§4.5), added after the initial 0.1.0 release, is the
  first real example of such a MINOR addition; `stations/lon`/`lat`
  (§3), added the same way afterwards, is the second; `model/resistivity_by_node`
  (§4.3), added for the generic AI/DL adapter, is the third — again
  without a literal version bump, following the same precedent.
- **MAJOR** (`0.x` -> `1.0`, or `1.x` -> `2.0`): any change that is not
  backward-compatible — renaming or removing a required field,
  changing an axis order, changing what "canonical" means. `0.1.0` is
  the initial, feature-complete-for-its-scope release; `1.0.0` is
  reserved for the point at which the schema is declared frozen for
  long-term archival compatibility, expected once the paper describing
  this format (see repository root's `PYCSAMT-PCSF-INVERSION-FORMAT-PLAN.md`,
  Phase 8) is accepted.

A reader must reject a file whose `pcsf_version` has a MAJOR component
it does not recognise, and may warn (but should not fail) on an
unrecognised MINOR component.

`read_pcsf` enforces this: it raises `ValueError` if `pcsf_version` is
missing or malformed, raises if the file's MAJOR component does not
match the reader's, and issues a `UserWarning` (without failing) if
only the MINOR component is newer. Since only `0.1.0` exists today,
every real file currently passes silently — the check exists so that
files written by a future MAJOR-incompatible writer fail loudly
instead of being silently misinterpreted.

## 6. Backend adapters

Each backend gets exactly one `to_pcsf()`-style adapter, added
alongside its existing result class rather than rewriting that class.

### 6.1 Smart topo attachment (`topo=`)

Every adapter below also accepts an optional `topo=` argument,
resolved by `pycsamt.format.topo_source.resolve_topo` into real
per-station `lon`/`lat`/elevation — the "smart" layer on top of the
lower-level `station_lonlat`/`station_elevations` mappings:

- **Accepted sources**: a `.bln` (Golden Software Surfer point list —
  no station identity, matched *positionally* in survey order,
  requiring an exact station-count match), a `.stn` (Zonge station
  file, parsed by the same reader `pycsamt.zonge` uses), a `.csv`/
  `.txt` with a recognisable header, an already-geo-located
  `Sites`/`MapData`/iterable-of-station-record object (e.g.
  `pycsamt.map.load_lines(edi_folder)`), or a plain `{name: (lon,
  lat[, elev])}` mapping. For `build_multiline_pcsf`: a
  `{line_id: <any of the above>}` mapping, or a sequence with one
  source per line, in line order.
- **Named vs. positional matching**: a source that carries station
  identity (`.stn`, a `.csv` with a station column, a dict, a
  Sites/MapData object) matches by id (exact, then
  `normalize_station_id`-normalized); a name-less source (a bare
  `.bln`) matches by position and therefore requires the topo point
  count to equal the expected station count — a mismatch raises
  (`on_mismatch="raise"`, the default) or warns and attributes only
  the overlapping prefix (`on_mismatch="warn"`).
- **Projected coordinates** (easting/northing, always true for
  `.stn`) need `epsg` or `utm_zone` to become lon/lat, via
  `pycsamt.gis.utils.to_ll` — the project's one existing UTM/EPSG
  utility; `pyproj` stays optional, imported lazily only when a
  conversion is actually requested.
- **Precedence**: when given, `topo` overrides any other lon/lat
  source (`station_lonlat`, or ModEM's own `.dat`-file
  `GG_Lat`/`GG_Lon`) per-station, with a `UserWarning` if both were
  supplied; a station `topo` has no data for keeps whatever the other
  source already gave it. Passing no `topo` leaves every adapter's
  pre-`topo` behaviour completely unchanged.
- **Multiline offsets ride the same mechanism already in place**:
  `build_multiline_pcsf` feeds `topo`'s resolved lon/lat into each
  line's own `sta_lat`/`sta_lon` *before* calling
  `line_offsets_from_stations` — the existing real-vs-synthetic
  offset computation (§6, `build_multiline_pcsf` bullet below) then
  runs unmodified. No separate offset step exists or is needed.

- `pycsamt.format.adapters.occam2d.occam2d_to_pcsf` — builds on
  `pycsamt.interp.ResistivityModel.from_occam2d`, which already
  recovers real station chainage from Occam2D's mesh-local coordinate
  frame. Produces `grid2d`. Occam2D itself has no real-world
  coordinate concept, so an optional `station_lonlat` mapping (mirrors
  `station_elevations`) is the only way this adapter's file becomes
  self-sufficiently geo-referenced (`stations/lon`/`lat`) — without
  it, placing a single Occam2D line needs a separate
  `known_stations` match at `load_pcsf_lines` time instead.
- `pycsamt.format.adapters.modem3d.modem3d_to_pcsf` — native `grid3d`,
  no synthesis. Converts `ln(rho)` to linear, carries `n_air` and any
  origin/rotation trailer the source `.rho` file has, and populates
  `stations/lon`/`lat` directly from the source `.dat` file's own
  `GG_Lat`/`GG_Lon` columns when present (`ModEmData.site_lonlat`) —
  no separate `known_stations` match needed either.
- `pycsamt.format.adapters.mare2dem.mare2dem_to_pcsf(result, mesh)` —
  produces `mesh_unstructured` from an explicit `TriMesh`, since
  `InversionResult` alone never loads mesh connectivity. Anisotropic
  models are rejected outright (`ValueError`) rather than silently
  dropping components — only isotropic support exists. Station
  identity/position is caller-supplied (no single reliable per-point
  name across MARE2DEM's MT/CSEM/DC variants); setting the supplied
  `StationTable`'s `lon`/`lat` is the caller's responsibility, the
  same field the other two adapters populate.
- `pycsamt.format.multiline.build_multiline_pcsf` — per-line
  `sta_lat`/`sta_lon` (when every line has them) already drove this
  builder's *real* cross-strike line-offset computation
  (`line_offsets_from_stations`, via `pycsamt.map.geometry.survey_uv`)
  before `stations/lon`/`lat` existed; the same values are now also
  persisted into the built file's `StationTable`, not just consumed
  transiently for the offset math.
- **DUHI (AI-inversion) has no separate adapter.** Its output only
  becomes a final resistivity model once folded back into an Occam2D
  run via `pycsamt.ai.inversion.mapping2d.map_ai_grid_to_occam`, so it
  rides the Occam2D adapter above; PCSF itself has no DUHI-specific
  code path.
- `pycsamt.format.adapters.generic.grid2d_to_pcsf` /
  `grid3d_to_pcsf` / `mesh_to_pcsf` — a solver-agnostic path for any
  *other* AI/DL inversion result (a UNet, a GCN, a ResNet, or a third
  party's own model) whose output is itself the final resistivity
  model, built from plain arrays with no dependency on any pycsamt
  result class. Same `topo=`/`station_lonlat`/`station_elevations`
  spatial mechanism as the adapters above; `origin`/`azimuth_deg`
  (`grid2d`) or `origin`/`rotation_deg` (`grid3d`) is this path's
  "offset" — real-world placement of an otherwise locally-referenced
  result. `mesh_to_pcsf` additionally accepts a per-node
  `resistivity_by_node` (the natural output shape of a graph-based
  model), projected onto per-triangle `resistivity` by averaging each
  triangle's 3 vertex values. Pairs with
  `pycsamt.format.provenance.ModelProvenance` for recording the
  model's architecture/framework/checkpoint/hyperparameters into
  `metadata['model_provenance']`, so a shared file is independently
  checkable, not just re-loadable — see `examples/ai_inversion_pcsf`.

## 7. Known limitations / deferred scope

These are deliberate, documented gaps, not oversights:

- **No embedded machine-readable schema (e.g. JSON Schema) is shipped
  in the file or the package yet.** This document is currently the
  only formal schema description; a self-describing, embeddable schema
  companion is tracked as pre-1.0 future work rather than included in
  0.1.0.
- **`stations/line_id` is optional even for `grid3d`/`grid2d`
  geometries**, required in practice only for `multiline` association.
  Whether to make it mandatory more broadly is an open question, not
  yet decided.
- **Uncertainty/sensitivity are opt-in, not mandatory-if-available.**
  An adapter may omit them even when the source `InversionResult`
  carries them; this is not currently enforced either way.

## 8. Conformance summary

A **conformant writer** must:

- set `pcsf_version`, `source_backend`, `resistivity_unit="ohm.m"` on
  the file root;
- store `model/resistivity` in linear ohm-m for every geometry kind
  except `multiline`, where it must be absent from the root and
  present per-line instead;
- declare `resistivity_native_encoding` whenever `resistivity_native`
  is present;
- declare `geometry`'s `kind` attribute before writing any
  kind-specific dataset.

A **conformant reader** must:

- read `geometry/kind` before interpreting any array beneath
  `geometry/`;
- never assume `model/resistivity`'s unit/scale beyond linear ohm-m;
- reject an unrecognised MAJOR `pcsf_version`;
- treat `multiline/derived_volume`, when present, as synthesized data,
  never as a native 3-D inversion result.

## 9. PCSM: the ASCII sibling encoding

**PCSF is canonical. PCSM is a derived, lossless, hand-editable text
projection of the same in-memory model — not a second schema.** The
relationship mirrors netCDF's own binary/text pair: netCDF is the
canonical binary container, and CDL (produced/consumed by `ncdump`/
`ncgen`) is its human-readable ASCII form. Every `.pcsm` file encodes
exactly the same `PCSFModel` a `.pcsf` file would; `pycsamt.format.text`
implements both directions.

### 9.1 Why a second encoding

PCSF (HDF5) is not something a user can open in a text editor or
hand-write a 20-line parser for. Occam2D's `.iter`, ModEM's `.rho`,
and MARE2DEM's own files, by contrast, are plain ASCII a user can read,
diff, hand-edit, and script against with nothing more than a text
editor and `str.split()`. PCSM exists so a PCSF-derived model has that
same property, without becoming a second canonical schema to keep in
sync with PCSF: the reference implementation reads/writes both
encodings of the identical `PCSFModel`/geometry dataclasses defined in
`schema.py` — no PCSM-specific data model exists.

### 9.2 Grammar

- `KEYWORD value` header lines encode scalars.
- `KEYWORD ... END_KEYWORD` blocks encode arrays. Values may be spread
  across any number of lines in any grouping; the `END_KEYWORD`
  terminator, not a declared count, ends the block, so a file stays
  easy to hand-edit — reflow lines, add blank lines, add or remove
  values freely. A reader compares the collected count against the
  count declared earlier (`NX`, `N_NODES`, ...) and raises a clear
  error on mismatch rather than silently misaligning subsequent data.
- `#` starts a comment — a whole line, or the remainder of a line
  after data — and is stripped before parsing. This is a deliberate
  improvement over the three native solver formats PCSM interoperates
  with, none of which support inline annotation.
- Six fields are free text and therefore exempt from comment-stripping
  — everything after the keyword to end of line is taken verbatim:
  `DESCRIPTION`, `CREATED_BY`, `CREATED_AT`, `CRS`, `SURVEY_JSON`,
  `METADATA_JSON`. A value in one of these that happens to contain `#`
  is not truncated.
- `multiline` geometry nests one `LINE_BEGIN <line_id> ... LINE_END
  <line_id>` block per line (itself containing a full `grid2d` block
  plus that line's own resistivity/offset fields) and an optional
  `DERIVED_VOLUME_BEGIN ... DERIVED_VOLUME_END` block wrapping a
  nested `grid3d`, mirroring `geometry/lines/<line_id>/` and
  `geometry/derived_volume/` in the HDF5 encoding (§4.4).
- Floats are written with Python's `repr()` — the shortest decimal
  string that round-trips to the exact same IEEE-754 value. A `.pcsm`
  file therefore round-trips a resistivity array **bit-exactly** back
  through `.pcsf`, unlike ModEM's own ASCII `.rho` format, which loses
  precision to ~5 significant figures.
- Within one block, every value is right-justified to that block's own
  widest formatted value, so columns line up visually — a writer-side
  cosmetic choice only; a reader parses by whitespace-splitting, so it
  is indifferent to padding, indentation, or line-wrapping width.
- Every resistivity block declares its own encoding as an inline
  comment on its header line, next to the data itself, not only in a
  separate header field — e.g. `RESISTIVITY  # linear ohm.m
  (canonical, ...)`, `RESISTIVITY_NATIVE  # source-native encoding:
  log10`. The canonical `RESISTIVITY` field is **always** linear
  ohm.m; this cannot be changed by any option (§2, principle 1) — a
  consumer relies on that invariant unconditionally, and reversing a
  `10**log10(x)` transform is not guaranteed bit-exact the way this
  format's `repr()`-based float round-trip already is. A writer may
  opt in (`write_pcsm(..., log10_view=True)`) to additionally emit a
  clearly-labelled `RESISTIVITY_LOG10` block (one per line, for
  `multiline`) — a convenience for reading resistivity in the log
  space many EM inversions (Occam2D among them) actually work in.
  This is write-only: `read_pcsm` discards it unconditionally, and it
  never becomes part of the returned `PCSFModel` or influences the
  canonical `RESISTIVITY` block in any way.

### 9.3 Versioning and conformance

PCSM shares PCSF's versioning policy exactly (§5) — a file's
`PCSM_VERSION` header is the same value a sibling `.pcsf` file's
`pcsf_version` root attribute would carry, checked by the identical
`pycsamt.format._version.check_pcsf_version` both encodings' readers
call. The conformance rules of §8 apply verbatim, substituting
`PCSM_VERSION`/`GEOMETRY_KIND` header lines for `pcsf_version`/
`geometry/kind` attributes.

### 9.4 Conversion

`pycsamt.format.text.pcsf_to_pcsm`/`pcsm_to_pcsf` convert between the
two encodings directly; both are thin compositions of the existing
read/write pairs (`write_pcsm(read_pcsf(path), out)` and vice versa) —
conversion is not a separate code path with its own failure modes.

### 9.5 Known limitations

- Large volumes (a native ModEM `grid3d`, a dense `mesh_unstructured`
  mesh) produce large text files — the same trade-off ModEM's own
  ASCII `.rho` format already accepts for the same reason, and not
  something a text encoding can eliminate while remaining plain,
  hand-editable text. A path ending in `.gz` (e.g. `model.pcsm.gz`) is
  written/read gzip-compressed, which mitigates this substantially in
  practice — measured on the real, bundled
  `willy_27freq_watex_line02_sample` 590,400-cell ModEM volume (see
  `examples/pcsm_conversion_demo`): 19.7 MB plain, 3.9 MB gzipped
  (5.0x smaller — real inverted resistivity compresses far better
  than synthetic/random test data would, and the column-aligned
  padding introduced for readability compresses away almost for free),
  vs. 4.1 MB for the
  equivalent `.pcsf` — but a gzipped file gives up casual text-editor
  inspection, so it stays opt-in rather than the default. PCSM's
  value proposition remains strongest for `grid2d`/`multiline` files
  and hand inspection generally; `.pcsf` remains preferable as a
  compact interchange form for large volumes.
- No embedded machine-readable grammar (e.g. an EBNF or PEG
  description) ships yet beyond this document and `text.py`'s module
  docstring.

## 10. Reference

- Full design rationale and phase-by-phase implementation history:
  `PYCSAMT-PCSF-INVERSION-FORMAT-PLAN.md` (repository root).
- User-facing usage guide with worked examples against real Occam2D,
  ModEM, and MARE2DEM datasets:
  `docs/source/user_guide/models/pcsf_format.rst`.
- Reference implementation: this package (`pycsamt.format`).
