# PCBH - pyCSAMT Common Borehole Format

**Version:** 0.1.0  
**Status:** stable pre-1.0 specification  
**Canonical encoding:** UTF-8 JSON  
**Canonical extension:** `.pcbh.json`

## 1. Purpose

PCBH is a self-contained exchange representation for one or more boreholes.
It stores borehole identity, collar position, measured trajectory, geological
logs, structural observations, shared vocabularies, units, and provenance.

PCBH is independent of an inversion model. PCSF stores subsurface models and
can embed or reference an unchanged PCBH document.

The keywords MUST, MUST NOT, SHOULD, SHOULD NOT, and MAY express normative
requirements.

This specification is citable as *PCBH 0.1 - pyCSAMT Common Borehole Format*,
pyCSAMT contributors, 2026. Cite the software release using the repository's
`CITATION.cff`; a release tag fixes the exact schema and specification text.

## 2. Encoding and version

- A PCBH 0.1 document MUST be a JSON object encoded as UTF-8.
- `pcbh_version` MUST equal `0.1.0`.
- JSON numbers MUST be finite. NaN and Infinity are invalid.
- Missing values MUST be omitted or encoded as JSON `null` where permitted.
- Unknown major versions MUST be rejected by future readers.
- `.pcbh.json` is the canonical extension. An extension is not proof of format.
- Explicit readers MAY open a file named `.pcbh`. Automatic detection MUST
  verify a JSON root with a compatible `pcbh_version`; `.pcbh` is only a hint.
- Extensions MUST appear in an `extensions` object whose top-level keys are
  namespaced, for example `organization:drilling_contract`. Readers and
  writers MUST preserve JSON-compatible extension values they do not
  interpret.

## 3. Root document

Required root members are:

| Member | Meaning |
|---|---|
| `pcbh_version` | PCBH semantic version |
| `document_id` | Stable document identifier |
| `created_at` | ISO-8601 creation timestamp |
| `created_by` | Creating application or organization |
| `crs` | Document coordinate reference |
| `units` | Document-wide units |
| `conventions` | Axis, depth, azimuth, and inclination rules |
| `boreholes` | Non-empty array of boreholes |

`document_id` and every borehole `id` MUST be non-empty. Borehole identifiers
MUST be unique and are case-sensitive. Producers SHOULD avoid identifiers that
differ only by case.

## 4. Coordinate reference system

All boreholes in a PCBH 0.1 document share one coordinate reference system.
`crs.horizontal` MUST be an EPSG identifier, another documented
pyproj-compatible value, or an explicit named local system such as
`LOCAL:mine-grid`. PCBH does not infer a missing CRS.

`collar.x` and `collar.y` use GIS XY order. `collar.z` is elevation positive
upward. Horizontal and vertical references are distinct. An unknown vertical
reference MUST be written explicitly as `unknown`; it MUST NOT be interpreted
as mean sea level.

Optional `longitude` and `latitude` are WGS84 decimal degrees and MUST occur
together. They do not replace projected X/Y.

## 5. Units and conventions

PCBH 0.1 declares units once for the document. Importers MUST convert mixed
source units before constructing the canonical document and SHOULD record the
conversion in provenance or an import report.

The following conventions are fixed:

- measured depth starts at the collar and increases along the borehole;
- collar/model elevation is positive upward;
- azimuth is clockwise from the declared north reference;
- borehole `inclination_deg` is measured from vertical downward;
- `inclination_deg = 0` is vertical down and `90` is horizontal;
- planar geological `dip_deg` is measured below horizontal;
- depth intervals are half-open: `[from_md, to_md)`.

## 6. Boreholes

Every borehole MUST provide:

- `id` and `name`;
- `kind` and `status`;
- a finite XYZ `collar`;
- `total_depth_md > 0`;
- a trajectory;
- zero or more interval log families and structural observations.

Standard kinds are `water`, `mining_exploration`, `mining_production`,
`geotechnical`, `environmental`, `petroleum`, `geothermal`, `scientific`,
`monitoring`, and `unknown`. Standard statuses are `planned`, `drilling`,
`completed`, `suspended`, `abandoned`, `decommissioned`, and `unknown`.
Extensions MUST be namespaced, for example `example:blast_hole`.

No survey station, interval, or structure may extend beyond
`total_depth_md`.

## 7. Trajectories

A vertical trajectory is represented by `method = vertical` and has no survey
stations. A surveyed trajectory uses `method = survey` and MUST contain at
least two stations sorted by strictly increasing measured depth. The first
station SHOULD be at `md = 0`; absence is a warning because a source adapter
may possess an explicit collar-to-first-station rule.

Survey station fields are:

- `md >= 0`;
- `azimuth_deg` in `[0, 360)`;
- `inclination_deg` in `[0, 180]`.

Values above 90 degrees describe upward/re-entry segments and SHOULD produce a
warning unless the workflow explicitly expects them. `north_reference` is
`true`, `grid`, `magnetic`, or `unknown`. Magnetic azimuth conversion requires
signed declination, evaluation date/epoch, source/model and version, evaluation
location, source and target north references, corrected value, and method.
Grid-north conversion additionally requires the target CRS and signed grid
convergence. These values belong in namespaced provenance. Without them the
source remains `magnetic`, validation warns, and readers MUST NOT invent a
correction.

Minimum curvature is the canonical PCBH 0.1 desurvey method. Derived XYZ paths
are not part of the source schema. They are generated deterministically with
X as easting, Y as northing, Z as elevation positive upward, and TVD positive
downward from the collar. A missing collar survey inherits the first measured
attitude, and the final measured attitude is extended to total depth. A
180-degree dogleg is rejected because it has no unique curvature plane.

Derived centerlines may be split at any measured-depth boundary. Boundary
points retain the requested MD exactly and are evaluated on the same
minimum-curvature arc. A SHA-256 source checksum covers only the collar, total
depth, and trajectory-defining fields so consumers can cache derived geometry
without treating display metadata as geometry.

### 7.1 Legacy geology compatibility

The legacy ``pycsamt.geology.Borehole`` is a vertical, profile-based view,
not an alternative PCBH representation. Promotion to PCBH therefore requires
an explicit absolute XYZ collar. Downgrade requires an explicit 2-D profile
distance; collar easting is never silently interpreted as that distance.

The compatibility view preserves half-open MD interval bounds, lithology
names, linear resistivity in ohm metres, and collar elevation. Embedded PCBH
lithology codes resolve through the document vocabulary. A ``RockDatabase``
may be supplied explicitly as a fallback classifier when a local code has no
name; the default conversion does not reinterpret recorded geology.

The legacy view cannot preserve the CRS and absolute XY collar, borehole ID
and classifications, deviated trajectory, additional log families,
structures, diameter, aliases, metadata, extensions, or detailed interval
metadata. ``legacy_conversion_losses`` reports the losses applicable to a
specific borehole. ``legacy_borehole_views`` supplies existing calibration
APIs with independent legacy objects and rejects documents whose depth and
resistivity units have not first been converted to metres and ohm metres.

## 7.2 Combined interval CSV

The combined CSV importer accepts one interval per row with repeated collar
fields. Its canonical mapping keys are dotted PCBH paths such as
``borehole.id``, ``collar.x``, ``interval.from_md``, and
``interval.lithology``. Callers may provide an explicit mapping and constants;
automatic recognition is restricted to the aliases documented by the public
mapping module. Every alias selection is recorded in the import report, and
ambiguous aliases require an explicit choice.

Delimiter detection considers only comma, semicolon, tab, and vertical bar.
Input must be UTF-8 with a header row and is subject to byte and row limits.
Missing tokens are the empty string, ``NA``, ``N/A``, ``NaN``, ``None``, and
``null``, compared case-insensitively. They become ``None`` and are never
converted into scientific labels.

Rows are grouped by borehole ID. Repeated X, Y, Z, CRS, kind, status, name, and
total depth must agree. Interval overlaps and collar conflicts reject the row;
there is no last-row-wins behavior. Missing total depth is inferred from the
deepest accepted interval and recorded. Lithology labels create reusable
document vocabulary entries, with deterministic codes when none are supplied.

Strict mode raises ``PCBHCSVImportError`` whenever an error is reported; its
``report`` attribute retains the complete diagnostics. Permissive mode returns
a valid document made only from accepted rows and leaves every rejected row,
conflict, mapping, default, inference, source checksum, and row count visible
in the accompanying ``ImportReport``.

## 7.3 Relational CSV projects

A relational project is controlled by a safe-loaded ``import.yaml`` manifest
with version ``0.1.0``, one horizontal CRS, canonical units, and a mapping of
logical table names to safe filenames in the same directory. ``collars`` is
required. ``surveys``, ``lithology``, ``structures``, ``samples``, and
``assays`` are optional. Absolute paths and parent traversal are rejected.

All child tables join through ``borehole_id``; assays join through
``sample_id``. Duplicate collars, duplicate sample IDs, duplicate
sample/analyte results, missing parents, invalid depths, and invalid PCBH
objects are reported rather than silently discarded. Export writes the same
manifest and only tables that contain records. File SHA-256 values are exposed
in the import report.

PCBH 0.1 stores normalized sample and assay records in the namespaced document
extensions ``pcbh:samples`` and ``pcbh:assays`` until first-class schema
objects are introduced. This representation round-trips through canonical
JSON and relational CSV without pretending that these records are interval
logs.

## 7.4 LAS subset

The LAS adapter supports unwrapped LAS 2.0 ASCII data. It preserves curve
mnemonics, descriptions, units, null samples, the shared MD index, and selected
well metadata under ``pcbh:continuous_curves`` and ``pcbh:las_metadata``.
Depth units ``M`` and ``FT`` are converted to canonical metres and conversions
are reported. Curve count, row width, finite depth, increasing depth, unique
mnemonics, and sample limits are validated.

When LITH is available it is grouped into PCBH lithology intervals. RESD is
used for representative interval resistivity only when its unit is explicitly
an ohm-metre spelling. LAS export writes inline curves when present and can
otherwise produce a small interval-derived DEPT/RESD/LITH subset. The returned
``LASExportReport`` lists written curves and unavoidable losses, including CRS,
trajectory, structures, non-lithology logs, samples/assays, and headers outside
the supported subset.

## 7.5 Viewer-neutral render model

The shared render builder produces typed Python primitives without importing
Plotly, Matplotlib, Dash, Qt, or a mesh backend. Its document-level contract
contains bounds, CRS/unit identity, per-borehole models, and material batches.
Each borehole contains a collar marker, sampled centerline, colored interval
segments, contact rings, and structure-glyph inputs. Trajectory points retain
MD, TVD, XYZ, azimuth, and inclination for hover and scientific inspection.

Geometry remains in the PCBH CRS. The builder does not transform coordinates
and rejects mismatched depth and coordinate units. Minimum-curvature paths are
sampled at a configurable MD step, with a per-hole vertex cap, and their
unsplit scientific trajectories are cached by source checksum. Display-only
Ramer-Douglas-Peucker simplification preserves interval boundaries exactly.

Vocabulary colors take precedence. Missing colors use a SHA-256-derived HSV
color, making fallbacks stable across processes and applications. Segments are
batched by log family and color for efficient backend rendering. Metadata on
collars, intervals, contacts, and structures contains selection state and
backend-neutral hover values; applications remain responsible for escaping
untrusted text before HTML display.

Display radius is explicitly view-only. ``auto`` combines model extent and
physical diameter, ``fixed`` uses a coordinate-unit radius, and
``exaggeration`` scales the physical radius or an automatic fallback. No mode
changes the borehole's scientific diameter.

## 8. Vocabularies and logs

Document vocabularies make files understandable without a local rock database.
Each entry has a unique `code`, a `name`, and optional `#RRGGBB` color,
description, external identifiers, and properties.

Standard interval families are:

- `lithology`;
- `formation`;
- `weathering`;
- `alteration`;
- `mineralization`;
- `oxidation`;
- `hydrostratigraphy`;
- `geotechnical`;
- `interpretation`.

Custom family names MUST be namespaced. Every interval MUST satisfy
`0 <= from_md < to_md` and MUST carry a code, label, or both. Codes in the
`lithology` and `formation` families MUST resolve in their document vocabulary.

Gaps are valid and mean unknown/unlogged material. Different families may
overlap. Intervals within one categorical family MUST NOT overlap in PCBH 0.1.
All nine standard interval-log families listed above are categorical and
exclusive within themselves. Continuous curves and structural/point
observations are not members of those interval families.

`data_nature` is `observed`, `interpreted`, `derived`, or `unknown`.
Resistivity-derived lithology MUST NOT be labeled `observed` automatically.

## 9. Structural observations

A structure is either a point (`at_md`) or a measured-depth zone
(`from_md`/`to_md`), never both. Its orientation representation is one of:

- `none`;
- `global_plane` using dip and dip direction, with optional strike;
- `global_line` using trend and plunge;
- `core_alpha_beta` using core-relative alpha and beta.

Global bearings lie in `[0, 360)`. Dip and plunge lie in `[0, 90]`. Alpha/beta
MUST NOT be interpreted as global orientation without the borehole orientation
at the observation depth and an explicit conversion method.

## 10. Validation

The reference implementation returns structured issues with `severity`,
`code`, `message`, `path`, and optional source information. Error-severity
issues make a document invalid. Warnings identify usable but ambiguous or
unusual states.

The JSON Schema validates representation-level constraints. Python semantic
validation additionally checks cross-object rules such as duplicate IDs,
vocabulary resolution, interval overlap, and total-depth bounds.

The reference reader applies configurable limits to file size, nesting,
borehole count, and interval count before returning an object. Duplicate JSON
object keys are invalid. The canonical writer validates by default, rejects
non-finite or non-JSON values, writes stable field order and UTF-8 text, and
atomically replaces the destination.

## 11. PCSF association

PCSF may embed an unchanged canonical PCBH object, reference one by URI and
SHA-256, or carry both. PCSF owns model geometry while PCBH retains its own CRS,
vertical reference, coordinates, and scientific observations. Alignment is a
derived view and MUST NOT mutate either source document. The normative HDF5
layout and coordinate conversion rules are specified in the PCSF specification.

## 12. Deferred features

Hydrogeology, construction, and PCBHZ attachments are specified by later
phases. Their absence from PCBH 0.1 MUST NOT be interpreted as removal from
the roadmap.

## 13. Visualization exports

PCBH JSON remains authoritative. GeoJSON, VTP, glTF, and GLB are derived
visualization views and return a ``PCBHExportReport`` enumerating omissions.

GeoJSON follows RFC 7946 and therefore transforms horizontal coordinates to
WGS84 longitude/latitude with XY order. Collar points and trajectory
LineStrings may carry a third coordinate. That Z value remains the PCBH
elevation in ``crs.coordinate_unit`` and ``crs.vertical``; RFC 7946 does not
standardize a vertical CRS, so consumers MUST inspect the accompanying
``z_reference`` and ``z_unit`` properties.

VTP exports triangulated interval tubes with RGB, measured depth, and material
identifier point arrays. glTF 2.0 and GLB export the same tube geometry as
indexed triangles with normals and vertex colors. Their coordinates remain in
the PCBH CRS, recorded in metadata/extras; neither format silently reprojects
or rebases local geometry. Log families other than the selected family,
structures, and non-geometric scientific fields are reported as losses rather
than discarded without notice.

## 14. Builder drafts and extension editors

The application builder uses a browser-safe flat draft containing project,
borehole/collar, survey, interval, and structure tables. Drafts are editing
state, not another exchange format. ``document_from_builder`` constructs and
validates the same ``PCBHDocument`` used by the Python API, while
``document_to_builder`` supports editing and session recovery.

Water, construction, sample, and assay tables are stored under namespaced
``pcbh:water``, ``pcbh:construction``, ``pcbh:samples``, and
``pcbh:assays`` borehole extensions until typed baseline classes are
standardized. CSV mapping profiles contain canonical-field/source-column
mappings and constants; they are browser-local configuration, not scientific
data.
