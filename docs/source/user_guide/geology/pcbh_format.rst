.. _geology_pcbh_format:

PCBH — Common Borehole Format
==============================

A field project rarely contains only one simple vertical log. Mining,
groundwater, geotechnical, and geothermal work may combine many collars,
deviation surveys, geological intervals, structures, samples, assays, water
observations, and construction records. A table can carry part of that
information, but it cannot reliably state how every depth, coordinate, code,
and relationship should be interpreted.

The pyCSAMT Common Borehole Format (PCBH) is the human-readable exchange
contract for that complete project. A canonical ``*.pcbh.json`` file contains
one or more boreholes and remains independent of an inversion model. It can be
read without pyCSAMT as ordinary UTF-8 JSON and checked against the public JSON
Schema at ``https://pycsamt.org/schemas/pcbh/0.1/schema.json``.

This page belongs in the geology guide because PCBH describes observed and
interpreted borehole evidence. :doc:`../models/pcsf_format` describes PCSF,
which stores a subsurface resistivity model. The two formats can be associated
for 3-D viewing, but neither is silently converted into the other.

.. important::

   :class:`pycsamt.geology.Borehole` is the established lightweight vertical
   log used by profile calibration. PCBH is the multi-hole spatial exchange
   format. Adapters connect them, so existing calibration workflows do not
   need to be rewritten. See :doc:`borehole` for the lightweight API.

What the document fixes explicitly
----------------------------------

PCBH removes several common sources of geological ambiguity:

* one horizontal CRS and vertical reference apply to the whole document;
* collar ``x``, ``y``, and ``z`` are spatial coordinates, with elevation
  positive upward;
* measured depth (MD) starts at the collar and increases along the hole;
* survey inclination is measured from vertical down, so 0 degrees is vertical,
  90 degrees is horizontal, and values above 90 describe upward/re-entry
  segments;
* intervals are ``[from_md, to_md)`` and may contain gaps, which mean
  unlogged or unknown material rather than zero thickness;
* observed, interpreted, derived, and unknown data are distinguished through
  ``data_nature``;
* geological codes are defined in document vocabularies, with optional color,
  description, physical properties, and external identifiers;
* unknown vertical datum is written as ``"unknown"`` rather than being
  treated as mean sea level.

The core document has the following logical shape. Extensions are preserved
when a reader does not understand them, provided their keys are namespaced.

.. code-block:: text

   project.pcbh.json
   ├── pcbh_version, document_id, title, provenance
   ├── crs                         # horizontal + vertical reference
   ├── units and conventions
   ├── dictionaries
   │   ├── lithologies
   │   └── formations
   ├── boreholes[]
   │   ├── id, name, kind, status
   │   ├── collar                  # x, y, z
   │   ├── total_depth_md, diameter
   │   ├── trajectory/stations[]   # md, azimuth, inclination
   │   ├── interval_logs/*[]
   │   ├── structures[]
   │   └── extensions
   └── extensions

Kinds such as ``water``, ``mining``, ``geotechnical``, ``geothermal``,
``petroleum``, ``monitoring``, and ``exploration`` share this core. Domain
records that have not yet been standardized as core classes use namespaced
extensions such as ``pcbh:water``, ``pcbh:construction``, ``pcbh:samples``,
and ``pcbh:assays``. An extension never changes the meaning of a core field.

Reading and validating a project
--------------------------------

The packaged minimal example is a vertical 80 m water borehole in UTM zone
29N with two lithology intervals. The following session uses only public APIs
and shows the actual object returned by the reference reader:

.. code-block:: pycon

   >>> from importlib import resources
   >>> from pycsamt.format import read_pcbh
   >>> example = resources.files("pycsamt.format.borehole").joinpath(
   ...     "examples/minimal-vertical.pcbh.json"
   ... )
   >>> document = read_pcbh(example)
   >>> document.document_id, len(document.boreholes)
   ('example:minimal-vertical', 1)
   >>> hole = document.boreholes[0]
   >>> hole.id, hole.kind, hole.total_depth_md
   ('BH-001', 'water', 80.0)
   >>> [(item.from_md, item.to_md, item.code)
   ...  for item in hole.interval_logs["lithology"]]
   [(0.0, 12.0, 'SOIL'), (12.0, 80.0, 'GRAN')]

``read_pcbh`` checks UTF-8 JSON, duplicate keys, compatible version, object
shape, finite numbers, resource limits, and semantic relationships. Semantic
validation catches problems that JSON Schema alone cannot express, including
duplicate borehole IDs, unknown vocabulary codes, overlapping categorical
logs, non-increasing survey stations, and observations below total depth.

.. warning::

   Do not call ``read_pcbh(..., validate=False)`` for routine exchange. That
   option supports diagnostics and migration of damaged historical data; it
   does not make an invalid document scientifically safe.

Writing canonical JSON is symmetric and atomic:

.. code-block:: pycon

   >>> from pathlib import Path
   >>> from tempfile import TemporaryDirectory
   >>> from pycsamt.format import write_pcbh
   >>> with TemporaryDirectory() as directory:
   ...     output = write_pcbh(document, Path(directory) / "project.pcbh.json")
   ...     print(output.name, output.exists())
   project.pcbh.json True

The writer validates before replacing the destination and emits indented,
stable UTF-8 JSON. ``.pcbh.json`` is canonical. An explicit reader may open a
file called ``.pcbh``, but software must confirm its JSON content and version
rather than trust the suffix.

Building a document without Python
----------------------------------

The web application's **Interpret → Borehole Builder** page exposes the same
contract through project/CRS fields and editors for collars, surveys,
intervals, structures, water, construction, samples, and assays. Validation
messages link back to the affected row. The 2-D log and 3-D trajectory previews
are derived from the current draft; they do not replace validation.

The builder can download canonical PCBH JSON or embed it in an uploaded PCSF
model. Drafts are recovered in browser-local storage, while exported scientific
files remain explicit user actions. CSV mapping profiles are also local editing
preferences and are not silently written into geological observations.

Importing repeated borehole rows from CSV
-----------------------------------------

The combined CSV importer is intended for the common layout in which collar
values repeat on each geological interval row. It recognizes conventional
headers, or accepts an explicit mapping from canonical fields to source
columns:

.. code-block:: pycon

   >>> from pycsamt.format import boreholes_from_csv, write_pcbh
   >>> document, report = boreholes_from_csv(
   ...     "drilling_intervals.csv",
   ...     columns={
   ...         "borehole.id": "HoleID",
   ...         "collar.x": "Easting",
   ...         "collar.y": "Northing",
   ...         "collar.z": "Elevation",
   ...         "interval.from_md": "From",
   ...         "interval.to_md": "To",
   ...         "interval.lithology": "Lithology",
   ...     },
   ...     constants={"crs.horizontal": "EPSG:32629"},
   ... )
   >>> write_pcbh(document, "drilling_intervals.pcbh.json")

The returned report records the source checksum, resolved mappings,
inferences, accepted and rejected rows, and localized warnings/errors. Repeated
collars and total depths must agree within a hole. Missing values remain
missing; an unavailable elevation is never silently replaced with zero.

Use ``strict=False`` only when a workflow deliberately accepts valid rows from
a partially damaged source and archives the complete rejection report. For a
relational exchange with separate collars, surveys, logs, and structures,
use :func:`~pycsamt.format.borehole.write_csv_directory` and
:func:`~pycsamt.format.borehole.boreholes_from_csv_directory`; their manifest
retains identity, CRS, units, checksums, and table relationships.

Trajectories, intervals, and structures
---------------------------------------

A vertical hole uses ``trajectory.method = "vertical"`` and needs no survey
stations. A deviated hole uses at least two stations with strictly increasing
MD. :func:`~pycsamt.format.borehole.desurvey` applies minimum curvature and
returns deterministic MD, true vertical depth, XYZ, azimuth, and inclination
points. Source azimuth remains tied to its declared ``true``, ``grid``,
``magnetic``, or ``unknown`` north reference; pyCSAMT never invents magnetic
declination or grid convergence.

The standard categorical log families are lithology, formation, weathering,
alteration, mineralization, oxidation, hydrostratigraphy, geotechnical, and
interpretation. Intervals cannot overlap within one family, while different
families may overlap because they describe different properties. Continuous
LAS curves are depth-indexed values rather than categorical intervals and are
preserved inline under ``pcbh:continuous_curves``.

Structural observations can represent a point or a depth zone. Global planes
use dip and dip direction; global lines use trend and plunge. Core-relative
alpha/beta values remain explicitly tagged and must not be treated as global
orientation until combined with the borehole attitude at the observation MD.

Using PCBH with legacy calibration
----------------------------------

Promotion from :class:`pycsamt.geology.Borehole` requires an absolute PCBH
collar. The legacy ``x`` value is profile chainage, not an easting, and is
retained only as ``metadata["legacy_profile_x"]``. Downgrade requires the
inverse choice explicitly:

.. code-block:: pycon

   >>> from pycsamt.format.borehole import legacy_borehole_views
   >>> views = legacy_borehole_views(
   ...     document,
   ...     profile_x={"BH-001": 425.0},
   ... )
   >>> views[0].name, views[0].x, len(views[0].intervals)
   ('BH-001', 425.0, 2)

This explicit mapping prevents an absolute UTM easting from being passed to a
2-D profile calibrator as chainage. A legacy view is necessarily lossy: it
cannot retain Y, CRS, deviation, multiple log families, structures, or most
provenance.

Associating boreholes with PCSF
-------------------------------

:func:`~pycsamt.format.borehole.embed_pcbh` places a validated PCBH document
inside a copied PCSF model. :func:`~pycsamt.format.borehole.reference_pcbh`
instead stores a URI and SHA-256 checksum, and the association may carry both.
The PCBH evidence and PCSF inversion remain independently identifiable.

Before a viewer inserts trajectories into a 3-D block,
:func:`~pycsamt.format.borehole.align_pcbh_to_pcsf` transforms horizontal
coordinates into the PCSF frame, applies model origin and rotation, and
reports whether each path is inside, outside, or intersects the block. If the
vertical references cannot be proven compatible, alignment stops until the
caller supplies an explicit ``vertical_offset``.

.. warning::

   A visually plausible collar is not evidence of correct vertical alignment.
   Record the vertical datum or a justified offset before interpreting contacts
   against model depth cells.

The shared render contract produced by
:func:`~pycsamt.format.borehole.build_render_model` supplies centerlines,
colored interval segments, contacts, structure glyphs, material batches, and
scientific hover values to both application front ends. Physical diameter and
minimum visible display radius remain separate so a narrow hole stays visible
without falsifying its true size.

Interchange and visualization exports
-------------------------------------

PCBH JSON remains authoritative. The other outputs serve narrower consumers:

.. list-table::
   :header-rows: 1
   :widths: 18 35 47

   * - Output
     - Best use
     - Important limitation
   * - LAS 2.0
     - One hole's depth-indexed logs
     - Cannot preserve the complete multi-hole project and rich CRS model
   * - relational CSV
     - Database and spreadsheet exchange
     - Requires its manifest to retain relationships and units
   * - GeoJSON
     - WGS84 collars and trajectory map features
     - Vertical CRS and detailed logs are not native GeoJSON semantics
   * - VTP
     - Scientific 3-D geometry and scalar arrays
     - Not a complete borehole archive
   * - glTF/GLB
     - Browser-ready colored 3-D presentation
     - Presentation geometry is not the source geological contract

Every visualization exporter returns a
:class:`~pycsamt.format.borehole.PCBHExportReport`. Applications should expose
its loss records rather than imply that a derived file is a lossless PCBH
replacement.

Versioning and independent use
------------------------------

PCBH is currently version ``0.1.0``: a stabilized pre-1.0 contract. Readers
must reject unknown major/minor versions rather than guess. The schema,
packaged reference fixture, semantic validator, governance rules, and reader
guide are distributed with pyCSAMT; the normative specification is available
as ``pycsamt/format/borehole/SPEC.md`` in the source and installed package.

The public API reference is under :doc:`../../api/format`. Extension authors
should also read ``GOVERNANCE.md`` in the PCBH package: public extension keys
use ``owner:name``, while the ``pcbh`` namespace is reserved for project-owned
fields.
