.. _models_pcsf_format:

PCSF — Common Subsurface Format
================================

Occam2D, ModEM, and MARE2DEM each write their own :term:`native file`\ s,
and each one disagrees with the others about almost everything that
matters for a downstream viewer: Occam2D stores resistivity as
:math:`\log_{10}(\rho)`, ModEM stores :math:`\ln(\rho)`, MARE2DEM stores
linear :math:`\rho` per region rather than per cell; Occam2D and ModEM
are rectilinear tensor grids with different axis orders, MARE2DEM is an
unstructured triangular mesh with no tensor grid at all; none of them
carry real-world station coordinates the way an EDI-derived survey does.
A tool that wants to plot "the resistivity model" therefore has to know
which backend produced it before it can draw a single pixel.

:mod:`pycsamt.format` is the answer: a single backend-neutral container,
:term:`PCSF`, that any of the three backends' results convert *to*.
Canonical resistivity in a PCSF file is always linear
:math:`\Omega\,\mathrm{m}`; the source backend's own encoding is kept
alongside it, explicitly labelled, never assumed. A PCSF file is a real
HDF5 container, so it can be opened with plain ``h5py`` or any HDF5
viewer as well as through pyCSAMT. Its ASCII sibling, :term:`PCSM`, is a
lossless, hand-editable projection of the same in-memory model. PCSF
remains the canonical interchange container; PCSM is useful for
inspection, version-control diffs, annotation, and small scripts that
should not depend on HDF5.

The adapters make the canonical resistivity invariant explicit. If
:math:`m` denotes a value in the backend's native model array, the
stored linear resistivity is

.. math::
   :label: eq-pcsf-canonical-resistivity

   \rho\;[\Omega\,\mathrm{m}] =
   \begin{cases}
      10^m, & \text{Occam2D native } \log_{10}(\rho),\\
      \exp(m), & \text{ModEM native } \ln(\rho),\\
      m, & \text{MARE2DEM native linear } \rho.
   \end{cases}

Consequently, a downstream reader never guesses a scale from the source
backend. ``model/resistivity`` is always linear
:math:`\Omega\,\mathrm{m}`; ``model/resistivity_native`` is optional and,
when present, is interpreted only together with its declared
``encoding`` attribute.

Four geometry kinds
--------------------

A PCSF file discriminates its geometry explicitly through
``geometry.kind`` rather than letting a reader infer it from array
shape.

.. list-table::
   :header-rows: 1
   :widths: 18 22 30 30

   * - Kind
     - Produced by
     - Resistivity shape
     - Notes
   * - ``grid2d``
     - Occam2D (and a DUHI-prepared result, folded back into Occam2D)
     - ``(n_z, n_x)``
     - Real station chainage, not the solver's mesh-local frame.
   * - ``grid3d``
     - ModEM 3-D
     - ``(n_z, n_y, n_x)``
     - Matches ModEM's own native axis order; carries a real-world grid
       origin and rotation when the source file has one.
   * - ``mesh_unstructured``
     - MARE2DEM
     - ``(n_triangles,)``
     - The real triangular mesh, never forced onto a rectilinear grid; a
       compact per-region table is kept alongside the per-cell
       expansion.
   * - ``multiline``
     - :func:`~pycsamt.format.multiline.build_multiline_pcsf`
     - one ``(n_z, n_x)`` array per line
     - A formalized fence/stack of independent ``grid2d`` lines, with
       real or synthetic cross-line offsets.

PCSF on disk
------------

The following ASCII tree maps the HDF5 hierarchy; it is not the content
of a text file. Items marked ``optional`` may be absent. A reader must
inspect the root version and ``geometry/kind`` attributes before it
interprets any array shape.

.. code-block:: text

   model.pcsf                         # one HDF5 file
   ├── attrs
   │   ├── pcsf_version = "0.1.0"
   │   ├── source_backend = "occam2d" | "modem3d" | "mare2dem" | ...
   │   └── resistivity_unit = "ohm.m"
   ├── geometry/
   │   ├── attrs: kind = "grid2d" | "grid3d" |
   │   │                 "mesh_unstructured" | "multiline"
   │   ├── x, z [, x_nodes, z_nodes]                 # grid2d
   │   ├── x, y, z [, x_nodes, y_nodes, z_nodes]     # grid3d
   │   ├── nodes, connectivity, region_ids           # mesh_unstructured
   │   └── lines/<line_id>/... [, derived_volume/]   # multiline
   ├── model/
   │   ├── resistivity                 # canonical linear Ω m
   │   ├── resistivity_native          # optional source encoding
   │   ├── resistivity_by_region       # optional mesh table
   │   ├── uncertainty                 # optional
   │   └── sensitivity                 # optional
   ├── stations/
   │   └── name, x, y, z [, line_id]                 # optional
   ├── topography/
   │   ├── attrs: kind = "per_station" | "raster"   # optional
   │   ├── per_station/{station_id, elevation}
   │   └── raster/{x, y, elevation}
   ├── history/<metric>                               # optional
   ├── survey_json                                    # optional UTF-8 JSON
   └── metadata_json                                  # optional UTF-8 JSON

For ``grid2d``, ``model/resistivity`` has shape ``(n_z, n_x)``. For
``grid3d``, its shape is ``(n_z, n_y, n_x)``. A multiline file is the
exception: it has no root ``model/resistivity`` dataset because each
``geometry/lines/<line_id>`` group owns its own ``(n_z, n_x)`` model.
An optional ``derived_volume`` is explicitly tagged as synthesized and
must not be interpreted as a native 3-D inversion.

Arrays containing at least 64 elements are currently written with gzip
level 4. Compression is an implementation default rather than a schema
requirement, so readers must accept compressed and uncompressed datasets.

.. figure:: /images/user_guide/models/pcsf_format_architecture.svg
   :align: center
   :width: 100%
   :alt: PCSF and PCSM architecture from solver-specific results through adapters, the shared PCSFModel, two encodings, and backend-neutral consumers.

   PCSF/PCSM architecture. Backend adapters normalize solver-specific
   geometry and resistivity into one validated ``PCSFModel``. PCSF and
   PCSM then serialize that same model as HDF5 or ASCII without changing
   its scientific meaning.

The central green block is the important boundary: backend knowledge ends
at the adapters. Everything to its right consumes explicit geometry and
linear resistivity rather than branching on Occam2D, ModEM, or MARE2DEM.
The dashed PCSF-to-PCSM connection denotes a lossless encoding conversion,
not another inversion or regridding step.

Converting an Occam2D result
------------------------------

:func:`pycsamt.format.adapters.occam2d.occam2d_to_pcsf` builds on
:meth:`pycsamt.interp.ResistivityModel.from_occam2d`, which already
recovers real station chainage from Occam2D's mesh-local coordinate
frame, and adds the cell-edge coordinates and iteration history PCSF
needs on top. The example below uses the bundled ``data/occam2D``
result (47 stations, Tongkeng CSAMT survey).

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.format import write_pcsf, read_pcsf
   >>> from pycsamt.format.adapters import occam2d_to_pcsf
   >>> from pycsamt.models.occam2d.results import InversionResult

   >>> result = InversionResult(workdir="data/occam2D")
   >>> model = occam2d_to_pcsf(result, description="Tongkeng CSAMT profile")
   >>> path = write_pcsf(model, "occam2d_run.pcsf")
   >>> restored = read_pcsf(path)
   >>> restored.kind
   'grid2d'
   >>> restored.resistivity.shape
   (31, 576)
   >>> print(f"{restored.resistivity.min():.1f}-{restored.resistivity.max():.1f} ohm.m")
   0.1-481684.9 ohm.m
   >>> len(restored.stations.name)
   47

``occam2d_to_pcsf`` is imported from :mod:`pycsamt.format.adapters`
directly in this example — that subpackage, and every other
second-level module (``schema``, ``io``, ``multiline``, ``topography``,
``pointcloud``), is importable right after ``import pycsamt.format``,
not only through an explicit ``from pycsamt.format.adapters.occam2d
import ...``.

The original log10 grid Occam2D actually wrote is not discarded; it
stays alongside the canonical linear array, explicitly labelled:

.. code-block:: pycon

   >>> restored.resistivity_native_encoding
   'log10'
   >>> import numpy as np
   >>> np.allclose(restored.resistivity, 10.0 ** restored.resistivity_native)
   True

Converting a ModEM 3-D result
--------------------------------

:func:`~pycsamt.format.adapters.modem3d.modem3d_to_pcsf` converts a
native 3-D ModEM volume the same way — this is the first genuinely
*persisted* 3-D resistivity volume in pyCSAMT; a web 3-D view built
purely from stacked 2-D sections only ever *synthesizes* one at render
time.

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.format.adapters import modem3d_to_pcsf
   >>> from pycsamt.models.modem.results import InversionResult

   >>> result = InversionResult(
   ...     workdir="data/modem/willy_27freq_watex_line02_sample", load_data=True
   ... )
   >>> model = modem3d_to_pcsf(result, description="Willy L18 line02, 27-freq")
   >>> model.kind
   'grid3d'
   >>> model.resistivity.shape
   (41, 50, 288)
   >>> model.geometry.origin
   array([-4509.828, -6752.725,     0.   ])
   >>> model.resistivity_native_encoding
   'ln'
   >>> np.allclose(model.resistivity, np.exp(model.resistivity_native))
   True

A real ModEM writer appends the grid's real-world centre and rotation
after the resistivity volume — ``model.geometry.origin`` above is that
real value, not a placeholder. Earlier versions of
:class:`~pycsamt.models.modem.model3d.ModEmModel3D` parsed that trailing
line and then silently discarded it; it is now exposed as
``ModEmModel3D.origin``/``.rotation``, matching the convention already
used by :func:`~pycsamt.models.modem.iotools.mackie.read_mackie3d` on
the same class.

Converting a MARE2DEM result
-------------------------------

MARE2DEM's own :class:`~pycsamt.models.mare2dem.results.InversionResult`
never loads a mesh — only the per-region ``.resistivity`` table — so
:func:`~pycsamt.format.adapters.mare2dem.mare2dem_to_pcsf` takes an
explicit :class:`~pycsamt.forward.maxwell.contracts_tri.TriMesh` second
argument rather than trying to resolve one on its own. The full
walkthrough, including rebuilding a real mesh in-process with the
``triangle`` Python package (already a hard pyCSAMT dependency — no
external Triangle binary required) from a run's own real ``.poly``
file, lives in ``examples/pcsf_conversion_demo/`` at the repository
root; running it reproduces the source run's real 6540-region partition
exactly.

Canonical resistivity there is the per-triangle array expanded from the
per-region table via each triangle's region id; the compact per-region
table itself is kept separately in ``resistivity_by_region`` for
provenance. MARE2DEM's own ``.resistivity`` file is already linear
:math:`\Omega\,\mathrm{m}`, so no native-encoding conversion applies —
``resistivity_native`` stays unset for this backend.

Stacking lines into a fence
------------------------------

:func:`~pycsamt.format.multiline.build_multiline_pcsf` formalizes what a
render-time-only 3-D fence view otherwise reconstructs from scratch on
every render: real per-line resistivity, plus a real-vs-synthetic
cross-line offset, computed from real station latitude/longitude when
available via :func:`pycsamt.map.geometry.survey_uv` (falling back to an
evenly-spaced synthetic stack otherwise).

.. code-block:: pycon
   :linenos:

   >>> import numpy as np
   >>> from pycsamt.format.multiline import build_multiline_pcsf

   >>> profiles = {
   ...     "L18": {"x": np.array([0.0, 250.0, 500.0]), "z": np.array([10.0, 100.0]),
   ...             "rho": np.array([[80.0, 120.0, 95.0], [40.0, 55.0, 45.0]])},
   ...     "L22": {"x": np.array([0.0, 250.0, 500.0]), "z": np.array([10.0, 100.0]),
   ...             "rho": np.array([[150.0, 160.0, 145.0], [70.0, 65.0, 72.0]])},
   ... }
   >>> model = build_multiline_pcsf(profiles, line_spacing=1.0)
   >>> model.kind
   'multiline'
   >>> [line.line_id for line in model.geometry.lines]
   ['L18', 'L22']
   >>> [line.offset_y for line in model.geometry.lines]
   [0.0, 1000.0]
   >>> [line.offset_kind for line in model.geometry.lines]
   ['synthetic', 'synthetic']

No station latitude/longitude was supplied here, so both lines fall
back to a synthetic 1000 m stack (``offset_kind="synthetic"``) rather
than mixing a real offset for one line with a synthetic one for the
other, which would misrepresent the survey's real geometry. An optional
cached ``derived_volume`` — each line resampled onto a shared grid, so a
large multiline file does not need to re-resample on every render — is
built by default (``cache_derived_volume=False`` to skip it), and is
explicitly tagged ``synthesized=True`` with its own
``derivation_method``, never presented as a native 3-D inversion.

Topography
------------

A PCSF ``topography`` group is scalar-per-station elevation, matching
the convention already used by :mod:`pycsamt.map.topo`.
:mod:`pycsamt.format.topography` wires the two together directly rather
than adding a third, independent parser:

.. code-block:: pycon

   >>> from pycsamt.map import load_lines
   >>> from pycsamt.format.topography import topography_from_map_data

   >>> data = load_lines("data/AMT/WILLY_DATA", detect="folder")
   >>> topo = topography_from_map_data(data)
   >>> len(topo.station_id)
   128
   >>> topo.station_id[0], float(topo.elevation[0])
   ('18-001A', 99.0)

The resulting :class:`~pycsamt.format.schema.TopographyPerStation` can
be passed straight into an adapter's own ``station_elevations``
parameter (available on both ``occam2d_to_pcsf`` and
``modem3d_to_pcsf``) via
:func:`~pycsamt.format.topography.topography_to_elev_map`. A station
without a known elevation is always recorded as ``nan``, never a
fabricated flat value.

A second topography kind, ``raster``, carries a standalone gridded
elevation surface instead of a per-station table.
:func:`pycsamt.format.topography.topography_from_grid` builds one from
plain ``x``/``y``/``elevation`` arrays — PCSF never parses a
georeferenced raster *file* itself (GeoTIFF, ASCII grid, ...), so this
introduces no GDAL/rasterio dependency; reading such a file into three
arrays is the caller's own responsibility, by whatever means it likes.
The example below grids the same real, EDI-derived ``WILLY_DATA``
station elevations used above onto a small lattice, using
:func:`pycsamt.map.geometry.equirect_xy` for the local metre projection
and :func:`scipy.interpolate.griddata` for the interpolation — both
pre-existing pycsamt/scipy functionality, not new logic:

.. code-block:: pycon
   :linenos:

   >>> import numpy as np
   >>> from pycsamt.map import load_lines
   >>> from pycsamt.map.geometry import equirect_xy
   >>> from pycsamt.format.topography import topography_from_grid
   >>> from scipy.interpolate import griddata

   >>> willy = load_lines("data/AMT/WILLY_DATA", detect="folder")
   >>> stations = [
   ...     s for s in willy.stations
   ...     if s.elevation is not None and s.latitude is not None
   ...     and s.longitude is not None
   ... ]
   >>> lat = [s.latitude for s in stations]
   >>> lon = [s.longitude for s in stations]
   >>> elev = np.asarray([s.elevation for s in stations])
   >>> x_pts, y_pts = equirect_xy(lat, lon)
   >>> xi = np.linspace(x_pts.min(), x_pts.max(), 10)
   >>> yi = np.linspace(y_pts.min(), y_pts.max(), 8)
   >>> grid_x, grid_y = np.meshgrid(xi, yi)
   >>> surface = griddata((x_pts, y_pts), elev, (grid_x, grid_y), method="linear")
   >>> surface = np.where(
   ...     np.isnan(surface),
   ...     griddata((x_pts, y_pts), elev, (grid_x, grid_y), method="nearest"),
   ...     surface,
   ... )
   >>> topo = topography_from_grid(xi, yi, surface)
   >>> topo.kind
   'raster'
   >>> topo.elevation.shape
   (8, 10)
   >>> print(f"{surface.min():.1f}-{surface.max():.1f} m")
   45.1-202.0 m

The result is an interpolated surface from 128 real station points, not
an independently-surveyed DEM — worth stating plainly wherever it is
used, the same way the illustrative Occam2D elevations earlier in this
page are labelled as illustrative rather than field-measured.

A terrain-draped Occam2D section
---------------------------------------

``occam2d_no_topo.pcsf`` carries no station elevation of its own -- the
Tongkeng ``data/occam2D`` run has no source EDI files to extract it from
(see :func:`~pycsamt.format.adapters.occam2d.occam2d_to_pcsf`'s
docstring). The same survey's original AVG station list
(``data/avg/k1.stn``) does carry real elevation, so it stands in as the
topography source, draped over the section with
:func:`pycsamt.topo.build_topo_section` and cropped to 1.5 km depth below
each station -- the same default ``OccamConfig.max_depth`` now used when
building a fresh Occam2D mesh:

.. figure:: /images/user_guide/models/pcsf_occam2d_topo_section.png
   :align: center
   :width: 100%

   The Occam2D ``grid2d`` PCSF section, cropped to 1.5 km depth and
   draped over real Tongkeng station topography, with every other
   station name labelled along the terrain (:mod:`pycsamt.topo`'s
   built-in label thinning keeps names legible without overlap). Colour
   uses a labelled 1st-99th percentile range so structure remains
   visible; compare values against the adjacent colour bar, not colour
   alone.

The section shows undulating terrain (roughly 400-575 m a.s.l.) above a
moderately resistive near-surface layer and a mixed conductive/resistive
sequence along the profile. The model is read through
:func:`~pycsamt.format.io.read_pcsf`, using its explicit grid nodes
(cropped and draped as described above).

One point cloud, every geometry kind
---------------------------------------

Two backends that store their geometry completely differently should
not need two different plotting code paths.
:func:`pycsamt.format.pointcloud.pcsf_to_point_cloud` flattens any of
the four geometry kinds into one ``(x, y, z, log10_rho)`` point cloud,
with a reproducible seeded subsample for a native ``grid3d``/
``mesh_unstructured`` file that can otherwise carry hundreds of
thousands of cells. It is the shared extraction behind pyCSAMT's first
desktop 3-D panel and is equally usable on its own -- a ModEM ``grid3d``
figure built from it will follow here once :mod:`pycsamt.app.mapview`
gains direct PCSF/PCSM support.

Consumers
-----------

A PCSF file is read directly by three different views, none of them
carrying backend-specific logic:

* :meth:`pycsamt.map.MapView.from_pcsf` builds real per-station 2-D
  curtains from a ``grid2d``/``multiline`` file the same way
  :meth:`~pycsamt.map.MapView.from_inversion_results` already does for
  a live ModEM folder, so the existing fence/depth-slice builders in
  :mod:`pycsamt.map.volume` render PCSF-sourced lines unmodified.
  ``grid3d``/``mesh_unstructured`` files raise ``NotImplementedError``
  there — a real per-station curtain has no way to slice a native
  volume without a station table walking a shared mesh.
* the web 3-D view's "PCSF file" data source currently accepts a
  ``multiline`` PCSF and reconstructs the profile dictionary used by its
  fence, block, and depth-slice views through
  :func:`~pycsamt.format.multiline.multiline_pcsf_to_profiles`.
* the desktop app's first 3-D/volume panel renders any ``.pcsf`` file
  through ``pcsf_to_point_cloud`` on a plain Matplotlib 3-D scatter,
  independent of the loaded EDI session.

.. important::

   The current web PCSF upload is not a general native-volume viewer. It
   rejects a single ``grid2d`` file and a native ModEM ``grid3d`` file
   because that callback is wired specifically to ``multiline`` profile
   reconstruction. The desktop PCSF viewer and
   ``pcsf_to_point_cloud`` accept both geometry kinds. Native ``grid2d``/
   ``grid3d`` web rendering is feasible, but requires a separate web
   point-cloud/section path rather than passing those files to the existing
   multiline callback.

PCSM: the browsable ASCII projection
-------------------------------------

PCSM stands for *pyCSAMT Common Subsurface Markup*. It serializes the
same :class:`~pycsamt.format.schema.PCSFModel` as PCSF, using keywords
and whitespace-delimited blocks instead of HDF5 groups and datasets.
There is no PCSM-specific scientific model and no additional conversion
of resistivity: the canonical ``RESISTIVITY`` block still follows
:eq:`eq-pcsf-canonical-resistivity`.

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Property
     - PCSF (``.pcsf``)
     - PCSM (``.pcsm`` / ``.pcsm.gz``)
   * - Representation
     - HDF5 groups, attributes, and typed datasets
     - UTF-8 keyword lines and delimited array blocks
   * - Best use
     - Compact interchange, applications, and large volumes
     - Inspection, diffs, comments, and hand editing
   * - Canonical resistivity
     - Linear :math:`\Omega\,\mathrm{m}`
     - Linear :math:`\Omega\,\mathrm{m}`
   * - Round trip
     - Reconstructs ``PCSFModel``
     - Reconstructs the same ``PCSFModel``
   * - Compression
     - Per-dataset HDF5 compression
     - Optional whole-file gzip with ``.gz``

The grammar is deliberately small. Scalar metadata uses
``KEYWORD value``. Arrays start with ``KEYWORD`` and terminate with
``END_KEYWORD``; values can be reflowed across lines because the
terminator, not line length, closes the block. ``#`` introduces a whole
line or trailing comment. The reader validates declared sizes such as
``NX`` and ``N_TRIANGLES`` against the values it collects, preventing a
hand edit from silently shifting the following block.

.. code-block:: text

   PCSM_VERSION 0.1.0
   SOURCE_BACKEND occam2d
   RESISTIVITY_UNIT ohm.m
   GEOMETRY_KIND grid2d
   NX 3
   NZ 2
   X_COORDS
   0.0 100.0 200.0
   END_X_COORDS
   Z_COORDS
   10.0 50.0
   END_Z_COORDS
   RESISTIVITY  # linear ohm.m (canonical)
   100.0 120.0 140.0
    50.0  60.0  70.0
   END_RESISTIVITY

Floats are written using a decimal representation that round-trips to
the same IEEE-754 value. ``write_pcsm(..., log10_view=True)`` can add a
clearly labelled ``RESISTIVITY_LOG10`` convenience block, but
``read_pcsm`` deliberately discards that derived view. It never replaces
or modifies canonical linear resistivity. ``RESISTIVITY_NATIVE``, by
contrast, preserves a real source-backend array and is accompanied by
``RESISTIVITY_NATIVE_ENCODING``.

Converting between the encodings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The conversion functions compose the normal public reader and writer;
there is no parallel conversion engine with different validation rules.
This example uses a real Occam2D-derived PCSF artifact from the bundled
conversion demo and verifies that the ASCII projection returns the same
geometry and canonical resistivity.

.. code-block:: pycon

   >>> import numpy as np
   >>> from pycsamt.format import pcsf_to_pcsm, read_pcsf, read_pcsm
   >>> source = "examples/pcsf_conversion_demo/output/occam2d_no_topo.pcsf"
   >>> text_path = pcsf_to_pcsm(source, "occam2d.pcsm")
   >>> binary_model = read_pcsf(source)
   >>> text_model = read_pcsm(text_path)
   >>> text_model.kind, text_model.resistivity.shape
   ('grid2d', (31, 576))
   >>> np.array_equal(text_model.resistivity, binary_model.resistivity)
   True

The reverse operation is symmetric:

.. code-block:: pycon

   >>> from pycsamt.format import pcsm_to_pcsf
   >>> restored_path = pcsm_to_pcsf("occam2d.pcsm", "occam2d_restored.pcsf")
   >>> restored_path.name
   'occam2d_restored.pcsf'

Real complete files
~~~~~~~~~~~~~~~~~~~

The expandable blocks below are generated from the real PCSF artifacts in
``examples/pcsf_conversion_demo/output`` by
``examples/pcsm_conversion/run_demo.py``. They are kept outside the page
so the narrative remains readable while users can inspect every value.

.. code-dropdown:: ../../../../examples/pcsm_conversion/output/occam2d.pcsm
   :language: text
   :linenos:
   :title: Browse the complete Occam2D grid2d PCSM file (602,093 bytes)

The Occam2D file contains the complete ``576 × 31`` model, its cell
centres and nodes, the source-native :math:`\log_{10}(\rho)` array, 47
stations, and iteration history. The linear and native blocks coexist so
consumers can use a stable canonical quantity without losing provenance.

.. code-dropdown:: ../../../../examples/pcsm_conversion/output/mare2dem.pcsm
   :language: text
   :linenos:
   :title: Browse the complete MARE2DEM mesh PCSM file (491,147 bytes)

The MARE2DEM file preserves 3,394 nodes, 6,774 triangular elements,
region identifiers, expanded per-element resistivity, and the compact
per-region table. It does not pretend that the unstructured mesh is a
rectilinear image; a consumer must respect its connectivity.

The real ModEM model contains ``41 × 50 × 288 = 590,400`` cells. Its
plain PCSM projection is 19,654,264 bytes, compared with about 4.1 MB for
PCSF, so embedding all values would make the documentation unnecessarily
heavy. This excerpt retains the actual headers and representative starts
and ends of every geometry, resistivity, and station block; each omission
is marked with ``...`` and the excerpt is intentionally not parseable.

.. code-dropdown:: ../../../../examples/pcsm_conversion/output/modem3d_excerpt.pcsm
   :language: text
   :linenos:
   :title: Browse the abbreviated real ModEM grid3d PCSM structure

Run the example with ``--full-modem`` to retain both the complete plain
file and ``modem3d.pcsm.gz``:

.. code-block:: console

   python examples/pcsm_conversion/run_demo.py --full-modem

For large native 3-D volumes, PCSF is normally the better interchange
choice. Gzip makes PCSM compact but removes its main advantage—casual
inspection in a text editor. Use plain PCSM when human readability is the
goal, compressed PCSM when text-schema portability matters, and PCSF for
routine application I/O.

Versioning and conformance
--------------------------

PCSF and PCSM share the same semantic version. A patch increment clarifies
documentation or implementation without changing representation; a minor
increment adds backward-compatible optional fields or geometry kinds; a
major increment permits incompatible schema changes. Readers reject a
missing, malformed, or unsupported major version and warn when a file's
minor version is newer than the reader understands.

A conforming writer declares the version, source backend, geometry kind,
and ``ohm.m`` unit, keeps canonical resistivity linear, and labels every
native encoding. A conforming reader branches on geometry kind before
interpreting arrays and treats a multiline derived volume as synthesized.
These rules are identical in both encodings; only their physical syntax
differs.

.. warning::

   PCSF and PCSM are at format version ``0.1.0`` (see
   :data:`pycsamt.format.schema.PCSF_VERSION`). The schema documented
   here has an explicit compatibility policy, but the pre-1.0 version
   still signals that the format is not yet frozen for long-term archival
   use. Preserve the original solver inputs and metadata alongside any
   PCSF/PCSM deliverable that must remain independently reproducible.
