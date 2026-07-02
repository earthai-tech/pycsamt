Map Commands
============

``pycsamt map`` is the CLI entry point for station-location and map-display
workflows.  The first command surface is intentionally small and stable:
it resolves the same survey inputs used by the rest of the CLI, lists or
exports station coordinates, and saves a static station map figure.

The command group is backed by the code-first ``pycsamt.map`` package.
For the static figure path, the CLI uses ``pycsamt.map.ensure_map_data``,
``pycsamt.map.StationMapOptions``, and ``pycsamt.map.build_station_map``
with the matplotlib backend.  That keeps terminal workflows aligned with
the public map API while leaving room for richer interactive/basemap/3-D
commands to be added gradually.

Command Map
-----------

.. list-table::
   :header-rows: 1
   :widths: 24 36 40

   * - Command
     - Purpose
     - Main output
   * - ``pycsamt map stations``
     - List or export station coordinate rows.
     - Text table, JSON, CSV, or an output file.
   * - ``pycsamt map plot``
     - Save a static station longitude/latitude figure.
     - PNG/PDF/SVG or another matplotlib-supported image.

Mental Model
------------

The map commands are deliberately station-first.  They do not perform
inversion, interpolation, or profile modelling by themselves.  They answer
two practical questions:

1. Which station coordinates did pyCSAMT read from this survey?
2. Can I quickly make a map-like station layout figure from those
   coordinates?

Both commands resolve input in the same order:

1. Positional ``EDI_DIR`` argument.
2. ``--survey EDI_DIR`` option.
3. Active survey context set with ``pycsamt survey set``.

That means a mapping workflow can start from an explicit directory:

.. code-block:: console

   pycsamt map stations data/AMT/WILLY_DATA/L18PLT
   pycsamt map plot data/AMT/WILLY_DATA/L18PLT --output l18_station_map.png

Or from the active survey:

.. code-block:: console

   pycsamt survey set data/AMT/WILLY_DATA/L18PLT
   pycsamt map stations
   pycsamt map plot --output l18_station_map.png

Input Resolution
----------------

Every current ``map`` subcommand accepts an optional ``EDI_DIR``:

.. code-block:: console

   pycsamt map stations EDI_DIR
   pycsamt map plot EDI_DIR

When ``EDI_DIR`` is omitted, pass ``--survey`` or rely on the active
survey:

.. code-block:: console

   pycsamt map stations --survey data/line01
   pycsamt map plot --survey data/line01 --output line01_map.png

Use ``--fresh`` when the EDI files changed on disk and you want pyCSAMT
to re-parse them instead of using cached active-survey state:

.. code-block:: console

   pycsamt map stations --fresh
   pycsamt map plot --fresh --output updated_map.png

If no positional path, no ``--survey`` path, and no active survey are
available, the command fails with the same active-survey error used by
the other survey-aware CLI commands.

Station Coordinates
-------------------

Usage:

.. code-block:: console

   pycsamt map stations [EDI_DIR] [OPTIONS]

``map stations`` extracts station coordinates into rows with the
following fields:

.. list-table::
   :header-rows: 1
   :widths: 22 22 56

   * - Field
     - Type
     - Meaning
   * - ``index``
     - integer
     - Station order in the resolved survey/map data.
   * - ``station``
     - string
     - Station identifier.
   * - ``lat``
     - float or null
     - Latitude in decimal degrees.
   * - ``lon``
     - float or null
     - Longitude in decimal degrees.
   * - ``elev``
     - float or null
     - Elevation in metres when available.
   * - ``nfreq``
     - integer or null
     - Number of frequencies when available from the site summary.

The command first tries the normalized ``pycsamt.map`` data path.  If the
input is a lightweight site-like object that the map API cannot normalize
fully, it falls back to tolerant station summaries and ``site.coords``.
For normal EDI directory workflows, the public map API is the expected
path.

Options
-------

.. list-table::
   :header-rows: 1
   :widths: 30 22 48

   * - Option
     - Default
     - Meaning
   * - ``--survey EDI_DIR``
     - none
     - Explicit survey source when positional ``EDI_DIR`` is omitted.
   * - ``--fresh``
     - false
     - Re-parse the survey from disk before extracting coordinates.
   * - ``--drop-missing``
     - false
     - Drop stations without finite latitude and longitude.
   * - ``--top INT``
     - all rows
     - Print or export only the first ``INT`` station rows.
   * - ``--output FILE``, ``-o FILE``
     - stdout
     - Write rendered output to a file instead of the terminal.
   * - ``--format {text,json,csv}``, ``-f``
     - ``text``
     - Output format.
   * - ``-v``, ``--verbose``
     - quiet
     - Increase logging verbosity.
   * - ``--no-color``
     - false
     - Disable colored logging/output.

Text Output
-----------

The default output is a compact table:

.. code-block:: console

   pycsamt map stations data/line01

Example shape:

.. code-block:: text

    idx  station                    lat           lon        elev    nfreq
   -----------------------------------------------------------------------
      1  S01                  5.123456    -3.123456      80.00       27
      2  S02                  5.124000    -3.122900      82.00       27

Text output is best for a quick terminal check.  Use JSON or CSV when
another tool will consume the result.

JSON Output
-----------

Use JSON when scripting:

.. code-block:: console

   pycsamt map stations data/line01 --format json

Example shape:

.. code-block:: json

   [
     {
       "index": 1,
       "station": "S01",
       "lat": 5.123456,
       "lon": -3.123456,
       "elev": 80.0,
       "nfreq": 27
     }
   ]

JSON preserves missing coordinates as ``null``.

CSV Output
----------

Use CSV for spreadsheets, GIS prep, or simple command-line joins:

.. code-block:: console

   pycsamt map stations data/line01 --format csv --output station_coords.csv

The CSV columns are:

.. code-block:: text

   index,station,lat,lon,elev,nfreq

When ``--output`` is supplied, the command writes the rendered table to
the requested path and prints a short confirmation.  Parent directories
are created automatically.

Filtering Rows
--------------

Drop stations that do not have usable map coordinates:

.. code-block:: console

   pycsamt map stations data/line01 --drop-missing

Limit output while inspecting a large survey:

.. code-block:: console

   pycsamt map stations data/large_survey --top 20

Both options apply before output is rendered, so they affect text, JSON,
CSV, and file output consistently.

Static Station Plot
-------------------

Usage:

.. code-block:: console

   pycsamt map plot [EDI_DIR] [OPTIONS]

``map plot`` saves a static station map figure using station longitude on
the x-axis and latitude on the y-axis.  It requires at least one station
with finite latitude and longitude.

Internally, the command:

1. Resolves the survey with the shared CLI survey resolver.
2. Extracts coordinate rows and fails early if all map coordinates are
   missing.
3. Normalizes the same input through ``pycsamt.map.ensure_map_data``.
4. Builds a figure with ``pycsamt.map.build_station_map`` and
   ``StationMapOptions(backend="matplotlib")``.
5. Saves the matplotlib figure to ``--output``.

Options
-------

.. list-table::
   :header-rows: 1
   :widths: 30 22 48

   * - Option
     - Default
     - Meaning
   * - ``--survey EDI_DIR``
     - none
     - Explicit survey source when positional ``EDI_DIR`` is omitted.
   * - ``--fresh``
     - false
     - Re-parse the survey from disk before plotting.
   * - ``--output FILE``, ``-o FILE``
     - ``station_map.png``
     - Figure path to write.
   * - ``--title TEXT``
     - ``Station Map``
     - Figure title passed to ``StationMapOptions``.
   * - ``--label`` / ``--no-label``
     - ``--label``
     - Show or hide station-name annotations.
   * - ``--dpi INT``
     - ``150``
     - Saved figure resolution.  Must be at least ``50``.
   * - ``--show``
     - false
     - Open a matplotlib window after saving.
   * - ``-v``, ``--verbose``
     - quiet
     - Increase logging verbosity.
   * - ``--no-color``
     - false
     - Disable colored logging/output.

Basic examples:

.. code-block:: console

   pycsamt map plot data/line01 --output line01_station_map.png
   pycsamt map plot data/line01 --title "Line 01 stations" --dpi 200
   pycsamt map plot data/line01 --no-label --output clean_station_map.svg

Use the active survey:

.. code-block:: console

   pycsamt survey set data/line01
   pycsamt map plot --output figures/line01_station_map.png

Display the figure after saving:

.. code-block:: console

   pycsamt map plot data/line01 --output line01_station_map.png --show

Output Formats
--------------

``map plot`` calls matplotlib ``savefig``.  The output format is inferred
from the filename extension.  Common choices include:

* ``.png`` for reports and notebooks;
* ``.pdf`` for vector-friendly documents;
* ``.svg`` for web or documentation figures.

Examples:

.. code-block:: console

   pycsamt map plot data/line01 --output station_map.png
   pycsamt map plot data/line01 --output station_map.pdf
   pycsamt map plot data/line01 --output station_map.svg

Current Scope And Limitations
-----------------------------

The current CLI surface focuses on station coordinates and a static
station layout.  It does not yet expose every option available in
``pycsamt.map.StationMapOptions``.  For example, the Python API already
contains options for overlays, themes, map backends, color maps, profile
lines, contours, and richer interactive figures.

Use Python directly when you need those lower-level controls:

.. code-block:: python

   from pycsamt.map import StationMapOptions, ensure_map_data, build_station_map

   data = ensure_map_data("data/line01")
   fig = build_station_map(
       data,
       StationMapOptions(
           backend="plotly",
           overlay="elevation",
           theme="light",
           show_labels=True,
       ),
   )
   fig.write_html("line01_station_map.html")

The CLI can grow into those workflows incrementally.  Keeping the first
commands narrow makes them easier to test and keeps terminal output
predictable.

Workflow Recipes
----------------

Check coordinates before plotting:

.. code-block:: console

   pycsamt map stations data/line01 --drop-missing
   pycsamt map plot data/line01 --output figures/line01_map.png

Export coordinates for GIS cleanup:

.. code-block:: console

   pycsamt map stations data/line01 --format csv --output gis/line01_coords.csv

Use cached active survey state, then force a refresh after editing EDIs:

.. code-block:: console

   pycsamt survey set data/line01
   pycsamt map stations
   pycsamt map stations --fresh --format csv --output updated_coords.csv

Generate a quick figure set for multiple lines:

.. code-block:: console

   pycsamt map plot data/line01 --output figures/line01_map.png
   pycsamt map plot data/line02 --output figures/line02_map.png
   pycsamt map plot data/line03 --output figures/line03_map.png

Common Failures
---------------

``No active survey``
    No positional ``EDI_DIR`` was supplied, no ``--survey`` was supplied,
    and no active survey was configured.  Provide a path or run
    ``pycsamt survey set EDI_DIR`` first.

``No stations with finite latitude and longitude were found``
    ``map plot`` found no usable coordinates.  Run
    ``pycsamt map stations --format csv`` to inspect the coordinate
    columns.  Use ``pycsamt site edit --set-coords`` or a coordinate table
    workflow to repair missing EDI header coordinates.

``matplotlib is required``
    Static plotting needs matplotlib.  Install the plotting dependencies
    in the active environment.

``Error writing FILE``
    The output path could not be written.  Check permissions, parent
    directory spelling, and whether the target file is open in another
    application.

Practical Notes
---------------

Run ``map stations`` before ``map plot`` on a new survey.  The coordinate
table is the fastest way to catch missing, swapped, or placeholder
coordinates before making figures.

Use ``--drop-missing`` when exporting to GIS tools that expect every row
to have valid longitude and latitude.

Keep station-coordinate exports near the figures they generated.  A CSV
beside the image makes it easier to audit exactly what was plotted later.

Prefer explicit input paths in automation:

.. code-block:: console

   pycsamt map stations data/line01 --format json
   pycsamt map plot data/line01 --output figures/line01_map.png

Active survey context is convenient interactively, but explicit paths are
clearer in scripts and batch jobs.
