Jones Commands
==============

``pycsamt jones`` works with A.G. Jones J-format MT files through the
native ``pycsamt.jones`` I/O layer. Use it to inspect station metadata,
validate file structure, review block contents, and export filtered
subsets of a J-file collection.

The Jones command group does not perform J-to-EDI conversion. For that,
use the transform layer:

.. code-block:: console

   pycsamt transform j ...

Command Map
-----------

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - Command
     - Purpose
   * - ``pycsamt jones info J_OR_DIR``
     - Summarize one J-file or a directory of J-files.
   * - ``pycsamt jones validate J_OR_DIR``
     - Validate file extension and, by default, data-block structure.
   * - ``pycsamt jones stations J_DIR``
     - Print station coordinates, azimuths, frequency counts, and data flags.
   * - ``pycsamt jones blocks J_FILE``
     - Inspect data blocks in one J-file.
   * - ``pycsamt jones select J_DIR``
     - Filter a collection and export selected J-files.

Global Conventions
------------------

Most Jones commands accept:

``--format {text,json,csv}``
    Output format for table-like commands. ``validate`` and ``select`` use
    only ``text`` and ``json``.

``-v`` / ``--verbose``
    Increase logging and parser detail. For ``info`` in single-file mode,
    verbose output also lists each parsed block.

``--no-color``
    Disable colored logs.

Internally, single files are loaded with ``pycsamt.jones.j.JFile`` and
directories are loaded with ``pycsamt.jones.collection.JCollection``.

J-Format Model
--------------

A Jones file is one station-oriented text file containing metadata and
homogeneous data blocks. The pyCSAMT parser expects the common Jones
layout:

* optional banner/comment lines;
* ``>KEY=VALUE`` metadata lines such as latitude, longitude, elevation,
  and azimuth;
* one or more station/data blocks, each with a station line, a data-type
  token, a row count, and numeric rows.

The data-type token is split into a block kind and component. Examples:
``RXY`` is rho/phase for component ``XY``; ``ZXY`` is complex impedance
for component ``XY``; ``TZX`` or ``TZY`` are tipper/GDS-style transfer
function blocks.

``jones info``
--------------

``pycsamt jones info`` accepts either one file or a directory.

.. code-block:: console

   pycsamt jones info data/j/S01.j
   pycsamt jones info data/j/
   pycsamt jones info data/j/ --top 10
   pycsamt jones info data/j/ --format csv
   pycsamt jones info data/j/S01.j --format json

Single-File Output
~~~~~~~~~~~~~~~~~~

Single-file mode reports:

``station``
    Parsed station name, or the file stem as fallback.

``path``
    Input file path.

``n_freq`` / ``freq_min`` / ``freq_max``
    Number of frequency points and frequency range in Hz.

``has_z`` / ``has_r`` / ``has_t``
    Whether parsed impedance, rho/phase, or tipper data are available.

``lat`` / ``lon`` / ``elev`` / ``azimuth``
    Parsed station metadata when present.

``software`` / ``date``
    Banner provenance fields when the banner supplies them.

``blocks``
    JSON-only detailed list of block ``kind``, ``comp``, and ``nrows``.
    In text mode this list is shown when ``-v`` is supplied.

Directory Output
~~~~~~~~~~~~~~~~

Directory mode loads a ``JCollection`` and prints the collection summary.
The default columns are selected from ``station``, ``n_freq``, ``has_z``,
``has_r``, ``has_t``, ``lat``, ``lon``, and ``az``.

``--top`` limits directory output to the first ``N`` rows after loading.
It does not apply special sorting; for sorted station tables use
``jones stations``.

Options
~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 28 20 52

   * - Option
     - Default
     - Meaning
   * - ``J_OR_DIR``
     - required
     - Existing Jones file or directory.
   * - ``--top``, ``-n``
     - none
     - Limit directory-mode output to the first ``N`` stations.
   * - ``--format``
     - ``text``
     - Print text, JSON, or CSV.
   * - ``-v``, ``--verbose``
     - ``0``
     - Increase parser/log detail.
   * - ``--no-color``
     - false
     - Disable colored logs.

``jones validate``
------------------

Validate one file or all J-like files below a directory.

.. code-block:: console

   pycsamt jones validate data/j/
   pycsamt jones validate data/j/S01.j
   pycsamt jones validate data/j/ --no-deep
   pycsamt jones validate data/j/ --format json

Validation Targets
~~~~~~~~~~~~~~~~~~

For a single file, the command validates exactly that file. For a
directory, it recursively scans files whose suffix is one of:
``.j``, ``.dat``, or ``.txt``.

Validation Modes
~~~~~~~~~~~~~~~~

``--deep``
    Default. Calls ``pycsamt.jones.validation.is_j_file(..., deep=True)``.
    The validator scans the text and requires at least one complete data
    block: station line, data-type line, count line, and a numeric row
    matching the expected layout.

``--no-deep``
    Extension-only mode. Accepts Jones-like extensions supported by the
    validator: ``.j``, ``.jones``, ``.txt``, and ``.dat``. In directory
    mode the CLI still only discovers ``.j``, ``.dat``, and ``.txt``.

Exit Codes
~~~~~~~~~~

``validate`` exits with status ``0`` when all files pass. It exits with
status ``1`` when no files are found or any file fails validation. JSON
output has this shape:

.. code-block:: json

   {
     "n_files": 2,
     "n_ok": 2,
     "n_fail": 0,
     "mode": "deep",
     "results": [
       {"path": "S01.j", "valid": true, "error": null}
     ]
   }

Options
~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 28 20 52

   * - Option
     - Default
     - Meaning
   * - ``J_OR_DIR``
     - required
     - Existing Jones file or directory.
   * - ``--deep`` / ``--no-deep``
     - ``--deep``
     - Scan structure deeply or check extension only.
   * - ``--format``
     - ``text``
     - Print text or JSON.
   * - ``-v``, ``--verbose``
     - ``0``
     - Increase logging verbosity.
   * - ``--no-color``
     - false
     - Disable colored logs.

``jones stations``
------------------

Print a station-level table for a directory of J-files.

.. code-block:: console

   pycsamt jones stations data/j/
   pycsamt jones stations data/j/ --sort-by lat --format csv
   pycsamt jones stations data/j/ --sort-by n_freq --top 5
   pycsamt jones stations data/j/ --format json

Columns are:

``station``
    Station/site name.

``lat`` / ``lon``
    Geographic coordinates when present.

``az``
    Azimuth in degrees.

``n_freq``
    Number of parsed frequency points.

``has_z`` / ``has_r`` / ``has_t``
    Data-family flags for impedance, rho/phase, and tipper.

Sort keys are ``station``, ``lat``, ``lon``, ``elev``, ``az``, and
``n_freq``. The current displayed table does not include ``elev`` even
though it can be used as a sort key when available in the collection
summary.

Options
~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 28 20 52

   * - Option
     - Default
     - Meaning
   * - ``J_DIR``
     - required
     - Existing directory of J-files.
   * - ``--sort-by``
     - ``station``
     - One of ``station``, ``lat``, ``lon``, ``elev``, ``az``, ``n_freq``.
   * - ``--top``, ``-n``
     - none
     - Limit output to the first ``N`` sorted rows.
   * - ``--format``
     - ``text``
     - Print text, JSON, or CSV.
   * - ``-v``, ``--verbose``
     - ``0``
     - Increase parser/log detail.
   * - ``--no-color``
     - false
     - Disable colored logs.

``jones blocks``
----------------

Inspect block-level contents of one J-file.

.. code-block:: console

   pycsamt jones blocks data/j/S01.j
   pycsamt jones blocks data/j/S01.j --kind R
   pycsamt jones blocks data/j/S01.j --comp XY --qa
   pycsamt jones blocks data/j/S01.j --format json

Block Summary
~~~~~~~~~~~~~

Each output row contains:

``token``
    Concatenated block kind and component, such as ``RXY``.

``kind``
    Data family. Common values are ``R`` for rho/phase, ``Z`` for
    impedance, ``T`` for tipper/GDS transfer functions, and ``C`` for
    Schmucker C. The parser may also encounter other Jones-family tokens.

``comp``
    Component such as ``XX``, ``XY``, ``YX``, ``YY``, ``ZX``, or ``ZY``.

``nrows``
    Number of numeric rows in the block.

``period_min`` / ``period_max``
    Period range in seconds, rounded for display.

QA Mode
~~~~~~~

``--qa`` calls each block's ``qa_summary()`` method. JSON output nests
the returned dictionary under ``qa``. Text output summarizes rejected rows
as ``rejected/n`` when that information is available.

Filtering
~~~~~~~~~

``--kind`` and ``--comp`` are case-insensitive and can be combined. If no
blocks match, the command prints an error and exits with status ``1``.

Options
~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 28 20 52

   * - Option
     - Default
     - Meaning
   * - ``J_FILE``
     - required
     - Existing Jones file.
   * - ``--kind``, ``-k``
     - none
     - Filter by block kind, for example ``R``, ``Z``, ``T``, or ``C``.
   * - ``--comp``, ``-c``
     - none
     - Filter by component, for example ``XY`` or ``ZY``.
   * - ``--qa``
     - false
     - Include block QA summaries.
   * - ``--format``
     - ``text``
     - Print text, JSON, or CSV.
   * - ``-v``, ``--verbose``
     - ``0``
     - Increase parser/log detail.
   * - ``--no-color``
     - false
     - Disable colored logs.

``jones select``
----------------

Filter a directory collection and export the selected station files.

.. code-block:: console

   pycsamt jones select data/j/ --has-z --output-dir subset/
   pycsamt jones select data/j/ --stations S01,S02 --output-dir subset/
   pycsamt jones select data/j/ --pattern "KB0-*" --min-nfreq 20 --dry-run
   pycsamt jones select data/j/ --has-r --has-t --output-dir subset/
   pycsamt jones select data/j/ --dry-run --format json

Filters are combined with logical AND. With no filters, all stations in
the collection are selected.

Filters
~~~~~~~

``--stations ID[,ID...]``
    Keep exact station IDs after trimming whitespace.

``--pattern GLOB``
    Keep station names matching a Python ``fnmatch`` glob pattern, such as
    ``S0*`` or ``KB0-*``.

``--has-z``
    Require complex impedance data.

``--has-r``
    Require rho/phase data.

``--has-t``
    Require tipper data.

``--min-nfreq N``
    Require at least ``N`` frequency points. Minimum accepted value is ``1``.

Dry Run
~~~~~~~

``--dry-run`` lists matching station names and writes no files. In dry-run
mode, ``--output-dir`` is not required.

Export Mode
~~~~~~~~~~~

Without ``--dry-run``, ``--output-dir`` is required. The CLI constructs a
new ``JCollection`` from the selected stations and calls
``JCollection.export(output_dir)``. The summary reports:

``n_total``
    Total stations loaded from the source directory.

``n_selected``
    Number that passed all filters.

``n_written`` / ``n_failed``
    Export result counts.

``selected`` / ``written`` / ``failed``
    Station names and failure details.

If no stations match, the command exits with status ``1``. If all selected
exports fail, it also exits with status ``1``.

Options
~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 28 20 52

   * - Option
     - Default
     - Meaning
   * - ``J_DIR``
     - required
     - Existing directory of J-files.
   * - ``--output-dir``, ``-o``
     - current directory
     - Required unless ``--dry-run`` is used.
   * - ``--stations``, ``-s``
     - none
     - Comma-separated exact station IDs.
   * - ``--pattern``, ``-p``
     - none
     - Glob pattern for station names.
   * - ``--has-z``
     - false
     - Keep stations with impedance data.
   * - ``--has-r``
     - false
     - Keep stations with rho/phase data.
   * - ``--has-t``
     - false
     - Keep stations with tipper data.
   * - ``--min-nfreq``
     - none
     - Minimum number of frequency points.
   * - ``--format``
     - ``text``
     - Print text or JSON summary.
   * - ``--dry-run``
     - false
     - Preview matching stations without writing.
   * - ``--overwrite``
     - false
     - Accepted by the CLI, but export replacement behavior is delegated to
       ``JCollection.export``.
   * - ``-v``, ``--verbose``
     - ``0``
     - Increase parser/export detail.
   * - ``--no-color``
     - false
     - Disable colored logs.

Typical Workflows
-----------------

Validate and inventory a raw J folder:

.. code-block:: console

   pycsamt jones validate raw_j/
   pycsamt jones info raw_j/ --format json > j_summary.json
   pycsamt jones stations raw_j/ --format csv > j_stations.csv

Inspect a station with suspicious data:

.. code-block:: console

   pycsamt jones info raw_j/S01.j -v
   pycsamt jones blocks raw_j/S01.j --qa
   pycsamt jones blocks raw_j/S01.j --kind R --comp XY --format json

Create a clean subset:

.. code-block:: console

   pycsamt jones select raw_j/ --has-r --min-nfreq 20 --dry-run
   pycsamt jones select raw_j/ --has-r --min-nfreq 20 --output-dir clean_j/

Troubleshooting
---------------

``No J-files (.j / .dat / .txt) found``
    ``validate`` directory mode found no target files with the extensions
    it scans. Check the path or file suffixes.

``Unrecognized Jones J file structure``
    Deep validation found no complete station/data/count/row block. Try
    ``--no-deep`` only when you intentionally want extension-level checks.

``No J-files found``
    ``info``, ``stations``, or ``select`` loaded an empty collection.

``No blocks match kind=... comp=...``
    The block filters were too narrow, or the file does not contain the
    requested data family/component.

``--output-dir is required``
    ``select`` was run without ``--dry-run`` and without an output
    directory.

``No stations match the filter``
    All filters are combined with logical AND. Relax one filter or preview
    with ``--dry-run --format json``.

API Equivalents
---------------

Load one file:

.. code-block:: python

   from pycsamt.jones.j import JFile

   jf = JFile.from_file("raw_j/S01.j")
   print(jf.station, jf.n_freq, jf.lat, jf.lon)

Load and summarize a collection:

.. code-block:: python

   from pycsamt.jones.collection import JCollection

   coll = JCollection.from_sources("raw_j")
   rows = coll.summary()

Validate a file:

.. code-block:: python

   from pycsamt.jones.validation import is_j_file

   is_j_file("raw_j/S01.j", deep=True)

Filter and export:

.. code-block:: python

   from pycsamt.jones.collection import JCollection

   coll = JCollection.from_sources("raw_j")
   selected = coll.where(lambda jf: jf.Res is not None and (jf.n_freq or 0) >= 20)
   selected.export("clean_j")

