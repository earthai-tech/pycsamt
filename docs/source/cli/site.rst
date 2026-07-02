Site Commands
=============

``pycsamt site`` operates on pyCSAMT ``Sites`` collections. It is the
higher-level survey command group for inspection, filtering, editing,
recomputing, exporting, and computing derived geophysical quantities from
EDI-backed stations.

Use ``pycsamt edi`` for low-level file checks. Use ``pycsamt site`` once
you want to work with a survey as stations.

The command group is a CLI wrapper around the ``pycsamt.site`` layer. It
resolves EDI files into a ``Sites`` collection, applies station-level
selection or editing helpers, and then either prints a report, writes EDI
files, or computes derived quantities from the underlying EDI objects.

Command map
-----------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Command
     - Purpose
   * - ``pycsamt site info [EDI_DIR]``
     - Report survey or single-station metadata.
   * - ``pycsamt site select [EDI_DIR]``
     - Filter stations and optionally write the subset.
   * - ``pycsamt site edit [EDI_DIR]``
     - Rotate, slice, fill, or re-coordinate sites.
   * - ``pycsamt site export [EDI_DIR]``
     - Write EDIs or a zip archive from a Sites collection.
   * - ``pycsamt site recompute [EDI_SOURCE]``
     - Normalize and rewrite foreign or raw EDI files.
   * - ``pycsamt site compute strike``
     - Estimate geoelectric strike.
   * - ``pycsamt site compute resistivity``
     - Extract apparent resistivity at a target frequency.
   * - ``pycsamt site compute phase-slope``
     - Compute a phase-slope depth proxy.
   * - ``pycsamt site compute tipper``
     - Compute tipper magnitude summary statistics.

Mental model
------------

``pycsamt site`` sits between raw EDI file utilities and full processing
pipelines:

``pycsamt edi``
    File-oriented inspection and validation. Use it when the question is
    about one file or EDI syntax.

``pycsamt site``
    Survey-oriented inspection, filtering, editing, export, and lightweight
    geophysical calculations. Use it when the question is about stations.

``pycsamt pipe``
    Reproducible multi-step processing workflows with QC plots, reports, and
    a run record.

Most ``site`` commands never modify the input folder directly. Commands that
write data require ``--output-dir`` and write new EDI files there. Use
``--overwrite`` only when replacing files in that output directory is
intended.

Shared options
--------------

The exact option set varies by sub-command, but most ``site`` commands share
these flags:

``-S, --survey EDI_DIR``
    Use this EDI directory or file as the survey source. If omitted, the
    active survey context is used when available.

``--fresh``
    Re-parse the survey from disk instead of using cached survey context.

``-f, --format text|json|csv``
    Select terminal output format for commands that emit tabular results.

``-o, --output-dir DIR``
    Directory where output files are written. Commands that transform or
    export EDIs generally require it unless ``--dry-run`` is used.

``--overwrite``
    Permit replacement of existing output files when the writer supports it.

``-v, --verbose`` and ``--no-color``
    Increase progress messages or disable ANSI color.

Input resolution
----------------

Most ``site`` commands resolve input in this priority order:

1. Positional ``EDI_DIR`` or ``EDI_SOURCE``.
2. ``--survey DIR``.
3. Active survey context set with ``pycsamt survey set``.

Examples:

.. code-block:: console

   pycsamt site info ./edis/
   pycsamt site info --survey ./edis/
   pycsamt survey set ./edis/
   pycsamt site info

Use ``--fresh`` on commands that expose it to reload from disk instead of
using cached survey state.

Inspect sites
-------------

.. code-block:: console

   pycsamt site info
   pycsamt site info survey/
   pycsamt site info --station HBH03
   pycsamt site info --top 20
   pycsamt site info --format json

Collection mode reports survey-level metadata and a station table. Single
station mode prints a detailed station report. CSV output for a station
contains its resistivity/phase table.

Behavior:

* ``--station`` matching is case-insensitive.
* When a requested station is not found, the command prints the available
  station names and exits nonzero.
* ``--top`` affects only the human-readable collection report. JSON and CSV
  still serialize the report data returned by ``SitesReport``.
* Collection CSV output is a station-level table. Single-station CSV output
  is the station's resistivity/phase table.

Use cases:

* confirm that the survey loaded as the expected number of stations;
* identify station names before using ``select`` or ``compute``;
* inspect one station's frequency and response table;
* generate a quick JSON inventory for automation.

Select stations
---------------

.. code-block:: console

   pycsamt site select --stations S01,S02,S03
   pycsamt site select --freq 0.1:1000 --output-dir filtered/
   pycsamt site select --bbox 27.0,101.0,29.0,103.0
   pycsamt site select --stations S01,S05 --freq 1:100 --drop-empty --output-dir subset/
   pycsamt site select --freq 0.1:1000 --dry-run

Filters:

``--stations S01,S02``
    Keep named stations. The parser accepts comma-separated identifiers. The
    selection helper supports flexible name matching, but the CLI path is
    usually used with explicit station IDs.

``--freq FMIN:FMAX``
    Keep stations that contain at least one finite frequency in this closed
    interval. This is a station selector, not a row slicer. Use
    ``site edit --freq`` when you want to trim each site's frequency rows.

``--bbox LAT_MIN,LON_MIN,LAT_MAX,LON_MAX``
    Keep stations inside a geographic box.

``--drop-empty``
    Remove stations with no data arrays.

``--keep-finite``
    Remove stations whose impedance tensor has no finite values.

``--phase-err-thresh DEGREES``
    Keep only stations whose maximum finite phase error is less than or equal
    to the threshold. Stations without a phase-error array are kept.

Without ``--output-dir``, the command prints a summary only. With
``--output-dir``, it writes filtered EDI files.

Output behavior:

* If no station remains after filtering, the command exits nonzero.
* ``--dry-run`` prints what would be selected and writes nothing.
* ``--format json`` emits ``n_input``, ``n_selected``, ``stations``, and, when
  files are written, ``written``.
* ``--format csv`` emits a one-column station table.
* The output directory is created when writing is requested.

Selection order is deterministic: filters are applied in this order:
station names, frequency membership, bounding box, empty-station removal,
finite-Z filter, and phase-error threshold.

Edit sites
----------

.. code-block:: console

   pycsamt site edit --rotate 30 --output-dir rotated/
   pycsamt site edit --freq 0.1:1000 --fill-missing --output-dir edited/
   pycsamt site edit --set-coords 28.5,102.3,200 --output-dir recoord/
   pycsamt site edit --coords-table coords.csv --output-dir recoord/
   pycsamt site edit --rotate 15 --freq 1:100 --dry-run

At least one edit operation is required. ``--output-dir`` is required
unless ``--dry-run`` is used. Available operations are tensor rotation,
frequency slicing, missing-value filling, rho/phase recomputation,
uniform coordinates, and per-station coordinate tables.

Edit operations:

``--rotate DEGREES``
    Rotate impedance tensors, and tipper vectors when present, in the
    horizontal plane.

``--freq FMIN:FMAX``
    Slice each site to the requested frequency interval and keep aligned
    frequency-indexed arrays together.

``--fill-missing``
    Fill missing/non-finite impedance values using the site editing helper.
    This is useful before computations that require complete arrays.

``--recompute-rho-phase``
    Recompute apparent resistivity and phase from the impedance tensor after
    earlier edits.

``--set-coords LAT,LON[,ELEV]``
    Apply one coordinate triple to every station.

``--coords-table FILE``
    Read per-station coordinates from a CSV/Excel-like table. The table should
    provide station, latitude, longitude, and optionally elevation columns.

Operational rules:

* At least one edit operation is required.
* ``--output-dir`` is required unless ``--dry-run`` is used.
* Edits are applied in CLI order: rotate, frequency slice, fill missing
  values, recompute rho/phase, set uniform coordinates, then table
  coordinates.
* The input survey is not overwritten. Edited EDIs are written to the output
  directory.

Export sites
------------

.. code-block:: console

   pycsamt site export --output-dir final_edis/
   pycsamt site export --template "{index:03d}_{station}.edi" --output-dir final_edis/
   pycsamt site export --zip --zip-name survey_export.zip --output-dir .
   pycsamt site export --manifest --output-dir final_edis/

For stable station-based filenames, use ``--template "{station}.edi"``.
The writer supports these placeholders:

``{station}``
    Resolved station name.

``{index}``
    Zero-based iteration index. Python format specs are supported, for
    example ``{index:03d}``.

``{lat}``, ``{lon}``, ``{elev}``
    Header coordinates when available, otherwise NaN.

``{chainage}``
    Optional profile chainage when available, otherwise NaN.

If the rendered filename does not end with ``.edi``, the extension is added.
Use ``--overwrite`` to replace existing files. ``--manifest`` writes a
``manifest.csv`` next to exported EDIs with
``index,station,lat,lon,elev,chainage,filename,path`` columns.

Zip behavior:

* ``--zip`` writes one archive instead of loose EDI files.
* ``--zip-name`` defaults to ``survey.zip``.
* Existing archives are protected unless ``--overwrite`` is supplied.

Recompute EDIs
--------------

``recompute`` normalizes foreign or raw EDI files and rewrites pyCSAMT
EDIs. It accepts a single file, a line folder, a survey folder, or the
active survey.

.. code-block:: console

   pycsamt site recompute WILLY_DATA/L18PLT --rotate 30
   pycsamt site recompute WILLY_DATA --flatten --template "{line}_{source_stem}.edi"
   pycsamt site recompute --survey WILLY_DATA --output-dir cleaned_edis
   pycsamt site recompute raw_edis/ --freq 0.1:1000 --fill-missing nan --dry-run

Useful options:

``--output-dir DIR``
    Root output directory. If omitted, a ``recomputed_edis`` folder is
    created near the input.

``--output-name NAME``
    Folder name used when ``--output-dir`` is omitted. The default is
    ``recomputed_edis``.

``--flatten``
    Write all EDIs directly under the output root instead of preserving
    line folders.

``--template TEMPLATE``
    Output filename template. Supported placeholders are ``{station}``,
    ``{index}``, ``{line}``, and ``{source_stem}``.

``--rotate DEGREES`` and ``--components Z,Tip``
    Rotate transfer functions before recomputing.

``--freq FMIN:FMAX``
    Keep only a frequency band.

``--fill-missing zero|nan``
    Fill missing impedance/tipper values before recomputing.

``--skip-rho-phase``
    Do not recompute apparent resistivity and phase from Z.

``--datatype mt|emap``
    Force writer datatype.

``--synthesize-spectra``
    Ask the EDI writer to synthesize spectra when possible.

``--no-manifest``
    Do not write ``manifest.csv``.

``--dry-run``
    Recompute in memory and report counts without writing files.

``--progress``
    Show a progress bar while recomputing.

Source and output rules:

* ``recompute`` accepts one EDI file, one line folder, a survey folder
  containing line folders, ``--survey``, or the active survey.
* If a line folder is passed and ``--output-dir`` is omitted, output is written
  under the parent survey folder as ``recomputed_edis/<line>/``.
* If a survey folder is passed and ``--flatten`` is not used, line folders are
  preserved under the recomputed output root.
* Any failed record causes a nonzero exit status after the summary is printed.

Compute quantities
------------------

Strike:

.. code-block:: console

   pycsamt site compute strike
   pycsamt site compute strike --method phase_diff --format json
   pycsamt site compute strike --format csv > strike.csv

Methods are ``swift`` and ``phase_diff``.

Output columns:

``station``
    Station name.

``method``
    Strike method used for the row.

``theta_deg``
    Estimated strike angle in degrees.

``swift`` performs a Swift-style off-diagonal power minimization on a
one-degree grid from 0 to 179 degrees. ``phase_diff`` is a simpler
phase-difference heuristic that returns 0 or 90 degrees for degraded or sparse
data.

Apparent resistivity at a frequency:

.. code-block:: console

   pycsamt site compute resistivity --freq 1.0
   pycsamt site compute resistivity --freq 10.0 --how interp
   pycsamt site compute resistivity --freq 0.1 --format csv > rho.csv

``--how`` is ``nearest`` or ``interp``.

``nearest``
    Uses the closest native frequency and reports it as ``f_used``.

``interp``
    Computes apparent resistivity at native frequencies and linearly
    interpolates to the requested frequency.

Output columns are ``station``, ``res_xy``, ``res_yx``, and ``f_used``.
``res_xy`` and ``res_yx`` are apparent resistivities for the off-diagonal
impedance components.

Phase slope:

.. code-block:: console

   pycsamt site compute phase-slope
   pycsamt site compute phase-slope --band 0.1:1000
   pycsamt site compute phase-slope --band 1:100 --format csv > slope.csv

``--band FMIN:FMAX`` restricts the regression interval. If omitted, the CLI
tries to infer the full frequency range from the loaded sites. Pass the band
explicitly when working with sparse or unusual data.

Output columns are ``station``, ``slope_xy``, and ``slope_yx``. Slopes are
phase-versus-log-frequency slopes in degrees per decade for the off-diagonal
components.

Tipper:

.. code-block:: console

   pycsamt site compute tipper
   pycsamt site compute tipper --freq 1.0
   pycsamt site compute tipper --format json

The current command emits summary statistics for all available tipper rows.
The accepted ``--freq`` option is reserved for frequency-specific evaluation,
but the CLI implementation currently calls the summary calculation. Output
columns are ``station``, ``mean``, ``median``, and ``max``.

Machine-readable formats
------------------------

All compute subcommands use the same emitter:

``--format text``
    Prints a pandas-style table.

``--format json``
    Prints records-oriented JSON.

``--format csv``
    Prints CSV with a header row.

The compute layer unwraps ``Site`` objects to their underlying EDI objects
before calling ``pycsamt.site.compute`` functions.

Typical workflow
----------------

.. code-block:: console

   pycsamt survey set data/AMT/WILLY_DATA/L18PLT
   pycsamt site info
   pycsamt site select --keep-finite --drop-empty --dry-run
   pycsamt site edit --rotate 30 --freq 0.1:1000 --output-dir cleaned/
   pycsamt site compute strike cleaned/ --format csv > strike.csv
   pycsamt site export cleaned/ --manifest --output-dir final_edis/

Troubleshooting
---------------

No active survey
    Pass an explicit ``EDI_DIR``, use ``--survey EDI_DIR``, or run
    ``pycsamt survey set EDI_DIR`` first.

Need to see new disk edits
    Add ``--fresh`` so the survey is re-parsed from disk instead of resolved
    from cached survey context.

No stations remain after ``select``
    Re-run with fewer filters and ``--dry-run``. Remember that
    ``select --freq`` keeps stations that touch a band; it does not slice
    frequency rows.

``edit`` refuses to run
    Provide at least one edit operation and, unless using ``--dry-run``, an
    ``--output-dir``.

Export reports existing files
    Choose a new output directory or pass ``--overwrite``.

Recompute exits nonzero
    At least one EDI record failed. Re-run with ``-v`` to print per-record
    status and inspect ``manifest.csv`` when it was written.

Missing tipper or impedance values
    Use ``site info`` first to confirm what arrays exist, then consider
    ``site edit --fill-missing`` or ``site recompute --fill-missing zero`` on
    a copy of the data.

Safety notes
------------

``site select``
    Filters station collections. Its frequency option is a station-membership
    test.

``site edit``
    Changes the data representation and writes edited EDIs to a new directory.

``site recompute``
    Normalizes and rewrites EDI files. Use ``--dry-run`` before a large survey
    rewrite and use ``--overwrite`` only for output directories you intend to
    replace.

Python equivalent
-----------------

.. code-block:: python

   from pycsamt.emtools._core import ensure_sites

   sites = ensure_sites("data/AMT/WILLY_DATA/L18PLT", recursive=True, verbose=0)

   # Inspect station names.
   names = [site.name for site in sites]

   # Use lower-level site helpers for scripted workflows.
   from pycsamt.site import selection

   subset = selection.keep_finite_z(sites)

