EDI Commands
============

``pycsamt edi`` is the raw EDI command group. It works directly with
``EDIFile`` and ``EDICollection`` objects from ``pycsamt.seg`` before the
data are wrapped into the higher-level ``Sites`` model.

Use this layer when you need to answer file-level questions:

* Are these files structurally valid EDI files?
* Which stations, coordinates, frequencies, and components are present?
* Does the line geometry look right before inversion or plotting?
* Do the tensors need a constant rotation before processing?
* Which files should be copied into a clean subset?

Once the files are valid and selected, move to ``pycsamt site`` for
survey-level operations or use ``ensure_sites`` from Python:

.. code-block:: python

   from pycsamt.emtools import ensure_sites

   sites = ensure_sites("clean_edis", recursive=True, verbose=1)

EDI versus site commands
------------------------

``pycsamt edi`` is intentionally close to the files on disk. It reads
EDI headers and transfer-function blocks, validates file structure, and
writes EDI files when rotating or exporting subsets.

``pycsamt site`` starts one level higher. It treats the survey as a
collection of stations and exposes operations such as filtering by data
quality, editing coordinates, recomputing output EDIs, exporting
manifests, and computing strike or resistivity tables.

In short:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Use
     - When
   * - ``pycsamt edi``
     - You are checking, rotating, or copying raw EDI files.
   * - ``pycsamt site``
     - You are working with a loaded survey as stations.
   * - ``pycsamt pipe``
     - You want a reproducible multi-step processing workflow.

Recommended workflow
--------------------

A careful EDI workflow usually starts with validation, then metadata,
then geometry:

.. code-block:: console

   pycsamt edi validate raw_edis/
   pycsamt edi info raw_edis/
   pycsamt edi stations raw_edis/ --format csv
   pycsamt edi profile raw_edis/ --distances

Then select or rotate only if needed:

.. code-block:: console

   pycsamt edi select raw_edis/ --min-nfreq 20 --dry-run
   pycsamt edi select raw_edis/ --min-nfreq 20 --output-dir clean_edis/

Finally, move to the higher-level survey API:

.. code-block:: console

   pycsamt site info clean_edis/

Command overview
----------------

.. list-table::
   :header-rows: 1
   :widths: 24 34 42

   * - Command
     - Main job
     - Typical output
   * - ``info``
     - Summarize one file or a directory.
     - Station, frequency range, tipper status, coordinates.
   * - ``validate``
     - Check EDI file structure.
     - Pass/fail counts and validation errors.
   * - ``stations``
     - Build a coordinate table.
     - Station, lat, lon, elevation, UTM fields, path.
   * - ``profile``
     - Estimate line geometry.
     - Bearing, station spacing, optional distances.
   * - ``rotate``
     - Rotate impedance and tipper tensors.
     - New EDI files in an output directory.
   * - ``select``
     - Filter and export a subset.
     - Copied EDI files and selection summary.

Validate EDI structure
----------------------

Validation should be the first command you run on a new directory. It
checks whether files are structurally recognizable as EDI files before
you spend time interpreting metadata or profile geometry.

.. code-block:: console

   pycsamt edi validate data/edis/
   pycsamt edi validate data/edis/S01.edi
   pycsamt edi validate data/edis/ --format json

By default validation is deep:

.. code-block:: console

   pycsamt edi validate data/edis/ --deep

Deep validation scans file content for required structural tags such as
``>HEAD``, ``>FREQ``, ``>=MTSECT``, and ``>END``. This is the mode you
want before processing or conversion.

Use extension-only validation only for a quick inventory:

.. code-block:: console

   pycsamt edi validate data/edis/ --no-deep

The command exits with status ``0`` when all files pass and status ``1``
when any file fails. That makes it suitable for scripts and CI checks.

Inspect files with ``info``
---------------------------

``info`` answers "what is inside this file or folder?" It can run on one
EDI file or on a directory.

Directory mode:

.. code-block:: console

   pycsamt edi info data/edis/
   pycsamt edi info data/edis/ --top 20
   pycsamt edi info data/edis/ --station S01
   pycsamt edi info data/edis/ --format csv
   pycsamt edi info data/edis/ --format json

Single-file mode:

.. code-block:: console

   pycsamt edi info data/edis/S01.edi

The reported fields include:

``station``
    Station name read from the EDI file, falling back to the filename
    stem when necessary.

``n_freq``
    Number of frequency samples in the impedance block.

``freq_min`` and ``freq_max``
    Frequency range in Hz.

``has_tipper``
    Whether the file includes tipper transfer functions.

``lat``, ``lon``, ``elev``
    Coordinates from the EDI header when available.

``dtype`` and ``channels``
    File data type and channel metadata where the reader exposes them.

Use ``info`` before ``select`` to learn station names and frequency
coverage. Use JSON output when feeding a station inventory into another
script.

Station coordinate tables
-------------------------

``stations`` focuses on coordinate metadata. It is the command to use
when checking survey geometry, looking for missing coordinates, or
creating a station inventory for GIS or spreadsheets.

.. code-block:: console

   pycsamt edi stations data/edis/
   pycsamt edi stations data/edis/ --pattern "S0*"
   pycsamt edi stations data/edis/ --sort-by lat
   pycsamt edi stations data/edis/ --top 10
   pycsamt edi stations data/edis/ --format csv

The table can include station name, latitude, longitude, elevation,
easting, northing, UTM zone, and path. Available fields depend on what is
stored in the EDI headers and what the coordinate conversion can infer.

Sorting options:

``station``
    Alphabetical station-name order.

``lat`` and ``lon``
    Geographic sorting. Useful for catching stations that are spatially
    out of sequence.

``elev``
    Elevation sorting. Useful for topographic checks.

Pattern filtering uses glob syntax:

.. code-block:: console

   pycsamt edi stations data/edis/ --pattern "18-02*"

This is useful when one directory contains multiple line names or station
families.

Profile geometry
----------------

``profile`` estimates survey-line geometry from the EDI collection. It
is especially useful before 2-D inversion, profile plotting, or any
workflow that assumes stations follow a coherent line.

.. code-block:: console

   pycsamt edi profile data/edis/
   pycsamt edi profile data/edis/ --bearing-method linear
   pycsamt edi profile data/edis/ --step-method median
   pycsamt edi profile data/edis/ --distances
   pycsamt edi profile data/edis/ --distances --format csv
   pycsamt edi profile data/edis/ --format json

Bearing methods:

``endpoints``
    Uses the vector from the first station to the last station. This is
    simple and works well for clean, ordered, nearly straight lines.

``linear``
    Fits the best profile axis from projected station coordinates. This
    is more robust when the line is slightly irregular or when endpoint
    coordinates are noisy.

Step methods:

``mean``
    Average inter-station spacing.

``median``
    Median inter-station spacing. This is often better when one station
    gap is unusually large or small.

Add ``--distances`` when you want a per-station along-profile distance
table. This is often the most useful output for checking station order
before preparing inversion files.

Rotate impedance and tipper tensors
-----------------------------------

``rotate`` applies a constant rotation angle to the impedance tensor and,
when present, the tipper vector. This is a file-writing operation, so use
``--dry-run`` first.

.. code-block:: console

   pycsamt edi rotate data/edis/ --angle 30 --dry-run

Then write the rotated files:

.. code-block:: console

   pycsamt edi rotate data/edis/ --angle 30 --output-dir rotated/
   pycsamt edi rotate data/edis/S01.edi --angle -15 --output-dir rotated/
   pycsamt edi rotate data/edis/ --angle 30 --output-dir rotated/ --format json

The angle is in degrees. Positive values follow the clockwise geophysical
convention used by the command implementation.

Rotation preserves filenames and writes the new files into
``--output-dir``. The command requires ``--output-dir`` unless
``--dry-run`` is used. If all writes fail, the command exits with a
non-zero status.

Use rotation when:

* all stations need a common coordinate-frame correction;
* you have a known survey or geoelectric strike convention to impose;
* downstream tools expect tensors in a consistent orientation.

Avoid rotating casually. Rotation changes the transfer functions that
later QC, static-shift, plotting, and inversion steps will use.

Select and export subsets
-------------------------

``select`` filters a directory of EDI files and exports the matching
subset to a new directory. Multiple filters are combined with logical
AND.

Preview the selection:

.. code-block:: console

   pycsamt edi select data/edis/ --has-tipper --dry-run
   pycsamt edi select data/edis/ --pattern "18-02*" --dry-run

Write the selected files:

.. code-block:: console

   pycsamt edi select data/edis/ \
       --stations S01,S02,S03 \
       --output-dir subset/

   pycsamt edi select data/edis/ \
       --has-tipper \
       --min-nfreq 40 \
       --output-dir subset/

   pycsamt edi select data/edis/ \
       --min-freq 0.1 \
       --max-freq 1000 \
       --min-nfreq 20 \
       --output-dir subset/

Filters:

``--stations ID[,ID...]``
    Keep exact station names. Names must match the station IDs stored in
    the EDI files.

``--pattern GLOB``
    Keep station names matching a glob pattern such as ``S0*`` or
    ``18-02*``.

``--has-tipper``
    Keep only stations with tipper data.

``--min-freq HZ``
    Keep stations whose minimum frequency is less than or equal to the
    requested value. This is useful when requiring low-frequency reach.

``--max-freq HZ``
    Keep stations whose maximum frequency is greater than or equal to
    the requested value. This is useful when requiring high-frequency
    coverage.

``--min-nfreq N``
    Keep stations with at least ``N`` frequency samples.

``--overwrite``
    Allow exported files to replace existing files in the output
    directory.

``--format json``
    Print a machine-readable summary of selected, written, and failed
    files.

Quality-control patterns
------------------------

Find unusable files early:

.. code-block:: console

   pycsamt edi validate raw_edis/ --format json

Inventory stations for a report:

.. code-block:: console

   pycsamt edi info raw_edis/ --format csv
   pycsamt edi stations raw_edis/ --format csv

Check profile assumptions:

.. code-block:: console

   pycsamt edi profile raw_edis/ --bearing-method linear --distances

Create a clean processing folder:

.. code-block:: console

   pycsamt edi select raw_edis/ \
       --min-nfreq 20 \
       --min-freq 0.1 \
       --max-freq 1000 \
       --dry-run

   pycsamt edi select raw_edis/ \
       --min-nfreq 20 \
       --min-freq 0.1 \
       --max-freq 1000 \
       --output-dir clean_edis/

Rotate after confirming geometry:

.. code-block:: console

   pycsamt edi rotate clean_edis/ --angle 30 --dry-run
   pycsamt edi rotate clean_edis/ --angle 30 --output-dir rotated_edis/

Move to Sites
-------------

After raw EDI checks are done, most pyCSAMT workflows should use the
``Sites`` abstraction:

.. code-block:: console

   pycsamt site info clean_edis/
   pycsamt site compute strike clean_edis/
   pycsamt pipe run clean_edis/ --preset basic_qc --out results/

The equivalent Python entry point is:

.. code-block:: python

   from pycsamt.emtools import ensure_sites

   sites = ensure_sites("clean_edis", recursive=True, verbose=1)

``ensure_sites`` is the canonical loader for EDI paths, directories, and
existing ``Sites`` objects. Files without valid impedance data are skipped
because downstream workflows require usable Z data.

Troubleshooting
---------------

``No .edi files found``
    The command searches for files ending in ``.edi`` when a directory is
    provided. Check the path and file extensions.

Validation fails in deep mode but passes with ``--no-deep``
    The file has an EDI extension but is missing required structural
    content. Use the deep-mode error message to identify the missing tag
    or malformed block.

Station names do not match your expected IDs
    ``select --stations`` uses station names read from the EDI content,
    not necessarily filenames. Run ``pycsamt edi info`` first to see the
    station IDs pyCSAMT will use.

Profile bearing looks wrong
    Check station coordinates with ``edi stations``. Then compare
    ``profile --bearing-method endpoints`` and
    ``profile --bearing-method linear``. A large difference often means
    the line is curved, unordered, or has a coordinate outlier.

Rotation writes nothing
    Run the same command with ``--dry-run`` to confirm files were found,
    then provide an explicit ``--output-dir``. Rotation does not write to
    the source directory by default.

Subset export selects no stations
    Remove filters one by one. Start with ``--dry-run`` and only
    ``--pattern`` or only ``--min-nfreq`` so you can see which condition
    is excluding the files.

