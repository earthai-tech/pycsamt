AVG Commands
============

``pycsamt avg`` is the command group for Zonge CSAMT/AMT processed-average
files. It helps you understand what is in an AVG file, check embedded QC
metrics, attach station coordinates from a companion ``.stn`` file, apply
common corrections, and re-export the result in a cleaner format.

AVG files are usually an intermediate survey product. In a practical
workflow you often inspect and clean AVG first, then convert to EDI for
the broader pyCSAMT processing stack.

For AVG-to-EDI conversion, use:

.. code-block:: console

   pycsamt transform avg

Use ``pycsamt avg`` when you are still working at the Zonge AVG layer.
Use ``pycsamt edi`` or ``pycsamt site`` after the data have been written
as EDI files.

What AVG commands are for
-------------------------

.. list-table::
   :header-rows: 1
   :widths: 26 34 40

   * - Command
     - Use it when
     - Main question answered
   * - ``info``
     - You receive a new AVG file.
     - What survey, stations, frequencies, and components are inside?
   * - ``validate``
     - You need to assess data quality before correction or conversion.
     - Which stations or components have high embedded errors?
   * - ``stations``
     - You need station positions or coordinates.
     - Do the AVG and ``.stn`` metadata line up?
   * - ``correct``
     - The line needs static-shift or capacitive-coupling correction.
     - What corrected AVG should be used downstream?
   * - ``export``
     - You need a normalized copy, legacy copy, or NetCDF dataset.
     - Which file should be handed to the next tool?

Supported AVG formats
---------------------

pyCSAMT can work with modern kind-2 AVG files through the full AVG object
model. Legacy kind-1 files are also supported for metadata inspection and
validation through a raw DataFrame loader.

This distinction matters mostly for optional dependencies:

* modern AVG inspection and correction use ``pycsamt.zonge.avg.AVG``;
* legacy/raw inspection and validation can still work without xarray;
* NetCDF export uses xarray and NetCDF support.

If an optional dependency is missing, the CLI prints the package that is
needed instead of silently producing partial output.

Recommended workflow
--------------------

Start by identifying the file:

.. code-block:: console

   pycsamt avg info raw/L18.avg

Then inspect quality metrics:

.. code-block:: console

   pycsamt avg validate raw/L18.avg --threshold 10 --top 20

Check station positions, preferably with the companion station file:

.. code-block:: console

   pycsamt avg stations raw/L18.avg --stn-file raw/L18.stn

Preview any correction before writing:

.. code-block:: console

   pycsamt avg correct raw/L18.avg --method both --filter tma --dry-run

Write the corrected AVG:

.. code-block:: console

   pycsamt avg correct raw/L18.avg \
       --method both \
       --filter tma \
       --output-dir corrected/

Finally, normalize the output format if needed:

.. code-block:: console

   pycsamt avg export corrected/L18_corrected.avg \
       --format modern \
       --output-dir final/

Inspect metadata with ``info``
------------------------------

``info`` gives the first high-level view of an AVG file. It reports the
data kind, project metadata, survey type, line name, operator when
available, station count, frequency count, frequency range, component
names, and row count.

.. code-block:: console

   pycsamt avg info data/avg/K2.AVG
   pycsamt avg info data/avg/K2.AVG --stn-file data/avg/K2.stn
   pycsamt avg info data/avg/K2.AVG --format json
   pycsamt avg info data/avg/K2.AVG --format csv

Use ``--stn-file`` when a companion Zonge station file is available. The
command attempts to attach topographic/station metadata before printing
the summary.

Text output is best for interactive checking. JSON is better when a
script needs to capture survey metadata. CSV gives a compact one-line
summary that can be added to a survey inventory.

What to look for:

``n_stations``
    Confirm that the file contains the number of stations expected for
    the line.

``n_frequencies`` and ``frequency_range``
    Confirm that the band is appropriate for the intended processing or
    inversion.

``components``
    Confirm that the expected impedance/resistivity components are
    present before converting or correcting.

``data_kind``
    Helps distinguish modern and legacy AVG structures.

Validate quality with ``validate``
----------------------------------

``validate`` summarizes quality-control columns already embedded in the
AVG table. It does not invent new QC metrics; it aggregates the error and
phase-spread fields present in the file.

.. code-block:: console

   pycsamt avg validate data/avg/K2.AVG
   pycsamt avg validate data/avg/K2.AVG --comp ZXY
   pycsamt avg validate data/avg/K2.AVG --threshold 8.0 --top 10
   pycsamt avg validate data/avg/K2.AVG --format csv
   pycsamt avg validate data/avg/K2.AVG --format json

QC columns recognized by the command include:

.. list-table::
   :header-rows: 1
   :widths: 24 48 28

   * - Column
     - Meaning
     - Default warning threshold
   * - ``pc_emag``
     - E-field magnitude percent error.
     - 5 percent
   * - ``pc_hmag``
     - H-field magnitude percent error.
     - 5 percent
   * - ``pc_rho``
     - Apparent-resistivity percent error.
     - 10 percent
   * - ``s_ephz``
     - E-field phase standard deviation.
     - 10 mrad
   * - ``s_hphz``
     - H-field phase standard deviation.
     - 10 mrad
   * - ``s_phz``
     - Phase standard deviation.
     - 10 mrad
   * - ``z.%err``
     - Composite impedance percent error.
     - 10 percent

The command groups rows by station and reports mean and maximum values.
A station is flagged when a mean QC value exceeds the relevant threshold.
Use ``--threshold`` to apply a single custom warning threshold when a
project has its own acceptance rule.

Component filtering
~~~~~~~~~~~~~~~~~~~

Use ``--comp`` when one component is responsible for most of the noise:

.. code-block:: console

   pycsamt avg validate data/avg/K2.AVG --comp ZXY
   pycsamt avg validate data/avg/K2.AVG --comp ZYX

Component names must match the values in the AVG ``comp`` column. If no
rows match, the command exits with an error rather than reporting an
empty table.

Worst-station review
~~~~~~~~~~~~~~~~~~~~

For a large line, start with the worst stations:

.. code-block:: console

   pycsamt avg validate data/avg/K2.AVG --top 15

This sorts flagged stations first, then by the largest mean QC value.
The result is a compact triage list for field notes or manual review.

Inspect stations with ``stations``
----------------------------------

``stations`` prints station positions stored in the AVG file. With a
companion ``.stn`` file, it can also include easting, northing,
elevation, latitude, and longitude.

.. code-block:: console

   pycsamt avg stations data/avg/K2.AVG
   pycsamt avg stations data/avg/K2.AVG --stn-file data/avg/K2.stn
   pycsamt avg stations data/avg/K2.AVG --stn-file data/avg/K2.stn --format csv
   pycsamt avg stations data/avg/K2.AVG --sort-by elev

Sorting options are:

``name``
    Sort alphabetically by station name.

``position``
    Sort by survey-line position. This is the default and is usually the
    best view for profile processing.

``elev``
    Sort by elevation when a station file has been attached.

When the ``.stn`` file is attached, pyCSAMT matches station coordinates
to AVG station positions. This is useful before export or conversion
because bad station metadata can cause profile geometry problems later.

Correct data with ``correct``
-----------------------------

``correct`` writes a new AVG file after applying one or both correction
families:

``static-shift``
    Estimates a smooth spatial reference at a selected frequency and
    applies station-wise shift factors.

``capacitive``
    Applies capacitive-coupling correction for high-frequency electric
    field distortion.

``both``
    Applies static-shift correction first, then capacitive-coupling
    correction.

Always preview first:

.. code-block:: console

   pycsamt avg correct data/avg/K2.AVG --dry-run

Then write the corrected file:

.. code-block:: console

   pycsamt avg correct data/avg/K2.AVG \
       --method static-shift \
       --filter tma \
       --output-dir corrected/

The output file is named from the source stem:

.. code-block:: text

   corrected/K2_corrected.avg

Choosing a static-shift filter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``--filter tma``
    Trimmed Moving Average. This works on apparent resistivity and uses
    the ``--window`` value. It is a conservative first choice for many
    survey lines.

``--filter flma``
    Fixed-Length Moving Average. This works on impedance and is useful
    when a fixed spatial window is preferred.

``--filter ama``
    Adaptive Moving Average. This also works on impedance and adapts to
    local behaviour along the profile.

Reference frequency
~~~~~~~~~~~~~~~~~~~

Static-shift correction needs a reference frequency:

.. code-block:: console

   pycsamt avg correct data/avg/K2.AVG --ref-freq 2048 --filter tma

If ``--ref-freq`` is omitted, the command uses the highest available
frequency. If the requested frequency is not present, pyCSAMT adjusts to
the nearest available frequency and reports that adjustment.

Window size
~~~~~~~~~~~

``--window`` controls the TMA window size in number of stations:

.. code-block:: console

   pycsamt avg correct data/avg/K2.AVG --filter tma --window 7

Small windows follow local variation more closely. Larger windows smooth
more aggressively. Choose a value that makes geological sense for the
station spacing and expected lateral structure.

Machine-readable correction summary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use JSON output when a processing script needs the correction summary:

.. code-block:: console

   pycsamt avg correct data/avg/K2.AVG \
       --method both \
       --output-dir corrected/ \
       --format json

For static shift, the JSON summary includes the selected filter,
reference frequency, number of stations, and minimum/maximum/mean shift
factor.

Export with ``export``
----------------------

``export`` rewrites an AVG file to a selected output format. This is
useful after correction, when standardizing file layout, or when handing
data to another tool.

.. code-block:: console

   pycsamt avg export data/avg/K2.AVG --output-dir clean/
   pycsamt avg export data/avg/K2.AVG --format legacy --output-dir legacy/
   pycsamt avg export data/avg/K2.AVG --format xarray --output-dir nc/
   pycsamt avg export data/avg/K2.AVG --stem L18_clean --output-dir clean/
   pycsamt avg export data/avg/K2.AVG --format modern --overwrite

Export formats:

``modern``
    Writes a clean CSV-based kind-2 AVG file. This is the default and is
    usually the best handoff format inside pyCSAMT.

``legacy``
    Writes fixed-width kind-1 AVG. Use this only when a downstream tool
    requires the older layout.

``xarray``
    Writes a NetCDF file through xarray. This is useful for analysis,
    archiving, and tools that prefer labelled multidimensional arrays.
    It requires optional ``xarray`` and ``netcdf4`` dependencies.

By default the output stem is the source stem plus ``_exported``. Use
``--stem`` when you need a stable project naming convention:

.. code-block:: console

   pycsamt avg export corrected/K2_corrected.avg \
       --format modern \
       --stem L18_after_static_shift \
       --output-dir final/

Use ``--format-out json`` to make the CLI response machine-readable. This
is separate from ``--format``, which controls the file format being
written:

.. code-block:: console

   pycsamt avg export corrected/K2_corrected.avg \
       --format modern \
       --format-out json \
       --output-dir final/

End-to-end examples
-------------------

Quick inventory
~~~~~~~~~~~~~~~

.. code-block:: console

   pycsamt avg info raw/L18.avg --format csv
   pycsamt avg stations raw/L18.avg --format csv
   pycsamt avg validate raw/L18.avg --top 10

This is the lightest pass when you only need to know whether an AVG file
is usable and what stations it contains.

QC-first correction
~~~~~~~~~~~~~~~~~~~

.. code-block:: console

   pycsamt avg info raw/L18.avg --stn-file raw/L18.stn
   pycsamt avg validate raw/L18.avg --threshold 10 --top 20
   pycsamt avg correct raw/L18.avg --method static-shift --filter tma --dry-run
   pycsamt avg correct raw/L18.avg --method static-shift --filter tma --output-dir corrected/
   pycsamt avg validate corrected/L18_corrected.avg --threshold 10 --top 20

This pattern validates before and after correction so the correction step
is auditable.

Correction plus normalized export
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: console

   pycsamt avg correct raw/L18.avg \
       --method both \
       --filter ama \
       --ref-freq 2048 \
       --output-dir corrected/

   pycsamt avg export corrected/L18_corrected.avg \
       --format modern \
       --stem L18_clean \
       --output-dir final/

This creates a clean AVG that can be archived or passed to conversion.

Troubleshooting
---------------

``No QC columns found``
    The AVG file does not contain the recognized QC columns
    ``pc_emag``, ``pc_hmag``, ``pc_rho``, ``s_ephz``, ``s_hphz``,
    ``s_phz``, or ``z.%err``. Use ``avg info`` to confirm the file type
    and inspect the source data.

``No rows match component``
    The component name passed to ``--comp`` is not present in the AVG
    ``comp`` column. Run ``avg info`` or ``avg validate`` without
    ``--comp`` to see available components.

Reference frequency was adjusted
    The requested ``--ref-freq`` is not exactly present in the file. The
    command uses the nearest available frequency and reports the value.

NetCDF export fails
    Install the optional dependencies needed for xarray export:

    .. code-block:: console

       pip install xarray netcdf4

Unexpected station coordinates
    Re-run ``avg stations`` with and without ``--stn-file``. If the
    positions differ, check that the station file belongs to the same
    line and acquisition.

