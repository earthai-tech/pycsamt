TDEM Commands
=============

``pycsamt tdem`` works with time-domain electromagnetic survey folders,
especially Zonge TEMAVG processed-average data. It can inspect a survey
folder, convert TEMAVG soundings into synthetic frequency-domain EDI files,
and create quick-look figures for decays, sections, maps, and dashboards.

Use this command group when the source is time-domain TEM data rather than
frequency-domain MT/AMT/CSAMT EDI data. After conversion, the written EDI
files can be inspected with ``pycsamt site`` and used by downstream pyCSAMT
workflows.

Command Map
-----------

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Command
     - Purpose
   * - ``pycsamt tdem info SURVEY_DIR``
     - Inspect a TEMAVG survey folder and report AVG, companion Z, LOG, and
       coordinate-table availability.
   * - ``pycsamt tdem convert SURVEY_DIR``
     - Convert TEMAVG soundings to frequency-domain EDI files.
   * - ``pycsamt tdem plot SURVEY_DIR``
     - Plot decay curves, pseudo-sections, station maps, elevation profiles,
       and dashboard figures.

Survey Folder Model
-------------------

``SURVEY_DIR`` must be a directory. The TDEM reader scans that directory for:

``*.AVG``
    Zonge TEMAVG processed-average files. These are the primary data source.

``*.Z``
    Optional Z-contour companion files used by ``plot --kind z-section`` and
    the Z panel in ``plot --kind dashboard``.

``*.LOG``
    Optional processing log files reported by ``tdem info``.

Coordinate tables
    Optional station/profile coordinate files. If ``--coord-file`` is not
    provided, pyCSAMT tries common names such as
    ``Coordinate of measuring point.xls``, ``Coordinate of measuring
    point.xlsx``, ``coordinates.csv``, ``Coordinates.csv``, and
    ``coordinate.csv``.

AVG files are keyed by file stem. For example, ``TEM100.AVG`` has stem
``TEM100``. The ``--stems`` option on ``convert`` and ``plot`` selects those
file stems.

Shared Options
--------------

``--pattern GLOB``
    Glob used to find processed-average files. The default is ``*.AVG``.
    Use ``--pattern "*.avg"`` for lowercase file extensions.

``--component COMP``
    TEMAVG component column used when extracting soundings. The default is
    ``Hz``.

``--stems STEM[,STEM...]``
    Comma-separated AVG file stems to process. When omitted, all discovered
    AVG files are used.

``-v, --verbose``
    Print more progress and companion-file detail.

``--no-color``
    Disable colored terminal output.

Inspect A TEMAVG Survey
-----------------------

``pycsamt tdem info`` parses the survey folder and reports available files and
coordinate metadata.

.. code-block:: console

   pycsamt tdem info data/TEMAVG/JIANGSU
   pycsamt tdem info data/TEMAVG/JIANGSU --format json
   pycsamt tdem info data/TEMAVG/JIANGSU --pattern "*.avg" -v

Options:

``--pattern GLOB``
    Select processed-average files. Defaults to ``*.AVG``.

``--format text|json|csv``
    Uses the shared format option. The current command implements text and
    JSON output directly; non-JSON formats fall back to text-style output.

``-v``
    In text mode, also lists companion ``.Z`` and ``.LOG`` stems.

JSON output includes:

``root``
    Survey directory path.

``n_avg_files``, ``n_z_files``, ``n_log_files``
    Counts of parsed primary and companion files.

``avg_stems``, ``z_stems``, ``log_stems``
    Sorted file stems.

``has_coordinates``
    Whether a coordinate table was loaded.

``coordinates``
    Coordinate-table summary with ``n_points`` and ``has_elevation`` when a
    table exists.

If no matching AVG files are found, the command exits with a usage error and
suggests checking the directory and pattern.

Convert TEMAVG To EDI
---------------------

``pycsamt tdem convert`` reads a TEMAVG survey folder, builds one
``TEMSounding`` per station, transforms each sounding to the frequency domain,
and writes synthetic impedance EDI files.

.. code-block:: console

   pycsamt tdem convert data/TEMAVG/JIANGSU --output-dir imp_edis/
   pycsamt tdem convert data/TEMAVG/JIANGSU --stems TEM100,TEM200 --output-dir imp_edis/
   pycsamt tdem convert data/TEMAVG/JIANGSU --method fourier --output-dir imp_edis/
   pycsamt tdem convert data/TEMAVG/JIANGSU --coord-file coordinates.csv --output-dir imp_edis/
   pycsamt tdem convert data/TEMAVG/JIANGSU --dry-run
   pycsamt tdem convert data/TEMAVG/JIANGSU --output-dir imp_edis/ --format json

Required output behavior:

* real conversion requires ``--output-dir``;
* ``--dry-run`` does not require ``--output-dir`` and writes nothing;
* the output directory is created when conversion runs.

Conversion options:

``--method late_time|fourier``
    Select the time-to-frequency transform. ``late_time`` is the default and
    is fast. ``fourier`` uses a Fourier transform path.

``--component COMP``
    TEMAVG component used for the sounding decay. The default is ``Hz``.

``--pattern GLOB``
    Select AVG files inside the survey directory.

``--stems STEM[,STEM...]``
    Convert only selected AVG file stems.

``--coord-file FILE``
    Coordinate table to attach to each sounding. When omitted, common
    coordinate filenames are discovered automatically.

``--freq-convention skin_depth|fourier_time``
    Pseudo-frequency convention used by the late-time transform. The default
    is ``skin_depth``.

``--phase-mode homogeneous|weidelt``
    Phase model used when synthesizing impedance from apparent resistivity.
    The default is ``homogeneous``.

``--no-geometry-correction``
    Disable the in-loop Biot-Savart geometry correction. By default, geometry
    correction is enabled.

``--format text|json``
    Select terminal summary format.

``--dry-run``
    Count and list soundings without writing EDI files.

Text summary:

.. code-block:: text

   Soundings : 12  | Written: 12  | Output : imp_edis/

JSON summary includes ``survey_dir``, ``output_dir``, ``method``,
``n_soundings``, ``n_written``, and ``written``.

Exit behavior:

* conversion exceptions produce a nonzero exit status;
* if soundings transform but no EDI files are written, the command prints a
  warning and exits nonzero;
* a successful conversion exits zero and reports written file paths when
  ``-v`` is used.

Plot TDEM Data
--------------

``pycsamt tdem plot`` visualizes data from a TEMAVG survey folder.

.. code-block:: console

   pycsamt tdem plot data/TEMAVG/JIANGSU --kind decay
   pycsamt tdem plot data/TEMAVG/JIANGSU --kind rho --output-dir figs/ --dpi 200
   pycsamt tdem plot data/TEMAVG/JIANGSU --kind dashboard --stems TEM100,TEM200 --output-dir figs/
   pycsamt tdem plot data/TEMAVG/JIANGSU --kind map --output-dir figs/

Output behavior:

* when ``--output-dir`` is omitted, the figure is shown interactively;
* when ``--output-dir`` is provided, Matplotlib uses the non-interactive
  ``Agg`` backend and saves a file;
* saved filenames use ``<survey_dir_name>_<kind>.<fmt>``;
* the output directory is created automatically.

Plot options:

``--kind, -k KIND``
    Plot kind. Choices are ``decay``, ``rho``, ``section``, ``z-section``,
    ``map``, ``overview``, ``gate-profile``, ``elevation``, and
    ``dashboard``.

``--stems STEM[,STEM...]``
    Include only selected AVG file stems. Missing stems are reported as a
    warning.

``--pattern GLOB``
    Select AVG files before plotting.

``--component COMP``
    Component used for sounding extraction in plot kinds that need soundings.

``--output-dir, -o DIR``
    Save the figure instead of showing it interactively.

``--fmt png|pdf|svg|jpg``
    Saved figure format. The default is ``png``.

``--dpi INT``
    Saved figure resolution. The default is ``150``.

Plot kinds:

``decay``
    Time-domain decay curves for each extracted sounding. Requires AVG data.

``rho``
    Transformed apparent-resistivity and phase curves. Requires AVG data and
    extracted soundings.

``section``
    TEMAVG gate-value pseudo-section using the first selected AVG file.

``z-section``
    Z-contour pseudo-section using companion ``.Z`` files. Exits nonzero if no
    matching ``.Z`` files are available.

``map``
    Station location map. Requires usable coordinate data in the survey.

``overview``
    Multi-panel survey overview with map and elevation profile.

``gate-profile``
    Gate values versus station position for representative windows.

``elevation``
    Station elevation profile from extracted sounding coordinates.

``dashboard``
    Four-panel dashboard combining selected decays, transformed resistivity,
    AVG pseudo-section, and Z-section when companion data are available.

If a plot kind cannot produce data, the command exits nonzero with a short
diagnostic such as ``No soundings extracted`` or ``No companion .Z files
found``.

Recommended Workflows
---------------------

Inspect first, then convert:

.. code-block:: console

   pycsamt tdem info data/TEMAVG/JIANGSU
   pycsamt tdem convert data/TEMAVG/JIANGSU --dry-run
   pycsamt tdem convert data/TEMAVG/JIANGSU --output-dir outputs/tem_edis -v
   pycsamt site info outputs/tem_edis

Convert only selected profiles:

.. code-block:: console

   pycsamt tdem info data/TEMAVG/JIANGSU --format json
   pycsamt tdem convert data/TEMAVG/JIANGSU \
       --stems TEM100,TEM200 \
       --method fourier \
       --output-dir outputs/tem_selected

Create a quick figure pack:

.. code-block:: console

   pycsamt tdem plot data/TEMAVG/JIANGSU --kind decay --output-dir figs/
   pycsamt tdem plot data/TEMAVG/JIANGSU --kind section --output-dir figs/
   pycsamt tdem plot data/TEMAVG/JIANGSU --kind map --output-dir figs/
   pycsamt tdem plot data/TEMAVG/JIANGSU --kind dashboard --output-dir figs/

Troubleshooting
---------------

No TEMAVG files found
    Check that the directory contains Zonge TEMAVG files and that the case of
    ``--pattern`` matches the file extension. The default is ``*.AVG``.

``--output-dir`` is required
    Real conversion requires an explicit output directory. Use ``--dry-run``
    when you only want to inspect soundings.

No soundings extracted
    Check ``--component``, ``--stems``, and the selected AVG files. A stem must
    match the AVG filename without its extension.

Missing coordinates
    Pass ``--coord-file`` to ``convert`` if automatic coordinate discovery does
    not find the survey table. Map and elevation plots need coordinate data.

``z-section`` fails
    That plot kind requires companion ``.Z`` contour files with matching
    stems.

Dashboard fails
    Dashboard needs enough AVG soundings and may need companion ``.Z`` data
    for the Z-section panel. Try ``decay`` or ``section`` first to isolate the
    missing input.

Conversion writes zero files
    Re-run with ``-v`` and inspect the dry-run sounding list. The command exits
    nonzero when soundings exist but no EDI file is written.

Safety Notes
------------

``tdem info`` and ``tdem plot`` read survey files and do not modify them.

``tdem convert`` writes new EDI files under ``--output-dir``. It does not
modify source ``.AVG``, ``.Z``, ``.LOG``, or coordinate files.

Use a fresh output directory when comparing ``late_time`` and ``fourier``
conversion results so files from separate runs are not mixed.

Python Equivalent
-----------------

Read a survey and inspect stems:

.. code-block:: python

   from pycsamt.tdem.survey import read_temavg_survey

   survey = read_temavg_survey("data/TEMAVG/JIANGSU", pattern="*.AVG")
   print(survey.stems)

Build soundings:

.. code-block:: python

   from pycsamt.tdem.workflow import read_temavg_soundings

   soundings = read_temavg_soundings(
       "data/TEMAVG/JIANGSU",
       stems=["TEM100"],
       component="Hz",
   )

Convert to EDI files:

.. code-block:: python

   from pycsamt.tdem.workflow import transform_temavg_survey

   result = transform_temavg_survey(
       "data/TEMAVG/JIANGSU",
       method="late_time",
       savepath="outputs/tem_edis",
   )

   print(result.n_soundings)
   print(result.written_paths)

Plot from Python:

.. code-block:: python

   from pycsamt.tdem import plot as tdem_plot

   ax = tdem_plot.plot_decay(soundings)
   fig = ax.figure
   fig.savefig("figs/decay.png", dpi=150, bbox_inches="tight")

