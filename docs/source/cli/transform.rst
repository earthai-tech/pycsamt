Transform Commands
==================

``pycsamt transform`` converts raw or intermediate MT data formats into
impedance EDI outputs. It is the conversion-oriented command group for turning
Spectra EDI, Zonge AVG, and Jones J-file sources into EDI files that can be
loaded by ``pycsamt site``, processed by ``pycsamt pipe``, or used as
inversion inputs.

Use this group when the goal is "make impedance EDIs from another transfer
function format." Use the dedicated ``pycsamt avg`` and ``pycsamt jones``
groups when the goal is inspection, validation, station listing, or format-
specific diagnostics before conversion.

Command Map
-----------

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Command
     - Purpose
   * - ``pycsamt transform spectra FILE_OR_DIR``
     - Convert Spectra-format EDI files to impedance EDI files, including
       tipper when a vertical magnetic channel is available.
   * - ``pycsamt transform avg FILE_OR_DIR``
     - Convert Zonge ``.avg`` files to impedance EDI files.
   * - ``pycsamt transform j FILE_OR_DIR``
     - Convert Jones ``.j`` files to impedance EDI files.

Mental Model
------------

The transform commands are thin CLI wrappers around
``pycsamt.transformers``:

``SpectraToEDI``
    Reads Spectra-format EDI files, estimates the impedance tensor from
    cross-spectral matrices, estimates tipper when possible, and writes
    impedance EDIs.

``AVGtoEDI``
    Reads Zonge AVG files, extracts transfer functions or resistivity/phase
    data, finalizes station metadata, and emits EDI files.

``JtoEDI``
    Reads Jones J-files, extracts impedance, tipper, or resistivity/phase
    payloads, finalizes metadata, and emits EDI files.

Every command accepts either a single file or a directory. Directory mode uses
only files with the expected suffix in the top-level directory:

* ``spectra`` uses ``*.edi``;
* ``avg`` uses ``*.avg``;
* ``j`` uses ``*.j``.

The commands process files in sorted order. They create the output directory
for real conversions, but require ``--output-dir`` unless ``--dry-run`` is
used.

Shared Options
--------------

``-o, --output-dir DIR``
    Destination directory for converted EDI files. Required for real writes.
    The default Click value is ``.``, but these commands treat that as
    "missing" and raise a usage error unless ``--dry-run`` is used.

``--station-name NAME``
    Override the station name for single-file inputs. For multi-file directory
    inputs, the converter keeps source-derived station names.

``--format text|json``
    Select terminal summary format. Transform commands support ``text`` and
    ``json`` summaries, not CSV.

``--dry-run``
    Validate source discovery and print what would be converted without
    writing files.

``-v, --verbose``
    Print more progress information. For spectra conversion, verbose output
    includes per-station frequency counts and tipper status.

``--no-color``
    Disable colored terminal output.

Spectra EDI Conversion
----------------------

``pycsamt transform spectra`` converts Spectra-format EDI files into
impedance EDIs.

.. code-block:: console

   pycsamt transform spectra HBH03.edi --output-dir imp_edis/
   pycsamt transform spectra spectra_dir/ --output-dir imp_edis/
   pycsamt transform spectra spectra_dir/ --station-suffix _IMP --output-dir imp_edis/
   pycsamt transform spectra spectra_dir/ --estimate-error --output-dir imp_edis/
   pycsamt transform spectra spectra_dir/ --dry-run
   pycsamt transform spectra spectra_dir/ --output-dir imp_edis/ --format json

The spectra conversion estimates, per frequency:

.. code-block:: text

   Z = S_EH * inv(S_HH)
   T = S_ZH * inv(S_HH)

where ``S_EH`` is the electric-magnetic cross-spectral block, ``S_HH`` is the
horizontal magnetic auto-spectral block, and ``S_ZH`` is used for tipper when a
vertical magnetic channel is present.

Spectra-specific options:

``--station-suffix STR``
    Append a suffix to every converted station name. A common convention is
    ``--station-suffix _IMP`` to distinguish impedance EDIs from source
    spectra EDIs.

``--station-name NAME``
    Override the station name for a single source file. In directory mode this
    is passed to the transformer but only applies when the resolved source is
    a single file.

``--estimate-error``
    Propagate first-order impedance errors into variance blocks such as
    ``>ZXX.VAR`` and ``>ZXY.VAR`` when enough spectra metadata are available.

``--ridge FLOAT``
    Apply Tikhonov regularization to the magnetic block inversion. This is
    useful when spectra are numerically ill-conditioned.

``--dof FLOAT``
    Override effective degrees of freedom used during error estimation. If
    omitted, the transformer tries to infer DoF from spectra metadata.

``--use-remote``
    Select the remote reference when both local and remote electric channels
    are present.

``--e-labels EX,EY``
    Electric channel labels used to identify the electric field channels.
    The default is ``EX,EY``.

``--h-labels HX,HY``
    Horizontal magnetic channel labels. The default is ``HX,HY``.

Spectra output summary:

``--format text``
    Prints ``Converted: N/M | Errors: K | Output: DIR``. With ``-v``, prints
    each converted station, number of frequencies, and whether tipper was
    present.

``--format json``
    Emits ``n_input``, ``n_ok``, ``n_fail``, ``output_dir``, ``converted``,
    and ``failures``. Each converted entry contains ``station``, ``n_freq``,
    and ``has_tipper``.

Failure behavior:

* the transformer runs with ``skip_errors=True``;
* per-file failures are collected and reported;
* the command exits nonzero only when all files fail;
* partial success exits zero while still listing failures.

Zonge AVG Conversion
--------------------

``pycsamt transform avg`` converts Zonge ``.avg`` files to impedance EDIs.

.. code-block:: console

   pycsamt transform avg survey.avg --output-dir imp_edis/
   pycsamt transform avg avg_dir/ --output-dir imp_edis/
   pycsamt transform avg survey.avg --station-name S001 --output-dir imp_edis/
   pycsamt transform avg survey.avg --dry-run
   pycsamt transform avg avg_dir/ --output-dir imp_edis/ --format json

Input behavior:

* a file input is converted directly;
* a directory input processes sorted ``*.avg`` files only;
* if no ``.avg`` files are found in a directory, the command exits nonzero;
* dry-run prints the files that would be converted and writes nothing.

Conversion behavior:

``AVGtoEDI`` extracts transfer functions from the AVG source. When impedance
is not directly present, it can fall back to resistivity and phase data through
the transformer finalization layer. If the AVG object carries topography or
coordinate metadata, the emitted EDI headers are enriched when possible.

Output summary:

``--format text``
    Prints ``Converted: N/M | Errors: K | Output: DIR`` and prints failures to
    stderr when present.

``--format json``
    Emits ``n_input``, ``n_ok``, ``n_fail``, ``converted``, and ``failures``.
    Converted entries include ``station`` and ``source``.

Failure behavior:

* conversion continues across files;
* write failures are recorded per source file;
* the command exits nonzero only if failures occurred and no EDI was written.

Jones J-File Conversion
-----------------------

``pycsamt transform j`` converts Jones ``.j`` files to impedance EDIs.

.. code-block:: console

   pycsamt transform j station.j --output-dir imp_edis/
   pycsamt transform j jfiles_dir/ --output-dir imp_edis/
   pycsamt transform j station.j --station-name J001 --output-dir imp_edis/
   pycsamt transform j station.j --dry-run
   pycsamt transform j jfiles_dir/ --output-dir imp_edis/ --format json

Input behavior:

* a file input is converted directly;
* a directory input processes sorted ``*.j`` files only;
* if no ``.j`` files are found in a directory, the command exits nonzero;
* dry-run prints the files that would be converted and writes nothing.

Conversion behavior:

``JtoEDI`` reads a Jones ``JFile`` and probes common transfer-function aliases:
impedance, tipper, and resistivity/phase. It finalizes station names and
frequency ordering, emits EDI objects, and copies available header coordinates
when present.

Output summary and failure behavior match ``transform avg``:

* text prints conversion counts;
* JSON emits ``n_input``, ``n_ok``, ``n_fail``, ``converted``, and
  ``failures``;
* converted entries include ``station`` and ``source``;
* the command exits nonzero only if every attempted conversion fails.

Dry Runs
--------

Use ``--dry-run`` before a large conversion to confirm source discovery:

.. code-block:: console

   pycsamt transform spectra spectra_dir/ --dry-run
   pycsamt transform avg avg_dir/ --dry-run
   pycsamt transform j jfiles_dir/ --dry-run

Dry-run mode:

* checks that the source path exists;
* expands the directory to the expected file suffix;
* prints the files that would be converted;
* does not require ``--output-dir``;
* writes no EDI files.

Output Files
------------

The concrete filename is chosen by the underlying EDI writer and station
metadata. In normal workflows, pass a clean output directory and inspect the
written ``*.edi`` files afterward:

.. code-block:: console

   pycsamt transform spectra data/AMT/SPECTRA --station-suffix _IMP --output-dir data/AMT/IMP_EDI
   pycsamt site info data/AMT/IMP_EDI

For conversion runs that feed downstream processing, a typical output layout
is:

.. code-block:: text

   converted_edis/
   |-- HBH01_IMP.edi
   |-- HBH02_IMP.edi
   `-- HBH03_IMP.edi

After conversion, validate the output as EDI-backed stations:

.. code-block:: console

   pycsamt site info converted_edis/
   pycsamt site select converted_edis/ --keep-finite --drop-empty --dry-run

Recommended Workflows
---------------------

Convert Spectra EDIs and run a first QC pipeline:

.. code-block:: console

   pycsamt transform spectra data/AMT/SPECTRA \
       --station-suffix _IMP \
       --estimate-error \
       --output-dir outputs/imp_edis \
       -v
   pycsamt site info outputs/imp_edis
   pycsamt survey set outputs/imp_edis
   pycsamt pipe run --preset basic_qc --out outputs/basic_qc

Convert AVG files and inspect the result:

.. code-block:: console

   pycsamt transform avg data/avg --output-dir outputs/avg_edis
   pycsamt site info outputs/avg_edis --format json

Convert Jones files and continue with site tools:

.. code-block:: console

   pycsamt transform j data/jfiles --output-dir outputs/j_edis
   pycsamt site select outputs/j_edis --keep-finite --dry-run

Troubleshooting
---------------

``--output-dir`` is required
    Real conversion writes require an explicit output directory. Add
    ``--output-dir DIR`` or use ``--dry-run``.

No files found
    Directory mode only scans the top-level directory for the expected suffix:
    ``*.edi`` for spectra, ``*.avg`` for AVG, and ``*.j`` for Jones files.

Spectra conversion has no tipper
    Tipper is emitted only when usable vertical magnetic-channel information is
    available in the spectra source.

Spectra conversion fails on ill-conditioned data
    Try ``--ridge`` with a small positive value, or first inspect the spectra
    source with lower-level tools.

Error estimation fails or gives empty variance blocks
    Use ``--dof`` when effective degrees of freedom cannot be inferred from
    spectra metadata, or omit ``--estimate-error`` for a basic conversion.

Partial conversion reports failures but exits zero
    This is expected when at least one file converts successfully. Inspect the
    ``failures`` array in JSON output or the stderr failure lines in text mode.

Station names are not what you expected
    Use ``--station-name`` for single-file conversions or ``--station-suffix``
    for spectra batches. For post-conversion renaming, use site editing/export
    workflows.

Safety Notes
------------

Transform commands write new EDI files under ``--output-dir``. They do not
modify the source Spectra, AVG, or Jones files.

Use separate output directories for separate conversion attempts. The spectra
CLI exposes ``--overwrite`` because it shares the common output option set,
but the current converter writes through the EDI backend; prefer a clean
directory when comparing experiments.

Python Equivalent
-----------------

Spectra conversion:

.. code-block:: python

   from pycsamt.transformers import SpectraToEDI

   result = SpectraToEDI(
       estimate_error=True,
       station_suffix="_IMP",
       verbose=1,
   ).transform_batch(
       "data/AMT/SPECTRA",
       output_dir="outputs/imp_edis",
   )

   print(result.n_ok, result.n_fail)

AVG conversion:

.. code-block:: python

   from pycsamt.transformers import AVGtoEDI

   collection = AVGtoEDI().transform("data/avg/survey.avg")
   for edi in collection:
       edi.write(savepath="outputs/avg_edis")

Jones conversion:

.. code-block:: python

   from pycsamt.transformers import JtoEDI

   edi = JtoEDI().transform("data/jfiles/station.j", name="J001")
   edi.write(savepath="outputs/j_edis")

