Inversion Commands
==================

``pycsamt invert`` is the command group for preparing, launching,
inspecting, summarising, and plotting inversion runs. It is designed
around a work directory: the build command writes solver-specific input
files into that directory, the run command launches the external solver
from that directory, and the reporting/plotting commands read the files
that the solver produced.

The command group currently supports two solver families:

* ``occam2d`` for 2-D Occam input, execution, result loading, and the
  Occam-focused diagnostic plots.
* ``modem`` for ModEM 2-D or 3-D input, execution, result loading, and
  ModEM-focused model/response plots.

The CLI does not hide the fact that these are different inversion
programs. It gives them a shared workflow shape, but the files, model
geometry, and some plots remain solver-specific.

Command Map
-----------

.. list-table::
   :header-rows: 1
   :widths: 24 34 42

   * - Command
     - Purpose
     - Main output
   * - ``pycsamt invert build``
     - Convert a survey or EDI directory into solver input files.
     - A populated inversion work directory.
   * - ``pycsamt invert run``
     - Launch the external Occam2D or ModEM executable.
     - Iteration/model/response files and solver logs.
   * - ``pycsamt invert status``
     - Inspect a work directory without loading a full result object.
     - Readiness flags, file counts, iteration counts, and RMS hints.
   * - ``pycsamt invert results``
     - Load a completed inversion result and print a compact summary.
     - Loaded iteration, final RMS, model shape, and resistivity stats.
   * - ``pycsamt invert plot model``
     - Plot a reconstructed model section.
     - Figure file or interactive matplotlib window.
   * - ``pycsamt invert plot misfit``
     - Plot convergence through iterations.
     - RMS curve, and for Occam2D optional roughness/Lagrange context.
   * - ``pycsamt invert plot response``
     - Compare observed and predicted responses.
     - Response panels for one or more stations.
   * - ``pycsamt invert plot pseudo``
     - Plot apparent-resistivity/phase pseudosections.
     - Pseudosection figure.
   * - ``pycsamt invert plot section``
     - Plot a ModEM depth/section view.
     - ModEM-only section figure.
   * - ``pycsamt invert plot 1d``
     - Extract station-centered 1-D profiles from an Occam2D model.
     - Occam2D-only profile figure.
   * - ``pycsamt invert plot per-site``
     - Plot per-site RMS and residual diagnostics.
     - Occam2D-only diagnostic figure.
   * - ``pycsamt invert plot grid``
     - Plot a compact response grid across stations.
     - Occam2D-only response-grid figure.

Mental Model
------------

An inversion workflow has four distinct stages.

1. Load a survey. ``invert build`` resolves the input as an explicit
   ``EDI_DIR``, an explicit ``--survey`` path, or the active survey set
   with ``pycsamt survey set``.

2. Write solver files. The builder creates the work directory if it does
   not exist, applies frequency/error/mesh options, and delegates to the
   relevant input builder:
   ``pycsamt.models.occam2d.builder.InputBuilder`` or
   ``pycsamt.models.modem.builder.InputBuilder``.

3. Run the solver. ``invert run`` does not rebuild the input files. It
   uses the existing work directory, detects or accepts the solver name,
   and calls the corresponding runner:
   ``pycsamt.models.occam2d.runner.OccamRunner`` or
   ``pycsamt.models.modem.runner.ModEmRunner``.

4. Read and plot outputs. ``status`` checks file signatures and simple
   counters. ``results`` and ``plot`` instantiate the solver-specific
   ``InversionResult`` class and therefore require completed result
   files.

This separation is deliberate. It lets you build once, adjust solver
files manually when needed, run on a workstation or cluster, and later
use the pyCSAMT CLI only for status, summaries, and figures.

Basic Workflow
--------------

Build an Occam2D run from an explicit EDI directory:

.. code-block:: console

   pycsamt invert build data/line01_edi --solver occam2d --workdir runs/line01_occam

Launch the solver:

.. code-block:: console

   pycsamt invert run runs/line01_occam --max-iter 100 --target-misfit 1.05

Inspect the directory:

.. code-block:: console

   pycsamt invert status runs/line01_occam

Summarise the loaded result:

.. code-block:: console

   pycsamt invert results runs/line01_occam

Save the standard figures:

.. code-block:: console

   pycsamt invert plot model    runs/line01_occam --save figures/line01_model.png
   pycsamt invert plot misfit   runs/line01_occam --save figures/line01_misfit.png
   pycsamt invert plot response runs/line01_occam --save figures/line01_response.png
   pycsamt invert plot pseudo   runs/line01_occam --save figures/line01_pseudo.png

For ModEM, use ``--solver modem`` at build time:

.. code-block:: console

   pycsamt invert build data/line01_edi --solver modem --modem-mode 3d \
       --initial-rho 100 --workdir runs/line01_modem
   pycsamt invert run runs/line01_modem --solver modem --async
   pycsamt invert status runs/line01_modem --solver modem

Input Resolution
----------------

``invert build`` accepts an optional positional ``EDI_DIR``:

.. code-block:: console

   pycsamt invert build EDI_DIR

If ``EDI_DIR`` is omitted, the command falls back to ``--survey`` and
then to the active survey context:

.. code-block:: console

   pycsamt survey set data/line01_edi
   pycsamt invert build --solver occam2d --workdir runs/line01_occam

Resolution priority is:

1. Positional ``EDI_DIR``.
2. ``--survey PATH``.
3. The active survey saved by ``pycsamt survey set``.

Use ``--fresh`` when you want the survey to be re-parsed from disk rather
than reused from cached context:

.. code-block:: console

   pycsamt invert build --fresh --workdir runs/line01_occam

If no explicit path and no active survey are available, the command fails
with a survey-context error instead of guessing.

Solver Detection
----------------

Most commands accept ``--solver occam2d`` or ``--solver modem``. For
commands that operate on an existing work directory, the option is
optional when the directory contains recognizable files.

Occam2D detection uses signatures such as:

* ``Occam2DMesh``
* ``Occam2DStartup``
* ``OccamData.dat``
* any ``*.iter`` or ``*.resp`` files

ModEM detection uses signatures such as:

* ``ModEMData.dat``
* ``ModEM.inv``
* ``Modular_NLCG.log``

If a directory is empty, partially copied, or uses non-standard names,
pass the solver explicitly:

.. code-block:: console

   pycsamt invert status runs/manual_occam --solver occam2d
   pycsamt invert results runs/manual_modem --solver modem

Build
-----

Usage:

.. code-block:: console

   pycsamt invert build [EDI_DIR] [OPTIONS]

The build command loads the survey, creates the work directory, prepares
solver configuration, and writes solver input files. The default solver
is ``occam2d`` and the default work directory is ``invert_run``.

General options:

.. list-table::
   :header-rows: 1
   :widths: 28 22 50

   * - Option
     - Default
     - Meaning
   * - ``--solver {occam2d,modem}``
     - ``occam2d``
     - Select the inversion input format and runner family.
   * - ``--workdir``, ``-w``
     - ``invert_run``
     - Directory where input files are written.
   * - ``--survey PATH``
     - none
     - Explicit survey path used when positional ``EDI_DIR`` is omitted.
   * - ``--fresh``
     - false
     - Rebuild survey state from disk before writing inputs.
   * - ``--freq FMIN:FMAX``
     - none
     - Keep only frequencies in this Hz range.
   * - ``--n-layers INT``
     - solver default
     - Occam2D earth layers or ModEM vertical cells.
   * - ``--cell-size METRES``
     - solver default
     - Target horizontal cell size near stations.
   * - ``-v``, ``--verbose``
     - quiet
     - Increase CLI logging.
   * - ``--no-color``
     - false
     - Disable colored CLI logging/output.

Occam2D options:

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Option
     - Default
     - Meaning
   * - ``--modes TE,TM[,DET]``
     - ``TE,TM``
     - Comma-separated Occam data modes. Values are upper-cased by the
       CLI before being passed to ``OccamConfig``.
   * - ``--error-floor-rho FLOAT``
     - config default
     - Relative apparent-resistivity error floor, for example ``0.05``.
   * - ``--error-floor-phase FLOAT``
     - config default
     - Absolute phase error floor in degrees.

ModEM options:

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Option
     - Default
     - Meaning
   * - ``--modem-mode {2d,3d}``
     - ``3d``
     - Select ModEM dimensionality.
   * - ``--initial-rho OHM_M``
     - config default
     - Initial half-space resistivity in ohm-m.

Examples:

.. code-block:: console

   pycsamt invert build survey/ --workdir run01
   pycsamt invert build survey/ --freq 0.1:1000 --n-layers 60 \
       --error-floor-rho 0.05 --workdir run01
   pycsamt invert build survey/ --solver modem --modem-mode 3d \
       --initial-rho 100 --workdir modem_run

Build is the only command in this group that resolves EDI/survey input.
All later commands expect a work directory that already exists.

Run
---

Usage:

.. code-block:: console

   pycsamt invert run WORKDIR [OPTIONS]

``run`` launches the external executable associated with the selected
solver. It uses ``OccamRunner`` for Occam2D and ``ModEmRunner`` for
ModEM. The solver can be detected from the work directory or supplied
explicitly.

Options:

.. list-table::
   :header-rows: 1
   :widths: 30 22 48

   * - Option
     - Default
     - Meaning
   * - ``--solver {occam2d,modem}``
     - auto
     - Solver family to launch.
   * - ``--max-iter INT``
     - file/config default
     - Override the maximum iteration count in the startup/control file.
   * - ``--target-misfit FLOAT``
     - file/config default
     - Override the target RMS misfit for this run.
   * - ``--async``
     - false
     - Start the solver and return immediately.
   * - ``-v``, ``--verbose``
     - quiet
     - Increase logging.
   * - ``--no-color``
     - false
     - Disable colored output.

Synchronous examples:

.. code-block:: console

   pycsamt invert run run01
   pycsamt invert run run01 --max-iter 100 --target-misfit 1.05

Asynchronous ModEM example:

.. code-block:: console

   pycsamt invert run modem_run --solver modem --async

When ``--async`` is used, the command prints the process ID and the
stdout log path:

* Occam2D: ``occam_stdout.log``
* ModEM: ``modem_stdout.log``

When a synchronous run exits with a non-zero code, the CLI reports the
solver-specific stderr log:

* Occam2D: ``occam_stderr.log``
* ModEM: ``modem_stderr.log``

The command assumes that the external solver binary can be discovered by
the configured runner. For ModEM this commonly means the configured
``Mod2DMT`` or ``Mod3DMT`` executable is on ``PATH`` or provided through
the model configuration. For Occam2D, the runner handles Occam binary
discovery according to the Occam2D runner implementation.

Status
------

Usage:

.. code-block:: console

   pycsamt invert status WORKDIR [OPTIONS]

``status`` is a lightweight directory inspection command. It does not
require a complete inversion result. It lists files, checks expected
solver signatures, and reports whether the directory appears ready to
run.

Options:

.. list-table::
   :header-rows: 1
   :widths: 30 22 48

   * - Option
     - Default
     - Meaning
   * - ``--solver {occam2d,modem}``
     - auto
     - Solver family. Required when detection fails.
   * - ``--format {text,json,csv}``
     - text
     - Human-readable table or JSON payload. The shared option accepts
       ``csv``, but this command currently treats non-JSON output as
       text.
   * - ``-v``, ``--verbose``
     - quiet
     - Increase logging.
   * - ``--no-color``
     - false
     - Disable colored table/log output.

Occam2D status reports:

* total file count and file names in JSON mode;
* whether data, mesh, model, and startup files are present;
* whether all required run inputs appear present;
* number of ``*.iter`` files;
* last ``*.iter`` file;
* number of ``*.resp`` files;
* last RMS parsed from an Occam log when available.

ModEM status reports:

* total file count and file names in JSON mode;
* whether data, control, and covariance files are present;
* whether data and control files appear sufficient to run;
* number of ``*.rho`` model files;
* last ``*.rho`` model file.

Examples:

.. code-block:: console

   pycsamt invert status run01
   pycsamt invert status run01 --format json
   pycsamt invert status modem_run --solver modem

JSON output is useful for automation:

.. code-block:: console

   pycsamt invert status run01 --format json

The JSON object includes keys such as ``workdir``, ``solver``,
``files_present``, ``n_files``, ``ready_to_run``, and solver-specific
iteration/model counters.

Results
-------

Usage:

.. code-block:: console

   pycsamt invert results WORKDIR [OPTIONS]

``results`` loads the solver-specific ``InversionResult`` object and
prints a compact summary. Unlike ``status``, it expects completed output
files that can be parsed by the result loader.

Options:

.. list-table::
   :header-rows: 1
   :widths: 30 22 48

   * - Option
     - Default
     - Meaning
   * - ``--solver {occam2d,modem}``
     - auto
     - Solver family. Required when detection fails.
   * - ``--iteration INT``
     - latest available
     - Load a specific iteration when the result reader supports it.
   * - ``--format {text,json,csv}``
     - text
     - Human-readable table or JSON payload. The shared option accepts
       ``csv``, but this command currently treats non-JSON output as
       text.
   * - ``-v``, ``--verbose``
     - quiet
     - Increase logging.
   * - ``--no-color``
     - false
     - Disable colored output.

The summary includes:

* work directory;
* solver;
* number of iteration files when available;
* loaded iteration;
* final RMS/misfit when available;
* model array shape when available;
* ``log10 rho`` minimum, maximum, and mean when ``rho_2d`` is available.

Examples:

.. code-block:: console

   pycsamt invert results run01
   pycsamt invert results run01 --iteration 50
   pycsamt invert results run01 --format json
   pycsamt invert results modem_run --solver modem

Use ``status`` first when you are not sure whether a directory contains
completed result files. ``results`` is intentionally stricter because it
constructs the full result object.

Plot Overview
-------------

Usage:

.. code-block:: console

   pycsamt invert plot COMMAND WORKDIR [OPTIONS]

All plot commands load an ``InversionResult`` first. They then construct
the corresponding plot class and return a matplotlib figure. A figure is
only written or displayed when you pass ``--save FILE`` or ``--show``.
Without either option the command creates the figure and prints a message
that it was not saved.

Shared plot options:

.. list-table::
   :header-rows: 1
   :widths: 30 22 48

   * - Option
     - Default
     - Meaning
   * - ``--solver {occam2d,modem}``
     - auto
     - Solver family. Required when detection fails.
   * - ``--iteration INT``
     - latest available
     - Load a specific result iteration.
   * - ``--save FILE``
     - none
     - Save the figure as PNG, PDF, SVG, or another matplotlib-supported
       output format.
   * - ``--show``
     - false
     - Display an interactive matplotlib window.
   * - ``--dpi INT``
     - ``100``
     - Figure resolution.
   * - ``--cmap NAME``
     - ``jet_r``
     - Matplotlib colormap. Available on model/pseudo/section/1d/per-site
       commands.
   * - ``-v``, ``--verbose``
     - quiet
     - Increase logging.
   * - ``--no-color``
     - false
     - Disable colored output.

The directory for ``--save`` is created automatically when needed:

.. code-block:: console

   pycsamt invert plot model run01 --save figures/inversion/model.png

Plot Model
----------

Usage:

.. code-block:: console

   pycsamt invert plot model WORKDIR [OPTIONS]

``model`` plots the reconstructed resistivity model. For Occam2D it uses
``pycsamt.models.occam2d.plot.PlotModel``. For ModEM it uses
``pycsamt.models.modem.plot.PlotModel2D`` through the current CLI.

Options in addition to the shared plot options:

.. list-table::
   :header-rows: 1
   :widths: 30 22 48

   * - Option
     - Default
     - Meaning
   * - ``--rho-min FLOAT``
     - ``1.0``
     - Lower resistivity color-scale bound in ohm-m.
   * - ``--rho-max FLOAT``
     - ``1000.0``
     - Upper resistivity color-scale bound in ohm-m.
   * - ``--depth-max METRES``
     - none
     - Maximum displayed depth.
   * - ``--no-stations``
     - false
     - Hide station markers where supported.

Examples:

.. code-block:: console

   pycsamt invert plot model run01 --save section.png
   pycsamt invert plot model run01 --rho-min 1 --rho-max 1000 --depth-max 5000
   pycsamt invert plot model modem_run --solver modem --show

Plot Misfit
-----------

Usage:

.. code-block:: console

   pycsamt invert plot misfit WORKDIR [OPTIONS]

``misfit`` plots convergence through iterations. Both solvers expose a
``PlotMisfit`` class. Occam2D additionally supports roughness and
Lagrange multiplier display controls through the CLI.

Options in addition to the shared plot options:

.. list-table::
   :header-rows: 1
   :widths: 30 22 48

   * - Option
     - Default
     - Meaning
   * - ``--no-roughness``
     - false
     - Occam2D: hide the roughness secondary axis.
   * - ``--lagrange``
     - false
     - Occam2D: add a lower panel for the Lagrange multiplier.

Examples:

.. code-block:: console

   pycsamt invert plot misfit run01 --show
   pycsamt invert plot misfit run01 --save convergence.png --lagrange

Plot Response
-------------

Usage:

.. code-block:: console

   pycsamt invert plot response WORKDIR [OPTIONS]

``response`` compares observed and predicted response curves. Use
``--station`` to focus on one station; omit it to let the plot class
select its default station set.

Options in addition to the shared plot options:

.. list-table::
   :header-rows: 1
   :widths: 30 22 48

   * - Option
     - Default
     - Meaning
   * - ``--station NAME``
     - all/default
     - Plot only the named station.

Examples:

.. code-block:: console

   pycsamt invert plot response run01 --show
   pycsamt invert plot response run01 --station S01 --save S01_resp.png

Plot Pseudo
-----------

Usage:

.. code-block:: console

   pycsamt invert plot pseudo WORKDIR [OPTIONS]

``pseudo`` plots apparent-resistivity and phase pseudosections from the
loaded inversion result. It is available for both Occam2D and ModEM.

Examples:

.. code-block:: console

   pycsamt invert plot pseudo run01 --show
   pycsamt invert plot pseudo run01 --cmap viridis --save pseudo.png

Plot Section
------------

Usage:

.. code-block:: console

   pycsamt invert plot section WORKDIR [OPTIONS]

``section`` is ModEM-only. It calls
``pycsamt.models.modem.plot.PlotSection``. When used on an Occam2D work
directory, the CLI raises a usage error and suggests ``plot 1d`` for
Occam2D depth profiles.

Options in addition to the shared plot options:

.. list-table::
   :header-rows: 1
   :widths: 30 22 48

   * - Option
     - Default
     - Meaning
   * - ``--depth METRES``
     - plot default
     - Depth of the requested slice/section control passed to the plotter.

Examples:

.. code-block:: console

   pycsamt invert plot section modem_run --depth 5000 --show
   pycsamt invert plot section modem_run --depth 10000 --save slice.png

Plot 1d
-------

Usage:

.. code-block:: console

   pycsamt invert plot 1d WORKDIR [OPTIONS]

``1d`` is Occam2D-only. It samples station-centered depth profiles from
the reconstructed 2-D model through
``pycsamt.models.occam2d.plot.PlotSounding1D``. When used on a ModEM
work directory, the CLI raises a usage error and suggests ``plot
section`` for ModEM depth slices.

Options in addition to the shared plot options:

.. list-table::
   :header-rows: 1
   :widths: 30 22 48

   * - Option
     - Default
     - Meaning
   * - ``--stations S01,S02``
     - all/default
     - Comma-separated station names to include.

Examples:

.. code-block:: console

   pycsamt invert plot 1d run01 --show
   pycsamt invert plot 1d run01 --stations S01,S05,S10 --save profiles.png

Plot Per-Site
-------------

Usage:

.. code-block:: console

   pycsamt invert plot per-site WORKDIR [OPTIONS]

``per-site`` is Occam2D-only. It calls
``pycsamt.models.occam2d.plot.PlotSiteMisfit`` and is intended for
checking which stations dominate residuals.

Examples:

.. code-block:: console

   pycsamt invert plot per-site run01 --show
   pycsamt invert plot per-site run01 --save site_rms.png

Plot Grid
---------

Usage:

.. code-block:: console

   pycsamt invert plot grid WORKDIR [OPTIONS]

``grid`` is Occam2D-only. It calls
``pycsamt.models.occam2d.plot.PlotResponseGrid`` and is useful when you
want a compact scan of response quality across many stations.

Examples:

.. code-block:: console

   pycsamt invert plot grid run01 --show
   pycsamt invert plot grid run01 --save response_grid.png

Choosing Plots
--------------

Use this quick guide when deciding what to generate after a run:

.. list-table::
   :header-rows: 1
   :widths: 32 34 34

   * - Question
     - Command
     - Notes
   * - Did the inversion converge?
     - ``invert plot misfit``
     - Start here before interpreting structure.
   * - What does the recovered model look like?
     - ``invert plot model``
     - Available for both solver families.
   * - Does the model fit individual stations?
     - ``invert plot response``
     - Add ``--station`` for focused checks.
   * - Are residuals concentrated by station?
     - ``invert plot per-site``
     - Occam2D only.
   * - Is there a broad data-pattern issue?
     - ``invert plot pseudo``
     - Useful before and after inversion.
   * - What vertical model is below each station?
     - ``invert plot 1d``
     - Occam2D only.
   * - What does a ModEM section/slice show?
     - ``invert plot section``
     - ModEM only.
   * - How do many station responses compare at once?
     - ``invert plot grid``
     - Occam2D only.

Batch Figure Export
-------------------

For a completed Occam2D run, a common export set is:

.. code-block:: console

   pycsamt invert plot model    run01 --save figures/run01/model.png
   pycsamt invert plot misfit   run01 --save figures/run01/misfit.png
   pycsamt invert plot response run01 --save figures/run01/response.png
   pycsamt invert plot pseudo   run01 --save figures/run01/pseudo.png
   pycsamt invert plot 1d       run01 --save figures/run01/profiles_1d.png
   pycsamt invert plot per-site run01 --save figures/run01/per_site.png
   pycsamt invert plot grid     run01 --save figures/run01/response_grid.png

For a completed ModEM run:

.. code-block:: console

   pycsamt invert plot model    modem_run --solver modem --save figures/modem/model.png
   pycsamt invert plot misfit   modem_run --solver modem --save figures/modem/misfit.png
   pycsamt invert plot response modem_run --solver modem --save figures/modem/response.png
   pycsamt invert plot pseudo   modem_run --solver modem --save figures/modem/pseudo.png
   pycsamt invert plot section  modem_run --solver modem --depth 10000 \
       --save figures/modem/section_10km.png

The ``--solver`` flag is optional if the run directory contains standard
solver signatures. It is included above for clarity and for manually
assembled directories.

Text Versus JSON Output
-----------------------

``status`` and ``results`` both support ``--format json``. The shared
CLI option also accepts ``csv``, but these inversion commands currently
emit JSON only for ``--format json`` and otherwise print the text table.
Use text for interactive work and JSON for scripts:

.. code-block:: console

   pycsamt invert status run01 --format json
   pycsamt invert results run01 --format json

Typical script checks are:

* ``ready_to_run`` from ``status`` before launching the solver;
* ``n_iterations_done`` for Occam2D progress tracking;
* ``n_model_files`` for ModEM progress tracking;
* ``rms`` and ``model_shape`` from ``results`` after completion.

Because ``results`` loads a full result object, it can fail on a work
directory that ``status`` can still inspect. That distinction is useful:
``status`` answers "what files are here?", while ``results`` answers
"can pyCSAMT parse this completed inversion?".

Common Failures
---------------

``Cannot detect solver``
    The work directory does not contain standard Occam2D or ModEM file
    signatures. Pass ``--solver occam2d`` or ``--solver modem``.

``No active survey``
    ``invert build`` was called without ``EDI_DIR`` or ``--survey`` and
    no active survey context has been set. Provide an input path or run
    ``pycsamt survey set PATH`` first.

``Figure created but not saved``
    A plot command ran without ``--save`` or ``--show``. Re-run with one
    of those options.

``section is a ModEM-only sub-command``
    Use ``pycsamt invert plot 1d`` for Occam2D depth profiles.

``1d is an Occam2D-only sub-command``
    Use ``pycsamt invert plot section`` for ModEM section/slice plots.

``Error loading results``
    The directory may not contain completed parseable result files, or a
    specific ``--iteration`` was requested but is unavailable. Check
    ``pycsamt invert status WORKDIR`` first.

Non-zero solver exit
    Check ``occam_stderr.log`` or ``modem_stderr.log`` in the work
    directory. Also confirm that the external executable is installed
    and discoverable by the runner.

Practical Notes
---------------

Keep build and run directories explicit. ``invert_run`` is convenient
for quick tests, but named directories such as ``runs/l18_occam_tm`` or
``runs/line03_modem_3d`` are much easier to audit later.

Record solver options in your command history or project notes. Build
options such as ``--freq``, ``--modes``, ``--error-floor-rho``,
``--n-layers``, and ``--cell-size`` materially change the inversion
input.

Use ``status`` before ``run`` when a work directory was copied from
another machine. It gives a quick check that the expected file family is
present before a solver process is launched.

Use ``results --iteration`` when comparing a specific iteration across
different parameter tests. Without it, the result loader uses the latest
available iteration.

Use explicit ``--solver`` in automated scripts. Auto-detection is
helpful interactively, but explicit solver names make scripts clearer and
reduce ambiguity for unusual file layouts.

Python Equivalents
------------------

The CLI is a thin orchestration layer over the model packages. For
Python workflows, use the same underlying modules directly.

Occam2D build/run/read pattern:

.. code-block:: python

   from pathlib import Path

   from pycsamt.emtools._core import ensure_sites
   from pycsamt.models.occam2d.config import OccamConfig
   from pycsamt.models.occam2d.builder import InputBuilder
   from pycsamt.models.occam2d.runner import OccamRunner
   from pycsamt.models.occam2d.results import InversionResult

   sites = ensure_sites("data/line01_edi", recursive=True, verbose=0)

   cfg = OccamConfig(modes=["TE", "TM"])
   cfg.freq_min = 0.1
   cfg.freq_max = 1000.0
   cfg.error_floor_rho = 0.05

   workdir = Path("runs/line01_occam")
   InputBuilder(source=sites, workdir=workdir, config=cfg, verbose=1).build()
   OccamRunner(workdir=workdir, verbose=1).run(max_iter=100, target_misfit=1.05)

   result = InversionResult(workdir=workdir)
   fig = result.plot_model(depth_max=5000)
   fig.savefig("figures/line01_model.png", dpi=150, bbox_inches="tight")

ModEM build/read pattern:

.. code-block:: python

   from pathlib import Path

   from pycsamt.emtools._core import ensure_sites
   from pycsamt.models.modem.config import ModEmConfig
   from pycsamt.models.modem.builder import InputBuilder
   from pycsamt.models.modem.results import InversionResult

   sites = ensure_sites("data/line01_edi", recursive=True, verbose=0)

   cfg = ModEmConfig(mode="3d")
   cfg.freq_min = 0.1
   cfg.freq_max = 1000.0
   cfg.initial_rho = 100.0

   workdir = Path("runs/line01_modem")
   InputBuilder(config=cfg).build(sites, workdir=workdir)

   result = InversionResult(workdir=workdir)

These Python snippets mirror the CLI modules. They are useful when a
project needs custom configuration that is not exposed as a command-line
option yet.

