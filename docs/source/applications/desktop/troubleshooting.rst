.. _applications-desktop-troubleshooting:

Troubleshooting
===============

This page lists common pyCSAMT desktop problems, what they usually mean, and
the first checks to make.  Start with the symptom you see in the GUI, status
bar, log dock, or terminal.

Quick Checks
------------

Before debugging deeply:

1. Confirm that the ``pycsamt-v2`` environment is active.
2. Launch the desktop from a terminal so startup errors are visible.
3. Keep the desktop log dock visible with **View > Log**.
4. Check the status bar after every load, correction, export, or inversion
   action.
5. Reload exported EDIs in a fresh session before using them downstream.

Useful launch commands:

.. code-block:: bash

   pycsamt-desktop
   pycsamt-gui
   python -m pycsamt.app.desktop

Missing PySide6
---------------

**Symptom**

The terminal prints:

.. code-block:: text

   PySide6 is required for the desktop app.
   Install it with:  pip install 'pycsamt[app]'

**Cause**

The desktop GUI entry point imports PySide6 at launch.  The core package can be
imported without Qt, but the desktop cannot open without the Qt bindings.

**Fix**

Install the application extras in the active environment:

.. code-block:: bash

   pip install "pycsamt[app]"

If you are working from the source checkout, reinstall the local project into
the environment:

.. code-block:: bash

   pip install -e ".[app]"

Then run:

.. code-block:: bash

   python -m pycsamt.app.desktop

Command Not Found
-----------------

**Symptom**

``pycsamt-desktop`` or ``pycsamt-gui`` is not recognized.

**Cause**

The package entry points were not installed into the active environment, or the
wrong Python environment is active.

**Fix**

Use the module entry point first:

.. code-block:: bash

   python -m pycsamt.app.desktop

If that works, reinstall the package so console scripts are regenerated:

.. code-block:: bash

   pip install -e ".[app]"

If the module entry point fails too, verify that ``python`` belongs to the
``pycsamt-v2`` environment.

Desktop Opens But No Data Load
------------------------------

**Symptom**

The status bar stays at **No data loaded**, or the load dialog does not enable
**Load Data**.

**Cause**

No compatible file was selected, the active filter does not match the files,
or a dropped folder did not contain matching files.

**Fix**

In the load dialog:

* choose the correct filter: EDI, AVG, J/ModEM, or all supported;
* use **Browse Files...** for a known set of station files;
* use **Browse Folder...** only when the folder really contains the target
  format;
* check that files end with the expected extension, such as ``.edi``;
* load one survey line at a time when troubleshooting.

If drag-and-drop reports that no files were found, switch the format filter and
try **Browse Folder...** so the selected folder is explicit.

Load Failed
-----------

**Symptom**

The status bar shows **Load failed: ...**, or the log dock records a loader
error.

**Cause**

The loader runs in a background thread and builds a ``Sites`` collection from
the selected paths.  Failures usually mean the selected files are not valid for
the chosen loader, the path no longer exists, or no usable stations with
impedance data were resolved.

**Fix**

1. Try loading one known-good EDI file.
2. Try loading one line folder instead of a parent project folder.
3. Remove mixed unsupported files from the selection.
4. Check whether the EDI file contains valid impedance ``Z`` data.
5. If you converted or recomputed EDIs, load the exported folder separately
   and compare station count with the original.

For script-level diagnosis, test the canonical loader:

.. code-block:: python

   from pycsamt.emtools import ensure_sites

   sites = ensure_sites(r"path\to\edi_folder", recursive=True, verbose=1)
   print(len(sites))

If strict loading is used in custom code, an empty input can raise:

.. code-block:: text

   ensure_sites(strict=True): no sites were resolved from the given input.

No Data Dialog
--------------

**Symptom**

A tool refuses to open and asks you to load data first.

**Cause**

Many scientific panels need the active ``Sites`` object.  The main window will
guard tools such as QC, corrections, profiles, maps, and advanced analysis
when no survey is loaded.

**Fix**

Load EDI data first with **Open / Load EDI...**.  After loading, check that the
status bar reports the expected number of stations.  Then reopen the tool.

Station Count Is Wrong
----------------------

**Symptom**

The station table loads, but the count is lower or higher than expected.

**Likely causes**

* a parent folder loaded several survey lines together;
* some files do not contain usable impedance data;
* duplicate station names were resolved by the loader;
* the wrong output folder was loaded after recomputation or conversion.

**Fix**

Load one line folder at a time.  Check the station metadata export or survey
overview, then compare with the raw folder contents.  If the station count is
lower after correction or recomputation, reload the exported EDI folder and
inspect the log/manifest for failed stations.

Blank Or Empty Plots
--------------------

**Symptom**

A map, profile, QC, or advanced figure opens but appears blank or says a plot
is unavailable.

**Likely causes**

* no station is selected for a station-specific profile;
* the active data do not contain the component required by the plot;
* the selected frequency/period band has no usable values;
* coordinates or elevations are missing for map/topography plots;
* an external result file was not loaded for inversion result plots.

**Fix**

1. Select a station in the main station table.
2. Return to the map/profile viewers and confirm the active survey is valid.
3. Use QC coverage to check whether the frequency band exists.
4. For map plots, confirm that latitude/longitude are present and plausible.
5. For inversion plots, confirm the result folder contains the expected solver
   outputs.

Figure Export Says No Figure
----------------------------

**Symptom**

The status bar says **No figure to export.**

**Cause**

The main **Export Figure** action searches visible panel windows for an active
matplotlib canvas.  No exportable canvas was visible or active.

**Fix**

Open the plot you want to save, make sure it is visible, then use either:

* the panel's own **Export** button, or
* **More > Export Figure** from the main toolbar.

For many open figures, use **Tools > Batch Export Plots...**.

Batch Export Saves Fewer Figures Than Expected
----------------------------------------------

**Symptom**

Batch export reports fewer figures than the number of windows you expected.

**Cause**

Only visible panel canvases are collected.  Hidden panels, tabs without
matplotlib canvases, and dialogs that were closed before export are not saved.

**Fix**

Open every panel and tab you want to capture before starting batch export.
After export, check the progress log and output folder.

Coordinate Or Map Problems
--------------------------

**Symptoms**

Stations plot in the wrong area, all stations stack on one point, contours look
unreasonable, or topography is missing.

**Likely causes**

* latitude/longitude are absent or swapped;
* projected coordinates were interpreted as geographic coordinates;
* elevation metadata are missing;
* multiple survey lines were loaded as one profile;
* an edited coordinate product was committed without rechecking the map.

**Fix**

Use the map viewer before correction or inversion.  Export station metadata to
CSV/JSON and inspect latitude, longitude, frequency count, and station names.
If coordinates were corrected, commit intentionally, then reopen the map and
profile viewers before processing.

Static Shift Correction Looks Too Smooth
----------------------------------------

**Symptom**

After correction, curves line up but geological contrast disappears.

**Cause**

The correction window can make a visually smoother profile even when the
method is over-applied.  A large station window or unjustified source/correction
assumption can suppress real lateral variation.

**Fix**

Return to QC and profile views:

* compare before/after curves and pseudosections;
* use **Overlay** and **Diff** views;
* reduce the smoothing or half-window parameter;
* check whether a strongly 3-D area makes static-shift assumptions unsafe;
* keep the raw survey available until the correction is defensible.

If AMA returns no useful factors, do not force the correction.  It may be
telling you that the assumptions are not satisfied.

Pipeline Stops Or Exports Nothing
---------------------------------

**Symptoms**

The pipeline stops on a step, reports an error, or the final export step writes
no EDIs.

**Likely causes**

* step 1 has no input survey;
* a required output folder was not selected;
* a correction or rotation step received data without required components;
* all stations were filtered out by an earlier step;
* the export folder is invalid or not writable.

**Fix**

Run the failing step alone.  Check the centre parameter panel, then read the
right-side log tab.  For export, choose a new empty output folder.  Save the
pipeline JSON only after the chain has produced acceptable intermediate
diagnostics.

Inversion Workdir Or Binary Errors
----------------------------------

**Symptoms**

The inversion window reports:

.. code-block:: text

   Working directory not found: ...
   Occam2D binary not found: ...

**Cause**

The selected solver working directory does not exist, or the configured
external binary path is wrong.

**Fix**

Create or select a run-specific working directory before running.  For external
solvers, set the executable path in **File > Preferences...**.  If the solver
can build input files without running the external binary, use that as a first
validation step and inspect the generated files.

Occam2D Or Solver Result Will Not Plot
--------------------------------------

**Symptoms**

Result plots say no response data, no model, no station offsets, or result
files are missing.

**Likely causes**

* the solver did not finish successfully;
* the result folder is not the same workdir used for the run;
* input files were generated but the external solver was never executed;
* the requested station name is not present in the result data;
* model, mesh, response, or log files are incomplete.

**Fix**

Open the inversion run folder and check that it contains solver output files,
not only generated inputs.  Read the solver log before plotting.  If the input
files exist but results do not, rerun the external solver or point the desktop
to the correct workdir.

AI Or Agent Runs Fail
---------------------

**Symptoms**

AI inversion or Agent Master stops, reports a worker error, or produces no
useful output.

**Likely causes**

* no active survey is loaded;
* the data cannot be converted to the expected observation array;
* required model dependencies are missing;
* an API key is missing for remote agent workflows;
* the request asks the agent to act without enough context.

**Fix**

Load data first and confirm QC/profile plots.  For agents, provide the survey
line, data state, and the question you want answered.  If an API key is needed,
set it in **File > Preferences...** and remember that it is stored in the
desktop session file.

Session Or Layout Problems
--------------------------

**Symptoms**

Windows reopen off-screen, theme or recent files look wrong, or old solver
paths keep returning.

**Cause**

Desktop state is stored in ``~/.pycsamt/session.json``.  It remembers window
geometry, recent files, selected station, theme, solver paths, API key setting,
and other preferences.

**Fix**

First try **File > Preferences...** and **File > Save Session**.  If geometry
is unusable, close the desktop and move or rename the session file:

.. code-block:: text

   ~/.pycsamt/session.json

The next launch will create a fresh default session.

Windows PowerShell Profile Warning
----------------------------------

**Symptom**

When running commands in PowerShell, you see a message similar to:

.. code-block:: text

   profile.ps1 cannot be loaded. The file is not digitally signed.

**Cause**

This is a PowerShell execution-policy warning for the user's shell profile.  It
is not a pyCSAMT desktop error.

**Fix**

Run the pyCSAMT command again in a shell where the profile warning is resolved,
or start PowerShell with a profile-free mode when debugging commands.  The GUI
itself is affected only if your profile was supposed to activate the
``pycsamt-v2`` environment.

When To Rebuild Or Reinstall
----------------------------

Reinstall the package when:

* console scripts such as ``pycsamt-desktop`` are missing;
* imports work in one shell but not another;
* PySide6 or optional app dependencies were added after installation;
* entry points were changed in the source checkout.

From the repository root:

.. code-block:: bash

   pip install -e ".[app]"

Then relaunch:

.. code-block:: bash

   pycsamt-desktop

What To Include In A Bug Report
-------------------------------

Include:

* the command used to launch the desktop;
* the active Python environment name;
* the exact status-bar or log-dock message;
* the file type and approximate number of stations loaded;
* whether the issue occurs with raw EDIs or exported/corrected EDIs;
* screenshots of the relevant panel if the issue is visual;
* the pipeline JSON, manifest CSV, or inversion log when the problem involves
  processing output.

This information usually identifies whether the problem is environment,
loading, data quality, configuration, or a real application bug.
