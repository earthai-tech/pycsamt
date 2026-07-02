.. _applications-desktop-workspace:

Workspace And Navigation
========================

The pyCSAMT desktop workspace is deliberately lean.  The main window keeps the
loaded survey visible, while maps, profiles, QC, corrections, modelling, and
agents open as independent panel windows.  This lets you keep the station list
anchored on one screen and move the scientific views around it.

.. figure:: ../../_static/applications/desktop/home.png
   :alt: pyCSAMT desktop main workspace after launch
   :class: pycsamt-screenshot

   The main workspace is the command centre: toolbar at the top, station list
   on the left, survey and station summaries on the right, and a compact log
   dock at the bottom.

Main Window Regions
-------------------

The main window has four persistent regions:

* **Toolbar** -- one-click access to loading, saving, viewers, processing, and
  export tools.
* **Station list** -- searchable table of loaded stations.
* **Survey overview and station detail** -- summary cards for the whole line
  and the currently selected station.
* **Log and status areas** -- loading messages, processing messages, progress,
  and the current data-state label.

The main window owns the active survey.  When data are loaded or replaced by a
committed correction, the station list and open panel windows are refreshed
from that active survey.

Station List
------------

The station list is the primary navigation surface after loading EDI files.
Use the search field above the table to filter station IDs while keeping the
loaded survey unchanged.  Selecting a station updates the detail card and also
notifies open profile/map windows.

Normal station-list workflow:

1. Load one survey line or a focused folder of EDI files.
2. Use the search field to isolate a station or station prefix.
3. Select a row to inspect its metadata in the detail card.
4. Open **Profile** or **Map** from the toolbar to inspect the selected station
   in context.
5. Clear the search text before checking line-wide coverage or profile trends.

Filtering is only a view operation.  It does not remove stations from the
active ``Sites`` object and it does not affect export or processing windows.

Survey Overview
---------------

The survey overview card summarizes the loaded data set.  Use it immediately
after loading to confirm that the number of stations, file set, and frequency
coverage look plausible before opening heavier processing tools.

If the overview does not match the line you intended to load, stop and reload
the data.  It is better to fix a folder-selection mistake before QC, maps, or
corrections start producing convincing figures from the wrong input.

Station Detail Card
-------------------

The station detail card shows metadata for the selected station.  It is the
fastest way to check whether a row has the expected station ID, coordinates,
elevation, and data availability.

Use the detail card before opening a station-specific profile:

* confirm the station name and location;
* check that coordinates and elevation are populated;
* use the profile action to jump directly into the profile viewer;
* compare neighbouring stations when a response curve looks unusual.

Panel Windows
-------------

Most scientific tools are independent floating windows.  Opening a tool from
the toolbar or **View** menu shows the same persistent window; closing a panel
normally hides it and preserves its geometry for the session.

Core panel windows:

.. list-table::
   :header-rows: 1
   :widths: 22 34 44

   * - Panel
     - Opens From
     - Use It For
   * - **Profile**
     - Toolbar, **View > Profile Viewer**, ``Ctrl+P``
     - Station curves, pseudosections, phase tensors, and profile sections.
   * - **Map**
     - Toolbar, **View > Map Viewer**, ``Ctrl+M``
     - Station geometry, contours, map editing, and spatial context.
   * - **QC**
     - Toolbar, **View > QC Dashboard**, ``Ctrl+Q``
     - Coverage, SNR, skew, dimensionality, static-shift, and distortion checks.
   * - **Corrections**
     - Toolbar, **View > Data Corrections**, ``Ctrl+R``
     - Previewing, stacking, and committing data corrections.
   * - **Forward**
     - Toolbar, **View > Forward Modelling**, ``Ctrl+F``
     - 1-D/2-D/3-D synthetic models and inversion starting models.
   * - **Inversion**
     - Toolbar, **View > Inversion Wizard**, ``Ctrl+I``
     - Solver setup, input-file generation, AI inversion, and run logs.
   * - **Interpret**
     - Toolbar, **View > Interpretation Studio**, ``Ctrl+Shift+I``
     - Interpreting inversion results and exporting interpretation products.
   * - **Pipeline**
     - Toolbar, **View > Processing Pipeline**, ``Ctrl+Shift+P``
     - Repeatable load/QC/edit/correct/rotate/export processing chains.
   * - **Agents**
     - Toolbar, **View > Agent Master**, ``Ctrl+Shift+A``
     - Agent-assisted review, explanations, and guided workflow checks.

The **More** toolbar menu contains secondary tools such as **TDEM Analysis**,
**Advanced Tools**, and **Export Figure**.  The same tools are also available
from the menu bar.

Toolbar And Menus
-----------------

The toolbar is optimized for the normal survey workflow:

1. **Open** and **Save** manage the survey session.
2. **Profile** and **Map** inspect loaded data.
3. **QC** and **Corrections** diagnose and condition the data.
4. **Forward**, **Inversion**, and **Interpret** move from modelling to
   interpretation.
5. **Pipeline** repeats a known processing chain.
6. **Agents** opens the conversational workflow layer.
7. **More** exposes secondary tools and figure export.
8. **Theme** toggles dark and light appearances.

The menu bar mirrors these actions and adds recent files, preferences, help,
documentation, GitHub, and about dialogs.  If a toolbar button is hidden by a
small screen, use the menu bar as the complete command list.

Log Dock And Status Bar
-----------------------

The log dock is the compact message area at the bottom of the main window.  It
records load events, session saves, errors, and messages emitted by workflows.
Use **View > Log** to hide or restore it.

The status bar provides the fastest state check:

* the left label reports whether data are loaded, corrected, or converted;
* the frequency label shows frequency-range context when available;
* the progress bar appears while background loading or processing is running;
* the readiness label changes when loading or errors occur;
* short status messages appear temporarily after actions such as saving.

When something looks wrong, read the status bar first and then the log dock.
Most load failures and missing-data messages are visible there before they are
visible in a plot.

Recent Files
------------

Use **File > Recent Files** to reload one of the last opened inputs.  Recent
files are stored in the desktop session and capped by the ``max_recent_files``
preference.  Choosing a recent entry starts the normal loading path, so open
panels refresh the same way they do after a manual file selection.

Theme And Preferences
---------------------

Use the theme button on the toolbar, or **View > Theme**, to switch between
dark and light mode.  The selected theme is saved in the desktop session and
is reapplied on the next launch.

Open **File > Preferences...** for persistent application settings such as
solver paths, default inversion working directory, logging level, map tile
provider, API key storage, and recent-file limits.

Saved Session State
-------------------

The desktop stores its session at ``~/.pycsamt/session.json``.  The file is
written when you save the session and when the application closes normally.

The session remembers:

* selected theme;
* recent files and last data directory;
* selected station;
* preferred frequency range and map overlay;
* main-window geometry and dock state;
* floating-panel geometries;
* solver paths and default inversion work directory;
* logging, map-tile, API-key, and recent-file preferences.

The session does not replace exported project files.  Save corrected EDIs,
pipeline JSON, figures, inversion inputs, and reports into the project output
folders described in :doc:`processing_workflows`.

Workspace Habits
----------------

For steady desktop work:

* keep the main window visible as the station and state reference;
* open profile and map viewers before processing windows;
* keep the log dock visible during loading, correction, inversion, and export;
* save the session after arranging panel windows on a multi-monitor setup;
* reload from exported EDIs when you need to verify a corrected data product;
* use one survey line per session when preparing inversion inputs.

This workspace style keeps inspection, processing, and export decisions tied
to the same active survey instead of scattering state across disconnected
windows.
