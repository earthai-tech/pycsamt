.. _applications-desktop-loading:

Loading Data And Sessions
=========================

The desktop app keeps one active survey session at a time.  Loading data
creates the in-memory :class:`pycsamt.site.Sites` collection that feeds the
station table, survey overview, profile viewer, map viewer, QC dashboard,
correction tools, forward modelling, inversion setup, and agent workflows.

Open Survey Data
----------------

Use **Open / Load EDI...** from the toolbar or **File > Open / Load EDI...**
from the menu.  The dialog is named **Open Survey Data** because it is the
front door for compatible survey files, not only a single EDI file picker.

.. figure:: ../../_static/applications/desktop/load-edis.png
   :alt: Selecting EDI files in the desktop loader
   :class: pycsamt-screenshot

   Selecting multiple EDI files from one survey-line folder before loading.

The loader offers four format filters:

.. list-table::
   :header-rows: 1
   :widths: 24 28 48

   * - Format
     - File pattern
     - Typical use
   * - ``EDI``
     - ``*.edi``
     - Standard MT, AMT, and CSAMT station transfer-function files.
   * - ``AVG``
     - ``*.avg``
     - Zonge AVG inputs that may be converted through desktop tools.
   * - ``J / ModEM``
     - ``*.j``
     - J-format or ModEM-style transfer-function sources.
   * - ``All supported``
     - ``*.edi``, ``*.avg``, ``*.j``
     - Mixed folders when you want the dialog to collect all compatible files.

Choose Files Or Folders
-----------------------

There are three normal ways to populate the file list:

* **Browse Files...** opens a multi-select file dialog for the active format.
* **Browse Folder...** scans the selected folder recursively and adds matching
  files.
* Drag files or folders onto the drop zone; folders are recursively searched
  with the same active format filter.

The dialog keeps the selection visible before any data are loaded.  Use
**Remove Selected** to drop individual paths or **Clear All** to restart the
selection.  **Load Data** stays disabled until at least one compatible file is
selected.

For line-based surveys, load one line folder at a time when you want the
desktop views to stay focused on a single profile.  In the WILLY example data,
that usually means opening one of these folders:

.. code-block:: text

   data/AMT/WILLY_DATA/L18PLT
   data/AMT/WILLY_DATA/L22PLT
   data/AMT/WILLY_DATA/L26PLT
   data/AMT/WILLY_DATA/L30PLT
   data/AMT/WILLY_DATA/L34PLT

For a quick inspection pass, loading a parent folder is acceptable; the dialog
will recurse and collect matching files.  For interpretation, correction, and
inversion preparation, a single survey line is usually easier to review and
less noisy in the station table.

What Happens During Loading
---------------------------

When you confirm the dialog, the desktop starts a background loader so the
window remains responsive.  The status bar and progress bar show loading
progress while pyCSAMT parses the files.

For EDI data, the loader builds a :class:`pycsamt.site.Sites` collection and a
station table with these fields:

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Column
     - Meaning
   * - ``ID``
     - Station identifier, normally derived from the EDI station name.
   * - ``Line``
     - Survey-line label when a line mapping is available.
   * - ``Latitude`` / ``Longitude``
     - Station coordinates read from file metadata.
   * - ``Elevation``
     - Station elevation, when available.
   * - ``N_freq``
     - Number of frequency samples in the loaded transfer function.
   * - ``Tipper``
     - Whether tipper information is present.

After a successful load, the same survey object is sent to the profile, map,
QC, correction, advanced tools, forward, inversion, interpretation, pipeline,
and agent windows.  This keeps the desktop session coherent: selecting a
station in the table updates the detail card, profile actions, map actions,
and downstream workflow panels.

Station Search And Selection
----------------------------

The station table search box filters visible stations by text.  Use it to
quickly isolate a station number, line label, or naming pattern before opening
profile and map views.

When a station is selected, the detail panel shows coordinates, elevation,
frequency coverage, period range, tipper availability, and a compact quality
indicator.  The **Open Profile** and **Show on Map** actions become available
once the station has enough data for those views.

Recent Files
------------

Every path accepted by the loader is added to **File > Recent Files**.  The
recent list is capped by the desktop preferences and is restored when the app
restarts.  Selecting an entry from the menu reloads that path directly.

Recent files are most useful for single-file checks and small station groups.
For full survey lines, prefer **Open / Load EDI...** and reload the folder so
new or removed files are picked up consistently.

Save And Restore The Desktop Session
------------------------------------

Use **Save Session** to persist the desktop UI state.  The session file is
stored at:

.. code-block:: text

   ~/.pycsamt/session.json

The session stores desktop preferences and layout information, including:

* theme;
* recent files;
* last data directory;
* selected station;
* frequency range settings;
* overlay choice;
* main-window geometry and dock state;
* per-window geometry for profile, map, QC, agent, inversion, and other tool
  windows;
* configured solver paths, default inversion working directory, logging level,
  tile provider, and API key if entered in preferences.

The session does not replace the survey files themselves.  Keep the original
EDI, AVG, or J files in a stable project folder; the desktop session stores
where you worked and how the interface was arranged, while the scientific data
remain on disk.

Loading Recomputed EDIs
-----------------------

If the desktop has just run an EDI recomputation workflow, the open-data dialog
can expose a **Load Recomputed EDIs** shortcut.  That button loads all
``*.edi`` files from the last recomputation output folder without manual
navigation.

Use this flow when you have rotated, filtered, filled, or rewritten stations
and want the corrected files to become the active desktop survey.  After the
load completes, the station table and downstream tools operate on the
recomputed EDI set.

Practical Loading Checklist
---------------------------

Before starting QC, correction, or inversion setup:

* load one coherent survey line or one intentionally mixed survey group;
* confirm the station count matches the expected line folder;
* scan the station table for missing coordinates or zero frequency counts;
* select a representative station and open its profile;
* open the map view to confirm station ordering and coordinate sanity;
* save the session once the workspace is arranged for the task.
