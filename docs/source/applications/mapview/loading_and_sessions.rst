.. _applications-mapview-loading:

Loading Data And Sessions
=========================

MapView reads the same inputs as the rest of pyCSAMT — folders of EDI station
files — plus ModEM inversion results that it can place in the same scene as the
stations that produced them.  Everything is loaded through one dialog with two
tabs.

Open the loader from the welcome screen's **Load Survey Lines - Start** button
or the **Load lines** button on the top bar.

EDI Stations
------------

.. figure:: ../../_static/applications/mapview/load-survey-data.png
   :alt: The Load Survey Lines dialog on the EDI stations tab
   :class: pycsamt-screenshot

   The **EDI stations** tab: replace or add lines, browse to a survey folder or
   drag files in, then **Load into view**.

On the **EDI stations** tab:

* choose **Replace** to load a fresh survey, or **Add lines** to merge more
  lines into the survey already loaded (existing stations are kept);
* use **Browse folder** to pick a survey directory, or drag EDI/AVG/J files (or
  whole folders) onto the drop zone;
* multi-line surveys load in one pass — each subfolder (or filename prefix)
  becomes a survey line;
* a progress bar reports parsing, and the footer confirms the result (for
  example *Loaded 128 station(s) from 5 line(s)*);
* click **Load into view** to place the survey on the map.

The same folder can be preloaded from the command line with
``pycsamt-mapview --data <folder>`` (see :doc:`installation`).

Once loaded, the **line picker** and **station list** in the left panel control
what is drawn.  The line chips (**All**, then one per line) toggle whole lines;
each station row has an eye toggle to show or hide it.  The selection applies to
**both** views — the 2-D map and the 3-D scene alike.

Inversion Results
-----------------

.. figure:: ../../_static/applications/mapview/inversion-load-data.png
   :alt: The Load Survey Lines dialog on the Inversion results tab
   :class: pycsamt-screenshot

   The **Inversion results** tab imports a ModEM 3-D result folder; a single run
   inverts one volume for the whole survey and is sliced into one panel per
   detected line.

The **Inversion results** tab imports model output so you can inspect
resistivity structure in the same scene as the survey geometry:

* **Browse folder** to a ModEM results directory;
* leave **Match coordinates from already-loaded EDI stations** enabled so the
  importer aligns the model laterally using the station positions already on the
  map — load the survey's EDIs first, then the results;
* click **Import inversion results**.

A single ModEM run inverts one volume for the whole survey; MapView slices it
into one panel per detected line, which is what you then see in the 3-D fence,
block, depth-slice, and iso-surface renderings (see :doc:`views`).

Sites
-----

The **Sites** button on the top bar opens station/sites settings — options that
control how stations are interpreted and displayed.  Use it when station
identifiers, grouping, or geometry need adjusting after loading.

Sessions
--------

Loaded data lives in a **per-browser-tab session cache** on the server:

* you can reload the page without re-uploading the survey;
* two browser tabs do not interfere with each other — a second tab starts
  empty by design;
* the **Session** button on the top bar saves and restores your working state;
* restarting the server clears the cache.

MapView never writes into your data folders, so loading, reloading, and
restarting are always safe.

.. seealso::

   :doc:`/user_guide/data_loading`
       How pyCSAMT reads EDI directories and builds survey containers.
