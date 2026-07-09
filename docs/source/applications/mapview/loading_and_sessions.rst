Loading Data And Sessions
=========================

MapView reads the same inputs as the rest of pyCSAMT — folders of EDI
station files — plus inversion results that it can place on the map next
to the stations that produced them.

EDI Stations
------------

Open the load dialog from the welcome panel or the **Load** button and
use the **EDI stations** tab:

* drag-and-drop EDI files into the drop zone, or pick a whole folder;
* multi-line surveys load in one go — each subfolder (or filename
  prefix) becomes a survey line you can toggle from the line picker;
* the same folder can be preloaded from the command line with
  ``pycsamt-mapview --data <folder>``.

Once loaded, the line picker in the sidebar controls which profiles are
visible in every view — station map, pseudosection, and 3-D scene alike.

Inversion Results
-----------------

The **Inversion results** tab of the load dialog imports model output so
you can inspect resistivity structure in the same scene as the survey
geometry. When the result file carries no coordinates of its own,
enable *match coordinates from already-loaded EDI stations* — the
importer aligns the model laterally using the station positions already
on the map.

Sessions
--------

Loaded data lives in a per-browser-tab session cache on the server: you
can reload the page without re-uploading, and two tabs do not interfere
with each other. Restarting the server clears the cache — MapView never
writes into your data folders.

.. seealso::

   :doc:`/user_guide/data_loading`
       How pyCSAMT reads EDI directories and builds survey containers.
