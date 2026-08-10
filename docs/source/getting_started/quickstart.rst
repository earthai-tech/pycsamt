.. _quickstart:

Quickstart
==========

This is the shortest complete pyCSAMT workflow: load a bundled EDI survey,
confirm what was read, normalize it for processing, and save one diagnostic
figure. It assumes pyCSAMT is already installed. See :doc:`installation` if
the import or command-line verification has not succeeded yet.

The example uses the WILLY L18PLT AMT line shipped with the source repository.
When working from an installed wheel, replace the bundled path with your own
EDI directory. For AVG, J-format, spectra, time-series, TEM/TDEM, or Stratagem
inputs, choose the correct route in :doc:`data_formats` first.

Load and inspect
----------------

Run these commands from the repository root:

.. code-block:: pycon

   >>> from pycsamt.api import read_edis
   >>> survey = read_edis(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     recursive=False,
   ...     strict=True,
   ...     on_dup="replace",
   ...     progress=False,
   ... )
   >>> survey.n_sites
   28
   >>> survey.stations[:3]
   ['23-18-001A', '23-18-002U', '23-18-003A']
   >>> len(survey.errors())
   0

The expected count is specific to this dataset. For field data, compare it
with the survey manifest. Zero parser errors means every discovered file in
this call was readable; it does not validate coordinates, units, uncertainty,
or transfer-function quality.

The survey inventory gives the fastest structural check:

.. code-block:: pycon

   >>> inventory = survey.summary()
   >>> print(inventory)
   APIFrame: edi_survey_summary
   kind: edi.summary
   shape: 28 rows x 6 columns
   columns: station, path, n_freq, tipper, spectra, ts
   numeric: 1 columns
   missing: 0.0%
   source: data/AMT/WILLY_DATA/L18PLT
   >>> table = inventory.to_pandas(copy=True)
   >>> sorted(table["n_freq"].unique().tolist())
   [53]

Every station has 53 impedance-frequency rows. The equal counts are useful,
but confirm the frequency values themselves before comparing arrays element by
element.

Normalize for processing
------------------------

The public survey view is convenient for inspection. Electromagnetic
processing functions normalize inputs to :class:`pycsamt.site.base.Sites` via
:func:`pycsamt.emtools.ensure_sites`:

.. code-block:: pycon

   >>> from pycsamt.emtools import ensure_sites
   >>> sites = ensure_sites(
   ...     survey.collection,
   ...     recursive=False,
   ...     order_by="input",
   ...     strict=True,
   ... )
   >>> type(sites).__name__, len(sites)
   ('Sites', 28)
   >>> sites[0].name, sites[0].summary()["nfreq"]
   ('18-001A', 53)

``order_by="input"`` preserves the supplied collection order. It is a safe
starting point only when the file sequence is already meaningful. Coordinate-
derived chainage and other ordering modes are explained in
:doc:`configuration` and :doc:`../user_guide/data_loading`.

Save one diagnostic figure
--------------------------

The survey fingerprint compares phase-sensitive quantities across station and
period. Save it after the inventory checks:

.. code-block:: pycon

   >>> from pathlib import Path
   >>> from pycsamt.emtools import plot_survey_fingerprint
   >>> output = Path("results/quickstart")
   >>> output.mkdir(parents=True, exist_ok=True)
   >>> fig = plot_survey_fingerprint(
   ...     survey.collection,
   ...     quantities=["skew", "ellipt", "s1"],
   ...     render="imshow",
   ...     plot_kws={"interpolation": "bilinear"},
   ...     station_grid=True,
   ...     period_range=(1e-4, 1.0),
   ...     recursive=False,
   ...     title="L18PLT quick survey fingerprint",
   ...     figsize=(11.2, 7.6),
   ... )
   >>> target = output / "survey_fingerprint.png"
   >>> fig.savefig(target, dpi=180, bbox_inches="tight")
   >>> target.as_posix()
   'results/quickstart/survey_fingerprint.png'

``pcolormesh`` is the function's default renderer; this quickstart selects
``imshow`` with bilinear interpolation to match the captured image. The detailed controls for
selecting one or more quantities and overriding panel colormaps are introduced
in :doc:`first_survey`.

.. figure:: /images/tutorials/read_edi_survey/survey_fingerprint.png
   :align: center
   :width: 100%
   :alt: Skew, ellipticity, and phase-tensor fingerprint for L18PLT.

   L18PLT phase-sensitive survey fingerprint over periods from
   :math:`10^{-4}` to :math:`1` s.

The latter half of the line contains coherent station-period changes across
several panels. That agreement makes the region worth further phase-tensor,
strike, noise, and source-effect checks, but it does not by itself identify a
geological boundary or justify removing stations.

Complete script
---------------

The complete executable version is kept outside the page so repeated setup and
file-writing code does not interrupt the quickstart narrative:

.. code-dropdown:: /../../scripts/quickstart_first_survey.py
   :language: python
   :linenos:
   :title: View the complete quickstart script

Run it from the repository root with:

.. code-block:: console

   python docs/scripts/quickstart_first_survey.py

The script prints:

.. code-block:: text

   Loaded 28 stations
         station  n_freq  tipper
   0  23-18-001A      53   False
   1  23-18-002U      53   False
   2  23-18-003A      53   False
   3  23-18-004A      53   False
   4  23-18-005U      53   False
   Saved diagnostic: results/quickstart/survey_fingerprint.png

Do not proceed solely because this script completes. Verify station identity,
coordinates, component availability, frequency alignment, and unexpected
diagnostic patterns before correction or inversion.

Continue with the appropriate depth
-----------------------------------

* :doc:`first_survey` explains every inventory check and interprets the figure
  in more detail.
* :doc:`../tutorials/read_edi_survey` provides a complete EDI-reading tutorial.
* :doc:`../tutorials/inspect_and_qc_survey` introduces richer quality-control
  diagnostics.
* :doc:`../user_guide/emtools/index` organizes processing methods and their
  scientific assumptions.
* :doc:`../user_guide/inversion/index` helps choose an inversion route only
  after the observations have been reviewed.
