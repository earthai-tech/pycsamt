.. _emtools_inspect:

First-Look Survey Inspection
============================

``pycsamt.emtools.inspect`` is the first module to use after loading a
survey. It answers practical questions before you run deeper
diagnostics, static-shift correction, dimensionality analysis, or
inversion:

* which stations loaded correctly?
* which stations have impedance and tipper data?
* do all stations share the same frequency grid?
* are the apparent-resistivity and phase curves sensible?
* where are the first obvious pseudo-section anomalies?
* which single station deserves a full response plot?

Full callable signatures live in the :doc:`API reference <../../api/emtools>`.
This page explains the workflow and gives concrete code patterns.

Why Inspect First
-----------------

Inspection is not interpretation yet. It is the quality gate between
"the files loaded" and "the data are ready for scientific decisions".
The inspection tools are intentionally plain: tables, coverage masks,
simple curves, pseudo-sections, tipper components, and one complete
station dashboard.

Run this stage before:

* deciding which frequency band to keep;
* comparing survey lines;
* trusting tipper-based induction arrows;
* building phase-tensor or dimensionality products;
* fitting a model response to observed data.

Load The Survey
---------------

The inspection functions all call ``ensure_sites`` internally, but it is
still useful to normalize once at the top of a script.

.. code-block:: python
   :linenos:

   from pathlib import Path

   from pycsamt.emtools import ensure_sites

   edi_dir = Path("data/AMT/WILLY_DATA/L18PLT")
   survey = ensure_sites(
       edi_dir,
       recursive=True,
       on_dup="replace",
       strict=True,
       verbose=1,
   )

Use ``strict=True`` when a missing or empty dataset should stop the
workflow. For exploratory notebooks, ``strict=False`` can be more
convenient because the plotting functions will draw "no data" messages
instead of failing immediately.

The Inspection Workflow
-----------------------

The module is easiest to use as a sequence.

.. list-table::
   :header-rows: 1
   :widths: 25 40 35

   * - Step
     - Question
     - Tool
   * - inventory
     - What stations, period ranges, coordinates, and tippers exist?
     - ``sites_summary``
   * - required sections
     - Which stations are missing ``mt`` or ``tipper``?
     - ``list_missing_sections``
   * - frequency grid
     - Do stations share the same frequency samples?
     - ``frequency_coverage`` and ``plot_coverage``
   * - quick curves
     - Are rho/phase curves plausible?
     - ``plot_rhoa_phi``
   * - survey image
     - Where are station-period anomalies?
     - ``pseudosection``
   * - tipper check
     - Is real tipper present and stable?
     - ``plot_tipper_components``
   * - station dashboard
     - What does one station look like in full?
     - ``plot_station_response``

Per-Site Summary
----------------

``sites_summary`` returns one row per station. By default it reports
station name, number of frequency samples, whether tipper data are
present, period range, and coordinates.

.. code-block:: python
   :linenos:

   import pandas as pd

   from pycsamt.emtools import ensure_sites, sites_summary

   survey = ensure_sites("data/AMT/WILLY_DATA/L18PLT", strict=True)

   summary = sites_summary(survey, api=False)
   print(summary.head())

   overview = {
       "n_sites": len(summary),
       "has_any_tipper": bool(summary["has_tipper"].any()),
       "n_freq_values": sorted(summary["n_freq"].unique()),
       "period_min": float(summary["period_min"].min()),
       "period_max": float(summary["period_max"].max()),
   }
   print(pd.Series(overview))

Read this table before plotting anything. A survey with mixed
``n_freq`` values needs frequency-grid attention. A survey with
``has_tipper=False`` everywhere should not be sent into tipper or
induction-arrow interpretation.

Choose Summary Columns
----------------------

The ``fields`` argument lets you keep the inventory narrow when you are
printing reports or comparing several lines.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import sites_summary

   compact = sites_summary(
       "data/AMT/WILLY_DATA/L18PLT",
       fields=(
           "station",
           "n_freq",
           "period_min",
           "period_max",
           "has_tipper",
       ),
       api=False,
   )

   print(compact.to_string(index=False))

The returned object may be an API-aware frame when the package API-view
mode is enabled. Passing ``api=False`` gives a plain pandas
``DataFrame`` for ordinary scripts.

Missing Sections
----------------

``list_missing_sections`` checks whether each station has required data
sections. The most common checks are ``"mt"`` for impedance and
``"tipper"`` for transfer-function data.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import ensure_sites, list_missing_sections

   survey = ensure_sites("data/AMT/WILLY_DATA/L18PLT", strict=True)

   missing = list_missing_sections(
       survey,
       require=("mt", "tipper"),
   )

   for station, sections in missing.items():
       print(f"{station}: missing {', '.join(sections)}")

This function uses the same internal extraction helpers as the plotting
tools. That matters because real EDI/Site objects can expose placeholder
attributes even when a section was not actually parsed.

Check Tipper Availability Explicitly
------------------------------------

AMT/CSAMT lines often have no tipper. MT surveys often do. Make that
difference explicit before writing code that assumes tipper exists.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import list_missing_sections

   amt_missing = list_missing_sections(
       "data/AMT/WILLY_DATA/L18PLT",
       require=("tipper",),
   )
   mt_missing = list_missing_sections(
       "data/MT/kap03lmt_edis",
       require=("tipper",),
   )

   print(f"L18PLT stations missing tipper: {len(amt_missing)}")
   print(f"KAP03 stations missing tipper: {len(mt_missing)}")

If every station is missing tipper, that is not necessarily a failure.
It simply means you should stay with impedance-based inspection and use
the transfer-function tools only on surveys that contain vertical-field
data.

Frequency Coverage Tables
-------------------------

``frequency_coverage`` has three modes.

.. code-block:: python
   :linenos:

   import numpy as np

   from pycsamt.emtools import ensure_sites, frequency_coverage

   survey = ensure_sites("data/MT/kap03lmt_edis", strict=True)

   per_site = frequency_coverage(survey, mode="per-site")
   union = frequency_coverage(survey, mode="union")
   intersection = frequency_coverage(survey, mode="intersection")

   print(f"stations: {len(per_site)}")
   print(f"union frequency count: {union.size}")
   print(f"common frequency count: {intersection.size}")

   for station, freq in per_site.items():
       missing_from_union = np.setdiff1d(union, freq)
       if missing_from_union.size:
           print(station, "missing", missing_from_union.size, "samples")

Use ``mode="per-site"`` when you need station names. Use ``"union"``
to know the full survey frequency grid. Use ``"intersection"`` to know
which frequencies are shared by every station.

Plot Frequency Coverage
-----------------------

``plot_coverage`` converts the frequency dictionary into a station by
frequency presence mask.

.. code-block:: python
   :linenos:

   import matplotlib.pyplot as plt

   from pycsamt.emtools import ensure_sites, plot_coverage

   survey = ensure_sites("data/MT/kap03lmt_edis", strict=True)

   fig, ax = plt.subplots(figsize=(8.0, 4.5))
   plot_coverage(
       survey,
       axis="period",
       ax=ax,
   )
   ax.set_title("KAP03 frequency coverage")

   fig.tight_layout()
   fig.savefig("kap03_frequency_coverage.png", dpi=200)
   plt.close(fig)

The colour value is presence, not data quality. A fully covered cell
only means the sample exists. Use QC, error, confidence, and frequency
editing tools to decide whether the sample is reliable.

Quick Rho And Phase Curves
--------------------------

``plot_rhoa_phi`` plots apparent resistivity and phase for one or more
stations. It accepts components such as ``"xy"``, ``"yx"``, ``"xx"``,
and ``"yy"`` when those columns exist in the station dataframe.

.. code-block:: python
   :linenos:

   import matplotlib.pyplot as plt

   from pycsamt.emtools import ensure_sites, plot_rhoa_phi

   survey = ensure_sites("data/AMT/WILLY_DATA/L18PLT", strict=True)

   station_names = survey.stations[:4]
   subset = [survey.get_site(name) for name in station_names]

   ax_rho, ax_phase = plot_rhoa_phi(
       subset,
       components=("xy", "yx"),
       axis="period",
       errorbar=True,
       figsize=(8.0, 6.0),
   )

   ax_rho.figure.savefig("l18plt_rho_phase_subset.png", dpi=200)
   plt.close(ax_rho.figure)

Do not plot every station at once unless the survey is tiny. The
function will draw the data, but the legend can become unreadable. Use
small station subsets for first inspection, then switch to
``plot_station_response`` for a station-level deep view.

Pseudo-Sections
---------------

``pseudosection`` creates a period by station image from a dataframe
quantity such as ``"rho_xy"``, ``"rho_yx"``, ``"phi_xy"``, or
``"phi_yx"``. Values are pivoted by station and period, with median
aggregation for duplicate cells.

.. code-block:: python
   :linenos:

   import matplotlib.pyplot as plt

   from pycsamt.emtools import ensure_sites, pseudosection

   survey = ensure_sites("data/AMT/WILLY_DATA/L18PLT", strict=True)

   fig, ax = plt.subplots(figsize=(10.0, 4.8))
   pseudosection(
       survey,
       quantity="rho_xy",
       period_range=(1e-4, 1.0),
       ax=ax,
       topo=False,
   )
   ax.set_title("L18PLT rho_xy pseudo-section")

   fig.tight_layout()
   fig.savefig("l18plt_rho_xy_pseudosection.png", dpi=200)
   plt.close(fig)

The x-axis is station order. The y-axis is period. Short periods are
drawn at the top because the image uses the common MT pseudo-section
convention: shallow-sensitive samples above deeper-sensitive samples.

Control The Pseudo-Section Scale
--------------------------------

Use fixed ``vmin`` and ``vmax`` when comparing two lines. Otherwise a
line with a narrow value range can look as dramatic as a line with a
much stronger anomaly.

.. code-block:: python
   :linenos:

   import matplotlib.pyplot as plt

   from pycsamt.emtools import ensure_sites, pseudosection

   line18 = ensure_sites("data/AMT/WILLY_DATA/L18PLT", strict=True)
   line22 = ensure_sites("data/AMT/WILLY_DATA/L22PLT", strict=True)

   fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.0), sharey=True)

   pseudosection(line18, quantity="rho_xy", vmin=10.0, vmax=5000.0, ax=axes[0])
   axes[0].set_title("L18PLT")

   pseudosection(line22, quantity="rho_xy", vmin=10.0, vmax=5000.0, ax=axes[1])
   axes[1].set_title("L22PLT")

   fig.tight_layout()
   fig.savefig("rho_xy_line_comparison.png", dpi=200)
   plt.close(fig)

If topography is configured globally, ``pseudosection`` can draw an
optional topography strip. Pass ``topo=False`` when you want a compact
data-only panel.

Tipper Components
-----------------

``plot_tipper_components`` draws real and imaginary parts of ``Tx`` and
``Ty`` versus period or frequency. Use it only after confirming that
the survey actually contains tipper data.

.. code-block:: python
   :linenos:

   import matplotlib.pyplot as plt

   from pycsamt.emtools import ensure_sites, plot_tipper_components

   survey = ensure_sites("data/MT/kap03lmt_edis", strict=True)

   station_names = ["kap103", "kap121", "kap142", "kap151"]
   subset = [survey.get_site(name) for name in station_names]

   fig, ax = plt.subplots(figsize=(8.5, 4.8))
   plot_tipper_components(
       subset,
       kind=("real", "imag"),
       axis="period",
       ax=ax,
   )
   ax.set_title("KAP03 selected tipper components")

   fig.tight_layout()
   fig.savefig("kap03_tipper_components.png", dpi=200)
   plt.close(fig)

The horizontal zero line is important. Sign changes, isolated spikes,
or one station separating strongly from the others are good reasons to
inspect induction arrows, tipper hodograms, and station metadata.

Full Station Response
---------------------

``plot_station_response`` is the richest first-look figure. It shows
apparent resistivity, phase, and, when available, the four tipper
sub-panels for one station.

.. code-block:: python
   :linenos:

   import matplotlib.pyplot as plt

   from pycsamt.emtools import ensure_sites, plot_station_response

   survey = ensure_sites("data/MT/kap03lmt_edis", strict=True)

   fig = plot_station_response(
       survey,
       station="kap151",
       components=("xx", "xy", "yx", "yy"),
       period_range=(1e-2, 2e4),
       show_tipper=True,
       show_error_bars=True,
       rho_lim=None,
       phase_lim=None,
       tipper_lim=(-2.5, 2.5),
       title="kap151 first-look response",
   )

   fig.savefig("kap151_station_response.png", dpi=200)
   plt.close(fig)

The first row is apparent resistivity on log-log axes. The second row
is phase on a log-period x-axis. The optional third row shows
``Re(Tx)``, ``Im(Tx)``, ``Re(Ty)``, and ``Im(Ty)``. If no tipper exists
or ``show_tipper=False``, the figure uses only the impedance rows.

Overlay A Model Response
------------------------

When ``sites_model`` is supplied, the station response overlays a
second dataset as dashed curves. If observed and model resistivity are
both available, the function appends an RMS value to component titles.
The RMS is computed in ``log10(rho)`` space.

.. code-block:: python
   :linenos:

   import matplotlib.pyplot as plt

   from pycsamt.emtools import ensure_sites, plot_station_response, smooth_mavg

   observed = ensure_sites("data/MT/kap03lmt_edis", strict=True)

   # For demonstration only: a smoothed copy behaves like a model response.
   # In production, pass forward-model or inversion-response EDI data here.
   model_like = smooth_mavg(observed, k=5)

   fig = plot_station_response(
       observed,
       station="kap151",
       sites_model=model_like,
       components=("xy", "yx"),
       period_range=(1e-2, 2e4),
       show_rms=True,
       title="kap151 observed vs model-like response",
   )

   fig.savefig("kap151_response_with_model_overlay.png", dpi=200)
   plt.close(fig)

Use this view after inversion or forward modelling to check whether the
model misses a whole component, a period band, or only local points. A
single RMS number is useful, but the curve shape tells you why the RMS
is high or low.

Build A First-Look Report Bundle
--------------------------------

The following script writes a compact first-look bundle for a survey:
summary table, missing-section table, frequency coverage, rho/phase
curves for a few stations, one pseudo-section, and one station
response.

.. code-block:: python
   :linenos:

   from pathlib import Path

   import matplotlib.pyplot as plt
   import pandas as pd

   from pycsamt.emtools import (
       ensure_sites,
       list_missing_sections,
       plot_coverage,
       plot_rhoa_phi,
       plot_station_response,
       pseudosection,
       sites_summary,
   )

   out = Path("inspect_report_l18plt")
   out.mkdir(parents=True, exist_ok=True)

   survey = ensure_sites("data/AMT/WILLY_DATA/L18PLT", strict=True)

   summary = sites_summary(survey, api=False)
   summary.to_csv(out / "sites_summary.csv", index=False)

   missing = list_missing_sections(survey, require=("mt", "tipper"))
   missing_rows = [
       {"station": station, "missing": ",".join(sections)}
       for station, sections in missing.items()
   ]
   pd.DataFrame(missing_rows).to_csv(out / "missing_sections.csv", index=False)

   fig, ax = plt.subplots(figsize=(8.0, 4.5))
   plot_coverage(survey, ax=ax)
   fig.tight_layout()
   fig.savefig(out / "frequency_coverage.png", dpi=200)
   plt.close(fig)

   subset_names = list(summary["station"].head(4))
   subset = [survey.get_site(name) for name in subset_names]

   ax_rho, ax_phase = plot_rhoa_phi(subset, components=("xy", "yx"))
   ax_rho.figure.savefig(out / "rho_phase_subset.png", dpi=200)
   plt.close(ax_rho.figure)

   fig, ax = plt.subplots(figsize=(10.0, 4.8))
   pseudosection(survey, quantity="rho_xy", ax=ax, topo=False)
   fig.tight_layout()
   fig.savefig(out / "rho_xy_pseudosection.png", dpi=200)
   plt.close(fig)

   first_station = str(summary["station"].iloc[0])
   fig = plot_station_response(
       survey,
       station=first_station,
       components=("xy", "yx"),
       show_tipper=False,
   )
   fig.savefig(out / f"station_response_{first_station}.png", dpi=200)
   plt.close(fig)

Reading The Inspection Results
------------------------------

Treat these outputs as a triage board:

``n_freq`` differs between stations
    Align or edit the frequency grid before making survey-wide
    pseudo-sections or station statistics that assume common samples.

All stations are missing tipper
    Skip tipper diagnostics for that survey. This can be normal for
    AMT/CSAMT lines.

Only some stations are missing tipper
    Keep those stations out of tipper maps or split the analysis into
    tipper-capable and impedance-only subsets.

Rho/phase curves are wildly separated
    Check station metadata, static shift, data quality, and whether a
    few stations dominate the line.

Pseudo-section anomalies appear only at one frequency
    Inspect frequency confidence and errors before interpreting that
    feature geologically.

Station response shows diagonal terms comparable to off-diagonals
    Follow up with impedance, tensor, dimensionality, and strike tools.

Worked Example
--------------

The gallery example below uses the bundled AMT and MT datasets to show
the same first-look workflow end to end.

.. literalinclude:: ../../../examples/emtools/plot_inspect.py
   :language: python
   :linenos:
