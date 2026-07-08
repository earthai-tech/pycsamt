.. _emtools_diag:

Polar Uncertainty Diagnostics
=============================

``pycsamt.emtools.diag`` evaluates predicted uncertainty intervals
against observed CSAMT/AMT apparent resistivity. It adapts the
``k-diagram`` polar uncertainty ideas to electromagnetic soundings:
coverage, relative interval width, and relative prediction error.

This module is different from most ``emtools`` pages. It cannot work
from observed EDI data alone. You must provide a prediction to evaluate:

- lower and upper apparent-resistivity bounds, ``q_lo`` and ``q_hi``;
- optionally, a point prediction ``model_rho`` for relative-error plots.

Full callable signatures live in the :doc:`API reference <../../api/emtools>`.
This user-guide page focuses on the workflow, data shapes, returned
tables, plots, and interpretation.

What The Diagnostics Measure
----------------------------

For each station and frequency, the module computes observed apparent
resistivity from one off-diagonal impedance component:

.. math::

   \rho_{a,obs} = 0.2 {|Z_{pq}|^2 \over f}

where ``pq`` is ``xy`` by default and ``f`` is frequency in hertz. This
is the same practical-unit EDI convention used by the other
``emtools`` resistivity diagnostics.

Given predicted bounds :math:`[L_j, U_j]`, coverage is:

.. math::

   c_j =
   \begin{cases}
   1, & L_j \le \rho_{a,obs,j} \le U_j \\
   0, & \text{otherwise}
   \end{cases}

The empirical coverage of a station is the mean of those binary values.
For a nominal 90 percent interval, a station with empirical coverage
above ``0.9`` is flagged as calibrated by the default rule.

The module also reports relative interval width:

.. math::

   width_j =
   100 {U_j - L_j \over \rho_{a,obs,j}}

and relative point-prediction error:

.. math::

   error_j =
   100 {\rho_{a,pred,j} - \rho_{a,obs,j} \over \rho_{a,obs,j}}

Coverage tells you whether observations fall inside predicted
intervals. Width tells you whether that coverage was useful or merely
overly cautious. Error tells you where a point prediction systematically
over- or under-predicts the observed sounding.

Inputs You Must Provide
-----------------------

The observed data input is the usual ``emtools`` ``sites`` argument:

- a directory containing EDI files,
- one EDI-like object,
- a ``Sites`` container,
- an iterable of site-like objects.

The prediction inputs can be shaped in three ways:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Input shape
     - Meaning
   * - scalar
     - Broadcast one value to every station and frequency.
   * - one array
     - Reuse the same per-frequency array for each station.
   * - ``dict[str, array]``
     - Use station-specific arrays keyed by station name.

For real work, the dictionary form is usually best. Each array must be
aligned with that station's frequency array. When a station key is
missing from the dictionary, that station is skipped.

Pure Coverage Score
-------------------

``coverage_score`` is the pure arithmetic helper. It does not load EDI
files. Use it when you already have observed values and interval bounds.

.. code-block:: python
   :linenos:

   import numpy as np

   from pycsamt.emtools.diag import coverage_score

   rho_obs = np.array([98.0, 105.0, 87.0, 130.0, 112.0])
   q_lo = np.array([90.0, 95.0, 90.0, 100.0, 100.0])
   q_hi = np.array([110.0, 115.0, 100.0, 120.0, 125.0])

   score = coverage_score(rho_obs, q_lo, q_hi)
   print(f"empirical coverage = {score:.2f}")

Line 9 computes the fraction of observations that fall inside their
interval. In this toy example, values below ``q_lo`` or above ``q_hi``
count as misses.

Building Example Bounds
-----------------------

The rest of the workflow needs prediction intervals. The example below
builds a simple baseline from real L18PLT observations: a rolling median
in log-resistivity space becomes the center line, and the interval width
grows toward longer periods.

This is not a forecasting model. It is a transparent way to demonstrate
the diagnostics using real observed EDI data.

.. code-block:: python
   :linenos:

   import numpy as np

   from pycsamt.emtools.diag import rho_coverage

   survey = "data/AMT/WILLY_DATA/L18PLT"

   raw = rho_coverage(
       survey,
       q_lo=0.0,
       q_hi=np.inf,
       rho_comp="xy",
   )

   q_lo = {}
   q_hi = {}
   model = {}

   log_period = np.log10(raw["period_s"])
   p_min = log_period.min()
   p_max = log_period.max()

   for station, group in raw.groupby("station", sort=False):
       group = group.reset_index(drop=True)
       smooth = (
           np.log10(group["rho_obs"])
           .rolling(5, center=True, min_periods=1)
           .median()
       )
       center = 10.0 ** smooth.to_numpy()
       t = (np.log10(group["period_s"]) - p_min) / (p_max - p_min + 1e-12)
       half_width = 0.15 + 0.30 * t

       q_lo[station] = center * (1.0 - half_width)
       q_hi[station] = center * (1.0 + half_width)
       model[station] = center

Lines 7-12 use ``rho_coverage`` with infinite bounds as a convenient
way to extract observed apparent resistivity. Lines 31-33 build
station-specific lower bounds, upper bounds, and point predictions.

Per-Frequency Coverage
----------------------

Use ``rho_coverage`` when you need one row per station and frequency.

.. code-block:: python
   :linenos:

   from pycsamt.emtools.diag import rho_coverage

   detail = rho_coverage(
       "data/AMT/WILLY_DATA/L18PLT",
       q_lo=q_lo,
       q_hi=q_hi,
       rho_comp="xy",
       recursive=True,
       on_dup="replace",
       strict=False,
       verbose=0,
   )

   print(detail.head())
   detail.to_csv("l18plt_coverage_detail.csv", index=False)

The output columns are:

- ``station``: station name.
- ``freq_hz`` and ``period_s``: frequency and inverse frequency.
- ``rho_obs``: observed apparent resistivity from ``Zxy`` or ``Zyx``.
- ``q_lo`` and ``q_hi``: prediction interval bounds.
- ``covered``: ``True`` when ``q_lo <= rho_obs <= q_hi``.
- ``width_pct``: interval width as a percentage of ``rho_obs``.

The most common mistake is misalignment. If your ``q_lo`` and ``q_hi``
arrays are not ordered the same way as the station's frequency array,
coverage will be meaningless even though the code can still run.

Single-Station Inspection
-------------------------

Before trusting summary statistics, inspect one station's observed
curve against its bounds.

.. code-block:: python
   :linenos:

   import matplotlib.pyplot as plt

   station = "18-001A"
   one = detail.loc[detail["station"] == station].sort_values("period_s")

   fig, ax = plt.subplots(figsize=(7, 4.5))
   ax.fill_between(
       one["period_s"],
       one["q_lo"],
       one["q_hi"],
       color="0.85",
       label="predicted interval",
   )
   ax.loglog(one["period_s"], one["rho_obs"], "o-", label="observed")
   ax.scatter(
       one.loc[~one["covered"], "period_s"],
       one.loc[~one["covered"], "rho_obs"],
       color="red",
       zorder=4,
       label="miss",
   )
   ax.set_xlabel("Period (s)")
   ax.set_ylabel("Apparent resistivity (ohm.m)")
   ax.set_title(f"{station} observed resistivity vs. prediction interval")
   ax.legend()
   fig.tight_layout()

Red points are misses. A few isolated misses may be acceptable. A whole
frequency band outside the interval usually means the model is biased or
the interval width is too narrow in that part of the spectrum.

Per-Station Coverage Table
--------------------------

Use ``coverage_table`` to summarize each station.

.. code-block:: python
   :linenos:

   from pycsamt.emtools.diag import coverage_table

   table = coverage_table(
       "data/AMT/WILLY_DATA/L18PLT",
       q_lo=q_lo,
       q_hi=q_hi,
       rho_comp="xy",
       nominal=0.9,
   )

   ranked = table.sort_values("empirical_cov")
   print(ranked.head(10))

The output columns are:

- ``station``: station name.
- ``n_freq``: number of evaluated frequencies.
- ``empirical_cov``: fraction of covered frequencies.
- ``mean_width_pct``: mean interval width as percent of observed
  resistivity.
- ``calibrated_flag``: ``True`` when ``empirical_cov >= nominal``.

A station can be calibrated because the model is genuinely good, or
because the intervals are very wide. Always read ``empirical_cov`` with
``mean_width_pct``.

Coverage Visualization
----------------------

``plot_polar_coverage`` maps frequency to polar angle and observed
resistivity to radius. Green points are covered. Red points are misses.
Thin radial segments show the prediction interval.

.. code-block:: python
   :linenos:

   import matplotlib.pyplot as plt

   from pycsamt.emtools.diag import plot_polar_coverage

   fig, ax = plt.subplots(
       subplot_kw={"projection": "polar"},
       figsize=(7, 7),
   )
   plot_polar_coverage(
       "data/AMT/WILLY_DATA/L18PLT",
       q_lo=q_lo,
       q_hi=q_hi,
       rho_comp="xy",
       n_freq_ticks=8,
       ax=ax,
   )
   fig.tight_layout()

This plot is useful when you want to know whether misses cluster in a
specific part of the frequency band. A red wedge suggests a systematic
frequency-dependent calibration problem. Scattered red points suggest
local noise or station-specific departures.

Width Drift
-----------

``plot_width_drift`` bins relative interval width by frequency band.
It answers a different question from coverage: how expensive was the
coverage in interval width?

.. code-block:: python
   :linenos:

   import matplotlib.pyplot as plt

   from pycsamt.emtools.diag import plot_width_drift

   fig, ax_cart = plt.subplots(figsize=(8, 4))
   plot_width_drift(
       "data/AMT/WILLY_DATA/L18PLT",
       q_lo=q_lo,
       q_hi=q_hi,
       n_bands=8,
       polar=False,
       ax=ax_cart,
   )

   fig2, ax2 = plt.subplots(
       subplot_kw={"projection": "polar"},
       figsize=(6, 6),
   )
   plot_width_drift(
       "data/AMT/WILLY_DATA/L18PLT",
       q_lo=q_lo,
       q_hi=q_hi,
       n_bands=8,
       polar=True,
       ax=ax2,
   )

If widths grow toward lower frequencies, uncertainty is increasing with
longer periods and, approximately, with greater investigation depth. If
widths are huge everywhere, high coverage may not be very informative.

Point-Prediction Error
----------------------

Use ``rho_error_stats`` and ``plot_polar_errors`` when you have a point
prediction, not only interval bounds.

.. code-block:: python
   :linenos:

   from pycsamt.emtools.diag import rho_error_stats, plot_polar_errors

   errors = rho_error_stats(
       "data/AMT/WILLY_DATA/L18PLT",
       model_rho=model,
       rho_comp="xy",
   )

   print(errors[["station", "freq_hz", "rel_err_pct", "abs_err_pct"]].head())

   ax = plot_polar_errors(
       "data/AMT/WILLY_DATA/L18PLT",
       model_rho=model,
       rho_comp="xy",
       n_bins=18,
   )

The output columns are:

- ``rho_obs``: observed apparent resistivity.
- ``rho_pred``: predicted apparent resistivity.
- ``rel_err_pct``: signed relative error.
- ``abs_err_pct``: absolute relative error.

The polar error plot uses red bars for over-prediction and blue bars for
under-prediction. Bar height is mean absolute relative error in each
frequency sector.

Comparing Calibration Scenarios
-------------------------------

A useful diagnostic exercise is to compare sensible, overconfident, and
underconfident intervals.

.. code-block:: python
   :linenos:

   import pandas as pd

   from pycsamt.emtools.diag import coverage_table

   scenarios = {
       "sensible": 1.0,
       "overconfident": 0.4,
       "underconfident": 3.0,
   }

   rows = []
   for name, multiplier in scenarios.items():
       lo_s = {}
       hi_s = {}
       for station in q_lo:
           center = model[station]
           lo_s[station] = center - (center - q_lo[station]) * multiplier
           hi_s[station] = center + (q_hi[station] - center) * multiplier
       t = coverage_table(
           "data/AMT/WILLY_DATA/L18PLT",
           q_lo=lo_s,
           q_hi=hi_s,
       )
       rows.append(
           {
               "scenario": name,
               "mean_coverage": t["empirical_cov"].mean(),
               "mean_width_pct": t["mean_width_pct"].mean(),
               "n_calibrated": int(t["calibrated_flag"].sum()),
           }
       )

   comparison = pd.DataFrame(rows)
   print(comparison)

Read this table as a trade-off. Overconfident intervals should have low
coverage and narrow width. Underconfident intervals should have high
coverage and wide intervals. A useful model is the one that reaches the
target coverage without making the intervals unnecessarily wide.

Lines 15-19 scale the interval half-width around the same point
prediction. That keeps the comparison fair: only the uncertainty width
changes, not the model center.

Reading The Results
-------------------

Use this interpretation order:

- Check ``coverage_table`` first for station-level calibration.
- Read ``mean_width_pct`` beside ``empirical_cov``.
- Use ``plot_polar_coverage`` to locate frequency bands where misses
  cluster.
- Use ``plot_width_drift`` to see whether the model becomes less
  certain at longer periods.
- Use ``plot_polar_errors`` to identify over- or under-prediction
  sectors when a point prediction is available.

Common Failure Modes
--------------------

Missing prediction keys
   If ``q_lo`` or ``q_hi`` is a dictionary and a station key is absent,
   that station is skipped. Check the station names in your prediction
   output.

Mismatched array length
   Prediction arrays must align with the loaded station frequency array.
   Build bounds from the same station order and frequency order used by
   pyCSAMT.

Intervals with high coverage but huge width
   This is underconfidence. The model is technically calibrated but not
   very useful.

Intervals with narrow width and low coverage
   This is overconfidence. The model misses too many observations for
   the claimed interval.

Scalar bounds
   Scalars are accepted for quick tests, but they are rarely meaningful
   for real apparent-resistivity uncertainty because resistivity varies
   strongly across frequency and station.

Saving A Reproducible Diagnostic Bundle
---------------------------------------

Save the detailed coverage table, station summary, error table, and
figures together.

.. code-block:: python
   :linenos:

   from pathlib import Path

   import matplotlib.pyplot as plt

   from pycsamt.emtools.diag import (
       coverage_table,
       plot_polar_coverage,
       plot_width_drift,
       rho_coverage,
       rho_error_stats,
   )

   survey = "data/AMT/WILLY_DATA/L18PLT"
   out = Path("outputs/diag_l18plt")
   out.mkdir(parents=True, exist_ok=True)

   detail = rho_coverage(survey, q_lo=q_lo, q_hi=q_hi)
   table = coverage_table(survey, q_lo=q_lo, q_hi=q_hi)
   errors = rho_error_stats(survey, model_rho=model)

   detail.to_csv(out / "coverage_detail.csv", index=False)
   table.to_csv(out / "coverage_table.csv", index=False)
   errors.to_csv(out / "relative_errors.csv", index=False)

   fig1, ax1 = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(7, 7))
   plot_polar_coverage(survey, q_lo=q_lo, q_hi=q_hi, ax=ax1)
   fig1.savefig(out / "polar_coverage.png", dpi=200)

   fig2, ax2 = plt.subplots(figsize=(8, 4))
   plot_width_drift(survey, q_lo=q_lo, q_hi=q_hi, ax=ax2)
   fig2.savefig(out / "width_drift.png", dpi=200)

Worked Example
--------------

The gallery example uses **L18PLT** from ``data/AMT/WILLY_DATA/`` and
builds synthetic prediction intervals around real observed apparent
resistivity. It demonstrates one-station inspection, per-station
coverage ranking, polar coverage, width drift, polar errors, and a
calibration scenario comparison.

Open the rendered example here:
:ref:`sphx_glr_examples_emtools_plot_diag.py`.

The source is included below so the page remains useful from the user
guide as well as from the Sphinx-Gallery page.

.. literalinclude:: ../../../examples/emtools/plot_diag.py
   :language: python
   :linenos:
