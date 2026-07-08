.. _emtools_gradient_imaging:

Gradient-Based Pseudo-Sections
==============================

``pycsamt.emtools.gradient_imaging`` builds CSAMT/AMT apparent
resistivity pseudo-sections from gradients rather than from
resistivity values alone. The goal is boundary emphasis: highlight
where apparent resistivity changes laterally, vertically, or in both
directions at once.

The implementation follows the gradient apparent-resistivity ideas of
Zhang, Farquharson, and Liu (2021). It is useful when a normal
``rho_a`` pseudo-section carries broad background variation that makes
target boundaries difficult to see.

Full callable signatures live in the :doc:`API reference <../../api/emtools>`.
This page explains the quantities, returned tables, plotting workflow,
component choices, and interpretation checks.

Why Use Gradients
-----------------

A plain apparent-resistivity pseudo-section answers:

   "How resistive is this station-frequency sample?"

A gradient pseudo-section answers a different question:

   "Where does apparent resistivity change sharply?"

That difference matters in CSAMT/AMT interpretation. Boundary-like
features can be clearer in derivatives than in raw values, especially
when the background has broad smooth variation.

The Three Quantities
--------------------

The module computes apparent resistivity from impedance first. By
default it uses the determinant-style geometric mean:

.. math::

   \rho_a =
   \sqrt{
   \left(0.2 {|Z_{xy}|^2 \over f}\right)
   \left(0.2 {|Z_{yx}|^2 \over f}\right)
   }

You can also use one component only with ``comp="xy"`` or ``comp="yx"``.

The three gradient products are:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Quantity
     - Meaning
   * - ``spatial``
     - Along-line station-to-station difference,
       :math:`\Delta\rho_a^x`.
   * - ``frequency``
     - Adjacent-frequency difference at each station,
       :math:`\Delta\rho_a^z`.
   * - ``joint``
     - Frequency difference of the spatial gradient,
       :math:`\Delta\rho_a^{zx}`.

The joint gradient is often the most useful image because it responds
where the apparent resistivity changes both laterally and with
frequency. Smooth background variation tends to be reduced.

Station Position And Spacing
----------------------------

The gradient tables need station order and station spacing. pyCSAMT uses
station position metadata when available. If no usable coordinates are
available, it falls back to regular spacing with ``spacing_m=200.0``.

Always report the spacing assumption when using fallback spacing. The
gradient values are apparent-resistivity differences, but the x-axis
and pair spacing still affect interpretation of lateral scale.

Spatial Gradient
----------------

``rho_spatial_gradient`` compares adjacent stations at the same
frequency.

.. math::

   \Delta\rho_a^x(j, f)
   \approx \rho_a(j, f) - \rho_a(j - 1, f)

.. code-block:: python
   :linenos:

   from pycsamt.emtools.gradient_imaging import rho_spatial_gradient

   spatial = rho_spatial_gradient(
       "data/AMT/WILLY_DATA/L18PLT",
       spacing_m=200.0,
       comp="det",
   )

   print(spatial.head())
   spatial.to_csv("l18plt_spatial_gradient.csv", index=False)

The output columns are:

- ``station_a`` and ``station_b``: adjacent station pair.
- ``x_m``: midpoint position of the pair.
- ``dx_m``: station spacing for the pair.
- ``freq_hz`` and ``period_s``: frequency and period.
- ``depth_m``: skin-depth-style apparent depth at the pair.
- ``rho_a_ohmm``: mean apparent resistivity of the station pair.
- ``delta_rho_x``: lateral apparent-resistivity difference.

Positive ``delta_rho_x`` means the right station has higher apparent
resistivity than the left station at that frequency. Negative values
mean the opposite.

Frequency Gradient
------------------

``rho_frequency_gradient`` compares adjacent frequencies at each
station.

.. math::

   \Delta\rho_a^z(j, f_k)
   \approx \rho_a(j, f_k) - \rho_a(j, f_{k-1})

.. code-block:: python
   :linenos:

   from pycsamt.emtools.gradient_imaging import rho_frequency_gradient

   vertical = rho_frequency_gradient(
       "data/AMT/WILLY_DATA/L18PLT",
       spacing_m=200.0,
       comp="det",
   )

   one = vertical.loc[vertical["station"] == "18-001A"].sort_values("period_s")
   print(one[["period_s", "rho_a_ohmm", "delta_rho_z"]].head())

The output columns are:

- ``station``: station name.
- ``x_m``: station position.
- ``freq_hz`` and ``period_s``: upper frequency of the adjacent pair.
- ``depth_m``: apparent depth at the mean resistivity of the pair.
- ``rho_a_ohmm``: mean apparent resistivity of the two frequencies.
- ``delta_rho_z``: frequency-direction apparent-resistivity difference.

This quantity is a proxy for vertical or depth-related changes because
different frequencies sample different investigation depths.

Joint Gradient
--------------

``rho_joint_gradient`` computes the frequency difference of the spatial
gradient.

.. math::

   \Delta\rho_a^{zx}(j, f_k)
   =
   \Delta\rho_a^x(j, f_k)
   -
   \Delta\rho_a^x(j, f_{k-1})

Expanded:

.. math::

   \Delta\rho_a^{zx}(j, f_k)
   =
   [\rho_a(j, f_k) - \rho_a(j-1, f_k)]
   -
   [\rho_a(j, f_{k-1}) - \rho_a(j-1, f_{k-1})]

.. code-block:: python
   :linenos:

   from pycsamt.emtools.gradient_imaging import rho_joint_gradient

   joint = rho_joint_gradient(
       "data/AMT/WILLY_DATA/L18PLT",
       spacing_m=200.0,
       comp="det",
   )

   print(joint.head())
   joint.to_csv("l18plt_joint_gradient.csv", index=False)

The output columns are:

- ``station_a`` and ``station_b``: adjacent station pair.
- ``x_m``: midpoint position of the pair.
- ``dx_m``: station spacing for the pair.
- ``freq_hz`` and ``period_s``: upper frequency of the adjacent
  frequency pair.
- ``depth_m``: apparent depth from the median resistivity of the four
  surrounding station-frequency cells.
- ``delta_rho_zx``: joint frequency-spatial gradient.

Joint gradients are strongest when lateral contrast changes with
frequency. That is why they are useful for boundary delineation.

Plot Gradient Sections
----------------------

``plot_gradient_section`` is the main visualization helper. It accepts
``quantity="spatial"``, ``"frequency"``, or ``"joint"``.

.. code-block:: python
   :linenos:

   import matplotlib.pyplot as plt

   from pycsamt.emtools.gradient_imaging import plot_gradient_section

   survey = "data/AMT/WILLY_DATA/L18PLT"

   fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
   plot_gradient_section(survey, quantity="spatial", ax=axes[0])
   axes[0].set_title("Spatial")

   plot_gradient_section(survey, quantity="frequency", ax=axes[1])
   axes[1].set_title("Frequency")

   plot_gradient_section(survey, quantity="joint", ax=axes[2])
   axes[2].set_title("Joint")

   fig.tight_layout()

The color map is centered at zero by default. Positive and negative
gradients are both meaningful: they indicate opposite directions of
apparent-resistivity change.

Use ``vlim`` when comparing several lines or components so the color
scale is consistent.

.. code-block:: python
   :linenos:

   from pycsamt.emtools.gradient_imaging import plot_gradient_section

   ax = plot_gradient_section(
       "data/AMT/WILLY_DATA/L18PLT",
       quantity="joint",
       comp="det",
       vlim=(-3000.0, 3000.0),
   )

Choosing The Impedance Component
--------------------------------

All gradient functions accept:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - ``comp``
     - Meaning
   * - ``"det"``
     - Geometric mean of ``xy`` and ``yx`` apparent resistivities.
   * - ``"xy"``
     - Use only ``Zxy``.
   * - ``"yx"``
     - Use only ``Zyx``.

The default ``"det"`` is usually more stable because it combines both
off-diagonal modes. Component-specific plots are still important when
the two modes disagree.

.. code-block:: python
   :linenos:

   import pandas as pd

   from pycsamt.emtools.gradient_imaging import rho_joint_gradient

   survey = "data/AMT/WILLY_DATA/L18PLT"

   rows = []
   for comp in ("det", "xy", "yx"):
       table = rho_joint_gradient(survey, comp=comp)
       rows.append(
           {
               "component": comp,
               "std": table["delta_rho_zx"].std(),
               "max_abs": table["delta_rho_zx"].abs().max(),
           }
       )

   component_sensitivity = pd.DataFrame(rows)
   print(component_sensitivity)

Report the component whenever you show a gradient section.

Single-Pair And Single-Station Curves
-------------------------------------

Before interpreting a full pseudo-section, inspect a station pair or
station curve.

.. code-block:: python
   :linenos:

   import matplotlib.pyplot as plt

   from pycsamt.emtools.gradient_imaging import (
       rho_frequency_gradient,
       rho_spatial_gradient,
   )

   survey = "data/AMT/WILLY_DATA/L18PLT"

   spatial = rho_spatial_gradient(survey)
   pair = spatial.loc[spatial["station_a"] == "18-001A"].sort_values("period_s")

   frequency = rho_frequency_gradient(survey)
   station = frequency.loc[frequency["station"] == "18-001A"].sort_values("period_s")

   fig, (ax_pair, ax_station) = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
   ax_pair.semilogx(pair["period_s"], pair["delta_rho_x"], "o-")
   ax_pair.axhline(0.0, color="0.4", linewidth=0.8)
   ax_pair.set_ylabel("Spatial gradient")

   ax_station.semilogx(station["period_s"], station["delta_rho_z"], "o-", color="C2")
   ax_station.axhline(0.0, color="0.4", linewidth=0.8)
   ax_station.set_xlabel("Period (s)")
   ax_station.set_ylabel("Frequency gradient")

   fig.tight_layout()

This makes it easier to tell whether a pseudo-section hotspot comes
from one extreme station pair, one frequency jump, or a coherent region.

Does The Joint Gradient Suppress Background?
--------------------------------------------

One practical check is to compare the spread of the spatial and joint
gradients. A joint gradient that is quieter in background areas should
often have a narrower distribution, while preserving strong localized
values.

.. code-block:: python
   :linenos:

   from pycsamt.emtools.gradient_imaging import (
       rho_joint_gradient,
       rho_spatial_gradient,
   )

   survey = "data/AMT/WILLY_DATA/L18PLT"

   spatial = rho_spatial_gradient(survey)
   joint = rho_joint_gradient(survey)

   spatial_std = spatial["delta_rho_x"].std()
   joint_std = joint["delta_rho_zx"].std()

   print(f"spatial std = {spatial_std:.1f} ohm.m")
   print(f"joint std   = {joint_std:.1f} ohm.m")
   print(f"ratio       = {joint_std / spatial_std:.2f}")

This is not a proof of geological correctness, but it is a useful
sanity check before relying on the joint image.

Compare Neighboring Lines
-------------------------

Use the same ``quantity``, ``comp``, and ``vlim`` when comparing lines.

.. code-block:: python
   :linenos:

   import matplotlib.pyplot as plt

   from pycsamt.emtools.gradient_imaging import plot_gradient_section

   l18 = "data/AMT/WILLY_DATA/L18PLT"
   l22 = "data/AMT/WILLY_DATA/L22PLT"

   fig, (ax18, ax22) = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
   plot_gradient_section(l18, quantity="joint", comp="det", vlim=(-3000, 3000), ax=ax18)
   ax18.set_title("L18PLT")

   plot_gradient_section(l22, quantity="joint", comp="det", vlim=(-3000, 3000), ax=ax22)
   ax22.set_title("L22PLT")

   fig.tight_layout()

Line-to-line similarity is a useful processing sanity check. Differences
can still be geological, but first confirm the same band, component, and
color scale were used.

Reading The Results
-------------------

Use this interpretation order:

- Start with ``comp="det"`` for a stable overview.
- Check ``xy`` and ``yx`` separately when off-diagonal modes disagree.
- Use the spatial gradient for lateral station-to-station contrasts.
- Use the frequency gradient for single-station vertical changes.
- Use the joint gradient for features that are both lateral and
  depth-dependent.
- Prefer coherent station-period regions over isolated extreme cells.
- Compare neighboring lines with the same ``vlim`` before making a
  structural interpretation.

Common Failure Modes
--------------------

Missing or irregular station positions
   The module falls back to ``spacing_m``. Report the assumed spacing
   when coordinates are unavailable.

Component-driven artifacts
   ``xy`` and ``yx`` can behave very differently. Always record
   ``comp`` and inspect determinant vs component-specific results.

Over-reading sign
   Positive and negative gradients indicate direction of change. A
   boundary can appear as adjacent positive and negative lobes.

Noisy frequency rows
   Gradients amplify row-to-row changes. Run frequency QC before
   interpreting subtle gradient features.

Color-scale mismatch
   Automatic color limits can make two lines look more different or
   more similar than they are. Use fixed ``vlim`` for comparisons.

Saving A Reproducible Bundle
----------------------------

Save the three gradient tables and the main pseudo-section.

.. code-block:: python
   :linenos:

   from pathlib import Path

   import matplotlib.pyplot as plt

   from pycsamt.emtools.gradient_imaging import (
       plot_gradient_section,
       rho_frequency_gradient,
       rho_joint_gradient,
       rho_spatial_gradient,
   )

   survey = "data/AMT/WILLY_DATA/L18PLT"
   out = Path("outputs/gradient_l18plt")
   out.mkdir(parents=True, exist_ok=True)

   spatial = rho_spatial_gradient(survey, comp="det")
   frequency = rho_frequency_gradient(survey, comp="det")
   joint = rho_joint_gradient(survey, comp="det")

   spatial.to_csv(out / "spatial_gradient.csv", index=False)
   frequency.to_csv(out / "frequency_gradient.csv", index=False)
   joint.to_csv(out / "joint_gradient.csv", index=False)

   fig, ax = plt.subplots(figsize=(10, 5))
   plot_gradient_section(survey, quantity="joint", comp="det", ax=ax)
   fig.tight_layout()
   fig.savefig(out / "joint_gradient_section.png", dpi=200)

Worked Example
--------------

The gallery example uses **L18PLT** and compares it with neighboring
**L22PLT** from ``data/AMT/WILLY_DATA/``. It demonstrates one station
pair, one station frequency-gradient curve, the joint gradient, all
three pseudo-sections, background-suppression checks, component
sensitivity, and line-to-line comparison.

Open the rendered example here:
:ref:`sphx_glr_examples_emtools_plot_gradient_imaging.py`.

The source is included below so the page remains useful from the user
guide as well as from the Sphinx-Gallery page.

.. literalinclude:: ../../../examples/emtools/plot_gradient_imaging.py
   :language: python
   :linenos:
