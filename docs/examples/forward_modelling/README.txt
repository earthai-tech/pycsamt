.. _forward_modelling:

Forward modelling
-----------------

Synthetic **forward modelling** with :mod:`pycsamt.forward` — computing
the electromagnetic response of a resistivity model you specify, with no
field data required. These examples build layered and gridded earth
models from scratch and run the 1-D, 2-D, and quasi-3-D solvers, which is
the natural way to design a survey, build intuition for what an anomaly
looks like in the data, or generate training targets for the
:ref:`AI-inversion <user_guide_ai_inversion>` workflows.

The gallery reads as a workflow, simplest first:

* **1-D layered models** — resistivity–depth profiles and geology priors;
* **1-D soundings** — MT, CSAMT, and TEM responses of those models;
* **2-D forward models** — finite-difference TE/TM pseudo-sections over
  conductive and resistive targets;
* **Plot styles and axis control** — the same figures in publication and
  dark styles, and switching the period axis;
* **Quasi-3-D forward modelling** — map views, pseudo-sections, and full
  impedance-tensor panels over a 3-D block anomaly.

Every model here is defined in the script itself, so each example is
self-contained and runnable on its own. See the :doc:`Forward modelling
guide </forward/index>` for the narrative reference.
