.. _getting_started_overview:

Start here
==========

Getting Started takes you from an unconfigured Python environment to a loaded
and inspected electromagnetic survey. It deliberately stops before correction,
inversion, or interpretation. Those stages require scientific choices that
cannot be made safely until the input format, station inventory, coordinates,
frequency coverage, and available response components are understood.

You can follow the pages in order or use :doc:`quickstart` when pyCSAMT is
already installed and the input is a directory of EDI files.

The recommended path
--------------------

.. list-table::
   :header-rows: 1
   :widths: 12 26 32 30

   * - Step
     - Page
     - Question it answers
     - Continue when
   * - 1
     - :doc:`installation`
     - Can this Python environment import and run pyCSAMT?
     - The package import and CLI checks succeed.
   * - 2
     - :doc:`data_formats`
     - What do the field files contain, and can they be loaded directly?
     - The input family and required conversion route are known.
   * - 3
     - :doc:`configuration`
     - Which station order, figure destination, and plot style should this
       session use?
     - The few settings that affect the first survey are explicit.
   * - 4
     - :doc:`first_survey`
     - Did the intended stations load, and what does the first diagnostic show?
     - Counts, identities, errors, coverage, and unusual patterns are reviewed.
   * - 5
     - :doc:`quickstart`
     - What does the complete minimal script look like?
     - The script reproduces its inventory and diagnostic figure.

The sequence separates different kinds of failure. Installation problems
should not be diagnosed as parser problems, unsupported field formats should
not be passed to an EDI reader, and a readable survey should not be mistaken
for an inversion-ready dataset.

Choose a shorter entry point
----------------------------

Start with :doc:`quickstart` if all of the following are true:

* pyCSAMT is already installed and importable;
* the data are frequency-domain EDI files;
* station metadata and units are known;
* you want a compact load, inventory, normalization, and plot example.

Start with :doc:`data_formats` if the delivery contains AVG/AMTAVG, Jones
J-format, spectral EDI, field time series, TEM/TDEM decay data, Stratagem raw
files, coordinate tables, or inversion files. Renaming an extension does not
convert the underlying measurements.

Start with :doc:`first_survey` if the data already load but the station count,
source-file mapping, frequency coverage, or diagnostic figure needs careful
review.

What Getting Started does not do
--------------------------------

This section does not prescribe filtering thresholds, static-shift factors,
dimensionality assumptions, inversion regularization, or geological meaning.
Those decisions depend on the acquisition, processing history, response
quality, survey geometry, and scientific objective.

After the first inspection, use:

* :doc:`../user_guide/data_loading` for reusable loading and normalization
  rules;
* :doc:`../user_guide/emtools/index` for electromagnetic QC and processing;
* :doc:`../user_guide/site/index` for station selection, coordinates, profiles,
  and metadata;
* :doc:`../user_guide/inversion/index` to choose an inversion path;
* :doc:`../tutorials/index` for complete worked workflows.

Before leaving this section
---------------------------

Confirm that:

* Python imports the expected pyCSAMT version and ``pycsamt --help`` works;
* the input family and original measurement units are known;
* raw inputs remain unchanged and generated outputs use a separate directory;
* the loaded station count matches the field manifest;
* parser failures and duplicate identities are understood;
* station order and coordinate reference information are plausible;
* frequency grids and required response components have been checked;
* the first diagnostic has been interpreted as a screening result, not a
  geological conclusion.

If any item is uncertain, return to the corresponding page before applying a
correction or preparing an inversion.
