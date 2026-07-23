.. _user_guide_prerequisites:

Prerequisites
=============

The user guide is the practical companion to the pyCSAMT science API,
applications, and command-line tools. It covers the full life of a survey:
turning raw instrument output or field-hardware exports into EDI, loading and
inspecting that data, running quality control and processing, correcting for
terrain and static shift, building maps and pseudosections, choosing an
inversion path, interpreting the result, and automating any of the above with
agents.

This section is organized as an operational workspace rather than a single
narrative. Each page documents one package and stands on its own: open it
directly when you already know what you need. The reading paths below exist
for the more common case, where the question is not "what does this function
do" but "where do I even start."

Reading Paths
-------------

Different users arrive here with different goals. Pick the path that matches
yours; they share pages on purpose, so do not be surprised to see the same
link twice.

- **New survey users**, starting from a folder of EDI files, should read
  :doc:`data_loading`, then :doc:`emtools/index` for quality control, then
  :doc:`map/index` to see stations and profiles on a map.
- **Users with non-EDI instrument output** (AVG, J-files, raw cross-power
  spectra, time series, or a transient TDEM decay) should start at
  :doc:`transformers` to reach EDI first, then continue as a new survey user
  above.
- **Geometrics/EMI Stratagem users** should read :doc:`stratagem/index`
  end to end; it covers the WinGLink hand-off, coordinate injection, and QC
  steps that are specific to that hardware before the data ever reaches
  :doc:`emtools/index`.
- **IoT and edge-acquisition users** should read :doc:`iot/index` for device
  configuration, edge QC, telemetry, and provenance, then bridge into the
  normal workflow through :doc:`data_loading`.
- **Processing and QC users** should work through :doc:`emtools/index` for
  static shift, frequency editing, noise removal, and dimensionality, pairing
  it with :doc:`site/index` for station-level selection, editing, and export.
- **Terrain and mapping users** should read :doc:`topo/index` before
  :doc:`map/index`; topography has to be extracted and draped before a
  section plot reflects real ground elevation instead of a flat datum.
- **Inversion users** should read :doc:`inversion/index` first -- it is the
  decision point between the classical engines in :doc:`models/index`, the
  learned workflows in :doc:`ai_inversion/index`, and, upstream of both,
  the synthetic-data and sensitivity work in :doc:`forward/index`.
- **Interpretation users**, once an inversion result exists, should read
  :doc:`interpretation/index` for hydrogeophysical context, uncertainty, and
  reporting, using :doc:`topo/index` and :doc:`map/index` for the section and
  overlay figures that go into that interpretation.
- **AI and agent users** who want automation rather than manual API calls
  should read :doc:`agents/index`; it links back into
  :doc:`ai_inversion/index` for the model-level detail behind the AI
  inversion agents specifically.
- **Reproducible / team workflows** should read :doc:`pipeline/index` to
  turn a manual sequence of the pages above into a named, versioned,
  re-runnable recipe.

Section Map
-----------

The table below lists every user guide section in reading order, the
package it documents, and what to expect on arrival.

.. list-table::
   :header-rows: 1
   :widths: 20 22 58

   * - Section
     - Package
     - What it covers
   * - :doc:`data_loading`
     - :mod:`pycsamt.emtools`, :mod:`pycsamt.site`
     - Turning files, directories, or parsed objects into a station-ordered
       ``Sites`` container.
   * - :doc:`transformers`
     - :mod:`pycsamt.transformers`, :mod:`pycsamt.tdem`
     - Converting AVG, J-file, raw spectra EDI, time-series, or TDEM decay
       input into EDI.
   * - :doc:`map/index`
     - :mod:`pycsamt.map`
     - Station maps, profile pseudosections, 3-D quick-look maps, and
       figure export from Python code.
   * - :doc:`topo/index`
     - :mod:`pycsamt.topo`
     - Embedding real station elevation into 2-D sections so plots drape
       over terrain.
   * - :doc:`inversion/index`
     - :mod:`pycsamt.inversion`
     - The decision point between classical, AI, and forward-modelling
       paths; the backend-neutral API and adapters.
   * - :doc:`interpretation/index`
     - --
     - Connecting reviewed inversion products to geology, uncertainty,
       and defensible reporting.
   * - :doc:`iot/index`
     - :mod:`pycsamt.iot`
     - Configuring and auditing the operational layer of an IoT-enabled
       field survey.
   * - :doc:`stratagem/index`
     - :mod:`pycsamt.stratagem`
     - The post-acquisition workflow for Geometrics/EMI Stratagem
       hardware, from raw files to processed EDI.
   * - :doc:`emtools/index`
     - :mod:`pycsamt.emtools`
     - Quality control, frequency editing, static shift, noise removal,
       tensor/strike analysis, and source diagnostics.
   * - :doc:`site/index`
     - :mod:`pycsamt.site`
     - Station-centric selection, editing, normalization, profiles, and
       export/reporting for site collections.
   * - :doc:`forward/index`
     - :mod:`pycsamt.forward`
     - Predicting responses from a known earth model: solvers, synthetic
       datasets, and noise models.
   * - :doc:`pipeline/index`
     - :mod:`pycsamt.pipeline`
     - Reproducible, named processing recipes connecting steps,
       configuration files, and CLI execution.

Two sections sit outside this table because they are not linear reading, but
they are just as central: :doc:`models/index` and :doc:`ai_inversion/index`
are the classical and learned engine layers nested under
:doc:`inversion/index`, and :doc:`agents/index` is documented as its own
top-level guide alongside this one, not a child page of it.

Relationship To Other Sections
------------------------------

* :doc:`../getting_started/index` covers installation, environment setup,
  and the first end-to-end run; read it first if pyCSAMT is not installed
  yet or a workflow has not been run once already.
* :doc:`../theory/index` gives the scientific background -- response
  functions, distortion effects, and inversion concepts -- behind the tools
  documented here.
* :doc:`../tutorials/index` gives longer, worked, end-to-end examples that
  combine several of the sections above.
* :doc:`../api/index` is the complete public class and function reference;
  every ``pycsamt.<subpackage>`` link on this page resolves there.
