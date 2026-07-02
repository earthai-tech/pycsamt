.. pyCSAMT documentation master file

.. container:: pycsamt-hero

   .. image:: _static/pycsamt-v2-symbol-logo.svg
      :alt: pyCSAMT
      :class: pycsamt-hero-logo

   .. rst-class:: pycsamt-kicker

      CSAMT / AMT / MT processing, modelling, and inversion

   .. rst-class:: pycsamt-hero-title

      Scientific electromagnetic workflows in Python

   .. rst-class:: pycsamt-hero-subtitle

      pyCSAMT v2 brings EDI loading, quality control, static-shift
      correction, phase-tensor analysis, forward modelling, AI-assisted
      inversion, mapping, and reports into one documented toolkit.

   .. grid:: 1 1 3 3
      :gutter: 2
      :class-container: pycsamt-hero-actions

      .. grid-item-card:: Install
         :link: getting_started/installation
         :link-type: doc

         Set up pyCSAMT and its optional scientific backends.

      .. grid-item-card:: First Survey
         :link: getting_started/first_survey
         :link-type: doc

         Load a survey line and inspect usable station data.

      .. grid-item-card:: API
         :link: api/index
         :link-type: doc

         Jump straight into the package reference.

.. rst-class:: pycsamt-status-note

   pyCSAMT v2 is under active development. APIs may change before the
   stable release.

.. rst-class:: pycsamt-section-heading

   Start Where You Work

.. grid:: 1 1 2 3
   :gutter: 3

   .. grid-item-card:: Getting started
      :link: getting_started/index
      :link-type: doc
      :img-top: _static/icons/getting-started-icon-rocket-takeoff.svg
      :class-card: pycsamt-card

      Installation, data formats, and the first survey workflow.

   .. grid-item-card:: User guide
      :link: user_guide/index
      :link-type: doc
      :img-top: _static/icons/user-guide-icon-book.svg
      :class-card: pycsamt-card

      Concepts and workflows for processing, inversion, and interpretation.

   .. grid-item-card:: Tutorials
      :link: tutorials/index
      :link-type: doc
      :img-top: _static/icons/quickstart-examples.svg
      :class-card: pycsamt-card

      End-to-end recipes for common survey tasks.

   .. grid-item-card:: Examples
      :link: examples/index
      :link-type: doc
      :img-top: _static/icons/gallery-icon-images.svg
      :class-card: pycsamt-card

      Runnable examples and gallery scripts.

   .. grid-item-card:: API reference
      :link: api/index
      :link-type: doc
      :img-top: _static/icons/api-reference-icon-braces.svg
      :class-card: pycsamt-card

      Complete module and class reference.

   .. grid-item-card:: CLI reference
      :link: cli/index
      :link-type: doc
      :img-top: _static/icons/cli-icon-terminal.svg
      :class-card: pycsamt-card

      Command-line workflows for survey processing and inversion.

   .. grid-item-card:: AI agents
      :link: agents/index
      :link-type: doc
      :img-top: _static/icons/applications.svg
      :class-card: pycsamt-card

      AI-assisted workflow automation and agent orchestration.

   .. grid-item-card:: Pipeline system
      :link: pipeline/index
      :link-type: doc
      :img-top: _static/icons/developer-notes-icon-tools.svg
      :class-card: pycsamt-card

      Reproducible workflow presets, steps, outputs, and reports.

   .. grid-item-card:: Forward modelling
      :link: forward/index
      :link-type: doc
      :img-top: _static/icons/scientific-foundations-icon-beaker.svg
      :class-card: pycsamt-card

      Synthetic models, forward responses, noise, datasets, and diagnostics.

   .. grid-item-card:: Model integrations
      :link: models/index
      :link-type: doc
      :img-top: _static/icons/bounds-vs-ridge-summary.svg
      :class-card: pycsamt-card

      Occam2D, ModEM, MARE2DEM, native files, and external runners.

   .. grid-item-card:: Site tools
      :link: site/index
      :link-type: doc
      :img-top: _static/icons/external-validation.svg
      :class-card: pycsamt-card

      Station containers, coordinate handling, selection, editing, export,
      and site reports.

   .. grid-item-card:: Map tools
      :link: map/index
      :link-type: doc
      :img-top: _static/icons/glossary-icon-book.svg
      :class-card: pycsamt-card

      Code-first station maps, pseudosections, 3-D quick-look views,
      overlays, multi-line loading, and export.

   .. grid-item-card:: Applications
      :link: applications/index
      :link-type: doc
      :img-top: _static/icons/applications.svg
      :class-card: pycsamt-card

      Desktop GUI, web app, and Agent Master interfaces for interactive work.

----

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Getting started

   getting_started/index
   installation
   quickstart
   api_view

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: User guide

   user_guide/index

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Tutorials and examples

   tutorials/index
   examples/index

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: API reference

   api/index

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Interfaces

   cli/index
   agents/index
   pipeline/index
   forward/index
   models/index
   site/index
   map/index
   applications/index

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Scientific background

   theory/index

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Development

   development/index
   release_notes/index
   changelog
   contributing
   references
