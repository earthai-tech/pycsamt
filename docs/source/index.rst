.. pyCSAMT documentation master file

.. image:: _static/pycsamt_logo.svg
   :alt: pyCSAMT
   :align: center
   :width: 380px

|

.. rst-class:: lead

   An open-source Python toolkit for audio-frequency magnetotelluric
   data processing, inversion, AI-based interpretation, and geological
   deliverable generation.

----

.. admonition:: Documentation status
   :class: warning

   pyCSAMT v2 is under active development.
   APIs are subject to change before the stable release.
   Use :ref:`genindex` to search existing content.

----

.. grid:: 2
   :gutter: 3

   .. grid-item-card:: Getting started
      :link: getting_started/index
      :link-type: doc

      Installation, data formats, and the first survey workflow.

   .. grid-item-card:: User guide
      :link: user_guide/index
      :link-type: doc

      Concepts and workflows for processing, inversion, and interpretation.

   .. grid-item-card:: Tutorials
      :link: tutorials/index
      :link-type: doc

      End-to-end recipes for common survey tasks.

   .. grid-item-card:: Examples
      :link: examples/index
      :link-type: doc

      Runnable examples and gallery scripts.

   .. grid-item-card:: API reference
      :link: api/index
      :link-type: doc

      Complete module and class reference.

   .. grid-item-card:: CLI reference
      :link: cli/index
      :link-type: doc

      Command-line workflows for survey processing and inversion.

   .. grid-item-card:: AI agents
      :link: agents/index
      :link-type: doc

      AI-assisted workflow automation and agent orchestration.

   .. grid-item-card:: Pipeline system
      :link: pipeline/index
      :link-type: doc

      Reproducible workflow presets, steps, outputs, and reports.

   .. grid-item-card:: Forward modelling
      :link: forward/index
      :link-type: doc

      Synthetic models, forward responses, noise, datasets, and diagnostics.

   .. grid-item-card:: Model integrations
      :link: models/index
      :link-type: doc

      Occam2D, ModEM, MARE2DEM, native files, and external runners.

   .. grid-item-card:: Site tools
      :link: site/index
      :link-type: doc

      Station containers, coordinate handling, selection, editing, export,
      and site reports.

   .. grid-item-card:: Map tools
      :link: map/index
      :link-type: doc

      Code-first station maps, pseudosections, 3-D quick-look views,
      overlays, multi-line loading, and export.

   .. grid-item-card:: Applications
      :link: applications/index
      :link-type: doc

      Desktop and web interfaces for interactive survey work.

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
