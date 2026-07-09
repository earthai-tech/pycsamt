.. _user_guide:

User Guide
==========

The user guide is the practical companion to the pyCSAMT science API,
applications, and command-line tools. It explains how users load survey data,
inspect quality, apply processing workflows, build maps, prepare inversions,
run AI-assisted tools, and move from station-level data to reproducible
project outputs.

This section is organized as an operational workspace. Start with the core
workflow pages when you want a guided path through common tasks. Open the
domain guides when you need deeper package-specific documentation for
``emtools``, site containers, forward modelling, external model integrations,
pipelines, agents, or maps.

Reading Paths
-------------

Different users usually arrive here with different goals:

- **New survey users** should begin with :doc:`data_loading`, then continue
  through :doc:`processing`, :doc:`mapping`, and :doc:`interpretation`.
- **Processing and QC users** should pair :doc:`processing` with
  :doc:`emtools/index` and :doc:`site/index`.
- **Inversion users** should read :doc:`inversion`, then move into the
  :doc:`modeling` guides and :doc:`pipeline/index`.
- **AI and agent users** should combine :doc:`ai_inversion` with
  :doc:`agents/index`.
- **Mapping users** can start with :doc:`mapping` for workflow guidance and
  then use :doc:`map/index` for the Python mapping layer.

Core Workflow Guides
--------------------

.. grid:: 1 1 2 3
   :gutter: 3
   :class-container: cta-tiles

   .. grid-item-card:: Data loading
      :link: data_loading
      :link-type: doc
      :img-top: ../_static/icons/user-guide-data-loading.svg
      :class-card: pycsamt-card sd-text-center

      Load EDI directories, explicit station files, API survey views, and
      existing EDI-like objects into the canonical ``Sites`` container.

   .. grid-item-card:: Processing
      :link: processing
      :link-type: doc
      :img-top: ../_static/icons/user-guide-processing.svg
      :class-card: pycsamt-card sd-text-center

      Move from loaded data to QC tables, cleaning, correction, static shift,
      tensor checks, and reproducible processing outputs.

   .. grid-item-card:: Mapping
      :link: mapping
      :link-type: doc
      :img-top: ../_static/icons/user-guide-mapping.svg
      :class-card: pycsamt-card sd-text-center

      Build station maps, survey-line views, pseudosections, and visual checks
      that make field data easier to inspect.

   .. grid-item-card:: Inversion
      :link: inversion
      :link-type: doc
      :img-top: ../_static/icons/user-guide-inversion.svg
      :class-card: pycsamt-card sd-text-center

      Prepare data, model settings, exports, and review steps for 1-D, 2-D,
      and external-engine inversion workflows.

   .. grid-item-card:: Interpretation
      :link: interpretation
      :link-type: doc
      :img-top: ../_static/icons/user-guide-interpretation.svg
      :class-card: pycsamt-card sd-text-center

      Connect processed responses, plots, inversion products, and geologic
      reasoning into reviewable interpretation workflows.

   .. grid-item-card:: IoT workflows
      :link: iot
      :link-type: doc
      :img-top: ../_static/icons/iot.svg
      :class-card: pycsamt-card sd-text-center

      Work with sensor-oriented or streaming workflows when survey data arrive
      from connected acquisition and monitoring contexts.

Domain And Package Guides
-------------------------

.. grid:: 1 1 2 3
   :gutter: 3
   :class-container: cta-tiles

   .. grid-item-card:: EM tools
      :link: emtools/index
      :link-type: doc
      :img-top: ../_static/icons/user-guide-emtools.svg
      :class-card: pycsamt-card sd-text-center

      Detailed narrative pages for ``pycsamt.emtools`` modules, including QC,
      static shift, tensor, strike, frequency, and plotting tools.

   .. grid-item-card:: Site tools
      :link: site/index
      :link-type: doc
      :img-top: ../_static/icons/user-guide-site.svg
      :class-card: pycsamt-card sd-text-center

      Understand station wrappers, ``Sites`` containers, selection, editing,
      coordinates, exports, reports, and recomputation helpers.

   .. grid-item-card:: Modeling
      :link: modeling
      :link-type: doc
      :img-top: ../_static/icons/user-guide-models.svg
      :class-card: pycsamt-card sd-text-center

      Forward modelling, external engine integrations (Occam2D, ModEM,
      MARE2DEM), and AI inversion — the whole model side in one section.

   .. grid-item-card:: Pipeline system
      :link: pipeline/index
      :link-type: doc
      :img-top: ../_static/icons/user-guide-pipeline.svg
      :class-card: pycsamt-card sd-text-center

      Run preset and configured processing chains, compare outputs, preserve
      manifests, and move between Python and CLI execution.

   .. grid-item-card:: Map tools
      :link: map/index
      :link-type: doc
      :img-top: ../_static/icons/user-guide-mapping.svg
      :class-card: pycsamt-card sd-text-center

      Work directly with the Python mapping layer for stations, profiles,
      overlays, 3-D quick looks, exports, and map APIs.

Reference And Next Steps
------------------------

.. grid:: 1 1 2 2
   :gutter: 3
   :class-container: cta-tiles

   .. grid-item-card:: Tutorials
      :link: ../tutorials/index
      :link-type: doc
      :img-top: ../_static/icons/quickstart-examples.svg
      :class-card: pycsamt-card sd-text-center

      Follow longer task-oriented walkthroughs after choosing the relevant
      user-guide path.

   .. grid-item-card:: Getting started
      :link: ../getting_started/index
      :link-type: doc
      :img-top: ../_static/icons/getting-started-icon-rocket-takeoff.svg
      :class-card: pycsamt-card sd-text-center

      Install pyCSAMT, verify the environment, and complete the first practical
      survey run before diving into advanced workflows.

Section Contents
----------------

.. toctree::
   :maxdepth: 1
   :hidden:

   data_loading
   processing
   mapping
   inversion
   interpretation
   iot

.. toctree::
   :maxdepth: 1
   :hidden:

   emtools/index
   site/index
   modeling
   pipeline/index
   map/index

.. toctree::
   :maxdepth: 1
   :hidden:

   ../tutorials/index
   ../getting_started/index
