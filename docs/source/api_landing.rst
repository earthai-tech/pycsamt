.. _api_landing:

API
===

Two entry points cover everything programmatic in pyCSAMT: the guide to
``pycsamt.api`` — the package front door of readers, tabular views, and
the ``configure_*`` / ``reset_*`` runtime settings — and the
auto-generated reference for every subpackage.

.. grid:: 1 1 2 2
   :gutter: 3
   :class-container: cta-tiles

   .. grid-item-card:: API configuration
      :link: api_guide/index
      :link-type: doc
      :img-top: _static/icons/developer-notes-icon-tools.svg
      :class-card: pycsamt-card sd-text-center

      The ``pycsamt.api`` front door: readers, APIFrame views, and the
      ``configure_*`` families that drive every pipeline stage.

   .. grid-item-card:: API reference
      :link: api/index
      :link-type: doc
      :img-top: _static/icons/api-reference-icon-braces.svg
      :class-card: pycsamt-card sd-text-center

      Full autodoc-generated reference for every public module, class,
      and function, organised by subpackage.

.. toctree::
   :maxdepth: 1
   :hidden:

   api_guide/index
   api/index
