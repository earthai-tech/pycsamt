.. _development:

Development
===========

This section documents contribution workflows, documentation rules, testing,
public API policy, and release practices.

.. grid:: 1 1 2 3
   :gutter: 3
   :class-container: cta-tiles

   .. grid-item-card:: Documentation build
      :link: documentation_build
      :link-type: doc
      :img-top: ../_static/icons/doc-build.svg
      :class-card: pycsamt-card sd-text-center

      Build the docs locally: environment, sphinx-gallery execution,
      incremental builds, and the output layout.

   .. grid-item-card:: Docstring style
      :link: docstring_style
      :link-type: doc
      :img-top: ../_static/icons/doc-style.svg
      :class-card: pycsamt-card sd-text-center

      The numpydoc conventions every public function, class, and module
      follows — with templates to copy from.

   .. grid-item-card:: Continuous integration
      :link: ci
      :link-type: doc
      :img-top: ../_static/icons/ci.svg
      :class-card: pycsamt-card sd-text-center

      What runs on every push and pull request: test matrix, linting,
      docs build, and release automation.

   .. grid-item-card:: API policy
      :link: api_policy
      :link-type: doc
      :img-top: ../_static/icons/policy.svg
      :class-card: pycsamt-card sd-text-center

      What counts as public API, deprecation rules, and the stability
      promises across the v2 line.

   .. grid-item-card:: Extending pyCSAMT
      :link: extending
      :link-type: doc
      :img-top: ../_static/icons/api-reference-icon-braces.svg
      :class-card: pycsamt-card sd-text-center

      Inherit from PyCSAMTObject, CoreObject, or MTBase to get repr,
      serialization, MT math, and EDI interop in your own classes.

.. toctree::
   :maxdepth: 1
   :hidden:

   documentation_build
   docstring_style
   ci
   api_policy
   extending

