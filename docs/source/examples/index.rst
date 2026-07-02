Examples
========

This section collects runnable examples for common pyCSAMT workflows.

The v2 example gallery is being organized around the real workflow APIs:
``ensure_sites`` for loading EDI surveys, ``build_qc_table`` for quality
control, ``estimate_ss_ama`` and ``correct_ss_ama`` for static-shift
analysis, and the agent wrappers for complete processing chains.

Available example sources
-------------------------

Until the gallery pages are expanded, start with the tutorial recipes:

.. grid:: 1 1 2 2
   :gutter: 3

   .. grid-item-card:: Read an EDI survey
      :link: ../tutorials/read_edi_survey
      :link-type: doc

      Load a survey line and inspect station impedance availability.

   .. grid-item-card:: Inspect and QC a survey
      :link: ../tutorials/inspect_and_qc_survey
      :link-type: doc

      Build quality-control tables and review flagged stations.

   .. grid-item-card:: Correct static shift
      :link: ../tutorials/correct_static_shift
      :link-type: doc

      Estimate AMA static-shift factors and apply finite positive factors.

   .. grid-item-card:: Run a pipeline from config
      :link: ../tutorials/run_pipeline_from_config
      :link-type: doc

      Execute reproducible survey-processing presets from configuration.
