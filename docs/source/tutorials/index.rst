.. _tutorials:

Tutorials
=========

Tutorials are practical, end-to-end recipes for common pyCSAMT survey tasks.
They are written in a workflow style: each page states the input assumptions,
shows runnable Python snippets, includes CLI equivalents where available,
describes expected outputs, and points to the related API and user-guide
sections.

The tutorial sequence follows the normal path from raw EDI files to reviewed,
processed data that can be used for modelling or inversion. If you are new to
pyCSAMT v2, read the tutorials in order below; the first two pages teach the
basic survey object and QC workflow, and the later pages build on that
foundation.

.. toctree::
   :maxdepth: 3
   :class: pycsamt-guide-toc

   overview
   read_edi_survey
   inspect_and_qc_survey
   compare_survey_lines_for_qc
   correct_static_shift
   condition_mt_line_with_tipper_and_rotation
   prepare_occam2d_inversion
   ai_inversion_from_corrected_edis
   essential_3d_ai_inversion
   run_pipeline_from_config
