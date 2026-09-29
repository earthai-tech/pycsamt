Generated-code validation and repair
====================================

Agent Master validates scripts before saving them. A script with a detected
failure is returned as an unsaved draft. Earlier script files are preserved;
successful publication uses a new filename and a matching
``.validation.json`` report containing the code hash and repair history.
The chat history retains the report for follow-up work.

Each check has a state: ``passed``, ``failed``, ``not_checked`` or
``unverifiable``. Reports distinguish:

* Syntax: parse and compile, without execution.
* Imports: static existence of supported pyCSAMT modules and symbols in source.
* Arguments: supported function/constructor/method signatures, positional and
  keyword arguments. Dynamic signatures, rebinding, inheritance and expanded
  ``*args`` / ``**kwargs`` can remain unverifiable.
* Scientific literals: finite positive layer resistivities/thicknesses with
  matching dimensions, and literal correction-factor mappings.
* Runtime imports: not checked; source inspection does not prove optional
  dependencies are installed or that module imports will succeed.
* Fixture execution and artifacts: not checked unless an isolated executor is
  available. No such executor is configured in this implementation.

External-library imports are explicitly not checked. Runtime array values,
dataframe contents, types and scientific results are not inferred from syntax.
The compatibility ``ok`` field means no detected static errors; it is not a
certificate that the script runs or produces valid science.

Bounded repair
--------------

When local generation fails a supported static or request-constraint check,
Agent Master can send the errors, source signatures, request contract and
current script back to the model. It performs at most two repair calls, each
limited to 1024 output tokens, within the original local request deadline and
a maximum 60-second generation/repair window. Model call limits and cancellation
still apply. Empty, unchanged or late responses stop repair; late responses are
discarded. Every changed candidate is validated again before publication.

Automatic repair is currently local-only. Cloud transports lack the shared
deadline required by this path, so cloud drafts still receive validation but
their reports explicitly state that automatic repair was skipped. Offline
templates and exact directory/DPI edits are validated without model repair.
Context overflow or exhausted time leaves an unsaved draft with its report.

Execution boundary
------------------

Generated scripts are never executed by the validator. Passing
``execute_fixture=True`` requests execution checking but returns
``unverifiable`` when no isolated executor is configured. It does not fall
back to a subprocess. Runtime artifact checks and scientific invariants remain
skipped in that case. Existing registered workflows are a separate execution
path and retain their current behavior.

.. code-block:: python

   from pycsamt.assistant.tools.validation_tools import validate_generated_code

   report = validate_generated_code(
       "from pycsamt.emtools._core import ensure_sites\n"
       "sites = ensure_sites('data/3edis')\n"
   )
   print(report["checks"])

``CodeGenerationAgent.execute`` accepts ``max_repairs`` (clamped to 0–2) and
``execute_fixture``. A saved artifact's full report is returned as
``validation`` and its sidecar path as ``validation_path``. General correctness
and constraint fulfillment still require review.
