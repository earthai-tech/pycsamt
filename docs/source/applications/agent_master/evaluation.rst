.. _agent-master-evaluation:

Evaluating answers and scripts
==============================

Agent Master ships a benchmark that scores what a user actually receives:
the presented answer, any returned script, its validation report, cited
sources, and latency. It complements the retrieval harness in
:mod:`pycsamt.assistant.evals.harness`, which scores routing and retrieval
only.

The benchmark runs every case through the same chat path as the application,
so offline rules, a local model and a cloud provider are compared on
identical requests, fixtures and retrieval settings.

The benchmark
-------------

``pycsamt/assistant/evals/answer_suites/agent_master.jsonl`` holds 48 cases in
eight categories: package questions, custom scripts, follow-up edits,
ambiguous requests, unsupported requests, scientific edge cases, developer
questions, and evidence and reporting. There are 32 **development** cases and
16 **held-out** cases, and the two splits share no scenario family.

Each case keeps its reviewed acceptance ``criteria`` and adds typed
``checks``. Every check refers to one criterion, and every criterion has at
least one check:

.. code-block:: json

   {"id": "CODE03", "split": "development", "fixture": "bundled_3edis",
    "query": "Write a script that estimates static-shift factors ...",
    "criteria": ["Verify estimation API", "Export factors",
                 "Do not apply correction", "Handle an empty factor result"],
    "checks": [
      {"id": "CODE03.6", "type": "code_none", "criterion": 2,
       "patterns": ["correct_ss_ama|apply_ss|StaticShiftAgent|\\.correct\\("],
       "stage": "generation", "constraint": true}]}

Held-out cases are for release comparison only. Do not use them as prompt
examples, repair rules or training data, and run them only after development
changes are frozen.

Outcomes
--------

Each check produces one of four outcomes:

``pass`` / ``fail``
   Decided automatically from the recorded job.
``not_assessed``
   The evidence was absent, for example a validation stage that did not run
   or a missing prior artifact, or the criterion needs a domain expert
   (``human`` checks).
``not_applicable``
   A conditional check whose ``when`` clause did not hold, such as a code
   check when the assistant correctly asked a question instead.

A case **passes** only when no check fails and none is unassessed. A case with
no failure but unassessed checks is **incomplete**; it is never counted as a
success. A runtime failure, such as a local-model timeout, is reported as
``environment_error`` rather than as an assistant failure, and a provider that
could not be run is reported as ``not_run`` with its reason.

Four checks run automatically on every completed case:

``api_claims``
   Every ``pycsamt.…`` reference and ``from pycsamt… import …`` statement in
   the answer or script is resolved statically against the checkout.
   References that do not exist fail; references the user named in the
   request are excluded, so explaining that a requested function does not
   exist is not penalized. Bare, unqualified function names are not checked.
``execution_claims``
   Prose outside code blocks must not claim that work was run, computed,
   plotted or saved unless the job recorded a workflow run or an isolated
   fixture execution.
``validation_honesty``
   Returned code must carry a validation report, the answer must not claim
   blanket validation when any check did not pass, and it must say that the
   script was not executed when it was not.
``artifact_exists``
   A reported saved script must exist on disk.

The pattern checks are screening heuristics, not semantic grading. They catch
omitted requirements and definite overclaims; ``human`` checks carry the
judgement that a regular expression cannot make.

Running the benchmark
---------------------

Use the environment in which Agent Master runs::

    python -m pycsamt.assistant.evals.runner --provider offline --split development --out eval-offline
    python -m pycsamt.assistant.evals.runner --provider ollama --model qwen2.5-coder:1.5b --split development --block-external --out eval-local
    python -m pycsamt.assistant.evals.runner --provider claude --allow-cloud --split development --out eval-claude

``--ids`` selects individual cases, and ``--timeout`` sets the local request
time limit. ``--rescore`` scores an earlier run's ``results.json`` again with
the current checks, without re-running any model.

Each run writes ``results.json`` (recorded jobs and presented answers),
``scored.json`` (per-check outcomes and metrics) and ``summary.txt``.

- **Cloud providers** send paid requests off the machine, so they run only
  with ``--allow-cloud``, even when a credential exists. Agent Master loads
  credentials from environment variables and from ``.env.local`` at the
  repository root. Without the opt-in, every case is recorded as
  ``not run: cloud requests not authorized``; with it but without a
  credential, as ``not run: no credentials``. Nothing is sent in either
  case, and credentials are never written to the results.
- **Local-only evidence.** ``--block-external`` refuses every connection to a
  non-loopback address and records the attempt in
  ``run.external_connection_attempts``. An empty list is the evidence that a
  local run made no cloud requests.
- **Execution.** The runner never executes generated scripts. Execution
  evidence comes only from the validation report's optional isolated
  fixture executor; see :doc:`code_validation`.

Fixtures are materialized per case in a fresh conversation. Data fixtures use
a disposable copy of ``data/3edis``, follow-up cases receive a reviewed prior
script as the previous assistant turn, and the injected-instruction case
appends its synthetic text to every retrieved context. Conditions that exist
only as a description in the request, such as a missing solver, are recorded
as ``prompt_only``.

Metrics
-------

``scored.json`` reports:

- case outcomes overall, by category and by split;
- the share of cases with no automated failure, and the share fully passed;
- constraint fulfilment over checks marked ``constraint``;
- unsupported API references and execution-claim violations;
- validation-state counts, fixture execution and scientific-check states;
- end-to-end latency and, for local runs, model time (p50, p90, maximum);
- failures by stage: routing, retrieval, generation, validation, execution,
  scientific, reporting, or environment.

Tests
-----

Routine tests in ``pycsamt/assistant/evals/tests/test_answers.py`` are
deterministic. They use recorded job records, a stub API resolver and a fake
chat runner. One integration test drives the real chat path with a mocked
local-model transport. Tests that need a real model are marked ``live``, and
CI deselects them with ``-m "not live"``.

Limitations
-----------

- A passing pattern check is not proof that an answer is scientifically
  correct. Cases stay incomplete until their ``human`` checks are reviewed.
- Static API resolution follows the checkout, not the installed package, and
  cannot follow dynamically constructed names.
- The edit-scope check compares non-comment lines. It rejects legitimate
  refactorings that the request did not ask for.
- Without an isolated executor, execution and artifact checks remain
  ``not_assessed``.
