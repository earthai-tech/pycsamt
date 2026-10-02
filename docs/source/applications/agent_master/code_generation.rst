.. _agent-master-code-generation:

Generating and editing scripts
==============================

Code requests carry your original wording, recent conversation, resolved
workflow and project context, retrieved evidence, and output requirements to
the code generator. Name the stations, quantities, units, filters, filenames,
and formats you need. Unusual requirements remain in the original request;
they are not discarded just because the workflow parser has no matching field.

For example::

    Write a Python script to load data/3edis, run QC, and save the table
    to qc_results.csv in an output folder.

Supported workflow templates remain available in Offline mode. A template
returned for a detailed request is labelled for review: Offline mode has not
applied arbitrary request-specific constraints. Select a local or cloud model
for custom compositions. An unspecified request such as “write me a script”
asks what the script should do, rather than defaulting to QC.

If a figure request specifies three panels but omits their quantities, the
assistant asks which quantities to show. Your answer is combined with the
pending request. Missing data directories are explicitly labelled as
placeholders. API evidence includes bounded static inspection of science
source signatures and documentation; this does not import or execute modules.

Follow-up edits
---------------

The conversation stores generated code and its artifact metadata. A recent
script can therefore be edited with requests such as::

    Change only the output directory to results/qc_review.
    Keep the same figure but save it at 300 dpi instead.
    Remove the correction step but keep the factor-estimation table.

The first two forms have a deterministic path for recognized literal output
paths and explicit ``savefig(..., dpi=...)`` calls. This changes only those
source locations. More complex expressions or edits go to the model with the
complete prior script. An explanation request does not regenerate the code.
Edited artifacts receive a new filename when the destination already exists,
preserving the earlier script.

Recent context is limited to six turns and 8000 characters of whole messages.
For a tight local context budget, optional templates are omitted first, then
retrieved evidence and API documentation are compacted with a recorded note.
The current request and prior script remain intact. Exact directory and DPI
edits bypass retrieval and model generation.
Older omitted turns are reported. The current request and previous script are
not silently shortened to fit a local model's context; an oversized prompt
can require a larger configured context or a smaller task.

Validation limits
-----------------

Generation does not execute the script. Static source/signature checks and
supported literal scientific checks are reported separately from skipped
runtime checks; see :doc:`code_validation`. They do not establish general
scientific validity or fulfillment of every requirement. Model edits still
need review. Small local models
may ignore instructions or invent APIs despite receiving accurate evidence.
Narrow checks flag absent requested filenames and recognized prohibited
static-shift correction calls. A draft that fails these checks is shown for
review without saving a script. Passing them is not a general correctness test.

Python interface
----------------

The shared request can also be supplied directly::

    from pycsamt.agents import CodeGenerationAgent, GenerationInput

    request = GenerationInput(
        original_request="Run QC and save the table as qc_results.csv.",
        workflow_config={"workflow": "qc", "data_path": "data/3edis"},
        output_requirements=["Export the QC table to qc_results.csv"],
    )
    result = CodeGenerationAgent(llm_provider="ollama").execute({
        "generation_input": request,
        "output_dir": "results/qc",
    })

The existing ``workflow_config``-only interface remains available for template
exports. The result records the generation mode, assumptions, script path,
and whether the edit used the narrow deterministic path. It never marks a
model-generated script as having all requirements verified.
