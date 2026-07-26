# Independent benchmark adjudication protocol

This protocol must be completed by a domain reviewer who did not construct the
benchmark or inspect system predictions. The benchmark and predictions are
frozen before review.

## Material supplied to the reviewer

- the 240 user-visible requests, with stable record IDs;
- definitions of the allowed intent and workflow labels;
- the package workflow registry and public API documentation;
- annotation guidance for line, symbol, unsupported, and review-required fields.

The reviewer must not receive `got_*` fields, retrieved rankings, benchmark
metrics, or the original labels during independent annotation.

## Required annotation

For every record:

1. assign the top-level intent;
2. assign the canonical workflow when applicable;
3. extract the survey line and other explicit slots;
4. mark whether the request is unsupported or requires clarification/review;
5. flag genuinely ambiguous records and explain why.

For retrieval evaluation, independently judge all expected symbols where
practical. If resources require sampling, use a pre-declared stratified random
sample covering every record category and report its size and seed.

## Analysis

Before adjudication, report:

- exact percentage agreement and Cohen's kappa for intent;
- exact percentage agreement and Cohen's kappa for workflow;
- micro-F1 agreement for extracted slots;
- agreement for unsupported and review-required decisions.

Resolve disagreements by documented discussion without inspecting system
predictions. Preserve the original, independent, and adjudicated labels.
Recompute all headline benchmark metrics on the adjudicated labels and report
both the original and adjudicated results. Do not alter queries or remove
difficult cases after predictions have been inspected.

## Archive additions

The completed archive should contain:

- reviewer role and relevant domain expertise;
- blinded annotation file;
- disagreement and adjudication log;
- agreement statistics;
- adjudicated benchmark version;
- regenerated metrics and consistency-audit output.
