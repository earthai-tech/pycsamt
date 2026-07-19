# pyCSAMT-Agent Paper Reproducibility Archive

This directory is the working reproducibility archive for the revised
Computers & Geosciences manuscript.

## Layout

- `field/l26plt/`
  Real L26PLT EDI field-workflow audit, static-shift factors, source arrays,
  field-result macros, and the three-station replay subset.

- `ai_inversion/ai1d_mt1d_resnet5/`
  Five-layer MT1D synthetic inversion benchmark, train/validation/test arrays,
  trained checkpoint, checksum, predictions, metrics, run configuration, and
  manuscript macros.

- `assistant_benchmark/development/`
  Development-only assistant/RAG benchmark over the current 47-record suite.
  This is not the final 240-request benchmark.

- `benchmark/v2.0.0/`
  Candidate 240-request assistant benchmark with stable JSONL records and
  category manifest.

- `assistant_benchmark/v2.0.0/`
  Results from running the candidate 240-request suite against the rebuilt
  local RAG index.

- `retrieval_ablation/v2.0.0/`
  Deterministic retrieval ablation over the frozen 240-request benchmark.

- `code_replay/v2.0.0/`
  Offline generated-code audit plus deterministic replay of the archived
  three-station L26PLT subset.

- `figures/`
  Regenerated paper figures backed by source arrays.

## Current AI-Inversion Status

The selected paper-scale AI run is the wide ResNet experiment in
`ai_inversion/ai1d_mt1d_resnet5_wide/`. It is real and reproducible, and it
improves the data-space result over the smaller baseline, but it is still not a
strong layer-thickness recovery result:

- Train/validation/test: 50000/5000/5000
- Noise: 3% Gaussian response perturbation
- Checkpoint SHA-256:
  `3f316ca0958d73f980e548e2eab4ff7c931465a91d76dc97927c8c030b91850c`
- Log-resistivity RMSE: 0.701
- Mean thickness relative error: 220.6%
- Apparent-resistivity misfit: 0.279 log10
- Phase RMSE: 7.06 degrees
- Median CPU inference time: 5.20 ms

The high thickness error is scientifically plausible for a freely parameterised
five-layer MT inverse problem: layer boundaries are weakly identifiable from a
single 1-D response, and many models can fit the same apparent resistivity and
phase. It should not be presented as highly accurate layer-boundary recovery.

## AI-Inversion Ablation Status

This archive now contains the small ablation ladder used to strengthen the
AI-inversion result before freezing the manuscript macros. The goal was not to
claim that a neural network uniquely recovers five-layer thicknesses from a
single MT1D response. The goal was to choose a reproducible configuration that
improves held-out data-space fit while reporting the remaining model-space
non-uniqueness honestly.

| Run | Samples train/val/test | Architecture | Dropout | Blocks | Checkpoint | Model-space result | Data-space result |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `ai1d_mt1d_resnet5` | 10000/2000/2000 | ResNet baseline | default | default | `0826c29c217b98ce9f10367b89884ca18bf4c74e6c39a6271481e32ff404b9f4` | log-resistivity RMSE 0.742; mean thickness relative error 226.4% | apparent-resistivity misfit 0.352 log10; phase RMSE 9.28 degrees |
| `ai1d_mt1d_resnet5_wide` | 50000/5000/5000 | ResNet channels `(64, 128, 256, 512)` | 0.1 | 3 | `3f316ca0958d73f980e548e2eab4ff7c931465a91d76dc97927c8c030b91850c` | log-resistivity RMSE 0.701; mean thickness relative error 220.6% | apparent-resistivity misfit 0.279 log10; phase RMSE 7.06 degrees |

The wide run is the selected paper configuration because it improves the
held-out re-forwarded response metrics and slightly improves model-space
resistivity error. It does not solve layer-boundary identifiability, so the
paper should emphasize response fit, resistivity-scale recovery, inference
speed, and limitations rather than exact layer-thickness accuracy.

The following strengthening steps have been completed in this archive:

1. **Fixed target contract.** The training, validation, test, prediction, and
   generated code paths retain `n_layers=5` metadata. The archive must not mix
   five-layer targets with merged arrays that have lost layer-count metadata.
2. **Larger held-out experiment.** The selected run uses 50000/5000/5000
   synthetic responses instead of the smaller 10000/2000/2000 baseline.
3. **Wider/deeper ResNet.** The selected run uses channels
   `(64, 128, 256, 512)`, 3 residual blocks, dropout 0.1, and 0.01
   augmentation noise.
4. **Data-space reporting.** The selected metrics include apparent-resistivity
   misfit and phase RMSE from re-forwarded responses, not only normalized
   parameter error.
5. **Frozen checkpoint and macros.** The selected checkpoint is checksummed and
   the paper macros in `ai_inversion/ai1d_mt1d_resnet5_wide/` are generated
   from that frozen result.

Remaining limitations to state explicitly:

- The current synthetic prior is still broad; a constrained geological prior
  with monotonic depth limits, minimum layer thickness, and scenario classes
  would make layer-boundary targets more interpretable.
- Exact per-layer thickness remains weakly identifiable. If the manuscript
  reports layer metrics, it should frame them as diagnostic rather than as the
  main scientific claim.
- Robust target metrics such as depth-to-first-conductor, basement
  resistivity, and response-domain fit are more defensible than exact
  individual layer thicknesses.
- The current frozen paper run is a single selected seed set. A future
  ensemble of 3-5 independently trained seeds would support uncertainty bands
  and expose training instability.

## Current Assistant-Benchmark Status

The candidate 240-request suite has been generated and run. It is structurally
complete but still marked as a candidate until human adjudication confirms the
labels and relevance judgements.

- Record count: 240
- Rebuilt RAG index: 10109 chunks
- Intent macro-F1: 80.0%
- Canonical workflow macro-F1: 81.7%
- Raw workflow macro-F1: 75.7%
- Slot micro-F1: 99.3%
- Symbol recall: 87.5%
- Retrieval non-empty: 99.6%
- Workflow in top-k: 88.2%
- Grounded hallucination rate: 0.0%
- Review recall: 100.0%
- Unsafe execution rate: 5.0%

The workflow score improved after adding domain-language coverage for phrases
such as bad frequencies, robust filtering, induction arrows, inversion export,
depth of investigation, and 1-D checkpoint validation. Remaining misses are
mostly concrete code-template gaps rather than syntax or import failures.

## Current Retrieval/Code/Replay Status

- Retrieval target records: 208
- Static package context Recall@1/Recall@5/MRR/nDCG@10:
  0.0% / 26.0% / 0.124 / 0.187
- BM25 Recall@1/Recall@5/MRR/nDCG@10: 70.2% / 80.8% / 0.743 / 0.760
- Hashed dense Recall@1/Recall@5/MRR/nDCG@10:
  39.9% / 67.8% / 0.519 / 0.572
- BM25+hashed dense RRF Recall@1/Recall@5/MRR/nDCG@10:
  63.0% / 78.4% / 0.700 / 0.722
- Full deterministic retrieval Recall@1/Recall@5/MRR/nDCG@10:
  75.5% / 83.7% / 0.787 / 0.800
- Full deterministic retrieval+local reranking Recall@1/Recall@5/MRR/nDCG@10:
  63.0% / 85.6% / 0.716 / 0.768
- Median selected retrieval latency: 5.0 ms
- The dense rows use deterministic hashed token vectors rather than an
  external embedding API; the reranker is local and deterministic rather than
  LLM-based.
- Plan validity: 91.2%
- Generated-code syntax/import success: 100.0% / 100.0%
- Generated-code scientific coverage: 100.0%
- Generated-script execution: 100.0%
- Replay plan agreement/output completeness/warning agreement:
  100.0% / 100.0% / 100.0%
- Replay maximum normalised difference: 0.00e+00 at tolerance 1.0e-09
- Replay overhead: 0.98x

## Frozen AI Result Policy

For the manuscript, use only the selected wide run unless a new ablation is
performed and the whole archive is regenerated. A replacement run must:

1. preserve the `n_layers=5` target contract from data generation through
   evaluation;
2. improve held-out data-space metrics, especially apparent-resistivity misfit
   and phase RMSE;
3. report model-space metrics and identifiable target metrics side by side;
4. include the run configuration, predictions, metrics, generated macros, and
   checkpoint SHA-256;
5. state whether the result is a single-seed run or an ensemble.
