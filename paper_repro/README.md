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

## How To Strengthen The AI Result

For a robust C&G paper result, run a small ablation ladder and report the best
frozen configuration only after it improves held-out data-space and model-space
metrics:

1. **Fix the target contract**
   Keep `n_layers=5` metadata through training and evaluation. This archive now
   does that; do not train on merged arrays without metadata.

2. **Use a constrained geological prior**
   Replace fully random five-layer thicknesses with a realistic prior:
   monotonic cumulative depth limits, minimum layer thickness, and scenario
   classes such as cover-conductor-basement. This makes thickness recovery a
   scientifically meaningful target instead of an underdetermined permutation.

3. **Train larger and regularise less aggressively**
   Test 50000/5000/5000 samples, dropout 0.05-0.15, and deeper/wider ResNet
   channels such as `(64, 128, 256, 512)` with 3 residual blocks. The current
   run overfits after about 12 epochs, suggesting the default model/data
   pairing is not optimal.

4. **Add a physics/data-space selection criterion**
   Select checkpoints by validation apparent-resistivity and phase misfit, not
   only normalised parameter MSE. Report both model-space and re-forwarded
   response metrics.

5. **Report identifiable targets**
   Add robust metrics for depth-to-first-conductor, basement resistivity, and
   data-space fit. These are more defensible than exact individual layer
   thicknesses for non-unique 1-D EM inversion.

6. **Use an ensemble or uncertainty bands**
   Train 3-5 independent seeds and report median plus bootstrap/seed intervals.
   This will expose instability instead of hiding it.

7. **Freeze only after ablation**
   Once the final configuration is chosen, regenerate this archive, checksum the
   checkpoint, and include only the frozen macros in the paper.
