# Presenting this project to a selection committee

## One-paragraph abstract (edit and add your numbers)
> Predicting how a single amino-acid substitution changes protein stability is routinely done with opaque tools whose results are hard to reproduce or compare. I built **PDB Mutant Impact Visualizer**, a containerised, DAG-orchestrated workflow that validates a mutation against a real PDB structure, computes interpretable structural and substitution descriptors, optionally adds a protein language-model score (ESM-2) and external predictors (DynaMut2/FoldX), and records full provenance (coordinate and configuration hashes, versions, per-stage status). Using [dataset], I evaluated whether descriptors alone, and descriptors plus ESM-2, rank experimental ΔΔG under **homology-aware grouped cross-validation** with cluster-bootstrap confidence intervals and mean/single-feature baselines. [Result in one sentence with CI.] The work emphasises reproducibility and honest evaluation over headline accuracy.

## What is genuinely yours / defensible
- A reproducible pipeline architecture (DAG with failure isolation and a self-generating diagram).
- A leakage-aware evaluation design, with ablations and uncertainty.
- Provenance and a `reproduce` command.
- Honest scoping: no composite score, sign conventions recorded, limitations documented.

## What you must **not** claim
- That it beats DynaMut2, FoldX, ESM-based or other published predictors (unless your held-out evaluation shows it, with CIs).
- That it is clinically validated or predicts pathogenicity.
- Any number you did not generate yourself from `pdbmi train` on a cited dataset.
- That contacts are hydrogen bonds or salt bridges (they are geometric candidates).

## 5-minute talk structure
1. Problem and why reproducibility matters (30 s) · 2. The DAG diagram (60 s) · 3. Demo: `pdbmi run` then `pdbmi reproduce` (60 s) · 4. Study design: hypotheses, data, grouped CV (90 s) · 5. Results with CIs, including what did *not* work (45 s) · 6. Limitations and next steps (15 s).

## Questions to expect — and honest answers
1. **Why not just use FoldX/DynaMut2?** They are integrated as optional backends; the contribution is reproducible orchestration and an evaluation design that makes tools comparable.
2. **How do you avoid data leakage?** Grouped CV by homology cluster; all variants of one cluster stay together (explain why PDB-ID grouping alone is insufficient).
3. **Why ridge regression?** Auditable and a strong honest baseline; coefficients are interpretable; complexity must earn its place against it.
4. **What does the ESM score mean?** log p(mutant) − log p(wild type) at a masked position: evolutionary plausibility, not ΔΔG.
5. **Sign conventions?** FoldX positive = destabilising; DynaMut2 the opposite; the model card records the training convention.
6. **Is the result reproducible?** Show the manifest, `pdbmi reproduce`, container digest, seed.
7. **What are the biggest limitations?** Wild-type structure only; dataset bias; no mutant structure; see `LIMITATIONS_AND_VALIDATION.md`.
8. **What would you do next with funding?** mmCIF/AlphaFold input, sequence-weighted conservation, mutant-structure modelling, external validation set, active-learning for experimental design.
9. **What did you learn when something failed?** Have one real example from your runs (e.g. a mapping failure or a leakage bug you caught).
10. **Which parts did you write and which are tools?** Be precise: Biopython, ESM-2, DynaMut2, FoldX are third-party; orchestration, evaluation and provenance are your work. If you used AI assistance for code, say so and be ready to explain every module.
