# Limitations and validation status

Read this before presenting results. Stating these plainly is a strength in a methods review.

## What this software is not
- Not a validated pathogenicity or clinical tool.
- Not a new ΔΔG predictor with proven accuracy. The ridge model is a transparent baseline; its accuracy is whatever your grouped-CV report shows on your dataset.
- Does not generate or relax a mutant structure (FoldX `BuildModel` does this internally; the structure is not exposed here).
- Single-point substitutions in the first model of a PDB-format file; no mmCIF, no multi-mutants, no insertions/deletions.

## Known scientific limits
- **Contacts** are distance-based candidates (no hydrogens, angles or energies).
- **SASA** is computed on the supplied chains only (no crystal-contact neighbours, no ligands).
- **Charge** model is approximate (His neutral); **B-factors** are context-dependent and not a mutation-induced-dynamics measure.
- **ESM-2 LLR** measures evolutionary plausibility in sequence context, not stability. Residues missing from the coordinates are absent from the sequence the model sees (reported as a warning). Sequences > 1022 residues are cropped around the site.
- **Conservation entropy** depends on the homolog set BLAST returns (no e-value cutoff, no sequence weighting).
- **AlphaFold models**: B-factor columns hold pLDDT; the tool does not yet detect this, so B-factor outputs are not meaningful for such files.
- **Dataset leakage**: grouping by PDB ID is *not* homology-aware. For publishable numbers use sequence-cluster IDs (≤ 25–30 % identity) in the `group` column.

## Validation status
| Item | Status |
|---|---|
| DAG engine, config, contacts maths, feature extraction, ESM scoring maths, ridge/CV/bootstrap, research workflows, pipeline wiring | Unit-tested offline with synthetic inputs and fakes. |
| Biopython-based parsing, SASA, neighbors, descriptors | Covered by the v1.x/2.0 tests that need Biopython; run `make test` in your environment. |
| Real-structure behaviour (1AKI etc.), remote services (DynaMut2, EBI), ESM weights download | **Needs a live run** (`make live`, `scripts/live_check.sh`). Not exercised during development of v2.0.0. |
| Type checking / lint | Run `make typecheck lint`; CI reports them non-blocking until you have a clean local run. |
| Scientific accuracy of the ML layer | **Unknown until you train and evaluate on a real dataset** (`docs/RESEARCH_PROTOCOL.md`). |

## Roadmap (not implemented)
mmCIF input; AlphaFold pLDDT detection; UniProt feature annotations; known-variant lookup (ClinVar/gnomAD); BLAST e-value cutoff + sequence weighting; mutant-structure overlay; resumable remote jobs; Zenodo/PyPI release automation.
