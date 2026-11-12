# Research protocol

A pre-specified plan makes the project a *study* rather than a tool demo. Fill in the bracketed parts with your own choices and results; keep this file under version control and date your changes.

## 1. Question
Do **interpretable structure-derived descriptors**, with and without a **zero-shot protein language-model score (ESM-2)**, explain experimental changes in protein stability for single substitutions **when evaluation is homology-aware**, and how much do they add over trivial baselines?

## 2. Hypotheses (state before running)
- **H1** Structure + substitution descriptors rank experimental ΔΔG better than a mean predictor and better than any single descriptor (grouped-CV Spearman, cluster-bootstrap CI excluding 0).
- **H2** Adding the ESM-2 LLR improves H1's model (CI of the paired improvement excludes 0).
- **H3** The ESM-2 LLR alone is weaker than the combined model for *stability* (it measures plausibility, not ΔΔG).
Negative results are reportable; say so in advance.

## 3. Data
| Need | Options (verify current access, licence and sign convention yourself) |
|---|---|
| Experimental stability changes | ProThermDB / ProTherm; FireProtDB; the S669 benchmark (Pancotti et al., *Brief. Bioinform.* 2022); the mega-scale cDNA-display proteolysis data (Tsuboyama et al., *Nature* 2023). |
| Structures | RCSB PDB entries named in the dataset (downloaded by the pipeline). |
| Homology clusters | Cluster protein sequences (e.g. MMseqs2 or CD-HIT at ≤ 25–30 % identity) and put the cluster ID in the `group` column. |

Required CSV columns: `pdb_id, chain, mutation` (+ `ddg`, `group`). Record dataset name, version, download date and licence in your notebook. **Fix one sign convention** (e.g. positive = destabilising) before training and pass it to `--target-convention`; the tool refuses to train without it.

## 4. Procedure (all commands are in the repo)
```bash
pdbmi features --input data/my_dataset.csv --out data/features.csv                 # offline descriptors
pip install -e ".[ml]" && pdbmi features --input data/my_dataset.csv --out data/features_esm.csv \
     --esm-model facebook/esm2_t12_35M_UR50D                                       # adds ESM-2 LLR
pdbmi train --features data/features_esm.csv --feature-groups substitution,structure,language_model \
     --target-convention "positive = destabilising" --k 5 --seed 0 --out data/results/model
```
Outputs: `model.json` (weights + model card with training-data SHA-256), `ablation.md` (substitution → +structure → +language model, and ESM-only), grouped-CV metrics against a mean baseline, and cluster-bootstrap 95 % CIs.

## 5. Evaluation design
- **Grouped k-fold** (all variants of a cluster in one fold) to avoid leakage.
- Metrics: Spearman, Pearson, RMSE, MAE, AUROC for "destabilising" (threshold 1 kcal/mol, an explicit, adjustable choice).
- Baselines: mean predictor; single-feature Spearman reported in the model card.
- Uncertainty: cluster bootstrap over proteins, `--bootstrap 1000`, fixed `--seed`.
- Sensitivity checks you should run and report: other seeds, other k, other ridge α, a different identity threshold, dropping the largest protein cluster.
- Hold-out: if you can, evaluate the final model once on an independent set with low similarity to training data; never tune on it.

## 6. Threats to validity (discuss in the write-up)
Experimental noise and dataset biases (stabilising mutations are rare; alanine-scanning-heavy sets); heterogeneous experimental conditions; wild-type structure ≠ mutant structure; remaining sequence similarity across folds; ridge linearity; ESM context truncation and coordinate gaps; sign-convention errors (the most common silent bug).

## 7. Reporting checklist
- [ ] Dataset name/version/date/licence and sign convention
- [ ] Homology clustering method and threshold
- [ ] Model card + `ablation.md` committed with the exact `features.csv` SHA-256
- [ ] Image digest or `constraints.txt`, git commit, seed
- [ ] Table with CIs for every comparison; negative results included
- [ ] Limitations section (copy from `LIMITATIONS_AND_VALIDATION.md` and adapt)
