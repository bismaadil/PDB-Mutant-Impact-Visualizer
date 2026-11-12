# Methods, assumptions, and limits

## Inputs and residue numbering

The application reads deposited PDB-format coordinates and uses the selected model's author chain and residue numbering. It rejects missing residues, non-standard amino acids, and mismatches between the requested wild-type amino acid and the structure. Insertion codes can be supplied for CLI analyses. The exact coordinate file is recorded by SHA-256 digest.

## Stability

- **DynaMut2** submits the supplied PDB file, chain, and single mutation to its public job API. The returned value and job identifier are preserved as service output. A service prediction is not an experimental measurement. Cite the DynaMut2 publication when reporting results.
- **FoldX** runs local `RepairPDB` followed by `BuildModel`, then reads the `total energy` field from the `Dif_*.fxout` difference output. FoldX is separate software with its own license terms; obtain and install it yourself. Check the raw output and FoldX version for research use.
- Stability prediction can be disabled. A failed or omitted backend is represented as missing; the program does not substitute zero.

## Conservation

The optional EMBL-EBI path resolves a PDB author residue through PDBe SIFTS, verifies that the resolved UniProt sequence has the expected wild-type residue, searches UniProtKB/Swiss-Prot using EBI NCBI BLAST+, retrieves hit sequences through Dbfetch, and aligns them with Clustal Omega. The output reports the observed amino-acid counts, non-gap coverage, and Shannon entropy at the aligned query position. It is not a ConSurf score, a phylogenetic conservation model, or a pathogenicity classifier. Homolog selection and database contents influence the result; inspect the returned alignment and identifiers.

EMBL-EBI requests are subject to service terms, rate limits, job retention, and availability. A valid contact email is required for programmatic submission. For a user-supplied PDB without a PDB ID, provide a curated UniProt accession and 1-based UniProt position; the program does not guess this mapping.

## Coordinate-derived observations

Neighbor residues are counted if any non-hydrogen atom is within the selected radius of the target residue's C-alpha atom. Site and protein mean B-factors are reported in Å² as present in the input. B-factors depend on experimental/model context and are not direct measurements of mutation-induced dynamics. The interactive view displays C-alpha positions and deposited B-factors; it does not display a generated mutant structure.

## Offline descriptors

- **Solvent accessibility:** Shrake-Rupley (probe 1.40 Å, 960 points) on the first model's protein atoms; relative SASA divides by the Tien et al. (2013) theoretical maximum for the wild-type residue.
- **Substitution properties:** BLOSUM62, Kyte-Doolittle hydropathy, Zamyatnin residue volume and an approximate pH 7 charge (His neutral). These describe the two amino acids, not the structural context.
- **Sign conventions:** FoldX positive ddG = destabilising; DynaMut2 positive = stabilising. Raw values are never comparable across backends.

## Interpretation

Measurements remain separate. The project deliberately does not combine stability, conservation, B-factor, or neighbor count into one impact score. These methods estimate different quantities and need task-specific validation before being combined. This application is a research workflow aid, not a clinical tool or substitute for experiments.

## Sources to cite

- Rodrigues CHM et al. DynaMut2: Assessing changes in stability and flexibility upon single and multiple point missense mutations. *Protein Science* (2021). https://doi.org/10.1002/pro.3942
- Delgado J et al. FoldX 5.0: working with RNA, small molecules and a new graphical interface. *Bioinformatics* (2019). https://doi.org/10.1093/bioinformatics/btz184
- Madeira F et al. The EMBL-EBI Job Dispatcher sequence analysis tools framework in 2024. *Nucleic Acids Research* (2024). https://doi.org/10.1093/nar/gkae241
- PDBe SIFTS mapping service: https://www.ebi.ac.uk/pdbe/coordinates/
- RCSB PDB file download service: https://www.rcsb.org/docs/programmatic-access/file-download-services

## Contacts (v2.0)

Computed from heavy-atom coordinates in the first model (standard residues only): side-chain contact residues within 4.5 Å (configurable); hydrogen-bond **candidates** (N/O to N/O within 3.5 Å, excluding trivial peptide-bond neighbours); salt-bridge **candidates** (Asp/Glu carboxylate O to Lys/Arg N within 4.0 Å); disulfides (Cys SG–SG within 2.5 Å); interface residue (any heavy atom within 5 Å of another chain); and the number of residues whose Cβ (Cα for Gly) lies within 10 Å of the site's Cβ. Hydrogen positions, angles and energies are not evaluated.

## Language-model score (v2.0, optional)

ESM-2 masked-marginal log-likelihood ratio: mask the site, read the model's distribution, report log p(mutant) − log p(wild type) (Meier et al., 2021; Lin et al., 2023). Negative = mutant less plausible. Sequences longer than 1022 residues are cropped around the site; residues missing from the coordinates are absent from the input sequence (warning issued). This is evolutionary plausibility, not ΔΔG.

## Research model (v2.0, optional)

Ridge regression on standardised descriptors, trained by the user (`pdbmi train`). Evaluation: protein-/cluster-grouped k-fold CV, cluster-bootstrap CIs, mean-predictor and single-feature baselines, feature-group ablations. The model card records the training-data SHA-256, target sign convention, CV metrics and limitations. No pretrained weights are distributed.

- Meier J et al. Language models enable zero-shot prediction of the effects of mutations on protein function. *NeurIPS* (2021). https://doi.org/10.1101/2021.07.09.450648
- Lin Z et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science* (2023). https://doi.org/10.1126/science.ade2574
