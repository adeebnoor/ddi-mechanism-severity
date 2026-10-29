# Mechanism class, attributed protein and the clinical severity of drug–drug interactions

Code and aggregate results for:

> Noor, A. Inferred Enzyme and Transporter Mechanisms and the Clinical Severity of Drug–Drug Interactions: A Cross-Resource Analysis with Drug-Level Sensitivity Testing. *Pharmaceuticals* (submitted).

## Contents

| Path | What it is |
|---|---|
| `code/class_level_analysis.py` | Chi-square with Monte Carlo p, Cramér's V, hypotheses H1–H5 (Fisher, BH, exact CI for H5), risk ratios and E-values — `python code/class_level_analysis.py results` |
| `code/protein_level_analysis.py` | Protein-level Fisher tests, BH adjustment, leave-one-drug-out over all drugs, leflunomide/methotrexate exclusion, conflicting-attribution exclusion, multiple-membership drug random-intercept model — `python code/protein_level_analysis.py <data_dir> results/protein_level_results.csv` |
| `code/external_checks.py` | DrugBank/KEGG confirmation by grade, Wilson CIs, Fisher and Cochran–Armitage tests — `python code/external_checks.py results` |
| `code/make_figures.py` | Figures 1–3 — `python code/make_figures.py results figures` |
| `data_pseudonymised/enzyme_pair_severity.csv` | Graded pair × protein records with drug identities replaced by codes (leflunomide and methotrexate kept by name, as the analysis needs them). Reproduces Table 3 and Table S2: `python code/protein_level_analysis.py data_pseudonymised results/protein_level_results.csv` |
| `results/class_level_contingency.csv` | Primary mechanism × severity counts (Table 1) |
| `results/protein_level_results.csv` | Output of the protein-level analysis (Table 3, Table S2) |
| `results/external_validation_by_grade.csv` | DrugBank and KEGG confirmation counts by grade (Table 5) |
| `results/attribution_audit_summary.csv` | Audit categories (Table 4) |
| `results/accession_label_corrections.csv` | UniProt accessions whose source labels were corrected (Table S1) |

## Input data

`protein_level_analysis.py` runs on the pseudonymised table supplied here, or on the original pair × protein table, with the columns
`CID_A, drug_A, CID_B, drug_B, enzyme_gene, uniprot, direction, severity`.

DDInter 2.0, KEGG, PubChem, UniProt and Wikidata are available from their providers. GoldD3R is available from the D3 dissertation repository (https://scholar.colorado.edu/csci_gradetds/106/). The identified (non-pseudonymised) protein-resolved attribution layer and the DrugBank-derived reference set contain content derived from third-party resources whose redistribution terms do not permit open release. An earlier version is archived at Zenodo (https://doi.org/10.5281/zenodo.23169653) with restricted access; the pair-level files are available from the author on reasonable request for non-commercial research, subject to the terms of the original providers.

## Licence

Code: MIT. Aggregate result tables and figures: CC BY 4.0. These licences cover the author's own work only and do not extend to third-party data.

## Contact

Adeeb Noor — arnoor@kau.edu.sa — ORCID 0000-0002-8251-1853
