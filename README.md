# Mechanism class, attributed protein and the clinical severity of drug–drug interactions

Code and aggregate results for:

> Noor, A. Inferred Enzyme and Transporter Mechanisms and the Clinical Severity of Drug–Drug Interactions: A Cross-Resource Analysis with Drug-Level Sensitivity Testing. *Pharmaceuticals* (submitted).

## Contents

| Path | What it is |
|---|---|
| `code/protein_level_analysis.py` | Protein-level Fisher tests, Benjamini–Hochberg adjustment, E-values (square-root scale for a common outcome), single-protein and most-frequent-drug sensitivity analyses |
| `code/make_figures.py` | Figures 1–3 of the article: `python code/make_figures.py results figures` |
| `results/protein_level_results.csv` | Output of the protein-level analysis (Table 3, Table S1) |
| `results/class_level_contingency.csv` | Primary mechanism × severity counts (Table 1); every class-level statistic in the article can be recomputed from it |
| `results/external_validation_by_grade.csv` | DrugBank and KEGG confirmation counts by grade (Table 5) |
| `results/attribution_audit_summary.csv` | Audit categories (Table 4) |
| `results/accession_label_corrections.csv` | UniProt accessions whose source labels were corrected |

## Input data

`protein_level_analysis.py` runs on a pair × protein table with the columns
`CID_A, drug_A, CID_B, drug_B, enzyme_gene, uniprot, direction, severity`.

DDInter 2.0, KEGG, PubChem, UniProt and Wikidata are available from their providers. GoldD3R is available from the D3 dissertation repository (https://scholar.colorado.edu/csci_gradetds/106/). The protein-resolved attribution layer and the DrugBank-derived reference set contain content derived from third-party resources whose redistribution terms do not permit open release; the pair-level files are available from the author on reasonable request for non-commercial research, subject to the terms of the original providers.

## Licence

Code: MIT. Aggregate result tables and figures: CC BY 4.0. These licences cover the author's own work only and do not extend to third-party data.

## Contact

Adeeb Noor — arnoor@kau.edu.sa — ORCID 0000-0002-8251-1853
