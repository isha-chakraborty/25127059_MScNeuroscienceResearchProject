# Exploring the contributions of brain DNA methylation changes to copper and iron metabolism dysregulation across neurodegenerative diseases

This repository contains the reproducible R workflow accompanying an MSc research project investigating copper- and iron-associated DNA methylation networks across neurodegenerative diseases. The analysis integrates published post-mortem brain co-methylation networks from Alzheimer’s disease (AD), frontotemporal lobar degeneration (FTLD), Parkinson’s disease (PD), progressive supranuclear palsy (PSP) and multiple system atrophy (MSA), curated metal-homeostasis gene sets, published cell-type enrichment results, functional enrichment and an independent FTLD RNA-sequencing dataset.

The final corrected workflow distinguishes:

- disease-associated modules selected using source-network-specific Bonferroni thresholds;
- modules enriched for copper, broad-iron or core-iron gene sets;
- individual metal-associated genes recurring across retained modules; and
- module-derived candidates showing differential expression in FTLDexp.

The workflow reproduces the final nine-module analysis used in the submitted thesis. It does not filter an obsolete 16-module recurrence analysis.

## Final analysis at a glance

- **9** metal-enriched disease-associated co-methylation modules
- **12** significant module–gene-set enrichment results
- **140** unique metal-associated candidates
- **89** candidates occurring in at least two modules
- **13** recurrent copper genes
- **77** recurrent broad-iron genes
- **57** recurrent core-iron genes
- **771** significant module–term functional-enrichment results across eight modules in the submitted g:Profiler run
- **0** FDR-significant terms for the three recurrent gene lists
- **40** significant FTLD module–gene expression occurrences representing **25** unique genes

## Repository structure

```text
.
├── R/
│   ├── utils.R
│   ├── 00_setup.R
│   ├── 01_preprocess_AD_methylation.R
│   ├── 02_curate_metal_gene_sets.R
│   ├── 03_select_disease_associated_modules.R
│   ├── 04_metal_gene_set_enrichment.R
│   ├── 05_integrate_published_EWCE.R
│   ├── 06_candidate_prioritisation_MM_GS.R
│   ├── 07_recurrence_analysis_9_modules.R
│   ├── 08_plot_occurrence_matrices.R
│   ├── 09_plot_circos_networks.R
│   ├── 10_functional_enrichment_gprofiler.R
│   ├── 11_integrate_FTLDexp.R
│   └── 12_validate_final_outputs.R
├── data/
│   ├── raw/                    # Not committed; optional AD preprocessing only
│   ├── gene_sets/              # AmiGO exports used to curate the three gene sets
│   ├── published_results/      # Neutral network annotations and FTLDexp input
│   └── derived/                # Documented non-sensitive analysis inputs
├── outputs/
│   ├── derived_data/
│   ├── figure_source_data/
│   ├── figures/
│   ├── logs/
│   ├── supplementary_tables/
│   ├── tables/
│   └── validation/
├── renv.lock
└── run_all.R
```

Raw methylation data, large intermediate objects and donor-level information should not be committed to the repository.

## Workflow

| Script | Purpose |
|---|---|
| `R/00_setup.R` | Defines packages, reproducible options and project-relative directories. |
| `R/01_preprocess_AD_methylation.R` | Imports, quality-controls, filters and BMIQ-normalises the AD methylation data. This resource-intensive step is optional when saved outputs or published networks are used. |
| `R/02_curate_metal_gene_sets.R` | Reconstructs the 47 copper, 199 broad-iron and 136 core-iron gene sets using human experimental-evidence annotations. |
| `R/03_select_disease_associated_modules.R` | Applies source-network-specific Bonferroni thresholds to published module–disease associations. |
| `R/04_metal_gene_set_enrichment.R` | Uses two-sided Fisher tests and cohort-specific annotated-gene backgrounds to identify metal-enriched modules. |
| `R/05_integrate_published_EWCE.R` | Integrates the cell-type enrichment outputs reported by the source methylation studies. |
| `R/06_candidate_prioritisation_MM_GS.R` | Calculates absolute module membership and gene significance, Spearman correlations and candidate rankings; generates Figure 6 and Supplementary Figure S2. |
| `R/07_recurrence_analysis_9_modules.R` | Builds recurrence directly from the corrected nine modules and produces Supplementary Table S5. |
| `R/08_plot_occurrence_matrices.R` | Generates Figures 7–8 and the complete 77-gene broad-iron Supplementary Figure S1. |
| `R/09_plot_circos_networks.R` | Generates Figures 9–10. Circos links represent repeated gene–module occurrence, not molecular interactions. |
| `R/10_functional_enrichment_gprofiler.R` | Runs module and recurrent-set enrichment using `gprofiler2::gost()`; generates Figure 11, Table 5 and Supplementary Table S6. |
| `R/11_integrate_FTLDexp.R` | Integrates four retained FTLD modules with differential-expression results; generates Figure 12, Table 6 and Supplementary Table S7. |
| `R/12_validate_final_outputs.R` | Checks reported totals, required deliverables, obsolete wording, package versions and output checksums. |

`R/utils.R` holds shared functions for symbol cleaning, module-name standardisation, Fisher testing, input checks and consistent figure export.

## Software environment

The submitted analysis used R 4.5.1. Package versions are recorded in `renv.lock`; the final validation script additionally writes the active session and installed package versions to `outputs/validation/`.

Install `renv` once, then restore the recorded environment:

```r
install.packages("renv")
renv::restore()
```

Bioconductor packages used for methylation preprocessing include ChAMP, minfi and limma. The downstream workflow uses packages including dplyr, tidyr, readr, ggplot2, ggrepel, patchwork, circlize and gprofiler2.

## Required inputs

The scripts use project-relative paths and contain no user-specific computer locations. The repository includes the non-sensitive inputs needed to reproduce the downstream nine-module analysis. Raw donor-level methylation files are excluded and are needed only for the optional AD preprocessing step.

Key public dataset accessions include:

- AD multi-region methylation data: **GSE125895**
- FTLD2 multi-omics data: **E-MTAB-12674**
- FTLD3/PSP methylation data: **GSE75704**
- FTLDexp RNA-sequencing data: **GSE153960**

The included neutral inputs are:

```text
data/gene_sets/copper/amigo/*_3.txt
data/gene_sets/iron/amigo/*_3.txt
data/published_results/network_gene_annotations/*.csv
data/published_results/FTLDexp_differential_expression.csv
data/derived/nine_module_MM_GS_annotated_records.csv
data/derived/submitted_nine_module_gprofiler_significant_results.rds
data/derived/submitted_recurrent_sets_gprofiler_significant_results.rds
```

The network annotation CSV files contain only gene symbols and cohort-specific module assignments. The MM–GS file contains the derived values used for the submitted plots, and the frozen g:Profiler objects preserve the results returned for the submitted analysis. Compatibility code can also locate the completed local project’s published-network RData objects by their contained object names, without embedding personal paths. Redistribution of original data remains subject to the source repositories and studies.

## Running the analysis

From the repository root:

```bash
Rscript run_all.R
```

The default workflow skips resource-intensive AD preprocessing and reuses the frozen submitted g:Profiler results. If those objects are absent, `auto` mode submits live queries.

Optional modes:

```bash
# Include raw AD preprocessing
Rscript run_all.R --run-ad-preprocessing

# Force new g:Profiler queries
Rscript run_all.R --run-gprofiler

# Require reuse of completed g:Profiler results
Rscript run_all.R --reuse-gprofiler
```

Individual scripts can also be run in numerical order from RStudio or with `Rscript`.

## Functional-enrichment settings

Module queries use:

- organism: *Homo sapiens* (`hsapiens`);
- unranked gene lists;
- GO Biological Process, Reactome and WikiPathways;
- the unique annotated genes in the corresponding published network as the custom background; and
- FDR-adjusted `P < 0.05`.

Recurrent copper, broad-iron and core-iron queries use their respective complete curated metal sets as custom backgrounds. g:Profiler is a live resource, so database updates can alter results from a new query. The submitted result objects should be retained to reproduce the thesis figures exactly; Script 10 reports any difference produced by a later release.

## FTLDexp integration

The expression analysis evaluates all curated metal-associated candidates assigned to FTLD1 red, FTLD2 darkmagenta, FTLD3 black and FTLD3 turquoise. Metal-set membership is retained separately, allowing genes directly contributing to module enrichment to be distinguished from other metal-associated genes in the same network.

Expression direction comes from the FTLDexp log2 fold change. The methylation analysis identifies the module in which each candidate occurred; it does not assign an expression direction.

## Validation

After all requested steps, `R/12_validate_final_outputs.R` confirms the reported totals, required figures and tables, and the absence of obsolete 16-module or six-FTLD-module wording. It writes:

```text
outputs/validation/final_validation_report.csv
outputs/validation/final_validation_summary.csv
outputs/validation/obsolete_wording_audit.csv
outputs/validation/package_versions.csv
outputs/validation/sessionInfo.txt
outputs/validation/final_output_manifest.csv
```

The workflow stops if a required validation check fails.

## Main outputs

The R workflow generates Figures 5–12, Tables 3A–6, Supplementary Figures S1–S2 and Supplementary Tables S1–S7, together with machine-readable source data. Conceptual or adapted figures created outside R are not regenerated by these scripts.

Supplementary Table S7 is exported in two clearly labelled forms: the complete 188-row candidate-integration audit and the 40 significant module–gene occurrences displayed in the submitted thesis. These are different views of the same integration analysis, not different results.

## Reproducibility and interpretation notes

- Module colours are cohort-specific labels. Identically coloured modules from separate datasets do not imply identical CpGs or genes.
- Recurrence counts repeated representation across distinct module–disease comparisons; it is not a measure of enrichment strength.
- Circos links connect occurrences of the same gene and do not represent molecular interactions.
- Module-level methylation networks cannot be labelled uniformly hypermethylated or hypomethylated without CpG-level genomic-context analysis.
- The analysis identifies associations and prioritised candidates; it does not establish causality.

## Source studies

The workflow integrates outputs or data described in the thesis and the associated source publications, including Fodder et al. (2023), Murthy et al. (2024), Fodder et al. (2026), Semick et al. (2019), Menden et al. (2023), Weber et al. (2018) and Hasan et al. (2022). Full citations are provided in the submitted thesis.
