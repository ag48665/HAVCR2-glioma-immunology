# HAVCR2 Is Associated With Immune and Macrophage Enrichment and Unfavourable Clinical Outcome in Glioma: An Integrated TCGA and CGGA Analysis

## Overview

This repository contains the analysis code, figures, and statistical results associated with an integrated analysis of HAVCR2 (TIM-3), immune infiltration, and myeloid-cell enrichment in diffuse glioma.

Transcriptomic and clinical data from The Cancer Genome Atlas (TCGA) lower-grade glioma (LGG) and glioblastoma (GBM) cohorts were analyzed, with key findings evaluated in the independent Chinese Glioma Genome Atlas (CGGA) mRNAseq_693 cohort.

The study evaluates relationships among immune infiltration, immune checkpoint expression, HAVCR2 expression, clinical outcome, myeloid-cell enrichment, macrophage- and microglia-associated transcriptional signatures, HAVCR2-centered co-expression, and immune-related pathway enrichment.

## Study Highlights

- High ESTIMATE immune scores were associated with significantly shorter overall survival in the combined TCGA glioma cohort.
- Higher immune score remained associated with worse overall survival after adjustment for age and tumor grade and stratification by IDH status.
- Immune-high tumors showed increased expression of multiple immune checkpoint genes, including PDCD1, CD274, CTLA4, LAG3, TIGIT, and HAVCR2.
- HAVCR2 showed a strong positive association with ESTIMATE-derived immune infiltration.
- Higher HAVCR2 expression was associated with shorter overall survival in univariable analyses in the independent CGGA cohort.
- HAVCR2 was not independently associated with overall survival after adjustment for age and tumor grade and stratification by IDH status in CGGA.
- HAVCR2 expression was strongly associated with computationally estimated myeloid-cell enrichment.
- Strong associations between HAVCR2 and macrophage-associated transcriptional enrichment were reproduced in both TCGA and CGGA.
- HAVCR2 was associated with macrophage- and microglia-related transcriptional signals in both LGG and GBM.
- HAVCR2-centered co-expression analyses identified broadly conserved immune- and myeloid-associated transcriptional programs.
- HAVCR2-associated transcriptional profiles were enriched for inflammatory and immune-related Hallmark pathways.

## Datasets

### TCGA

- TCGA-LGG (Lower Grade Glioma)
- TCGA-GBM (Glioblastoma)

### CGGA

- CGGA mRNAseq_693 cohort

All analyses were performed using publicly available, de-identified datasets.

## Repository Structure

figures/ - Manuscript figures

results/ - Statistical results and summary tables

scripts/ - R analysis scripts

## Main Figures

### Figure 1
Immune score survival analysis (TCGA)

### Figure 2
Immune checkpoint expression in immune-high and immune-low tumors

### Figure 3
Correlations between immune score and checkpoint genes

### Figure 4
Association between HAVCR2 expression and immune score

### Figure 5
CGGA survival analysis according to HAVCR2 expression

### Figure 6
HAVCR2 expression according to grade and IDH status

### Figure 7
HAVCR2 and macrophage marker score in TCGA

### Figure 8
CGGA validation of the HAVCR2-macrophage association

## Key Results

- High immune score and overall survival: log-rank p < 0.0001.

- TCGA multivariable immune-score analysis: HR per 100-unit increase = 1.032 (95% CI 1.015-1.050), p = 2.54 x 10^-4.

- HAVCR2-immune score association: Spearman rho = 0.887 (n = 601).

- HAVCR2-macrophage marker score in TCGA: Spearman rho = 0.905 (n = 800).

- HAVCR2-macrophage marker score in CGGA: Spearman rho = 0.860.

- CGGA HAVCR2 survival analysis: log-rank p < 0.0001.

- HAVCR2 univariable Cox analysis in CGGA: HR = 1.035 (95% CI 1.022-1.047), p = 2.61 x 10^-8.

- HAVCR2 multivariable Cox analysis in CGGA: HR per 1-SD increase = 1.052 (95% CI 0.961-1.151), p = 0.274.

## Reproducibility

Analyses were performed using R version 4.5.2.

Major tools and packages used in the analysis include:

- TCGAbiolinks
- SummarizedExperiment
- survival
- survminer
- ggplot2
- dplyr
- xCell2
- fgsea
- msigdbr

Statistical analyses included Kaplan-Meier survival analysis, log-rank tests, Cox proportional hazards regression, Wilcoxon rank-sum tests, Spearman rank correlations, cell-type deconvolution, marker-based transcriptional scoring, co-expression analysis, and preranked gene set enrichment analysis.

## Data Availability

TCGA data are publicly available through the Genomic Data Commons Data Portal.

CGGA data are publicly available through the Chinese Glioma Genome Atlas data portal.

No new patient-level datasets were generated during this study.

The repository provides the analysis workflow, code, statistical results, and supporting materials required to reproduce the analyses from the publicly available source datasets.

## Author

Agata Gabara

## License

This repository is provided for academic and research purposes.

## Citation

If you use this repository, please cite the associated manuscript:

Gabara A. HAVCR2 Is Associated With Immune and Macrophage Enrichment and Unfavourable Clinical Outcome in Glioma: An Integrated TCGA and CGGA Analysis. Cancer Informatics. Under review.
