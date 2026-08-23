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

# Figures
# Figures

This directory contains the main figures associated with the manuscript:

**Gabara A. HAVCR2 Is Associated With Immune and Macrophage Enrichment and Unfavourable Clinical Outcome in Glioma: An Integrated TCGA and CGGA Analysis.**

## Figure 1 – Immune Score and Overall Survival in TCGA Glioma

![Figure 1 – Immune Score and Overall Survival in TCGA Glioma](Figure1_FINAL_TCGA_immune_score_KM_publication.png)

**File:** `Figure1_FINAL_TCGA_immune_score_KM_publication.png`

Kaplan–Meier analysis of overall survival according to ESTIMATE-derived immune score in the combined TCGA glioma cohort.

High immune score was associated with significantly shorter overall survival (log-rank p < 0.0001).

---

## Figure 2 – Immune Checkpoint Expression According to Immune Score

![Figure 2 – Immune Checkpoint Expression](Figure2_Checkpoint_Boxplots.png)

**File:** `Figure2_Checkpoint_Boxplots.png`

Comparison of immune checkpoint gene expression between immune-high and immune-low tumors.

The analysis includes PDCD1, CD274, CTLA4, LAG3, TIGIT, and HAVCR2.

---

## Figure 3 – Correlations Between Immune Score and Immune Checkpoint Genes

![Figure 3 – Immune Score and Checkpoint Correlations](Figure3_Checkpoint_Correlations.png)

**File:** `Figure3_Checkpoint_Correlations.png`

Spearman correlations between ESTIMATE-derived immune score and immune checkpoint gene expression in TCGA glioma.

---

## Figure 4 – Association Between HAVCR2 Expression and Immune Score

![Figure 4 – HAVCR2 and ESTIMATE Immune Score](Figure4_REVISED_HAVCR2_ESTIMATE.png)

**File:** `Figure4_REVISED_HAVCR2_ESTIMATE.png`

Association between HAVCR2 expression and ESTIMATE-derived immune infiltration in TCGA glioma.

HAVCR2 showed a strong positive association with immune score (Spearman rho = 0.887; n = 601).

---

## Figure 5 – HAVCR2 Expression and Overall Survival in CGGA

![Figure 5 – CGGA HAVCR2 Survival Analysis](Figure5_CGGA_HAVCR2_KM.png)

**File:** `Figure5_CGGA_HAVCR2_KM.png`

Kaplan–Meier analysis of overall survival according to HAVCR2 expression in the independent CGGA mRNAseq_693 cohort.

Higher HAVCR2 expression was associated with shorter overall survival in univariable survival analysis (log-rank p < 0.0001).

---

## Figure 6 – HAVCR2 Expression According to Grade and IDH Status

![Figure 6 – HAVCR2 Expression According to Grade and IDH Status](Figure6_HAVCR2_Grade_IDH.png)

**File:** `Figure6_HAVCR2_Grade_IDH.png`

HAVCR2 expression according to tumor grade and IDH status in the CGGA validation cohort.

---

## Figure 7 – HAVCR2 and Macrophage-Associated Transcriptional Enrichment in TCGA

![Figure 7 – HAVCR2 and Macrophage Score](Figure7_HAVCR2_MacrophageScore.png)

**File:** `Figure7_HAVCR2_MacrophageScore.png`

Association between HAVCR2 expression and macrophage-associated transcriptional enrichment in TCGA glioma.

HAVCR2 showed a strong positive correlation with the macrophage marker score (Spearman rho = 0.905; n = 800).

---

## Figure 8 – CGGA Validation of the HAVCR2–Macrophage Association

![Figure 8 – CGGA HAVCR2 and Macrophage Association](Figure8_CGGA_HAVCR2_Macrophages.png)

**File:** `Figure8_CGGA_HAVCR2_Macrophages.png`

Independent validation of the association between HAVCR2 expression and macrophage-associated transcriptional enrichment in the CGGA mRNAseq_693 cohort.

HAVCR2 showed a strong positive correlation with the macrophage marker score (Spearman rho = 0.860).

---

## Reproducibility

The figures were generated using the R analysis scripts provided in the [`scripts/`](../scripts/) directory.

Statistical results and supporting output files are available in the [`results/`](../results/) directory.

Manuscript tables are available in the [`tables/`](../tables/) directory.

## Notes

Figures are provided as PNG files for direct viewing through GitHub.

All analyses were performed using publicly available, de-identified TCGA and CGGA datasets.
---

## Reproducibility

The figures were generated using the R analysis scripts provided in the [`scripts/`](../scripts/) directory.

Statistical results and supporting output files are available in the [`results/`](../results/) directory.

Manuscript tables are available in the [`tables/`](../tables/) directory.

## Notes

Figures are provided as PNG files for direct viewing through GitHub.

All analyses were performed using publicly available, de-identified TCGA and CGGA datasets.
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
