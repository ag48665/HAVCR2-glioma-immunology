# HAVCR2 Is Associated With Immune and Macrophage Enrichment and Unfavourable Clinical Outcome in Glioma: An Integrated TCGA and CGGA Analysis

## Overview

This repository contains the analysis code, figures, and statistical results associated with an integrated analysis of HAVCR2 (TIM-3), immune infiltration, and myeloid-cell enrichment in diffuse glioma.

Transcriptomic and clinical data from The Cancer Genome Atlas (TCGA) lower-grade glioma (LGG) and glioblastoma (GBM) cohorts were analyzed, with key findings evaluated in the independent Chinese Glioma Genome Atlas (CGGA) mRNAseq_693 cohort.

The study evaluates relationships among immune infiltration, immune checkpoint expression, HAVCR2 expression, clinical outcome, myeloid-cell enrichment, macrophage- and microglia-associated transcriptional signatures, HAVCR2-centered co-expression, and immune-related pathway enrichment.

## Study Highlights

* High ESTIMATE immune scores were associated with significantly shorter overall survival in the combined TCGA glioma cohort.
* Higher immune score remained associated with worse overall survival after adjustment for age and tumor grade and stratification by IDH status.
* Immune-high tumors showed increased expression of multiple immune checkpoint genes, including PDCD1, CD274, CTLA4, LAG3, TIGIT, and HAVCR2.
* HAVCR2 showed a strong positive association with ESTIMATE-derived immune infiltration.
* Higher HAVCR2 expression was associated with shorter overall survival in univariable analyses in the independent CGGA cohort.
* HAVCR2 was not independently associated with overall survival after adjustment for age and tumor grade and stratification by IDH status in CGGA.
* HAVCR2 expression was strongly associated with computationally estimated myeloid-cell enrichment.
* Strong associations between HAVCR2 and macrophage-associated transcriptional enrichment were reproduced in both TCGA and CGGA.
* HAVCR2 was associated with macrophage- and microglia-related transcriptional signals in both LGG and GBM.
* HAVCR2-centered co-expression analyses identified broadly conserved immune- and myeloid-associated transcriptional programs.
* HAVCR2-associated transcriptional profiles were enriched for inflammatory and immune-related Hallmark pathways.

## Datasets

### TCGA

* TCGA-LGG (Lower Grade Glioma)
* TCGA-GBM (Glioblastoma)

### CGGA

* CGGA mRNAseq_693 cohort

All analyses were performed using publicly available, de-identified datasets.

## Repository Structure

```text
figures/      Manuscript figures
results/      Statistical results and summary tables
scripts/      R analysis scripts
