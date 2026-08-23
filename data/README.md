# Data

This directory documents the datasets used in the study and provides information required to reproduce the analyses.

## TCGA

Transcriptomic and clinical data from the TCGA Lower Grade Glioma (TCGA-LGG) and Glioblastoma Multiforme (TCGA-GBM) cohorts were obtained from the Genomic Data Commons (GDC).

The analyses use publicly available, de-identified data. Large downloaded TCGA files and intermediate R objects are not stored in this repository.

The TCGA download and preprocessing workflow is provided in:

`../scripts/01_download_tcga.R`

## CGGA

Independent validation was performed using the CGGA mRNAseq_693 cohort obtained from the Chinese Glioma Genome Atlas (CGGA).

CGGA source data are not redistributed in this repository. They should be obtained directly from the CGGA data portal in accordance with the database access conditions.

## Intermediate Data

Intermediate objects, including processed expression matrices, STAR count matrices, and RDS/RData files, were generated locally during the analyses.

These large derived or intermediate files are not included in this repository. They can be regenerated from the publicly available source datasets using the analysis scripts provided in the `scripts/` directory.

## Data Availability

All source datasets used in this study are publicly available.

TCGA data are available through the Genomic Data Commons Data Portal.

CGGA data are available through the Chinese Glioma Genome Atlas data portal.

No new patient-level datasets were generated during this study, and no individual-level source data are redistributed through this repository.
