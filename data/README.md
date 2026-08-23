\# Data



This directory documents the datasets used in the study.



\## TCGA



Transcriptomic and clinical data from the TCGA-LGG and TCGA-GBM cohorts were obtained from the Genomic Data Commons (GDC).



The analyses use publicly available, de-identified data. Large downloaded TCGA files and intermediate R objects are not stored in this repository.



The TCGA download and preprocessing workflow is provided in:



`../scripts/01\_download\_tcga.R`



\## CGGA



Independent validation was performed using the CGGA mRNAseq\_693 cohort obtained from the Chinese Glioma Genome Atlas (CGGA).



CGGA source data are not redistributed in this repository. They should be obtained directly from the CGGA data portal in accordance with the database access conditions.



\## Intermediate Data



Intermediate objects, including processed expression matrices, STAR count matrices, and RDS/RData files, were generated locally during the analyses.



These large derived or intermediate files are not included in the repository and can be regenerated from the publicly available source datasets using the analysis scripts.



No new patient-level datasets were generated in this study.

