# ============================================================
# Revision analysis:
# HAVCR2 association with macrophage and microglia marker scores
# stratified by TCGA-LGG and TCGA-GBM
#
# Addresses Reviewer 1 Comment 2 and Reviewer 2 Comments 2-4
#
# Input:
#   C:/GDC/TCGA_LGG_GBM_STAR_counts.rds
#
# Outputs:
#   results/HAVCR2_macrophage_microglia_LGG_GBM.csv
#   results/macrophage_microglia_LGG_GBM_patient_level.csv
# ============================================================


# ------------------------------------------------------------
# 1. Required package
# ------------------------------------------------------------

if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  stop("Install SummarizedExperiment first.")
}

library(SummarizedExperiment)


# ------------------------------------------------------------
# 2. Load prepared TCGA object
# ------------------------------------------------------------

se <- readRDS(
  "C:/GDC/TCGA_LGG_GBM_STAR_counts.rds"
)

cat("\nSummarizedExperiment dimensions:\n")
print(dim(se))


# ------------------------------------------------------------
# 3. Extract TPM expression
# ------------------------------------------------------------

expr <- SummarizedExperiment::assay(
  se,
  "tpm_unstrand"
)

rd <- as.data.frame(
  SummarizedExperiment::rowData(se)
)

if (!"gene_name" %in% colnames(rd)) {
  stop("gene_name not found in rowData(se).")
}

gene_symbols <- as.character(
  rd$gene_name
)


# ------------------------------------------------------------
# 4. Clean expression matrix
# ------------------------------------------------------------

keep <- (
  !is.na(gene_symbols) &
    gene_symbols != "" &
    rowSums(expr, na.rm = TRUE) > 0
)

expr <- expr[
  keep,
  ,
  drop = FALSE
]

gene_symbols <- gene_symbols[
  keep
]

rownames(expr) <- gene_symbols

expr <- rowsum(
  expr,
  group = rownames(expr),
  reorder = FALSE
)

expr_log <- log2(
  expr + 1
)

cat(
  "\nGenes after filtering/collapsing:",
  nrow(expr_log),
  "\n"
)

cat(
  "Samples:",
  ncol(expr_log),
  "\n"
)


# ------------------------------------------------------------
# 5. Check HAVCR2
# ------------------------------------------------------------

if (!"HAVCR2" %in% rownames(expr_log)) {
  stop("HAVCR2 not found.")
}


# ------------------------------------------------------------
# 6. Define marker sets
# ------------------------------------------------------------

macrophage_markers <- c(
  "AIF1",
  "CD68",
  "CD163",
  "MSR1",
  "MRC1",
  "CSF1R",
  "ITGAM",
  "LST1",
  "TYROBP"
)

microglia_markers <- c(
  "TMEM119",
  "P2RY12",
  "CX3CR1",
  "SALL1"
)

available_macrophage <- intersect(
  macrophage_markers,
  rownames(expr_log)
)

available_microglia <- intersect(
  microglia_markers,
  rownames(expr_log)
)

cat("\nAvailable macrophage markers:\n")
print(available_macrophage)

cat("\nAvailable microglia markers:\n")
print(available_microglia)

if (length(available_macrophage) < 2) {
  stop("Too few macrophage markers available.")
}

if (length(available_microglia) < 2) {
  stop("Too few microglia markers available.")
}


# ------------------------------------------------------------
# 7. Sample-level scores
# ------------------------------------------------------------

havcr2_sample <- as.numeric(
  expr_log[
    "HAVCR2",
  ]
)

names(havcr2_sample) <- colnames(expr_log)

macrophage_score_sample <- colMeans(
  expr_log[
    available_macrophage,
    ,
    drop = FALSE
  ],
  na.rm = TRUE
)

microglia_score_sample <- colMeans(
  expr_log[
    available_microglia,
    ,
    drop = FALSE
  ],
  na.rm = TRUE
)


# ------------------------------------------------------------
# 8. Obtain project information
# ------------------------------------------------------------

cd <- as.data.frame(
  SummarizedExperiment::colData(se)
)

sample_ids <- colnames(expr_log)

cd <- cd[
  match(
    sample_ids,
    rownames(cd)
  ),
  ,
  drop = FALSE
]

project_candidates <- c(
  "project_id",
  "project.project_id",
  "project"
)

project_column <- project_candidates[
  project_candidates %in% colnames(cd)
]

if (length(project_column) == 0) {
  
  project_column <- grep(
    "project.*id|project_id",
    colnames(cd),
    ignore.case = TRUE,
    value = TRUE
  )
}

if (length(project_column) == 0) {
  stop("No project ID column found.")
}

project_column <- project_column[1]

project_id <- as.character(
  cd[[project_column]]
)

cat(
  "\nUsing project column:",
  project_column,
  "\n"
)

print(
  table(
    project_id,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# 9. Construct patient-level dataset
# ------------------------------------------------------------

patient_id <- substr(
  sample_ids,
  1,
  12
)

sample_data <- data.frame(
  sample_id = sample_ids,
  patient_id = patient_id,
  project_id = project_id,
  HAVCR2 = havcr2_sample,
  macrophage_score = macrophage_score_sample,
  microglia_score = microglia_score_sample,
  stringsAsFactors = FALSE
)

sample_data <- sample_data[
  sample_data$project_id %in%
    c(
      "TCGA-LGG",
      "TCGA-GBM"
    ),
  ,
  drop = FALSE
]

cat(
  "\nSample counts:\n"
)

print(
  table(
    sample_data$project_id
  )
)


# ------------------------------------------------------------
# 10. Collapse to patient level using median
# ------------------------------------------------------------

patients <- unique(
  sample_data$patient_id
)

patient_data <- do.call(
  rbind,
  lapply(
    patients,
    function(p) {
      
      tmp <- sample_data[
        sample_data$patient_id == p,
        ,
        drop = FALSE
      ]
      
      data.frame(
        patient_id = p,
        project_id = unique(tmp$project_id)[1],
        HAVCR2 = median(
          tmp$HAVCR2,
          na.rm = TRUE
        ),
        macrophage_score = median(
          tmp$macrophage_score,
          na.rm = TRUE
        ),
        microglia_score = median(
          tmp$microglia_score,
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

cat(
  "\nPatient-level counts:\n"
)

print(
  table(
    patient_data$project_id
  )
)

write.csv(
  patient_data,
  file =
    "results/macrophage_microglia_LGG_GBM_patient_level.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Helper function for Spearman correlations
# ------------------------------------------------------------

run_test <- function(
    cohort_name,
    x_name,
    y_name,
    data
) {
  
  d <- data[
    data$project_id == cohort_name,
    ,
    drop = FALSE
  ]
  
  x <- d[[x_name]]
  y <- d[[y_name]]
  
  ok <- complete.cases(
    x,
    y
  )
  
  test <- suppressWarnings(
    cor.test(
      x[ok],
      y[ok],
      method = "spearman",
      exact = FALSE
    )
  )
  
  data.frame(
    cohort = cohort_name,
    comparison = paste(
      x_name,
      "vs",
      y_name
    ),
    rho = unname(
      test$estimate
    ),
    p_value = test$p.value,
    n = sum(ok),
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 12. Run LGG correlations
# ------------------------------------------------------------

lgg_havcr2_macro <- run_test(
  "TCGA-LGG",
  "HAVCR2",
  "macrophage_score",
  patient_data
)

lgg_havcr2_micro <- run_test(
  "TCGA-LGG",
  "HAVCR2",
  "microglia_score",
  patient_data
)

lgg_macro_micro <- run_test(
  "TCGA-LGG",
  "macrophage_score",
  "microglia_score",
  patient_data
)


# ------------------------------------------------------------
# 13. Run GBM correlations
# ------------------------------------------------------------

gbm_havcr2_macro <- run_test(
  "TCGA-GBM",
  "HAVCR2",
  "macrophage_score",
  patient_data
)

gbm_havcr2_micro <- run_test(
  "TCGA-GBM",
  "HAVCR2",
  "microglia_score",
  patient_data
)

gbm_macro_micro <- run_test(
  "TCGA-GBM",
  "macrophage_score",
  "microglia_score",
  patient_data
)


# ------------------------------------------------------------
# 14. Combine results
# ------------------------------------------------------------

results <- rbind(
  lgg_havcr2_macro,
  lgg_havcr2_micro,
  lgg_macro_micro,
  gbm_havcr2_macro,
  gbm_havcr2_micro,
  gbm_macro_micro
)

results$FDR <- p.adjust(
  results$p_value,
  method = "BH"
)

write.csv(
  results,
  file =
    "results/HAVCR2_macrophage_microglia_LGG_GBM.csv",
  row.names = FALSE
)

cat(
  "\n====================================\n"
)

cat(
  "STRATIFIED MACROPHAGE / MICROGLIA RESULTS\n"
)

cat(
  "====================================\n"
)

print(
  results
)


# ------------------------------------------------------------
# 15. Direct comparison table for manuscript
# ------------------------------------------------------------

manuscript_table <- data.frame(
  comparison = c(
    "HAVCR2 vs macrophage score",
    "HAVCR2 vs microglia score",
    "Macrophage score vs microglia score"
  ),
  
  rho_LGG = c(
    lgg_havcr2_macro$rho,
    lgg_havcr2_micro$rho,
    lgg_macro_micro$rho
  ),
  
  p_LGG = c(
    lgg_havcr2_macro$p_value,
    lgg_havcr2_micro$p_value,
    lgg_macro_micro$p_value
  ),
  
  n_LGG = c(
    lgg_havcr2_macro$n,
    lgg_havcr2_micro$n,
    lgg_macro_micro$n
  ),
  
  rho_GBM = c(
    gbm_havcr2_macro$rho,
    gbm_havcr2_micro$rho,
    gbm_macro_micro$rho
  ),
  
  p_GBM = c(
    gbm_havcr2_macro$p_value,
    gbm_havcr2_micro$p_value,
    gbm_macro_micro$p_value
  ),
  
  n_GBM = c(
    gbm_havcr2_macro$n,
    gbm_havcr2_micro$n,
    gbm_macro_micro$n
  )
)

write.csv(
  manuscript_table,
  file =
    "results/Table_revision_HAVCR2_macrophage_microglia_LGG_GBM.csv",
  row.names = FALSE
)

cat(
  "\nDirect LGG vs GBM comparison:\n"
)

print(
  manuscript_table
)


# ------------------------------------------------------------
# 16. Save session information
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file =
    "results/sessionInfo_macrophage_microglia_LGG_GBM.txt"
)


cat(
  "\n============================================\n"
)

cat(
  "MACROPHAGE/MICROGLIA ANALYSIS COMPLETED\n"
)

cat(
  "============================================\n"
)

cat(
  "Results written to results/ folder.\n"
)
