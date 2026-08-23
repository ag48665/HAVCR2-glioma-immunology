# ============================================================
# Revision analysis:
# xCell2 HAVCR2 correlations stratified by TCGA-LGG / TCGA-GBM
#
# Addresses concern that pooled LGG+GBM associations may be
# driven by differences between tumor entities.
#
# Inputs:
#   C:/GDC/TCGA_LGG_GBM_STAR_counts.rds
#   results/xCell2_TCGA_LGG_GBM_scores.rds
#
# Outputs:
#   results/xCell2_HAVCR2_LGG_GBM_stratified_all.csv
#   results/xCell2_HAVCR2_LGG_stratified.csv
#   results/xCell2_HAVCR2_GBM_stratified.csv
#   results/xCell2_HAVCR2_LGG_GBM_myeloid.csv
# ============================================================


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

required <- c(
  "SummarizedExperiment"
)

missing <- required[
  !vapply(
    required,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing)) {
  stop(
    "Install required package(s): ",
    paste(missing, collapse = ", ")
  )
}

library(SummarizedExperiment)


# ------------------------------------------------------------
# 2. Load data
# ------------------------------------------------------------

se <- readRDS(
  "C:/GDC/TCGA_LGG_GBM_STAR_counts.rds"
)

xcell_scores <- readRDS(
  "results/xCell2_TCGA_LGG_GBM_scores.rds"
)

cat("\nSummarizedExperiment dimensions:\n")
print(dim(se))

cat("\nxCell2 dimensions:\n")
print(dim(xcell_scores))


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
  stop(
    "gene_name not found in rowData(se)."
  )
}

gene_symbols <- as.character(
  rd$gene_name
)

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

if (!"HAVCR2" %in% rownames(expr_log)) {
  stop(
    "HAVCR2 not found in expression matrix."
  )
}


# ------------------------------------------------------------
# 4. Confirm sample alignment
# ------------------------------------------------------------

common_samples <- intersect(
  colnames(expr_log),
  colnames(xcell_scores)
)

cat(
  "\nSamples shared between expression and xCell2:",
  length(common_samples),
  "\n"
)

if (length(common_samples) < 100) {
  stop(
    "Too few matching samples between expression and xCell2."
  )
}

expr_log <- expr_log[
  ,
  common_samples,
  drop = FALSE
]

xcell_scores <- xcell_scores[
  ,
  common_samples,
  drop = FALSE
]


# ------------------------------------------------------------
# 5. Obtain project information
# ------------------------------------------------------------

cd <- as.data.frame(
  SummarizedExperiment::colData(se)
)

cat(
  "\nAvailable colData columns:\n"
)

print(
  colnames(cd)
)


# Match colData rows to expression sample columns
cd <- cd[
  match(
    common_samples,
    rownames(cd)
  ),
  ,
  drop = FALSE
]

if (any(is.na(rownames(cd)))) {
  stop(
    "Could not match all samples to colData."
  )
}


# Look for a project-ID column.
project_candidates <- c(
  "project_id",
  "project.project_id",
  "project"
)

project_column <- project_candidates[
  project_candidates %in% colnames(cd)
]

if (length(project_column) == 0) {
  
  # Try automatic pattern detection
  project_column <- grep(
    "project.*id|project_id",
    colnames(cd),
    ignore.case = TRUE,
    value = TRUE
  )
}

if (length(project_column) == 0) {
  stop(
    paste(
      "No project ID column found.",
      "Inspect colnames(colData(se))."
    )
  )
}

project_column <- project_column[1]

cat(
  "\nUsing project column:",
  project_column,
  "\n"
)

project_id <- as.character(
  cd[[project_column]]
)

cat(
  "\nProject IDs found:\n"
)

print(
  table(
    project_id,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# 6. Construct patient IDs from TCGA barcode
# ------------------------------------------------------------

sample_ids <- common_samples

# TCGA patient barcode = first 12 characters:
# TCGA-XX-YYYY
patient_id <- substr(
  sample_ids,
  1,
  12
)

metadata <- data.frame(
  sample_id = sample_ids,
  patient_id = patient_id,
  project_id = project_id,
  stringsAsFactors = FALSE
)

metadata <- metadata[
  metadata$project_id %in%
    c(
      "TCGA-LGG",
      "TCGA-GBM"
    ),
  ,
  drop = FALSE
]

cat(
  "\nSamples retained by project:\n"
)

print(
  table(
    metadata$project_id
  )
)

cat(
  "\nUnique patients retained by project:\n"
)

print(
  table(
    metadata$project_id[
      !duplicated(metadata$patient_id)
    ]
  )
)


# ------------------------------------------------------------
# 7. Restrict expression and xCell to retained samples
# ------------------------------------------------------------

keep_samples <- metadata$sample_id

havcr2_sample <- as.numeric(
  expr_log[
    "HAVCR2",
    keep_samples
  ]
)

names(havcr2_sample) <- keep_samples

xcell_sample <- xcell_scores[
  ,
  keep_samples,
  drop = FALSE
]


# ------------------------------------------------------------
# 8. Collapse multiple samples from same patient
#
# Median is used so each patient contributes only once.
# ------------------------------------------------------------

patients <- unique(
  metadata$patient_id
)

havcr2_patient <- numeric(
  length(patients)
)

names(havcr2_patient) <- patients

xcell_patient <- matrix(
  NA_real_,
  nrow = nrow(xcell_sample),
  ncol = length(patients),
  dimnames = list(
    rownames(xcell_sample),
    patients
  )
)

patient_project <- character(
  length(patients)
)

names(patient_project) <- patients


for (p in patients) {
  
  samp <- metadata$sample_id[
    metadata$patient_id == p
  ]
  
  havcr2_patient[p] <- median(
    havcr2_sample[samp],
    na.rm = TRUE
  )
  
  xcell_patient[
    ,
    p
  ] <- apply(
    xcell_sample[
      ,
      samp,
      drop = FALSE
    ],
    1,
    median,
    na.rm = TRUE
  )
  
  patient_project[p] <- unique(
    metadata$project_id[
      metadata$patient_id == p
    ]
  )[1]
}


cat(
  "\nPatient-level dataset:\n"
)

print(
  table(
    patient_project
  )
)


# ------------------------------------------------------------
# 9. Function for stratified correlations
# ------------------------------------------------------------

run_correlations <- function(
    project_name
) {
  
  patient_keep <- names(
    patient_project
  )[
    patient_project == project_name
  ]
  
  cat(
    "\nRunning correlations for",
    project_name,
    "n =",
    length(patient_keep),
    "\n"
  )
  
  h <- havcr2_patient[
    patient_keep
  ]
  
  xs <- xcell_patient[
    ,
    patient_keep,
    drop = FALSE
  ]
  
  results <- lapply(
    rownames(xs),
    function(celltype) {
      
      score <- as.numeric(
        xs[
          celltype,
        ]
      )
      
      ok <- complete.cases(
        h,
        score
      )
      
      if (sum(ok) < 10) {
        return(NULL)
      }
      
      test <- suppressWarnings(
        cor.test(
          h[ok],
          score[ok],
          method = "spearman",
          exact = FALSE
        )
      )
      
      data.frame(
        cohort = project_name,
        cell_type = celltype,
        rho = unname(
          test$estimate
        ),
        p_value = test$p.value,
        n = sum(ok),
        stringsAsFactors = FALSE
      )
    }
  )
  
  results <- do.call(
    rbind,
    results
  )
  
  results$FDR <- p.adjust(
    results$p_value,
    method = "BH"
  )
  
  results <- results[
    order(
      -abs(results$rho)
    ),
    ,
    drop = FALSE
  ]
  
  results
}


# ------------------------------------------------------------
# 10. LGG analysis
# ------------------------------------------------------------

res_lgg_xcell <- run_correlations(
  "TCGA-LGG"
)

cat(
  "\nTop LGG associations:\n"
)

print(
  head(
    res_lgg_xcell,
    15
  )
)


# ------------------------------------------------------------
# 11. GBM analysis
# ------------------------------------------------------------

res_gbm_xcell <- run_correlations(
  "TCGA-GBM"
)

cat(
  "\nTop GBM associations:\n"
)

print(
  head(
    res_gbm_xcell,
    15
  )
)


# ------------------------------------------------------------
# 12. Combine
# ------------------------------------------------------------

stratified_xcell <- rbind(
  res_lgg_xcell,
  res_gbm_xcell
)

write.csv(
  stratified_xcell,
  file =
    "results/xCell2_HAVCR2_LGG_GBM_stratified_all.csv",
  row.names = FALSE
)

write.csv(
  res_lgg_xcell,
  file =
    "results/xCell2_HAVCR2_LGG_stratified.csv",
  row.names = FALSE
)

write.csv(
  res_gbm_xcell,
  file =
    "results/xCell2_HAVCR2_GBM_stratified.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. Extract myeloid/macrophage-related results
# ------------------------------------------------------------

myeloid_pattern <- paste(
  c(
    "myeloid",
    "macroph",
    "monocyte",
    "microgl",
    "dendritic"
  ),
  collapse = "|"
)

myeloid_stratified <- stratified_xcell[
  grepl(
    myeloid_pattern,
    stratified_xcell$cell_type,
    ignore.case = TRUE
  ),
  ,
  drop = FALSE
]

write.csv(
  myeloid_stratified,
  file =
    "results/xCell2_HAVCR2_LGG_GBM_myeloid.csv",
  row.names = FALSE
)

cat(
  "\nStratified myeloid results:\n"
)

print(
  myeloid_stratified
)


# ------------------------------------------------------------
# 14. Direct LGG vs GBM summary table
# ------------------------------------------------------------

lgg_summary <- res_lgg_xcell[
  grepl(
    myeloid_pattern,
    res_lgg_xcell$cell_type,
    ignore.case = TRUE
  ),
  c(
    "cell_type",
    "rho",
    "p_value",
    "FDR",
    "n"
  )
]

gbm_summary <- res_gbm_xcell[
  grepl(
    myeloid_pattern,
    res_gbm_xcell$cell_type,
    ignore.case = TRUE
  ),
  c(
    "cell_type",
    "rho",
    "p_value",
    "FDR",
    "n"
  )
]

colnames(lgg_summary)[
  -1
] <- paste0(
  colnames(lgg_summary)[-1],
  "_LGG"
)

colnames(gbm_summary)[
  -1
] <- paste0(
  colnames(gbm_summary)[-1],
  "_GBM"
)

comparison_table <- merge(
  lgg_summary,
  gbm_summary,
  by = "cell_type",
  all = TRUE
)

write.csv(
  comparison_table,
  file =
    "results/xCell2_HAVCR2_LGG_vs_GBM_myeloid_comparison.csv",
  row.names = FALSE
)

cat(
  "\nLGG vs GBM myeloid comparison:\n"
)

print(
  comparison_table
)


# ------------------------------------------------------------
# 15. Save patient-level data for reproducibility
# ------------------------------------------------------------

saveRDS(
  list(
    HAVCR2 = havcr2_patient,
    xCell2 = xcell_patient,
    project = patient_project
  ),
  file =
    "results/xCell2_patient_level_revision_data.rds"
)


# ------------------------------------------------------------
# 16. Session information
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file =
    "results/sessionInfo_xCell2_LGG_GBM_stratified.txt"
)


cat(
  "\n=============================================\n"
)

cat(
  "xCell2 LGG/GBM STRATIFIED ANALYSIS COMPLETED\n"
)

cat(
  "=============================================\n"
)

cat(
  "Results written to results/ folder.\n"
)
