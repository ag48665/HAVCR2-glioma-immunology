# ============================================================
# 10_FINAL_xCell2_HAVCR2_validation.R
#
# FINAL VALIDATION
# HAVCR2 vs xCell2 cell-type enrichment
#
# Analyses:
#   1. Combined TCGA-LGG + TCGA-GBM
#   2. TCGA-LGG separately
#   3. TCGA-GBM separately
#
# Patient-level analysis
# Spearman correlations
# Benjamini-Hochberg FDR
#
# INPUTS:
#   C:/GDC/TCGA_LGG_GBM_STAR_counts.rds
#   results/xCell2_TCGA_LGG_GBM_scores.rds
#
# ============================================================


# ------------------------------------------------------------
# 1. Required package
# ------------------------------------------------------------

if (!requireNamespace(
  "SummarizedExperiment",
  quietly = TRUE
)) {
  stop(
    "Install SummarizedExperiment first."
  )
}

library(SummarizedExperiment)


# ------------------------------------------------------------
# 2. Input files
# ------------------------------------------------------------

se_file <- "C:/GDC/TCGA_LGG_GBM_STAR_counts.rds"

xcell_file <- file.path(
  "results",
  "xCell2_TCGA_LGG_GBM_scores.rds"
)


if (!file.exists(se_file)) {
  stop(
    "STAR-counts RDS not found: ",
    se_file
  )
}

if (!file.exists(xcell_file)) {
  stop(
    "xCell2 score file not found: ",
    xcell_file
  )
}


cat("\n=============================================\n")
cat("FINAL xCell2 HAVCR2 VALIDATION\n")
cat("=============================================\n")


# ------------------------------------------------------------
# 3. Load data
# ------------------------------------------------------------

se <- readRDS(
  se_file
)

xcell_scores <- readRDS(
  xcell_file
)


cat(
  "\nSummarizedExperiment dimensions:\n"
)

print(
  dim(se)
)


cat(
  "\nxCell2 score dimensions:\n"
)

print(
  dim(xcell_scores)
)


# ------------------------------------------------------------
# 4. Extract TPM expression
# ------------------------------------------------------------

if (!"tpm_unstrand" %in% assayNames(se)) {
  stop(
    "tpm_unstrand assay not found."
  )
}

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


# ------------------------------------------------------------
# 5. Clean genes
# ------------------------------------------------------------

keep <- (
  !is.na(gene_symbols) &
    gene_symbols != "" &
    rowSums(
      expr,
      na.rm = TRUE
    ) > 0
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


# ------------------------------------------------------------
# 6. Collapse duplicate symbols
# ------------------------------------------------------------

expr <- rowsum(
  expr,
  group = rownames(expr),
  reorder = FALSE
)


# ------------------------------------------------------------
# 7. Transform expression
# ------------------------------------------------------------

expr_log <- log2(
  expr + 1
)

if (!"HAVCR2" %in% rownames(expr_log)) {
  stop(
    "HAVCR2 not found."
  )
}


# ------------------------------------------------------------
# 8. Align expression and xCell2 samples
# ------------------------------------------------------------

common_samples <- intersect(
  colnames(expr_log),
  colnames(xcell_scores)
)

cat(
  "\nShared expression/xCell2 samples:",
  length(common_samples),
  "\n"
)

if (length(common_samples) < 100) {
  stop(
    "Too few matching samples."
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
# 9. Metadata
# ------------------------------------------------------------

cd <- as.data.frame(
  SummarizedExperiment::colData(se)
)

matched <- match(
  common_samples,
  rownames(cd)
)

if (any(is.na(matched))) {
  stop(
    "Some samples could not be matched to colData."
  )
}

cd <- cd[
  matched,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 10. Find project column
# ------------------------------------------------------------

project_candidates <- c(
  "project_id",
  "project.project_id",
  "project"
)

project_column <- project_candidates[
  project_candidates %in%
    colnames(cd)
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
  stop(
    "Project ID column not found."
  )
}

project_column <- project_column[1]

project_id <- as.character(
  cd[[project_column]]
)


cat(
  "\nProject IDs:\n"
)

print(
  table(
    project_id,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# 11. Patient IDs
# ------------------------------------------------------------

sample_ids <- common_samples

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


# ------------------------------------------------------------
# 12. Keep only LGG / GBM
# ------------------------------------------------------------

metadata <- metadata[
  metadata$project_id %in%
    c(
      "TCGA-LGG",
      "TCGA-GBM"
    ),
  ,
  drop = FALSE
]


keep_samples <- metadata$sample_id


cat(
  "\nSamples retained:\n"
)

print(
  table(
    metadata$project_id
  )
)


cat(
  "\nUnique patients retained:\n"
)

unique_patient_table <- metadata[
  !duplicated(metadata$patient_id),
  ,
  drop = FALSE
]

print(
  table(
    unique_patient_table$project_id
  )
)


# ------------------------------------------------------------
# 13. HAVCR2 sample-level values
# ------------------------------------------------------------

havcr2_sample <- as.numeric(
  expr_log[
    "HAVCR2",
    keep_samples
  ]
)

names(havcr2_sample) <- keep_samples


# ------------------------------------------------------------
# 14. xCell2 sample-level matrix
# ------------------------------------------------------------

xcell_sample <- xcell_scores[
  ,
  keep_samples,
  drop = FALSE
]


# ------------------------------------------------------------
# 15. Collapse to patient level
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


cat(
  "\nCollapsing to patient level...\n"
)


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
  
  projects_p <- unique(
    metadata$project_id[
      metadata$patient_id == p
    ]
  )
  
  if (length(projects_p) != 1) {
    stop(
      paste(
        "Inconsistent project for patient",
        p
      )
    )
  }
  
  patient_project[p] <- projects_p
}


# ------------------------------------------------------------
# 16. Patient cohort sizes
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("PATIENT-LEVEL COHORTS\n")
cat("=============================================\n")

print(
  table(
    patient_project
  )
)

cat(
  "\nTotal patients:",
  length(patients),
  "\n"
)


# ------------------------------------------------------------
# 17. Correlation function
# ------------------------------------------------------------

run_correlations <- function(
    patient_keep,
    cohort_name
) {
  
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
          ,
          drop = TRUE
        ]
      )
      
      ok <- complete.cases(
        h,
        score
      )
      
      if (
        sum(ok) < 10 ||
        length(unique(score[ok])) < 2
      ) {
        return(NULL)
      }
      
      ct <- suppressWarnings(
        cor.test(
          h[ok],
          score[ok],
          method = "spearman",
          exact = FALSE
        )
      )
      
      data.frame(
        cohort = cohort_name,
        cell_type = celltype,
        rho = unname(
          ct$estimate
        ),
        p_value = ct$p.value,
        n = sum(ok),
        stringsAsFactors = FALSE
      )
    }
  )
  
  results <- Filter(
    Negate(is.null),
    results
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
      -results$rho
    ),
    ,
    drop = FALSE
  ]
  
  rownames(results) <- NULL
  
  results
}


# ------------------------------------------------------------
# 18. POOLED analysis
# ------------------------------------------------------------

pooled_patients <- names(
  patient_project
)

res_pooled <- run_correlations(
  pooled_patients,
  "TCGA-LGG+GBM"
)


# ------------------------------------------------------------
# 19. LGG analysis
# ------------------------------------------------------------

lgg_patients <- names(
  patient_project
)[
  patient_project == "TCGA-LGG"
]

res_lgg <- run_correlations(
  lgg_patients,
  "TCGA-LGG"
)


# ------------------------------------------------------------
# 20. GBM analysis
# ------------------------------------------------------------

gbm_patients <- names(
  patient_project
)[
  patient_project == "TCGA-GBM"
]

res_gbm <- run_correlations(
  gbm_patients,
  "TCGA-GBM"
)


# ------------------------------------------------------------
# 21. Save complete results
# ------------------------------------------------------------

dir.create(
  "results",
  showWarnings = FALSE,
  recursive = TRUE
)


write.csv(
  res_pooled,
  "results/xCell2_HAVCR2_POOLED_FINAL.csv",
  row.names = FALSE
)

write.csv(
  res_lgg,
  "results/xCell2_HAVCR2_LGG_FINAL.csv",
  row.names = FALSE
)

write.csv(
  res_gbm,
  "results/xCell2_HAVCR2_GBM_FINAL.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 22. Print strongest associations
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("TOP xCell2 ASSOCIATIONS - POOLED\n")
cat("=============================================\n\n")

print(
  head(
    res_pooled,
    20
  ),
  digits = 10
)


cat("\n=============================================\n")
cat("TOP xCell2 ASSOCIATIONS - LGG\n")
cat("=============================================\n\n")

print(
  head(
    res_lgg,
    20
  ),
  digits = 10
)


cat("\n=============================================\n")
cat("TOP xCell2 ASSOCIATIONS - GBM\n")
cat("=============================================\n\n")

print(
  head(
    res_gbm,
    20
  ),
  digits = 10
)


# ------------------------------------------------------------
# 23. Extract myeloid-related cell types
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


pooled_myeloid <- res_pooled[
  grepl(
    myeloid_pattern,
    res_pooled$cell_type,
    ignore.case = TRUE
  ),
  ,
  drop = FALSE
]


lgg_myeloid <- res_lgg[
  grepl(
    myeloid_pattern,
    res_lgg$cell_type,
    ignore.case = TRUE
  ),
  ,
  drop = FALSE
]


gbm_myeloid <- res_gbm[
  grepl(
    myeloid_pattern,
    res_gbm$cell_type,
    ignore.case = TRUE
  ),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 24. Print MYELOID results
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("MYELOID xCell2 RESULTS - POOLED\n")
cat("=============================================\n\n")

print(
  pooled_myeloid,
  digits = 10
)


cat("\n=============================================\n")
cat("MYELOID xCell2 RESULTS - LGG\n")
cat("=============================================\n\n")

print(
  lgg_myeloid,
  digits = 10
)


cat("\n=============================================\n")
cat("MYELOID xCell2 RESULTS - GBM\n")
cat("=============================================\n\n")

print(
  gbm_myeloid,
  digits = 10
)


# ------------------------------------------------------------
# 25. Save myeloid results
# ------------------------------------------------------------

write.csv(
  pooled_myeloid,
  "results/xCell2_HAVCR2_myeloid_POOLED_FINAL.csv",
  row.names = FALSE
)

write.csv(
  lgg_myeloid,
  "results/xCell2_HAVCR2_myeloid_LGG_FINAL.csv",
  row.names = FALSE
)

write.csv(
  gbm_myeloid,
  "results/xCell2_HAVCR2_myeloid_GBM_FINAL.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 26. Combined comparison table
# ------------------------------------------------------------

prepare_compare <- function(
    x,
    suffix
) {
  
  out <- x[
    ,
    c(
      "cell_type",
      "rho",
      "p_value",
      "FDR",
      "n"
    )
  ]
  
  colnames(out)[-1] <- paste0(
    colnames(out)[-1],
    "_",
    suffix
  )
  
  out
}


pooled_comp <- prepare_compare(
  pooled_myeloid,
  "POOLED"
)

lgg_comp <- prepare_compare(
  lgg_myeloid,
  "LGG"
)

gbm_comp <- prepare_compare(
  gbm_myeloid,
  "GBM"
)


comparison <- merge(
  pooled_comp,
  lgg_comp,
  by = "cell_type",
  all = TRUE
)

comparison <- merge(
  comparison,
  gbm_comp,
  by = "cell_type",
  all = TRUE
)


cat("\n=============================================\n")
cat("POOLED vs LGG vs GBM MYELOID COMPARISON\n")
cat("=============================================\n\n")

print(
  comparison,
  digits = 10
)


write.csv(
  comparison,
  "results/xCell2_HAVCR2_myeloid_POOLED_LGG_GBM_FINAL.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 27. Check values currently reported in manuscript
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("MANUSCRIPT SECTION 3.5 CHECK\n")
cat("=============================================\n\n")

cat(
  "Current manuscript pooled values to compare:\n",
  "Myeloid: rho = 0.872\n",
  "Monocytes: rho = 0.820\n",
  "Macrophages: rho = 0.803\n",
  "Alternatively activated macrophages: rho = 0.782\n",
  "Inflammatory macrophages: rho = 0.696\n\n"
)

print(
  pooled_myeloid,
  digits = 10
)


# ------------------------------------------------------------
# 28. Save exact patient-level data
# ------------------------------------------------------------

saveRDS(
  list(
    HAVCR2 = havcr2_patient,
    xCell2 = xcell_patient,
    project = patient_project
  ),
  "results/xCell2_HAVCR2_patient_level_FINAL.rds"
)


# ------------------------------------------------------------
# 29. Session information
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file =
    "results/sessionInfo_xCell2_HAVCR2_FINAL.txt"
)


# ------------------------------------------------------------
# 30. Final summary
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("FINAL xCell2 VALIDATION COMPLETED\n")
cat("=============================================\n")

cat(
  "\nLGG patients:",
  length(lgg_patients),
  "\n"
)

cat(
  "GBM patients:",
  length(gbm_patients),
  "\n"
)

cat(
  "Total patients:",
  length(pooled_patients),
  "\n"
)

cat(
  "\nResults saved in results/ folder.\n"
)

cat(
  "\n=============================================\n"
)
