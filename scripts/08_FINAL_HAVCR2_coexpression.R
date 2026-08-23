# ============================================================
# 08_FINAL_HAVCR2_coexpression.R
#
# FINAL REVISION ANALYSIS
# HAVCR2-centered co-expression analysis
#
# TCGA-LGG and TCGA-GBM analyzed separately
#
# Strategy:
#   1. Load TCGA STAR-Counts SummarizedExperiment
#   2. Extract TPM
#   3. Collapse duplicated gene symbols
#   4. Transform expression as log2(TPM + 1)
#   5. Collapse multiple samples to patient level using median
#   6. Analyze TCGA-LGG and TCGA-GBM separately
#   7. Calculate gene-wise Spearman correlations with HAVCR2
#   8. Compare rho values between LGG and GBM
#
# Outputs:
#   results/HAVCR2_coexpression_LGG_FINAL.csv
#   results/HAVCR2_coexpression_GBM_FINAL.csv
#   results/HAVCR2_coexpression_LGG_GBM_comparison_FINAL.csv
#   results/HAVCR2_coexpression_TOP20_FINAL.csv
#   results/sessionInfo_HAVCR2_coexpression_FINAL.txt
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
# 2. Find TCGA STAR-Counts object
# ------------------------------------------------------------

possible_files <- c(
  "C:/GDC/TCGA_LGG_GBM_STAR_counts.rds",
  "data/TCGA_LGG_GBM_STAR_counts.rds",
  "TCGA_LGG_GBM_STAR_counts.rds"
)

existing_files <- possible_files[
  file.exists(possible_files)
]

if (length(existing_files) == 0) {
  
  stop(
    paste(
      "TCGA_LGG_GBM_STAR_counts.rds not found.",
      "Checked:",
      paste(
        possible_files,
        collapse = " | "
      )
    )
  )
}

infile <- existing_files[1]

cat("\n=============================================\n")
cat("FINAL HAVCR2 CO-EXPRESSION ANALYSIS\n")
cat("=============================================\n")

cat(
  "\nUsing input file:\n",
  infile,
  "\n"
)


# ------------------------------------------------------------
# 3. Load SummarizedExperiment
# ------------------------------------------------------------

se <- readRDS(
  infile
)

cat(
  "\nSummarizedExperiment dimensions:\n"
)

print(
  dim(se)
)

cat(
  "\nAvailable assays:\n"
)

print(
  assayNames(se)
)


# ------------------------------------------------------------
# 4. Check TPM assay
# ------------------------------------------------------------

if (!"tpm_unstrand" %in% assayNames(se)) {
  
  stop(
    "Assay 'tpm_unstrand' not found."
  )
}


# ------------------------------------------------------------
# 5. Extract TPM
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


# ------------------------------------------------------------
# 6. Clean genes
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
# 7. Collapse duplicated gene symbols
# ------------------------------------------------------------
#
# TPM values belonging to duplicated gene symbols
# are summed before log transformation.
# ------------------------------------------------------------

expr <- rowsum(
  expr,
  group = rownames(expr),
  reorder = FALSE
)

cat(
  "\nGenes after filtering/collapsing:",
  nrow(expr),
  "\n"
)

cat(
  "Expression samples:",
  ncol(expr),
  "\n"
)


# ------------------------------------------------------------
# 8. Log2 transformation
# ------------------------------------------------------------

expr_log <- log2(
  expr + 1
)


# ------------------------------------------------------------
# 9. HAVCR2 check
# ------------------------------------------------------------

if (!"HAVCR2" %in% rownames(expr_log)) {
  
  stop(
    "HAVCR2 not found in expression matrix."
  )
}

cat(
  "\nHAVCR2 successfully found.\n"
)

cat(
  "\nHAVCR2 sample-level summary:\n"
)

print(
  summary(
    as.numeric(
      expr_log[
        "HAVCR2",
      ]
    )
  )
)


# ------------------------------------------------------------
# 10. Obtain sample metadata
# ------------------------------------------------------------

meta <- as.data.frame(
  SummarizedExperiment::colData(se)
)

sample_ids <- colnames(
  expr_log
)

if (is.null(rownames(meta))) {
  
  stop(
    "colData(se) has no row names."
  )
}

matched_rows <- match(
  sample_ids,
  rownames(meta)
)

if (any(is.na(matched_rows))) {
  
  stop(
    paste(
      "Could not match",
      sum(is.na(matched_rows)),
      "expression samples to colData."
    )
  )
}

meta <- meta[
  matched_rows,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 11. Find project ID column
# ------------------------------------------------------------

project_candidates <- c(
  "project_id",
  "project.project_id",
  "project"
)

project_column <- project_candidates[
  project_candidates %in% colnames(meta)
]

if (length(project_column) == 0) {
  
  project_column <- grep(
    "project.*id|project_id",
    colnames(meta),
    ignore.case = TRUE,
    value = TRUE
  )
}

if (length(project_column) == 0) {
  
  cat(
    "\nAvailable metadata columns:\n"
  )
  
  print(
    colnames(meta)
  )
  
  stop(
    "No project ID column found in colData(se)."
  )
}

project_column <- project_column[1]

cat(
  "\nUsing project column:",
  project_column,
  "\n"
)

project_id <- as.character(
  meta[[project_column]]
)

cat(
  "\nProject IDs at sample level:\n"
)

print(
  table(
    project_id,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# 12. Construct TCGA patient IDs
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
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 13. Keep TCGA-LGG and TCGA-GBM only
# ------------------------------------------------------------

keep_samples <- (
  sample_data$project_id %in%
    c(
      "TCGA-LGG",
      "TCGA-GBM"
    )
)

sample_data <- sample_data[
  keep_samples,
  ,
  drop = FALSE
]

expr_log <- expr_log[
  ,
  sample_data$sample_id,
  drop = FALSE
]

cat(
  "\nSamples retained by project:\n"
)

print(
  table(
    sample_data$project_id
  )
)

cat(
  "\nUnique patients before collapse:\n"
)

print(
  table(
    sample_data$project_id[
      !duplicated(
        sample_data$patient_id
      )
    ]
  )
)


# ------------------------------------------------------------
# 14. Collapse expression to PATIENT level
# ------------------------------------------------------------
#
# If more than one sample is available for one patient,
# use the median log2(TPM + 1) expression.
# ------------------------------------------------------------

patients <- unique(
  sample_data$patient_id
)

expr_patient <- matrix(
  NA_real_,
  nrow = nrow(expr_log),
  ncol = length(patients),
  dimnames = list(
    rownames(expr_log),
    patients
  )
)

patient_project <- character(
  length(patients)
)

names(patient_project) <- patients


cat(
  "\nCollapsing expression to patient level...\n"
)

for (p in patients) {
  
  samples_p <- sample_data$sample_id[
    sample_data$patient_id == p
  ]
  
  projects_p <- unique(
    sample_data$project_id[
      sample_data$patient_id == p
    ]
  )
  
  if (length(projects_p) != 1) {
    
    stop(
      paste(
        "Patient",
        p,
        "has inconsistent project assignments."
      )
    )
  }
  
  patient_project[p] <- projects_p[1]
  
  if (length(samples_p) == 1) {
    
    expr_patient[
      ,
      p
    ] <- expr_log[
      ,
      samples_p
    ]
    
  } else {
    
    expr_patient[
      ,
      p
    ] <- apply(
      expr_log[
        ,
        samples_p,
        drop = FALSE
      ],
      1,
      median,
      na.rm = TRUE
    )
  }
}


# ------------------------------------------------------------
# 15. Patient-level cohort counts
# ------------------------------------------------------------

cat(
  "\n=============================================\n"
)

cat(
  "PATIENT-LEVEL COHORTS\n"
)

cat(
  "=============================================\n"
)

print(
  table(
    patient_project
  )
)

cat(
  "\nTotal unique patients:",
  ncol(expr_patient),
  "\n"
)


# ------------------------------------------------------------
# 16. Final HAVCR2 check at patient level
# ------------------------------------------------------------

havcr2_patient <- as.numeric(
  expr_patient[
    "HAVCR2",
  ]
)

names(havcr2_patient) <- colnames(
  expr_patient
)

if (all(is.na(havcr2_patient))) {
  
  stop(
    "Patient-level HAVCR2 expression is entirely missing."
  )
}


# ------------------------------------------------------------
# 17. Spearman co-expression helper function
# ------------------------------------------------------------

run_coexpression <- function(
    expression_matrix,
    cohort_name
) {
  
  cat(
    "\n=============================================\n"
  )
  
  cat(
    "RUNNING CO-EXPRESSION:",
    cohort_name,
    "\n"
  )
  
  cat(
    "=============================================\n"
  )
  
  n_patients <- ncol(
    expression_matrix
  )
  
  cat(
    "Patients:",
    n_patients,
    "\n"
  )
  
  hav <- as.numeric(
    expression_matrix[
      "HAVCR2",
    ]
  )
  
  genes <- setdiff(
    rownames(expression_matrix),
    "HAVCR2"
  )
  
  results <- vector(
    "list",
    length(genes)
  )
  
  for (i in seq_along(genes)) {
    
    g <- genes[i]
    
    x <- as.numeric(
      expression_matrix[
        g,
      ]
    )
    
    ok <- complete.cases(
      hav,
      x
    )
    
    n_ok <- sum(ok)
    
    if (
      n_ok < 20 ||
      length(unique(x[ok])) < 2 ||
      length(unique(hav[ok])) < 2
    ) {
      
      results[[i]] <- data.frame(
        gene = g,
        rho = NA_real_,
        p_value = NA_real_,
        n = n_ok,
        stringsAsFactors = FALSE
      )
      
      next
    }
    
    test <- suppressWarnings(
      cor.test(
        hav[ok],
        x[ok],
        method = "spearman",
        exact = FALSE
      )
    )
    
    results[[i]] <- data.frame(
      gene = g,
      rho = unname(
        test$estimate
      ),
      p_value = test$p.value,
      n = n_ok,
      stringsAsFactors = FALSE
    )
  }
  
  res <- do.call(
    rbind,
    results
  )
  
  res$FDR <- p.adjust(
    res$p_value,
    method = "BH"
  )
  
  res <- res[
    is.finite(res$rho),
    ,
    drop = FALSE
  ]
  
  # Rank positive correlations from strongest to weakest
  res <- res[
    order(
      -res$rho
    ),
    ,
    drop = FALSE
  ]
  
  res$positive_rank <- seq_len(
    nrow(res)
  )
  
  rownames(res) <- NULL
  
  res
}


# ------------------------------------------------------------
# 18. LGG patient-level expression
# ------------------------------------------------------------

lgg_patients <- names(
  patient_project
)[
  patient_project == "TCGA-LGG"
]

expr_lgg <- expr_patient[
  ,
  lgg_patients,
  drop = FALSE
]


# ------------------------------------------------------------
# 19. GBM patient-level expression
# ------------------------------------------------------------

gbm_patients <- names(
  patient_project
)[
  patient_project == "TCGA-GBM"
]

expr_gbm <- expr_patient[
  ,
  gbm_patients,
  drop = FALSE
]


# ------------------------------------------------------------
# 20. Run LGG co-expression
# ------------------------------------------------------------

res_lgg <- run_coexpression(
  expr_lgg,
  "TCGA-LGG"
)

res_lgg$cohort <- "TCGA-LGG"


# ------------------------------------------------------------
# 21. Run GBM co-expression
# ------------------------------------------------------------

res_gbm <- run_coexpression(
  expr_gbm,
  "TCGA-GBM"
)

res_gbm$cohort <- "TCGA-GBM"


# ------------------------------------------------------------
# 22. Save full LGG and GBM results
# ------------------------------------------------------------

dir.create(
  "results",
  showWarnings = FALSE,
  recursive = TRUE
)

write.csv(
  res_lgg,
  file =
    "results/HAVCR2_coexpression_LGG_FINAL.csv",
  row.names = FALSE
)

write.csv(
  res_gbm,
  file =
    "results/HAVCR2_coexpression_GBM_FINAL.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 23. Print strongest positive correlations
# ------------------------------------------------------------

cat(
  "\n=============================================\n"
)

cat(
  "TOP 20 POSITIVE CORRELATIONS - LGG\n"
)

cat(
  "=============================================\n\n"
)

print(
  head(
    res_lgg[
      ,
      c(
        "gene",
        "rho",
        "p_value",
        "FDR",
        "n"
      )
    ],
    20
  ),
  digits = 10
)


cat(
  "\n=============================================\n"
)

cat(
  "TOP 20 POSITIVE CORRELATIONS - GBM\n"
)

cat(
  "=============================================\n\n"
)

print(
  head(
    res_gbm[
      ,
      c(
        "gene",
        "rho",
        "p_value",
        "FDR",
        "n"
      )
    ],
    20
  ),
  digits = 10
)


# ------------------------------------------------------------
# 24. Create top-20 combined table
# ------------------------------------------------------------

top_lgg <- head(
  res_lgg[
    ,
    c(
      "gene",
      "rho",
      "p_value",
      "FDR",
      "n"
    )
  ],
  20
)

top_lgg$cohort <- "TCGA-LGG"

top_gbm <- head(
  res_gbm[
    ,
    c(
      "gene",
      "rho",
      "p_value",
      "FDR",
      "n"
    )
  ],
  20
)

top_gbm$cohort <- "TCGA-GBM"

top20 <- rbind(
  top_lgg,
  top_gbm
)

write.csv(
  top20,
  file =
    "results/HAVCR2_coexpression_TOP20_FINAL.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 25. Compare LGG and GBM correlations
# ------------------------------------------------------------

lgg_compare <- res_lgg[
  ,
  c(
    "gene",
    "rho",
    "p_value",
    "FDR",
    "n"
  )
]

gbm_compare <- res_gbm[
  ,
  c(
    "gene",
    "rho",
    "p_value",
    "FDR",
    "n"
  )
]

colnames(
  lgg_compare
)[-1] <- paste0(
  colnames(lgg_compare)[-1],
  "_LGG"
)

colnames(
  gbm_compare
)[-1] <- paste0(
  colnames(gbm_compare)[-1],
  "_GBM"
)

comparison <- merge(
  lgg_compare,
  gbm_compare,
  by = "gene",
  all = FALSE
)

comparison$delta_rho_LGG_minus_GBM <- (
  comparison$rho_LGG -
    comparison$rho_GBM
)

comparison$abs_delta_rho <- abs(
  comparison$delta_rho_LGG_minus_GBM
)

comparison <- comparison[
  order(
    -comparison$abs_delta_rho
  ),
  ,
  drop = FALSE
]

rownames(comparison) <- NULL


# ------------------------------------------------------------
# 26. Save LGG vs GBM comparison
# ------------------------------------------------------------

write.csv(
  comparison,
  file =
    "results/HAVCR2_coexpression_LGG_GBM_comparison_FINAL.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 27. Print genes with largest LGG/GBM differences
# ------------------------------------------------------------

cat(
  "\n=============================================\n"
)

cat(
  "LARGEST LGG vs GBM DIFFERENCES IN RHO\n"
)

cat(
  "=============================================\n\n"
)

print(
  head(
    comparison[
      ,
      c(
        "gene",
        "rho_LGG",
        "rho_GBM",
        "delta_rho_LGG_minus_GBM",
        "abs_delta_rho"
      )
    ],
    20
  ),
  digits = 10
)


# ------------------------------------------------------------
# 28. Specifically check genes currently reported
#     in the manuscript
# ------------------------------------------------------------

manuscript_genes <- c(
  "LAPTM5",
  "CD53",
  "NCKAP1L",
  "ITGB2",
  "LILRB4",
  "BTK",
  "SASH3",
  "SAMSN1",
  "CD86"
)

manuscript_check <- comparison[
  comparison$gene %in%
    manuscript_genes,
  ,
  drop = FALSE
]

manuscript_check <- manuscript_check[
  match(
    manuscript_genes,
    manuscript_check$gene
  ),
  ,
  drop = FALSE
]

manuscript_check <- manuscript_check[
  !is.na(manuscript_check$gene),
  ,
  drop = FALSE
]


cat(
  "\n=============================================\n"
)

cat(
  "GENES CURRENTLY REPORTED IN MANUSCRIPT\n"
)

cat(
  "=============================================\n\n"
)

print(
  manuscript_check[
    ,
    c(
      "gene",
      "rho_LGG",
      "rho_GBM",
      "n_LGG",
      "n_GBM"
    )
  ],
  digits = 10
)


# ------------------------------------------------------------
# 29. Save manuscript-gene check
# ------------------------------------------------------------

write.csv(
  manuscript_check,
  file =
    "results/HAVCR2_coexpression_manuscript_genes_CHECK.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 30. Save exact patient-level data
# ------------------------------------------------------------

saveRDS(
  list(
    expression = expr_patient,
    project = patient_project
  ),
  file =
    "results/HAVCR2_coexpression_patient_level_FINAL.rds"
)


# ------------------------------------------------------------
# 31. Save session information
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file =
    "results/sessionInfo_HAVCR2_coexpression_FINAL.txt"
)


# ------------------------------------------------------------
# 32. Final summary
# ------------------------------------------------------------

cat(
  "\n=============================================\n"
)

cat(
  "FINAL CO-EXPRESSION ANALYSIS COMPLETED\n"
)

cat(
  "=============================================\n"
)

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
  length(patients),
  "\n"
)

cat(
  "\nFiles saved in results/ folder.\n"
)

cat(
  "\n=============================================\n"
)