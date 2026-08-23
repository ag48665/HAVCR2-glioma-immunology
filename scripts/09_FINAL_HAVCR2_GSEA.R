# ============================================================
# 09_FINAL_HAVCR2_GSEA.R
#
# FINAL HAVCR2-ASSOCIATED GSEA
#
# Cohorts:
#   TCGA-LGG + TCGA-GBM pooled
#   TCGA-LGG
#   TCGA-GBM
#
# Expression:
#   TPM -> log2(TPM + 1)
#   patient-level
#
# Ranking:
#   Spearman correlation with HAVCR2
#
# Gene sets:
#   MSigDB Hallmark
#
# GSEA:
#   fgseaMultilevel
#   minSize = 15
#   maxSize = 500
# ============================================================


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

required <- c(
  "SummarizedExperiment",
  "fgsea",
  "msigdbr"
)

missing <- required[
  !vapply(
    required,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing) > 0) {
  stop(
    "Install required package(s): ",
    paste(missing, collapse = ", ")
  )
}

library(SummarizedExperiment)


# ------------------------------------------------------------
# 2. Input file
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
    "TCGA_LGG_GBM_STAR_counts.rds not found."
  )
}

infile <- existing_files[1]

cat("\n=============================================\n")
cat("FINAL HAVCR2 GSEA ANALYSIS\n")
cat("=============================================\n")

cat(
  "\nUsing input file:",
  infile,
  "\n"
)


# ------------------------------------------------------------
# 3. Load TCGA STAR-Counts
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
# 4. Extract TPM
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
# 6. Collapse duplicated gene symbols
# ------------------------------------------------------------

expr <- rowsum(
  expr,
  group = rownames(expr),
  reorder = FALSE
)

cat(
  "\nGenes after filtering/collapse:",
  nrow(expr),
  "\n"
)

cat(
  "Samples:",
  ncol(expr),
  "\n"
)


# ------------------------------------------------------------
# 7. log2(TPM + 1)
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
# 8. Obtain metadata
# ------------------------------------------------------------

meta <- as.data.frame(
  SummarizedExperiment::colData(se)
)

sample_ids <- colnames(
  expr_log
)

matched_rows <- match(
  sample_ids,
  rownames(meta)
)

if (any(is.na(matched_rows))) {
  stop(
    "Could not align expression samples with metadata."
  )
}

meta <- meta[
  matched_rows,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 9. Find project column
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
  stop(
    "Project ID column not found."
  )
}

project_column <- project_column[1]

project_id <- as.character(
  meta[[project_column]]
)

cat(
  "\nProject IDs:\n"
)

print(
  table(project_id)
)


# ------------------------------------------------------------
# 10. Patient IDs
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
# 11. Keep LGG and GBM
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


# ------------------------------------------------------------
# 12. Collapse to patient level
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
  "\nCollapsing to patient level...\n"
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
        "Inconsistent project assignment for",
        p
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
# 13. Cohort counts
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
  "\nTotal patients:",
  ncol(expr_patient),
  "\n"
)


# ------------------------------------------------------------
# 14. Remove genes with minimal variability
# ------------------------------------------------------------

gene_sd <- apply(
  expr_patient,
  1,
  sd,
  na.rm = TRUE
)

keep_variable <- (
  is.finite(gene_sd) &
    gene_sd > 0
)

expr_patient <- expr_patient[
  keep_variable,
  ,
  drop = FALSE
]

cat(
  "\nVariable genes retained:",
  nrow(expr_patient),
  "\n"
)

if (!"HAVCR2" %in% rownames(expr_patient)) {
  stop(
    "HAVCR2 removed during variability filtering."
  )
}


# ------------------------------------------------------------
# 15. Function to create HAVCR2 Spearman ranking
# ------------------------------------------------------------

make_ranking <- function(
    expression_matrix
) {
  
  hav <- as.numeric(
    expression_matrix[
      "HAVCR2",
    ]
  )
  
  genes <- setdiff(
    rownames(expression_matrix),
    "HAVCR2"
  )
  
  rho <- vapply(
    genes,
    function(g) {
      
      x <- as.numeric(
        expression_matrix[
          g,
        ]
      )
      
      ok <- complete.cases(
        hav,
        x
      )
      
      if (
        sum(ok) < 20 ||
        length(unique(x[ok])) < 2
      ) {
        return(
          NA_real_
        )
      }
      
      suppressWarnings(
        cor(
          hav[ok],
          x[ok],
          method = "spearman"
        )
      )
    },
    numeric(1)
  )
  
  rho <- rho[
    is.finite(rho)
  ]
  
  rho <- sort(
    rho,
    decreasing = TRUE
  )
  
  rho
}


# ------------------------------------------------------------
# 16. Prepare Hallmark gene sets
# ------------------------------------------------------------

hallmark_df <- msigdbr::msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

cat(
  "\nMSigDB Hallmark rows:",
  nrow(hallmark_df),
  "\n"
)

cat(
  "Hallmark gene sets:",
  length(
    unique(
      hallmark_df$gs_name
    )
  ),
  "\n"
)

hallmark_pathways <- split(
  hallmark_df$gene_symbol,
  hallmark_df$gs_name
)


# ------------------------------------------------------------
# 17. GSEA function
# ------------------------------------------------------------

run_gsea <- function(
    ranking,
    cohort_name
) {
  
  cat(
    "\n=============================================\n"
  )
  
  cat(
    "RUNNING GSEA:",
    cohort_name,
    "\n"
  )
  
  cat(
    "=============================================\n"
  )
  
  cat(
    "Ranking genes:",
    length(ranking),
    "\n"
  )
  
  gsea <- fgsea::fgseaMultilevel(
    pathways = hallmark_pathways,
    stats = ranking,
    minSize = 15,
    maxSize = 500,
    eps = 0
  )
  
  gsea <- as.data.frame(
    gsea
  )
  
  gsea$cohort <- cohort_name
  
  gsea <- gsea[
    order(
      gsea$padj,
      -abs(gsea$NES)
    ),
    ,
    drop = FALSE
  ]
  
  rownames(gsea) <- NULL
  
  gsea
}


# ------------------------------------------------------------
# 18. Safe CSV export
# ------------------------------------------------------------

prepare_gsea_for_csv <- function(x) {
  
  list_cols <- vapply(
    x,
    is.list,
    logical(1)
  )
  
  x[list_cols] <- lapply(
    x[list_cols],
    function(col) {
      
      vapply(
        col,
        function(z) {
          paste(
            z,
            collapse = ";"
          )
        },
        character(1)
      )
    }
  )
  
  x
}


# ------------------------------------------------------------
# 19. Pooled cohort ranking
# ------------------------------------------------------------

rank_all <- make_ranking(
  expr_patient
)

cat(
  "\nRanking length - ALL:",
  length(rank_all),
  "\n"
)


# ------------------------------------------------------------
# 20. LGG ranking
# ------------------------------------------------------------

lgg_patients <- names(
  patient_project
)[
  patient_project == "TCGA-LGG"
]

rank_lgg <- make_ranking(
  expr_patient[
    ,
    lgg_patients,
    drop = FALSE
  ]
)

cat(
  "Ranking length - LGG:",
  length(rank_lgg),
  "\n"
)


# ------------------------------------------------------------
# 21. GBM ranking
# ------------------------------------------------------------

gbm_patients <- names(
  patient_project
)[
  patient_project == "TCGA-GBM"
]

rank_gbm <- make_ranking(
  expr_patient[
    ,
    gbm_patients,
    drop = FALSE
  ]
)

cat(
  "Ranking length - GBM:",
  length(rank_gbm),
  "\n"
)


# ------------------------------------------------------------
# 22. Run GSEA
# ------------------------------------------------------------

gsea_all <- run_gsea(
  rank_all,
  "TCGA-LGG+GBM"
)

gsea_lgg <- run_gsea(
  rank_lgg,
  "TCGA-LGG"
)

gsea_gbm <- run_gsea(
  rank_gbm,
  "TCGA-GBM"
)


# ------------------------------------------------------------
# 23. Create results folder
# ------------------------------------------------------------

dir.create(
  "results",
  showWarnings = FALSE,
  recursive = TRUE
)


# ------------------------------------------------------------
# 24. Save full results
# ------------------------------------------------------------

write.csv(
  prepare_gsea_for_csv(
    gsea_all
  ),
  "results/GSEA_HAVCR2_Hallmark_ALL_FINAL.csv",
  row.names = FALSE
)

write.csv(
  prepare_gsea_for_csv(
    gsea_lgg
  ),
  "results/GSEA_HAVCR2_Hallmark_LGG_FINAL.csv",
  row.names = FALSE
)

write.csv(
  prepare_gsea_for_csv(
    gsea_gbm
  ),
  "results/GSEA_HAVCR2_Hallmark_GBM_FINAL.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 25. Print TOP GSEA results
# ------------------------------------------------------------

cat(
  "\n=============================================\n"
)

cat(
  "TOP GSEA RESULTS - POOLED\n"
)

cat(
  "=============================================\n\n"
)

print(
  head(
    gsea_all[
      ,
      c(
        "pathway",
        "NES",
        "pval",
        "padj",
        "size"
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
  "TOP GSEA RESULTS - LGG\n"
)

cat(
  "=============================================\n\n"
)

print(
  head(
    gsea_lgg[
      ,
      c(
        "pathway",
        "NES",
        "pval",
        "padj",
        "size"
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
  "TOP GSEA RESULTS - GBM\n"
)

cat(
  "=============================================\n\n"
)

print(
  head(
    gsea_gbm[
      ,
      c(
        "pathway",
        "NES",
        "pval",
        "padj",
        "size"
      )
    ],
    20
  ),
  digits = 10
)


# ------------------------------------------------------------
# 26. Immune pathways currently relevant to manuscript
# ------------------------------------------------------------

manuscript_pathways <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"
)


extract_manuscript_pathways <- function(
    x,
    cohort_label
) {
  
  out <- x[
    x$pathway %in%
      manuscript_pathways,
    ,
    drop = FALSE
  ]
  
  out$cohort <- cohort_label
  
  out[
    ,
    c(
      "cohort",
      "pathway",
      "NES",
      "pval",
      "padj",
      "size"
    )
  ]
}


manuscript_check <- rbind(
  extract_manuscript_pathways(
    gsea_all,
    "TCGA-LGG+GBM"
  ),
  extract_manuscript_pathways(
    gsea_lgg,
    "TCGA-LGG"
  ),
  extract_manuscript_pathways(
    gsea_gbm,
    "TCGA-GBM"
  )
)


# ------------------------------------------------------------
# 27. Print manuscript pathway check
# ------------------------------------------------------------

cat(
  "\n=============================================\n"
)

cat(
  "PATHWAYS CURRENTLY DESCRIBED IN MANUSCRIPT\n"
)

cat(
  "=============================================\n\n"
)

print(
  manuscript_check,
  digits = 10
)


# ------------------------------------------------------------
# 28. Save manuscript pathway check
# ------------------------------------------------------------

write.csv(
  manuscript_check,
  "results/GSEA_HAVCR2_manuscript_pathways_CHECK_FINAL.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 29. Immune-related pathway summary
# ------------------------------------------------------------

immune_pattern <- paste(
  c(
    "INTERFERON",
    "INFLAMMATORY",
    "IL6",
    "IL2",
    "TNFA",
    "COMPLEMENT",
    "ALLOGRAFT"
  ),
  collapse = "|"
)

immune_all <- gsea_all[
  grepl(
    immune_pattern,
    gsea_all$pathway,
    ignore.case = TRUE
  ),
  ,
  drop = FALSE
]

immune_lgg <- gsea_lgg[
  grepl(
    immune_pattern,
    gsea_lgg$pathway,
    ignore.case = TRUE
  ),
  ,
  drop = FALSE
]

immune_gbm <- gsea_gbm[
  grepl(
    immune_pattern,
    gsea_gbm$pathway,
    ignore.case = TRUE
  ),
  ,
  drop = FALSE
]

immune_summary <- rbind(
  immune_all,
  immune_lgg,
  immune_gbm
)

write.csv(
  prepare_gsea_for_csv(
    immune_summary
  ),
  "results/GSEA_HAVCR2_immune_pathways_summary_FINAL.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 30. Save gene rankings
# ------------------------------------------------------------

write.csv(
  data.frame(
    gene = names(rank_all),
    rho = as.numeric(rank_all)
  ),
  "results/HAVCR2_gene_ranking_ALL_FINAL.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    gene = names(rank_lgg),
    rho = as.numeric(rank_lgg)
  ),
  "results/HAVCR2_gene_ranking_LGG_FINAL.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    gene = names(rank_gbm),
    rho = as.numeric(rank_gbm)
  ),
  "results/HAVCR2_gene_ranking_GBM_FINAL.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 31. Save R objects
# ------------------------------------------------------------

saveRDS(
  list(
    ALL = gsea_all,
    LGG = gsea_lgg,
    GBM = gsea_gbm
  ),
  "results/GSEA_HAVCR2_Hallmark_results_FINAL.rds"
)


# ------------------------------------------------------------
# 32. Save session info
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file =
    "results/sessionInfo_GSEA_HAVCR2_FINAL.txt"
)


# ------------------------------------------------------------
# 33. Final summary
# ------------------------------------------------------------

cat(
  "\n=============================================\n"
)

cat(
  "FINAL HAVCR2 GSEA ANALYSIS COMPLETED\n"
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
  "\nResults saved in results/ folder.\n"
)

cat(
  "\n=============================================\n"
)