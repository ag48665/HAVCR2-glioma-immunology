# ============================================================
# Revision analysis:
# HAVCR2-associated gene set enrichment analysis (GSEA)
#
# Addresses Reviewer 2 Comment 6
#
# Strategy:
# 1. Use TCGA-LGG + TCGA-GBM TPM expression
# 2. Collapse to patient level
# 3. Rank genes by Spearman correlation with HAVCR2
# 4. Perform preranked GSEA using MSigDB Hallmark gene sets
# 5. Repeat separately in LGG and GBM
#
# Input:
#   C:/GDC/TCGA_LGG_GBM_STAR_counts.rds
#
# Outputs:
#   results/GSEA_HAVCR2_*.csv
# ============================================================


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("fgsea", quietly = TRUE)) {
  BiocManager::install(
    "fgsea",
    ask = FALSE,
    update = FALSE
  )
}

if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  BiocManager::install(
    "SummarizedExperiment",
    ask = FALSE,
    update = FALSE
  )
}

if (!requireNamespace("msigdbr", quietly = TRUE)) {
  install.packages("msigdbr")
}

library(SummarizedExperiment)
library(fgsea)
library(msigdbr)


# ------------------------------------------------------------
# 2. Load TCGA expression data
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

# Collapse duplicated gene symbols
expr <- rowsum(
  expr,
  group = rownames(expr),
  reorder = FALSE
)

# log2 TPM
expr_log <- log2(
  expr + 1
)

cat(
  "\nGenes after filtering:",
  nrow(expr_log),
  "\n"
)

cat(
  "Samples:",
  ncol(expr_log),
  "\n"
)


# ------------------------------------------------------------
# 4. Obtain sample/project metadata
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

if (!"project_id" %in% colnames(cd)) {
  stop(
    "project_id not found in colData(se)."
  )
}

project_id <- as.character(
  cd$project_id
)

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
  "\nSamples by project:\n"
)

print(
  table(
    metadata$project_id
  )
)


# ------------------------------------------------------------
# 5. Collapse expression to patient level
# ------------------------------------------------------------

patients <- unique(
  metadata$patient_id
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

for (p in patients) {
  
  samp <- metadata$sample_id[
    metadata$patient_id == p
  ]
  
  if (length(samp) == 1) {
    
    expr_patient[
      ,
      p
    ] <- expr_log[
      ,
      samp
    ]
    
  } else {
    
    expr_patient[
      ,
      p
    ] <- apply(
      expr_log[
        ,
        samp,
        drop = FALSE
      ],
      1,
      median,
      na.rm = TRUE
    )
  }
  
  patient_project[p] <- unique(
    metadata$project_id[
      metadata$patient_id == p
    ]
  )[1]
}

cat(
  "\nPatient-level cohort:\n"
)

print(
  table(
    patient_project
  )
)


# ------------------------------------------------------------
# 6. Remove genes with almost no variability
# ------------------------------------------------------------

gene_sd <- apply(
  expr_patient,
  1,
  sd,
  na.rm = TRUE
)

keep_var <- (
  is.finite(gene_sd) &
    gene_sd > 0.05
)

expr_patient <- expr_patient[
  keep_var,
  ,
  drop = FALSE
]

cat(
  "\nVariable genes retained:",
  nrow(expr_patient),
  "\n"
)


# ------------------------------------------------------------
# 7. Check HAVCR2
# ------------------------------------------------------------

if (!"HAVCR2" %in% rownames(expr_patient)) {
  stop(
    "HAVCR2 not found after filtering."
  )
}


# ------------------------------------------------------------
# 8. Function to create HAVCR2-correlated gene ranking
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
      
      if (sum(ok) < 20) {
        return(NA_real_)
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
# 9. Prepare MSigDB Hallmark gene sets
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
  length(unique(hallmark_df$gs_name)),
  "\n"
)

hallmark_pathways <- split(
  hallmark_df$gene_symbol,
  hallmark_df$gs_name
)


# ------------------------------------------------------------
# 10. Function to run fgsea
# ------------------------------------------------------------

run_gsea <- function(
    ranking,
    cohort_name
) {
  
  cat(
    "\nRunning GSEA for",
    cohort_name,
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
  
  gsea
}


# ------------------------------------------------------------
# Helper function for safe CSV export
# Converts list columns such as leadingEdge to text
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
        function(z) paste(z, collapse = ";"),
        character(1)
      )
    }
  )
  
  x
}


# ------------------------------------------------------------
# 11. Pooled LGG + GBM GSEA
# ------------------------------------------------------------

rank_all <- make_ranking(
  expr_patient
)

cat(
  "\nRanking length - ALL:",
  length(rank_all),
  "\n"
)

gsea_all <- run_gsea(
  rank_all,
  "TCGA-LGG+GBM"
)

write.csv(
  prepare_gsea_for_csv(gsea_all),
  file = "results/GSEA_HAVCR2_Hallmark_ALL.csv",
  row.names = FALSE
)

cat(
  "\nTop pooled GSEA results:\n"
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
    15
  )
)


# ------------------------------------------------------------
# 12. LGG-specific GSEA
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
  "\nRanking length - LGG:",
  length(rank_lgg),
  "\n"
)

gsea_lgg <- run_gsea(
  rank_lgg,
  "TCGA-LGG"
)

write.csv(
  prepare_gsea_for_csv(gsea_lgg),
  file = "results/GSEA_HAVCR2_Hallmark_LGG.csv",
  row.names = FALSE
)

cat(
  "\nTop LGG GSEA results:\n"
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
    15
  )
)


# ------------------------------------------------------------
# 13. GBM-specific GSEA
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
  "\nRanking length - GBM:",
  length(rank_gbm),
  "\n"
)

gsea_gbm <- run_gsea(
  rank_gbm,
  "TCGA-GBM"
)

write.csv(
  prepare_gsea_for_csv(gsea_gbm),
  file = "results/GSEA_HAVCR2_Hallmark_GBM.csv",
  row.names = FALSE
)

cat(
  "\nTop GBM GSEA results:\n"
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
    15
  )
)


# ------------------------------------------------------------
# 14. Save gene rankings
# ------------------------------------------------------------

write.csv(
  data.frame(
    gene = names(rank_all),
    rho = as.numeric(rank_all)
  ),
  file = "results/HAVCR2_gene_ranking_ALL.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    gene = names(rank_lgg),
    rho = as.numeric(rank_lgg)
  ),
  file = "results/HAVCR2_gene_ranking_LGG.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    gene = names(rank_gbm),
    rho = as.numeric(rank_gbm)
  ),
  file = "results/HAVCR2_gene_ranking_GBM.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 15. Extract immune-related pathways
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
  prepare_gsea_for_csv(immune_summary),
  file = "results/GSEA_HAVCR2_immune_pathways_summary.csv",
  row.names = FALSE
)

cat(
  "\nImmune-related Hallmark pathways:\n"
)

print(
  immune_summary[
    ,
    c(
      "cohort",
      "pathway",
      "NES",
      "pval",
      "padj"
    )
  ]
)


# ------------------------------------------------------------
# 16. Save R objects
# ------------------------------------------------------------

saveRDS(
  list(
    ALL = gsea_all,
    LGG = gsea_lgg,
    GBM = gsea_gbm
  ),
  file = "results/GSEA_HAVCR2_Hallmark_results.rds"
)


# ------------------------------------------------------------
# 17. Session information
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = "results/sessionInfo_GSEA_HAVCR2.txt"
)


cat(
  "\n============================================\n"
)

cat(
  "HAVCR2 GSEA ANALYSIS COMPLETED\n"
)

cat(
  "============================================\n"
)

cat(
  "Results written to results/ folder.\n"
)