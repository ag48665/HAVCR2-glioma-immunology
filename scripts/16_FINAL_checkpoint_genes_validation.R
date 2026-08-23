# ============================================================
# 16_FINAL_checkpoint_genes_validation.R
#
# FINAL VALIDATION
# Immune checkpoint expression vs ESTIMATE immune infiltration
#
# Genes:
#   PDCD1, CD274, CTLA4, LAG3, TIGIT, HAVCR2
#
# Outputs:
#   Table3_FINAL_checkpoint_expression.csv
#   Figure3_FINAL_checkpoint_correlations.csv
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

clinical_file <- file.path(
  "results",
  "tcga_clinical_clean.csv"
)

if (!file.exists(se_file)) {
  stop(
    "STAR-counts file not found: ",
    se_file
  )
}

if (!file.exists(clinical_file)) {
  stop(
    "Clinical file not found: ",
    clinical_file
  )
}


cat("\n=============================================\n")
cat("FINAL CHECKPOINT GENE VALIDATION\n")
cat("=============================================\n")


# ------------------------------------------------------------
# 3. Load TCGA expression object
# ------------------------------------------------------------

se <- readRDS(
  se_file
)

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
# 4. Clean expression matrix
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
# 5. Collapse duplicate gene symbols
# ------------------------------------------------------------

expr <- rowsum(
  expr,
  group = rownames(expr),
  reorder = FALSE
)


# ------------------------------------------------------------
# 6. Check checkpoint genes
# ------------------------------------------------------------

checkpoint_genes <- c(
  "PDCD1",
  "CD274",
  "CTLA4",
  "LAG3",
  "TIGIT",
  "HAVCR2"
)

missing_genes <- setdiff(
  checkpoint_genes,
  rownames(expr)
)

if (length(missing_genes) > 0) {
  stop(
    "Missing checkpoint genes: ",
    paste(
      missing_genes,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 7. Sample metadata
# ------------------------------------------------------------

meta <- as.data.frame(
  SummarizedExperiment::colData(se)
)

sample_ids <- colnames(
  expr
)

matched_rows <- match(
  sample_ids,
  rownames(meta)
)

if (any(is.na(matched_rows))) {
  stop(
    "Could not align expression samples to colData."
  )
}

meta <- meta[
  matched_rows,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 8. Find project column
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
    "Project column not found."
  )
}

project_column <- project_column[1]

project_id <- as.character(
  meta[[project_column]]
)


# ------------------------------------------------------------
# 9. Build sample-level dataset
# ------------------------------------------------------------

patient_id <- substr(
  sample_ids,
  1,
  12
)

sample_data <- data.frame(
  sample_id = sample_ids,
  patient = patient_id,
  project_id = project_id,
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

expr <- expr[
  ,
  sample_data$sample_id,
  drop = FALSE
]


# ------------------------------------------------------------
# 10. Collapse checkpoint expression to patient level
# ------------------------------------------------------------
#
# IMPORTANT:
# Use TPM scale here because Table 3 values in manuscript
# are on the original TPM scale, not log2(TPM + 1).
# ------------------------------------------------------------

patients <- unique(
  sample_data$patient
)

checkpoint_patient <- matrix(
  NA_real_,
  nrow = length(checkpoint_genes),
  ncol = length(patients),
  dimnames = list(
    checkpoint_genes,
    patients
  )
)

for (p in patients) {
  
  samples_p <- sample_data$sample_id[
    sample_data$patient == p
  ]
  
  if (length(samples_p) == 1) {
    
    checkpoint_patient[
      ,
      p
    ] <- expr[
      checkpoint_genes,
      samples_p
    ]
    
  } else {
    
    checkpoint_patient[
      ,
      p
    ] <- apply(
      expr[
        checkpoint_genes,
        samples_p,
        drop = FALSE
      ],
      1,
      median,
      na.rm = TRUE
    )
  }
}


cat(
  "\nPatient-level expression patients:",
  ncol(checkpoint_patient),
  "\n"
)


# ------------------------------------------------------------
# 11. Load ESTIMATE immune scores
# ------------------------------------------------------------

clin <- read.csv(
  clinical_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_clin <- c(
  "patient",
  "paper_ESTIMATE.immune.score"
)

missing_clin <- setdiff(
  required_clin,
  colnames(clin)
)

if (length(missing_clin) > 0) {
  stop(
    "Missing clinical column(s): ",
    paste(
      missing_clin,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 12. Collapse clinical data to patient level
# ------------------------------------------------------------

clin$.has_immune <- !is.na(
  clin[["paper_ESTIMATE.immune.score"]]
)

clin$.row_order <- seq_len(
  nrow(clin)
)

clin <- clin[
  order(
    clin$patient,
    -as.integer(clin$.has_immune),
    clin$.row_order
  ),
]

clin <- clin[
  !duplicated(clin$patient),
  ,
  drop = FALSE
]

clin$immune_score <- as.numeric(
  clin[["paper_ESTIMATE.immune.score"]]
)

clin2 <- clin[
  ,
  c(
    "patient",
    "immune_score"
  ),
  drop = FALSE
]


# ------------------------------------------------------------
# 13. Build merged patient-level dataset
# ------------------------------------------------------------

expr_df <- data.frame(
  patient = colnames(checkpoint_patient),
  t(checkpoint_patient),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

merged <- merge(
  expr_df,
  clin2,
  by = "patient",
  all = FALSE
)

merged <- merged[
  !is.na(merged$immune_score),
  ,
  drop = FALSE
]


cat("\n=============================================\n")
cat("CHECKPOINT ANALYSIS COHORT\n")
cat("=============================================\n")

cat(
  "Patients with checkpoint expression + immune score:",
  nrow(merged),
  "\n"
)


# ------------------------------------------------------------
# 14. Median immune-score cutoff
# ------------------------------------------------------------

immune_median <- median(
  merged$immune_score,
  na.rm = TRUE
)

merged$immune_group <- ifelse(
  merged$immune_score > immune_median,
  "High",
  "Low"
)

merged$immune_group <- factor(
  merged$immune_group,
  levels = c(
    "Low",
    "High"
  )
)

cat(
  "Median immune score:",
  immune_median,
  "\n"
)

cat(
  "\nImmune groups:\n"
)

print(
  table(
    merged$immune_group
  )
)


# ------------------------------------------------------------
# 15. Wilcoxon comparisons
# ------------------------------------------------------------

wilcox_results <- lapply(
  checkpoint_genes,
  function(g) {
    
    high <- merged[
      merged$immune_group == "High",
      g
    ]
    
    low <- merged[
      merged$immune_group == "Low",
      g
    ]
    
    high <- as.numeric(high)
    low <- as.numeric(low)
    
    test <- suppressWarnings(
      wilcox.test(
        high,
        low,
        exact = FALSE
      )
    )
    
    data.frame(
      Gene = g,
      Median_High = median(
        high,
        na.rm = TRUE
      ),
      Median_Low = median(
        low,
        na.rm = TRUE
      ),
      P_value = test$p.value,
      stringsAsFactors = FALSE
    )
  }
)

table3 <- do.call(
  rbind,
  wilcox_results
)

table3$FDR <- p.adjust(
  table3$P_value,
  method = "BH"
)


# ------------------------------------------------------------
# 16. Order Table 3 like manuscript
# ------------------------------------------------------------

table3$Gene <- factor(
  table3$Gene,
  levels = c(
    "HAVCR2",
    "PDCD1",
    "CD274",
    "CTLA4",
    "LAG3",
    "TIGIT"
  )
)

table3 <- table3[
  order(table3$Gene),
  ,
  drop = FALSE
]

table3$Gene <- as.character(
  table3$Gene
)


# ------------------------------------------------------------
# 17. Print Table 3 results
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("TABLE 3 - FINAL VALIDATION\n")
cat("=============================================\n\n")

print(
  table3,
  digits = 12
)


# ------------------------------------------------------------
# 18. Spearman correlations with immune score
# ------------------------------------------------------------

cor_results <- lapply(
  checkpoint_genes,
  function(g) {
    
    x <- as.numeric(
      merged[[g]]
    )
    
    y <- merged$immune_score
    
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
      Gene = g,
      rho = unname(
        test$estimate
      ),
      p_value = test$p.value,
      n = sum(ok),
      stringsAsFactors = FALSE
    )
  }
)

figure3 <- do.call(
  rbind,
  cor_results
)

figure3$FDR <- p.adjust(
  figure3$p_value,
  method = "BH"
)

figure3 <- figure3[
  order(
    -figure3$rho
  ),
  ,
  drop = FALSE
]

rownames(figure3) <- NULL


# ------------------------------------------------------------
# 19. Print Figure 3 correlation results
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("FIGURE 3 - SPEARMAN CORRELATIONS\n")
cat("=============================================\n\n")

print(
  figure3,
  digits = 12
)


# ------------------------------------------------------------
# 20. Explicit manuscript checks
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("MANUSCRIPT CHECK\n")
cat("=============================================\n")

hav_row <- table3[
  table3$Gene == "HAVCR2",
  ,
  drop = FALSE
]

cat(
  "\nHAVCR2 Median High =",
  hav_row$Median_High,
  "\n"
)

cat(
  "HAVCR2 Median Low =",
  hav_row$Median_Low,
  "\n"
)

cat(
  "HAVCR2 FDR =",
  format(
    hav_row$FDR,
    scientific = TRUE,
    digits = 12
  ),
  "\n"
)

cat(
  "\nStrongest Spearman gene:",
  figure3$Gene[1],
  "\n"
)

cat(
  "Strongest rho:",
  figure3$rho[1],
  "\n"
)


# ------------------------------------------------------------
# 21. Save outputs
# ------------------------------------------------------------

dir.create(
  "results",
  showWarnings = FALSE,
  recursive = TRUE
)

write.csv(
  table3,
  "results/Table3_FINAL_checkpoint_expression.csv",
  row.names = FALSE
)

write.csv(
  figure3,
  "results/Figure3_FINAL_checkpoint_correlations.csv",
  row.names = FALSE
)

write.csv(
  merged,
  "results/checkpoint_analysis_patient_level_FINAL.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 22. Session info
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file =
    "results/sessionInfo_FINAL_checkpoint_validation.txt"
)


# ------------------------------------------------------------
# 23. Final message
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("CHECKPOINT VALIDATION COMPLETED\n")
cat("=============================================\n")

cat(
  "\nFiles saved:\n",
  "results/Table3_FINAL_checkpoint_expression.csv\n",
  "results/Figure3_FINAL_checkpoint_correlations.csv\n",
  "results/checkpoint_analysis_patient_level_FINAL.csv\n"
)

cat("\n=============================================\n")
