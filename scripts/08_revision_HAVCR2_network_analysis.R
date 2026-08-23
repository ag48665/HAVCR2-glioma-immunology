# ============================================================
# Revision analysis: HAVCR2 co-expression network
# Reviewer 2, Comment 6
# TCGA-LGG and TCGA-GBM
# ============================================================

library(SummarizedExperiment)

# ------------------------------------------------------------
# 1. Load TCGA STAR-Counts object
# ------------------------------------------------------------

se <- readRDS(
  "C:/GDC/TCGA_LGG_GBM_STAR_counts.rds"
)

cat("\n=== SUMMARIZED EXPERIMENT ===\n")
print(dim(se))

cat("\nAvailable assays:\n")
print(assayNames(se))

# ------------------------------------------------------------
# 2. Extract TPM
# ------------------------------------------------------------

expr <- assay(
  se,
  "tpm_unstrand"
)

rd <- as.data.frame(
  rowData(se)
)

if (!"gene_name" %in% colnames(rd)) {
  stop("gene_name not found in rowData(se)")
}

gene_symbols <- as.character(
  rd$gene_name
)

# ------------------------------------------------------------
# 3. Clean genes
# ------------------------------------------------------------

keep <- !is.na(gene_symbols) &
  gene_symbols != "" &
  rowSums(expr, na.rm = TRUE) > 0

expr <- expr[
  keep,
  ,
  drop = FALSE
]

gene_symbols <- gene_symbols[
  keep
]

# Collapse duplicate symbols
expr <- rowsum(
  expr,
  group = gene_symbols,
  reorder = FALSE
)

expr_log <- log2(
  expr + 1
)

cat(
  "\nGenes:",
  nrow(expr_log),
  "\n"
)

cat(
  "Samples:",
  ncol(expr_log),
  "\n"
)

# ------------------------------------------------------------
# 4. HAVCR2 check
# ------------------------------------------------------------

if (!"HAVCR2" %in% rownames(expr_log)) {
  stop("HAVCR2 not found.")
}

cat("\nHAVCR2 found successfully.\n")

print(
  summary(
    as.numeric(
      expr_log["HAVCR2", ]
    )
  )
)

# ------------------------------------------------------------
# 5. Inspect sample metadata
# ------------------------------------------------------------

meta <- as.data.frame(
  colData(se)
)

cat("\n=== METADATA COLUMNS ===\n")
print(
  colnames(meta)
)

cat("\n=== PROJECT INFORMATION ===\n")

project_candidates <- grep(
  "project",
  colnames(meta),
  ignore.case = TRUE,
  value = TRUE
)

print(project_candidates)

for (x in project_candidates) {
  
  cat(
    "\nCOLUMN:",
    x,
    "\n"
  )
  
  print(
    table(
      meta[[x]],
      useNA = "ifany"
    )
  )
}

# ------------------------------------------------------------
# 6. Sample identifiers
# ------------------------------------------------------------

cat("\n=== FIRST SAMPLE IDS ===\n")

print(
  head(
    colnames(expr_log),
    10
  )
)

cat("\n=== FIRST METADATA ROW NAMES ===\n")

print(
  head(
    rownames(meta),
    10
  )
)

cat("\n===========================================\n")
cat("NETWORK ANALYSIS — PREPARATION COMPLETED\n")
cat("===========================================\n")

