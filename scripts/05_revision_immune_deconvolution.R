# ============================================================
# Revision analysis: immune cell-type enrichment with xCell2
# Addresses Reviewer 1 Comments 1-2 and Reviewer 2 Comments 1-3
# Input: TCGA-LGG + TCGA-GBM STAR-Counts prepared by TCGAbiolinks
# ============================================================

# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_bioc <- c(
  "SummarizedExperiment",
  "xCell2",
  "BiocParallel"
)

for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(SummarizedExperiment)
library(xCell2)
library(BiocParallel)

# ------------------------------------------------------------
# 2. Load prepared TCGA expression object
# ------------------------------------------------------------

se <- readRDS(
  "C:/GDC/TCGA_LGG_GBM_STAR_counts.rds"
)

cat("Dimensions of SummarizedExperiment:\n")
print(dim(se))

cat("\nAvailable assays:\n")
print(SummarizedExperiment::assayNames(se))

# ------------------------------------------------------------
# 3. Extract TPM expression
# ------------------------------------------------------------

expr <- SummarizedExperiment::assay(
  se,
  "tpm_unstrand"
)

cat("\nExpression matrix dimensions:\n")
print(dim(expr))

# ------------------------------------------------------------
# 4. Obtain gene symbols
# ------------------------------------------------------------

rd <- as.data.frame(
  SummarizedExperiment::rowData(se)
)

cat("\nAvailable rowData columns:\n")
print(colnames(rd))

if (!"gene_name" %in% colnames(rd)) {
  stop(
    "Column 'gene_name' not found in rowData(se)."
  )
}

gene_symbols <- as.character(
  rd$gene_name
)

# ------------------------------------------------------------
# 5. Clean expression matrix
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

# Collapse duplicate gene symbols by SUM.
# TPM values assigned to duplicated gene symbols are combined.
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
  "Samples:",
  ncol(expr),
  "\n"
)

# ------------------------------------------------------------
# 6. Expression transformation
# ------------------------------------------------------------

expr_log <- log2(
  expr + 1
)

# ------------------------------------------------------------
# 7. Basic HAVCR2 check
# ------------------------------------------------------------

if (!"HAVCR2" %in% rownames(expr_log)) {
  stop(
    "HAVCR2 not found in expression matrix."
  )
}

havcr2 <- as.numeric(
  expr_log[
    "HAVCR2",
  ]
)

cat(
  "\nHAVCR2 expression summary:\n"
)

print(
  summary(havcr2)
)

# ------------------------------------------------------------
# 8. Load pre-trained xCell2 reference
# Tumor Microenvironment Compendium
# ------------------------------------------------------------

data(
  "TMECompendium.xCell2Ref",
  package = "xCell2"
)

cat(
  "\nLoaded xCell2 reference:\n"
)

print(
  class(TMECompendium.xCell2Ref)
)

# ------------------------------------------------------------
# 9. Check gene overlap with reference
# ------------------------------------------------------------

ref_genes <- xCell2::getGenesUsed(
  TMECompendium.xCell2Ref
)

shared_genes <- intersect(
  rownames(expr_log),
  ref_genes
)

cat(
  "\nGenes in TCGA expression matrix:",
  nrow(expr_log),
  "\n"
)

cat(
  "Genes in xCell2 reference:",
  length(ref_genes),
  "\n"
)

cat(
  "Shared genes:",
  length(shared_genes),
  "\n"
)

cat(
  "Fraction of reference genes shared:",
  round(
    length(shared_genes) /
      length(ref_genes),
    4
  ),
  "\n"
)

# ------------------------------------------------------------
# 10. Run xCell2
# ------------------------------------------------------------

cat(
  "\nStarting xCell2 analysis...\n"
)

# SerialParam is safest on Windows.
xcell_scores <- xCell2::xCell2Analysis(
  mix = expr_log,
  xcell2object = TMECompendium.xCell2Ref,
  minSharedGenes = 0.80,
  rawScores = FALSE,
  spillover = TRUE,
  spilloverAlpha = 0.5,
  BPPARAM = BiocParallel::SerialParam()
)

cat(
  "\nxCell2 completed.\n"
)

cat(
  "Cell types:",
  nrow(xcell_scores),
  "\n"
)

cat(
  "Samples:",
  ncol(xcell_scores),
  "\n"
)

# ------------------------------------------------------------
# 11. Save full xCell2 output
# ------------------------------------------------------------

write.csv(
  xcell_scores,
  file =
    "results/xCell2_TCGA_LGG_GBM_all_celltypes.csv",
  row.names = TRUE
)

# ------------------------------------------------------------
# 12. Correlate every xCell2 cell-type score with HAVCR2
# ------------------------------------------------------------

cor_results <- lapply(
  rownames(xcell_scores),
  function(celltype) {
    
    score <- as.numeric(
      xcell_scores[
        celltype,
      ]
    )
    
    ok <- complete.cases(
      havcr2,
      score
    )
    
    if (sum(ok) < 10) {
      return(NULL)
    }
    
    test <- suppressWarnings(
      cor.test(
        havcr2[ok],
        score[ok],
        method = "spearman",
        exact = FALSE
      )
    )
    
    data.frame(
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

cor_results <- do.call(
  rbind,
  cor_results
)

cor_results$FDR <- p.adjust(
  cor_results$p_value,
  method = "BH"
)

cor_results <- cor_results[
  order(
    -abs(cor_results$rho)
  ),
]

write.csv(
  cor_results,
  file =
    "results/xCell2_HAVCR2_correlations_all_celltypes.csv",
  row.names = FALSE
)

cat(
  "\nTop HAVCR2-cell type correlations:\n"
)

print(
  head(
    cor_results,
    20
  )
)

# ------------------------------------------------------------
# 13. Myeloid / macrophage / microglia-related results
# ------------------------------------------------------------

myeloid_terms <- grep(
  paste(
    c(
      "macroph",
      "microgl",
      "monocyte",
      "dendritic",
      "myeloid"
    ),
    collapse = "|"
  ),
  cor_results$cell_type,
  ignore.case = TRUE,
  value = TRUE
)

cat(
  "\nRelevant myeloid / macrophage / microglial cell types:\n"
)

print(
  myeloid_terms
)

myeloid_results <- cor_results[
  cor_results$cell_type %in%
    myeloid_terms,
  ,
  drop = FALSE
]

write.csv(
  myeloid_results,
  file =
    "results/xCell2_HAVCR2_myeloid_microglia_correlations.csv",
  row.names = FALSE
)

cat(
  "\nHAVCR2 correlations with myeloid-related cell types:\n"
)

print(
  myeloid_results
)

# ------------------------------------------------------------
# 14. Independent microglia marker score
# Reviewer explicitly mentioned TMEM119 and P2RY12
# ------------------------------------------------------------

microglia_markers <- c(
  "TMEM119",
  "P2RY12",
  "CX3CR1",
  "SALL1"
)

available_microglia <- intersect(
  microglia_markers,
  rownames(expr_log)
)

cat(
  "\nAvailable microglia markers:\n"
)

print(
  available_microglia
)

if (
  length(available_microglia) >= 2
) {
  
  microglia_score <- colMeans(
    expr_log[
      available_microglia,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )
  
  test_microglia <- cor.test(
    havcr2,
    microglia_score,
    method = "spearman",
    exact = FALSE
  )
  
  microglia_result <- data.frame(
    signature =
      "microglia_marker_score",
    markers =
      paste(
        available_microglia,
        collapse = ";"
      ),
    rho =
      unname(
        test_microglia$estimate
      ),
    p_value =
      test_microglia$p.value,
    n =
      sum(
        complete.cases(
          havcr2,
          microglia_score
        )
      )
  )
  
  write.csv(
    microglia_result,
    file =
      "results/HAVCR2_microglia_marker_score_correlation.csv",
    row.names = FALSE
  )
  
  cat(
    "\nHAVCR2 vs microglia marker score:\n"
  )
  
  print(
    microglia_result
  )
  
} else {
  
  warning(
    paste(
      "Fewer than 2 microglia markers",
      "were present;",
      "microglia marker score",
      "was not calculated."
    )
  )
}

# ------------------------------------------------------------
# 15. Original macrophage marker score
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

available_macrophage <- intersect(
  macrophage_markers,
  rownames(expr_log)
)

cat(
  "\nAvailable macrophage markers:\n"
)

print(
  available_macrophage
)

if (
  length(available_macrophage) >= 2
) {
  
  macrophage_score <- colMeans(
    expr_log[
      available_macrophage,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )
  
  test_macrophage <- cor.test(
    havcr2,
    macrophage_score,
    method = "spearman",
    exact = FALSE
  )
  
  macrophage_result <- data.frame(
    signature =
      "original_macrophage_marker_score",
    markers =
      paste(
        available_macrophage,
        collapse = ";"
      ),
    rho =
      unname(
        test_macrophage$estimate
      ),
    p_value =
      test_macrophage$p.value,
    n =
      sum(
        complete.cases(
          havcr2,
          macrophage_score
        )
      )
  )
  
  write.csv(
    macrophage_result,
    file =
      "results/HAVCR2_original_macrophage_score_correlation.csv",
    row.names = FALSE
  )
  
  cat(
    "\nHAVCR2 vs original macrophage score:\n"
  )
  
  print(
    macrophage_result
  )
}

# ------------------------------------------------------------
# 16. Compare macrophage and microglia marker scores directly
# ------------------------------------------------------------

if (
  exists("microglia_score") &&
  exists("macrophage_score")
) {
  
  test_macro_micro <- cor.test(
    macrophage_score,
    microglia_score,
    method = "spearman",
    exact = FALSE
  )
  
  macro_micro_result <- data.frame(
    comparison =
      "macrophage_score_vs_microglia_score",
    rho =
      unname(
        test_macro_micro$estimate
      ),
    p_value =
      test_macro_micro$p.value,
    n =
      sum(
        complete.cases(
          macrophage_score,
          microglia_score
        )
      )
  )
  
  write.csv(
    macro_micro_result,
    file =
      "results/macrophage_vs_microglia_score_correlation.csv",
    row.names = FALSE
  )
  
  cat(
    "\nMacrophage vs microglia score:\n"
  )
  
  print(
    macro_micro_result
  )
}

# ------------------------------------------------------------
# 17. Save R objects for later figures / revision
# ------------------------------------------------------------

saveRDS(
  xcell_scores,
  file =
    "results/xCell2_TCGA_LGG_GBM_scores.rds"
)

saveRDS(
  cor_results,
  file =
    "results/xCell2_HAVCR2_correlations.rds"
)

# ------------------------------------------------------------
# 18. Session information
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file =
    "results/sessionInfo_xCell2_revision.txt"
)

cat(
  "\n===========================================\n"
)

cat(
  "xCell2 REVISION ANALYSIS COMPLETED\n"
)

cat(
  "===========================================\n"
)

cat(
  "Results written to results/ folder.\n"
)
