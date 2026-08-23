# ============================================================
# CHECK: HAVCR2 vs macrophage score
# TCGA + locate CGGA files
# ============================================================

cat("\n=============================================\n")
cat("TCGA HAVCR2-MACROPHAGE CHECK\n")
cat("=============================================\n")

tcga_file <- "results/macrophage_microglia_LGG_GBM_patient_level.csv"

if (!file.exists(tcga_file)) {
  stop("File not found: ", tcga_file)
}

tcga <- read.csv(
  tcga_file,
  stringsAsFactors = FALSE
)

cat("\nTCGA rows:", nrow(tcga), "\n")

cat("\nTCGA projects:\n")
print(table(tcga$project_id))

ok <- complete.cases(
  tcga$HAVCR2,
  tcga$macrophage_score
)

test_tcga <- cor.test(
  tcga$HAVCR2[ok],
  tcga$macrophage_score[ok],
  method = "spearman",
  exact = FALSE
)

cat("\n=============================================\n")
cat("TCGA RESULT\n")
cat("=============================================\n")

cat("n =", sum(ok), "\n")
cat("Spearman rho =", unname(test_tcga$estimate), "\n")
cat("p =", format(test_tcga$p.value, scientific = TRUE), "\n")


# ------------------------------------------------------------
# Find possible CGGA files
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("CGGA FILES FOUND\n")
cat("=============================================\n")

cgga_files <- list.files(
  path = ".",
  pattern = "CGGA|cgga|693",
  recursive = TRUE,
  full.names = TRUE
)

print(cgga_files)

cat("\n=============================================\n")
cat("CHECK COMPLETED\n")
cat("=============================================\n")
