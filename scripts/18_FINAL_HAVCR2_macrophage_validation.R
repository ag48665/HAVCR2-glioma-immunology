# ============================================================
# 18_FINAL_HAVCR2_macrophage_validation.R
#
# FINAL VALIDATION
# HAVCR2 vs macrophage marker score
#
# TCGA + CGGA
#
# Expected manuscript values:
#
# TCGA:
#   rho = 0.905
#   p < 2.2e-16
#   n = 800
#
# CGGA:
#   rho = 0.860
#   p < 2.2e-16
# ============================================================


# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------

tcga_file <- file.path(
  "results",
  "macrophage_microglia_LGG_GBM_patient_level.csv"
)

cgga_file <- file.path(
  "results",
  "Table4_CGGA_Macrophage_Validation.csv"
)


# ------------------------------------------------------------
# 2. Check files
# ------------------------------------------------------------

if (!file.exists(tcga_file)) {
  stop(
    "TCGA patient-level file not found: ",
    tcga_file
  )
}

if (!file.exists(cgga_file)) {
  stop(
    "CGGA validation file not found: ",
    cgga_file
  )
}


cat("\n=============================================\n")
cat("FINAL HAVCR2-MACROPHAGE VALIDATION\n")
cat("=============================================\n")


# ------------------------------------------------------------
# 3. Load TCGA patient-level data
# ------------------------------------------------------------

tcga <- read.csv(
  tcga_file,
  stringsAsFactors = FALSE
)

cat("\n=============================================\n")
cat("TCGA DATA CHECK\n")
cat("=============================================\n")

cat(
  "Rows:",
  nrow(tcga),
  "\n"
)

cat(
  "\nColumns:\n"
)

print(
  colnames(tcga)
)


required_tcga <- c(
  "patient_id",
  "HAVCR2",
  "macrophage_score"
)

missing_tcga <- setdiff(
  required_tcga,
  colnames(tcga)
)

if (length(missing_tcga) > 0) {
  stop(
    "Missing TCGA column(s): ",
    paste(
      missing_tcga,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 4. TCGA complete cases
# ------------------------------------------------------------

ok_tcga <- complete.cases(
  tcga$HAVCR2,
  tcga$macrophage_score
)

tcga_cor <- tcga[
  ok_tcga,
  ,
  drop = FALSE
]


cat(
  "\nTCGA complete HAVCR2 + macrophage score:",
  nrow(tcga_cor),
  "\n"
)


# ------------------------------------------------------------
# 5. TCGA Spearman correlation
# ------------------------------------------------------------

test_tcga <- suppressWarnings(
  cor.test(
    tcga_cor$HAVCR2,
    tcga_cor$macrophage_score,
    method = "spearman",
    exact = FALSE
  )
)

rho_tcga <- unname(
  test_tcga$estimate
)

p_tcga <- test_tcga$p.value


cat("\n=============================================\n")
cat("TCGA RESULT\n")
cat("=============================================\n")

cat(
  "n =",
  nrow(tcga_cor),
  "\n"
)

cat(
  "Spearman rho =",
  format(
    rho_tcga,
    digits = 12
  ),
  "\n"
)

cat(
  "Rounded rho =",
  format(
    round(
      rho_tcga,
      3
    ),
    nsmall = 3
  ),
  "\n"
)

cat(
  "p =",
  format(
    p_tcga,
    scientific = TRUE,
    digits = 12
  ),
  "\n"
)

cat(
  "p < 2.2e-16 ?",
  p_tcga < 2.2e-16,
  "\n"
)


# ------------------------------------------------------------
# 6. TCGA manuscript check
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("TCGA MANUSCRIPT CHECK\n")
cat("=============================================\n")

cat(
  "n = 800 ?",
  nrow(tcga_cor) == 800,
  "\n"
)

cat(
  "rho rounds to 0.905 ?",
  round(
    rho_tcga,
    3
  ) == 0.905,
  "\n"
)


# ------------------------------------------------------------
# 7. Load CGGA validation result
# ------------------------------------------------------------

cgga <- read.csv(
  cgga_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n=============================================\n")
cat("CGGA VALIDATION FILE\n")
cat("=============================================\n\n")

print(
  cgga
)

cat(
  "\nColumns:\n"
)

print(
  colnames(cgga)
)


# ------------------------------------------------------------
# 8. Detect CGGA columns
# ------------------------------------------------------------

cn <- colnames(
  cgga
)

rho_col <- grep(
  "correlation|rho",
  cn,
  ignore.case = TRUE,
  value = TRUE
)

p_col <- grep(
  "^p$|pvalue|p.value|p_value",
  cn,
  ignore.case = TRUE,
  value = TRUE
)

n_col <- grep(
  "^n$|sample|patient",
  cn,
  ignore.case = TRUE,
  value = TRUE
)


cat("\n=============================================\n")
cat("CGGA COLUMN DETECTION\n")
cat("=============================================\n")

cat(
  "Correlation column:",
  paste(
    rho_col,
    collapse = " | "
  ),
  "\n"
)

cat(
  "P-value column:",
  paste(
    p_col,
    collapse = " | "
  ),
  "\n"
)

cat(
  "n column:",
  paste(
    n_col,
    collapse = " | "
  ),
  "\n"
)


# ------------------------------------------------------------
# 9. Extract CGGA values
# ------------------------------------------------------------

if (length(rho_col) < 1) {
  stop(
    "Could not identify CGGA correlation column."
  )
}

if (length(p_col) < 1) {
  stop(
    "Could not identify CGGA p-value column."
  )
}

rho_cgga <- as.numeric(
  cgga[[rho_col[1]]][1]
)

p_cgga_raw <- as.character(
  cgga[[p_col[1]]][1]
)


if (length(n_col) >= 1) {
  
  n_cgga <- as.numeric(
    cgga[[n_col[1]]][1]
  )
  
} else {
  
  n_cgga <- NA_real_
}


cat("\n=============================================\n")
cat("CGGA RESULT\n")
cat("=============================================\n")

cat(
  "Spearman rho =",
  rho_cgga,
  "\n"
)

cat(
  "Rounded rho =",
  format(
    round(
      rho_cgga,
      3
    ),
    nsmall = 3
  ),
  "\n"
)

cat(
  "p =",
  p_cgga_raw,
  "\n"
)

cat(
  "n =",
  n_cgga,
  "\n"
)


# ------------------------------------------------------------
# 10. CGGA manuscript check
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("CGGA MANUSCRIPT CHECK\n")
cat("=============================================\n")

cat(
  "rho rounds to 0.860 ?",
  round(
    rho_cgga,
    3
  ) == 0.860,
  "\n"
)

cat(
  "Reported p compatible with < 2.2e-16 ?",
  grepl(
    "2.2e-16|2.2.*10",
    p_cgga_raw,
    ignore.case = TRUE
  ),
  "\n"
)


# ------------------------------------------------------------
# 11. Combined final summary
# ------------------------------------------------------------

summary_table <- data.frame(
  cohort = c(
    "TCGA",
    "CGGA"
  ),
  n = c(
    nrow(tcga_cor),
    n_cgga
  ),
  rho = c(
    rho_tcga,
    rho_cgga
  ),
  rho_rounded = c(
    round(
      rho_tcga,
      3
    ),
    round(
      rho_cgga,
      3
    )
  ),
  p_value = c(
    format(
      p_tcga,
      scientific = TRUE,
      digits = 12
    ),
    p_cgga_raw
  ),
  stringsAsFactors = FALSE
)


cat("\n=============================================\n")
cat("FINAL TCGA + CGGA SUMMARY\n")
cat("=============================================\n\n")

print(
  summary_table
)


# ------------------------------------------------------------
# 12. Save outputs
# ------------------------------------------------------------

dir.create(
  "results",
  showWarnings = FALSE,
  recursive = TRUE
)

write.csv(
  summary_table,
  "results/HAVCR2_macrophage_TCGA_CGGA_FINAL_validation.csv",
  row.names = FALSE
)

write.csv(
  tcga_cor[
    ,
    c(
      "patient_id",
      "project_id",
      "HAVCR2",
      "macrophage_score"
    )
  ],
  "results/HAVCR2_macrophage_TCGA_patient_level_FINAL.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. Session information
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file =
    "results/sessionInfo_FINAL_HAVCR2_macrophage_validation.txt"
)


# ------------------------------------------------------------
# 14. Final message
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("HAVCR2-MACROPHAGE VALIDATION COMPLETED\n")
cat("=============================================\n")

cat(
  "\nSaved:\n",
  "results/HAVCR2_macrophage_TCGA_CGGA_FINAL_validation.csv\n",
  "results/HAVCR2_macrophage_TCGA_patient_level_FINAL.csv\n",
  "results/sessionInfo_FINAL_HAVCR2_macrophage_validation.txt\n"
)

cat("\n=============================================\n")
