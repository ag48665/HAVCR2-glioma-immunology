# ============================================================
# 15b_FINAL_HAVCR2_ESTIMATE_merge_check.R
#
# FINAL VALIDATION
# HAVCR2 expression vs ESTIMATE immune score
#
# Merge:
#   HAVCR2 patient-level expression
#   + TCGA ESTIMATE immune score
#
# Expected manuscript result:
#   rho ~ 0.887
#   n = 601
# ============================================================


# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------

havcr2_file <- file.path(
  "results",
  "macrophage_microglia_LGG_GBM_patient_level.csv"
)

clinical_file <- file.path(
  "results",
  "tcga_clinical_clean.csv"
)

if (!file.exists(havcr2_file)) {
  stop(
    "HAVCR2 patient-level file not found: ",
    havcr2_file
  )
}

if (!file.exists(clinical_file)) {
  stop(
    "Clinical file not found: ",
    clinical_file
  )
}


# ------------------------------------------------------------
# 2. Load HAVCR2 patient-level data
# ------------------------------------------------------------

hav <- read.csv(
  havcr2_file,
  stringsAsFactors = FALSE
)

cat("\n=============================================\n")
cat("HAVCR2 PATIENT-LEVEL DATA\n")
cat("=============================================\n")

cat(
  "Rows:",
  nrow(hav),
  "\n"
)

cat(
  "\nColumns:\n"
)

print(
  colnames(hav)
)

required_hav <- c(
  "patient_id",
  "HAVCR2"
)

missing_hav <- setdiff(
  required_hav,
  colnames(hav)
)

if (length(missing_hav) > 0) {
  stop(
    "Missing HAVCR2 file column(s): ",
    paste(
      missing_hav,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 3. Load clinical/ESTIMATE data
# ------------------------------------------------------------

clin <- read.csv(
  clinical_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n=============================================\n")
cat("CLINICAL / ESTIMATE DATA\n")
cat("=============================================\n")

cat(
  "Rows:",
  nrow(clin),
  "\n"
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
# 4. Prepare clinical immune-score data
# ------------------------------------------------------------
#
# Prefer row with available immune score,
# then collapse to one row per patient.
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


cat(
  "\nUnique clinical patients:",
  nrow(clin2),
  "\n"
)

cat(
  "Patients with immune score:",
  sum(!is.na(clin2$immune_score)),
  "\n"
)


# ------------------------------------------------------------
# 5. Prepare HAVCR2 data
# ------------------------------------------------------------

hav2 <- hav[
  ,
  c(
    "patient_id",
    "HAVCR2"
  ),
  drop = FALSE
]

colnames(hav2)[1] <- "patient"

hav2$HAVCR2 <- as.numeric(
  hav2$HAVCR2
)


# ------------------------------------------------------------
# 6. Merge by patient ID
# ------------------------------------------------------------

merged <- merge(
  hav2,
  clin2,
  by = "patient",
  all = FALSE
)


cat("\n=============================================\n")
cat("MERGED COHORT\n")
cat("=============================================\n")

cat(
  "Matched patients:",
  nrow(merged),
  "\n"
)


# ------------------------------------------------------------
# 7. Complete cases
# ------------------------------------------------------------

ok <- complete.cases(
  merged$HAVCR2,
  merged$immune_score
)

cor_dat <- merged[
  ok,
  ,
  drop = FALSE
]


cat(
  "Complete HAVCR2 + immune score:",
  nrow(cor_dat),
  "\n"
)


# ------------------------------------------------------------
# 8. Spearman correlation
# ------------------------------------------------------------

test <- suppressWarnings(
  cor.test(
    cor_dat$HAVCR2,
    cor_dat$immune_score,
    method = "spearman",
    exact = FALSE
  )
)

rho <- unname(
  test$estimate
)

p_value <- test$p.value


cat("\n=============================================\n")
cat("FINAL SPEARMAN RESULT\n")
cat("=============================================\n")

cat(
  "n =",
  nrow(cor_dat),
  "\n"
)

cat(
  "Spearman rho =",
  format(
    rho,
    digits = 12
  ),
  "\n"
)

cat(
  "Rounded rho =",
  format(
    round(rho, 3),
    nsmall = 3
  ),
  "\n"
)

cat(
  "p-value =",
  format(
    p_value,
    scientific = TRUE,
    digits = 12
  ),
  "\n"
)

cat(
  "p < 0.001 ?",
  p_value < 0.001,
  "\n"
)


# ------------------------------------------------------------
# 9. Save outputs
# ------------------------------------------------------------

write.csv(
  cor_dat,
  "results/HAVCR2_ESTIMATE_patient_level_FINAL.csv",
  row.names = FALSE
)

summary_table <- data.frame(
  n = nrow(cor_dat),
  rho = rho,
  p_value = p_value,
  stringsAsFactors = FALSE
)

write.csv(
  summary_table,
  "results/HAVCR2_ESTIMATE_summary_FINAL.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Final message
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("HAVCR2 vs ESTIMATE FINAL CHECK COMPLETED\n")
cat("=============================================\n")
