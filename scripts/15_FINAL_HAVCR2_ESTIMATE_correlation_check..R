# ============================================================
# 15_FINAL_HAVCR2_ESTIMATE_correlation_check.R
#
# FINAL VALIDATION
# HAVCR2 expression vs ESTIMATE immune score
#
# Goal:
#   confirm Spearman rho, p-value, and n
#   for the TCGA glioma cohort
#
# Expected manuscript values:
#   rho = 0.887
#   p < 0.001
#   n = 601
# ============================================================


# ------------------------------------------------------------
# 1. Input file
# ------------------------------------------------------------

infile <- file.path(
  "results",
  "tcga_clinical_clean.csv"
)

if (!file.exists(infile)) {
  stop(
    "Input file not found: ",
    infile
  )
}


# ------------------------------------------------------------
# 2. Load data
# ------------------------------------------------------------

d <- read.csv(
  infile,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n=============================================\n")
cat("FINAL HAVCR2 vs ESTIMATE CHECK\n")
cat("=============================================\n")

cat(
  "\nRows in source file:",
  nrow(d),
  "\n"
)


# ------------------------------------------------------------
# 3. Check required columns
# ------------------------------------------------------------

required_columns <- c(
  "patient",
  "paper_ESTIMATE.immune.score"
)

missing_columns <- setdiff(
  required_columns,
  colnames(d)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required column(s): ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 4. Find HAVCR2 column
# ------------------------------------------------------------

havcr2_candidates <- grep(
  "HAVCR2",
  colnames(d),
  ignore.case = TRUE,
  value = TRUE
)

cat(
  "\nHAVCR2 candidate columns:\n"
)

print(
  havcr2_candidates
)

if (length(havcr2_candidates) == 0) {
  stop(
    "No HAVCR2 column found in tcga_clinical_clean.csv."
  )
}


# ------------------------------------------------------------
# 5. Choose HAVCR2 column
# ------------------------------------------------------------
#
# If there is only one HAVCR2 column, use it.
# If there are several, stop and print them so that
# we do not guess which one was used in the manuscript.
# ------------------------------------------------------------

if (length(havcr2_candidates) > 1) {
  stop(
    paste(
      "More than one HAVCR2 column found:",
      paste(
        havcr2_candidates,
        collapse = " | "
      )
    )
  )
}

havcr2_col <- havcr2_candidates[1]

cat(
  "\nUsing HAVCR2 column:",
  havcr2_col,
  "\n"
)


# ------------------------------------------------------------
# 6. Collapse to patient level
# ------------------------------------------------------------
#
# Prefer rows with both immune score and HAVCR2 available.
# ------------------------------------------------------------

d$.has_both <- (
  !is.na(
    d[["paper_ESTIMATE.immune.score"]]
  ) &
    !is.na(
      d[[havcr2_col]]
    )
)

d$.row_order <- seq_len(
  nrow(d)
)

d <- d[
  order(
    d$patient,
    -as.integer(d$.has_both),
    d$.row_order
  ),
]

d <- d[
  !duplicated(d$patient),
]


cat(
  "Unique patients after patient-level collapse:",
  nrow(d),
  "\n"
)


# ------------------------------------------------------------
# 7. Prepare variables
# ------------------------------------------------------------

d$immune_score <- as.numeric(
  d[["paper_ESTIMATE.immune.score"]]
)

d$HAVCR2 <- as.numeric(
  d[[havcr2_col]]
)


# ------------------------------------------------------------
# 8. Complete cases
# ------------------------------------------------------------

ok <- complete.cases(
  d$immune_score,
  d$HAVCR2
)

cor_dat <- d[
  ok,
  c(
    "patient",
    "immune_score",
    "HAVCR2"
  ),
  drop = FALSE
]


cat("\n=============================================\n")
cat("FINAL CORRELATION COHORT\n")
cat("=============================================\n")

cat(
  "Patients with both HAVCR2 and immune score:",
  nrow(cor_dat),
  "\n"
)


# ------------------------------------------------------------
# 9. Spearman correlation
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
cat("SPEARMAN RESULT\n")
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
  "p-value =",
  format(
    p_value,
    scientific = TRUE,
    digits = 12
  ),
  "\n"
)

cat(
  "Rounded rho (3 decimals) =",
  format(
    round(
      rho,
      3
    ),
    nsmall = 3
  ),
  "\n"
)

cat(
  "p < 0.001 ?",
  p_value < 0.001,
  "\n"
)


# ------------------------------------------------------------
# 10. Save exact cohort
# ------------------------------------------------------------

dir.create(
  "results",
  showWarnings = FALSE,
  recursive = TRUE
)

write.csv(
  cor_dat,
  file =
    "results/HAVCR2_ESTIMATE_correlation_cohort_FINAL.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Save summary
# ------------------------------------------------------------

summary_table <- data.frame(
  analysis =
    "HAVCR2 vs ESTIMATE immune score",
  n = nrow(cor_dat),
  rho = rho,
  p_value = p_value,
  stringsAsFactors = FALSE
)

write.csv(
  summary_table,
  file =
    "results/HAVCR2_ESTIMATE_correlation_summary_FINAL.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 12. Final message
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("HAVCR2 vs ESTIMATE CHECK COMPLETED\n")
cat("=============================================\n")

cat(
  "\nSaved:\n",
  "results/HAVCR2_ESTIMATE_correlation_cohort_FINAL.csv\n",
  "results/HAVCR2_ESTIMATE_correlation_summary_FINAL.csv\n"
)

cat("\n=============================================\n")

