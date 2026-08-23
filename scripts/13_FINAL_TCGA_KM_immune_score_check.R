# ============================================================
# 13_FINAL_TCGA_KM_immune_score_check.R
#
# FINAL VALIDATION
# Kaplan-Meier overall survival by ESTIMATE immune score
#
# TCGA-LGG + TCGA-GBM combined cohort
#
# Rules:
#   - patient-level analysis
#   - prefer patient row with available immune score
#   - require survival time, event and immune score
#   - exclude negative recorded survival times
#   - median immune score cutoff
#   - High: immune score > median
#   - Low:  immune score <= median
# ============================================================


# ------------------------------------------------------------
# 1. Required package
# ------------------------------------------------------------

if (!requireNamespace(
  "survival",
  quietly = TRUE
)) {
  stop(
    "Install the survival package first."
  )
}

library(survival)


# ------------------------------------------------------------
# 2. Input/output files
# ------------------------------------------------------------

infile <- file.path(
  "results",
  "tcga_clinical_clean.csv"
)

outfile_cohort <- file.path(
  "results",
  "TCGA_FINAL_KM_immune_score_cohort.csv"
)

outfile_summary <- file.path(
  "results",
  "Figure1_FINAL_KM_summary.csv"
)

outfile_figure <- file.path(
  "figures",
  "Figure1_FINAL_TCGA_immune_score_KM.png"
)

outfile_session <- file.path(
  "results",
  "sessionInfo_FINAL_TCGA_KM.txt"
)


# ------------------------------------------------------------
# 3. Check input
# ------------------------------------------------------------

if (!file.exists(infile)) {
  stop(
    "Input file not found: ",
    infile
  )
}


# ------------------------------------------------------------
# 4. Load data
# ------------------------------------------------------------

d <- read.csv(
  infile,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n=============================================\n")
cat("FINAL TCGA KAPLAN-MEIER VALIDATION\n")
cat("=============================================\n")

cat(
  "\nRows in source file:",
  nrow(d),
  "\n"
)


# ------------------------------------------------------------
# 5. Check required columns
# ------------------------------------------------------------

required_columns <- c(
  "patient",
  "paper_Survival..months.",
  "paper_Vital.status..1.dead.",
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
# 6. Collapse to patient level
# ------------------------------------------------------------
#
# Prefer a row with available ESTIMATE immune score.
# ------------------------------------------------------------

d$.has_immune <- !is.na(
  d[["paper_ESTIMATE.immune.score"]]
)

d$.row_order <- seq_len(
  nrow(d)
)

d <- d[
  order(
    d$patient,
    -as.integer(d$.has_immune),
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

d$os_months <- as.numeric(
  d[["paper_Survival..months."]]
)

d$event <- as.numeric(
  d[["paper_Vital.status..1.dead."]]
)

d$immune_score <- as.numeric(
  d[["paper_ESTIMATE.immune.score"]]
)


# ------------------------------------------------------------
# 8. Complete cases required for KM
# ------------------------------------------------------------
#
# IMPORTANT:
# KM does NOT require age, grade or IDH.
# Therefore its final n may differ from the Cox n = 598.
# ------------------------------------------------------------

cc <- complete.cases(
  d[
    ,
    c(
      "os_months",
      "event",
      "immune_score"
    )
  ]
)

kmdat <- d[
  cc,
  ,
  drop = FALSE
]


cat(
  "\nComplete cases for KM before survival-time validation:",
  nrow(kmdat),
  "\n"
)


# ------------------------------------------------------------
# 9. Validate event variable
# ------------------------------------------------------------

if (!all(
  kmdat$event %in% c(0, 1)
)) {
  stop(
    "Event variable contains values other than 0 and 1."
  )
}


# ------------------------------------------------------------
# 10. Identify negative survival times
# ------------------------------------------------------------

negative_survival <- (
  !is.na(kmdat$os_months) &
    kmdat$os_months < 0
)


cat(
  "Patients with negative survival time:",
  sum(negative_survival),
  "\n"
)


if (any(negative_survival)) {
  
  cat(
    "\nExcluded negative survival records:\n"
  )
  
  print(
    kmdat[
      negative_survival,
      c(
        "patient",
        "os_months",
        "event",
        "immune_score"
      ),
      drop = FALSE
    ]
  )
}


# ------------------------------------------------------------
# 11. Exclude negative survival times
# ------------------------------------------------------------

kmdat <- kmdat[
  !negative_survival,
  ,
  drop = FALSE
]


cat(
  "\nPatients remaining for final KM:",
  nrow(kmdat),
  "\n"
)

cat(
  "Deaths:",
  sum(kmdat$event == 1),
  "\n"
)

cat(
  "Censored:",
  sum(kmdat$event == 0),
  "\n"
)


# ------------------------------------------------------------
# 12. Median immune-score cutoff
# ------------------------------------------------------------

immune_median <- median(
  kmdat$immune_score,
  na.rm = TRUE
)


kmdat$immune_group <- ifelse(
  kmdat$immune_score > immune_median,
  "High",
  "Low"
)

kmdat$immune_group <- factor(
  kmdat$immune_group,
  levels = c(
    "Low",
    "High"
  )
)


cat("\n=============================================\n")
cat("FINAL KM COHORT\n")
cat("=============================================\n")

cat(
  "Patients:",
  nrow(kmdat),
  "\n"
)

cat(
  "Deaths:",
  sum(kmdat$event == 1),
  "\n"
)

cat(
  "Median immune score cutoff:",
  format(
    immune_median,
    digits = 12
  ),
  "\n"
)

cat(
  "\nImmune groups:\n"
)

print(
  table(
    kmdat$immune_group
  )
)


cat(
  "\nDeaths by immune group:\n"
)

print(
  with(
    kmdat,
    table(
      immune_group,
      event
    )
  )
)


# ------------------------------------------------------------
# 13. Kaplan-Meier model
# ------------------------------------------------------------

km_fit <- survfit(
  Surv(
    os_months,
    event
  ) ~ immune_group,
  data = kmdat
)


# ------------------------------------------------------------
# 14. Log-rank test
# ------------------------------------------------------------

lr <- survdiff(
  Surv(
    os_months,
    event
  ) ~ immune_group,
  data = kmdat
)

logrank_p <- pchisq(
  lr$chisq,
  df = length(lr$n) - 1,
  lower.tail = FALSE
)


cat("\n=============================================\n")
cat("LOG-RANK RESULT\n")
cat("=============================================\n")

cat(
  "Chi-square:",
  format(
    lr$chisq,
    digits = 12
  ),
  "\n"
)

cat(
  "Exact log-rank p-value:",
  format(
    logrank_p,
    scientific = TRUE,
    digits = 12
  ),
  "\n"
)

cat(
  "p < 0.0001 ?",
  logrank_p < 0.0001,
  "\n"
)


# ------------------------------------------------------------
# 15. Median survival by group
# ------------------------------------------------------------

km_table <- summary(
  km_fit
)$table

cat("\n=============================================\n")
cat("KAPLAN-MEIER SUMMARY\n")
cat("=============================================\n\n")

print(
  km_table
)


# ------------------------------------------------------------
# 16. Save exact KM cohort
# ------------------------------------------------------------

dir.create(
  "results",
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  "figures",
  showWarnings = FALSE,
  recursive = TRUE
)


write.csv(
  kmdat,
  outfile_cohort,
  row.names = FALSE
)


# ------------------------------------------------------------
# 17. Save summary table
# ------------------------------------------------------------

group_counts <- table(
  kmdat$immune_group
)

group_deaths <- tapply(
  kmdat$event,
  kmdat$immune_group,
  sum
)

summary_table <- data.frame(
  analysis = "TCGA-LGG+GBM Kaplan-Meier",
  total_n = nrow(kmdat),
  total_deaths = sum(kmdat$event == 1),
  immune_score_median = immune_median,
  low_n = unname(group_counts["Low"]),
  high_n = unname(group_counts["High"]),
  low_deaths = unname(group_deaths["Low"]),
  high_deaths = unname(group_deaths["High"]),
  logrank_chisq = unname(lr$chisq),
  logrank_p = logrank_p,
  excluded_negative_OS = sum(negative_survival),
  stringsAsFactors = FALSE
)

write.csv(
  summary_table,
  outfile_summary,
  row.names = FALSE
)


# ------------------------------------------------------------
# 18. Create final Figure 1
# ------------------------------------------------------------
#
# Base R survival plot.
# This is primarily for validation of the curve.
# ------------------------------------------------------------

png(
  filename = outfile_figure,
  width = 1800,
  height = 1400,
  res = 200
)

plot(
  km_fit,
  lty = c(1, 1),
  lwd = 3,
  mark.time = TRUE,
  xlab = "Overall survival (months)",
  ylab = "Survival probability",
  main = "Overall survival by ESTIMATE immune score",
  xlim = c(
    0,
    max(
      kmdat$os_months,
      na.rm = TRUE
    )
  )
)

legend(
  "topright",
  legend = c(
    paste0(
      "Low (n = ",
      unname(group_counts["Low"]),
      ")"
    ),
    paste0(
      "High (n = ",
      unname(group_counts["High"]),
      ")"
    )
  ),
  lty = 1,
  lwd = 3,
  bty = "n"
)

text(
  x = max(kmdat$os_months, na.rm = TRUE) * 0.08,
  y = 0.20,
  labels = if (
    logrank_p < 0.0001
  ) {
    "Log-rank p < 0.0001"
  } else {
    paste0(
      "Log-rank p = ",
      format(
        logrank_p,
        digits = 3
      )
    )
  },
  adj = 0
)

dev.off()


# ------------------------------------------------------------
# 19. Save session info
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = outfile_session
)


# ------------------------------------------------------------
# 20. Final report
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("FINAL KM VALIDATION COMPLETED\n")
cat("=============================================\n")

cat(
  "\nFinal KM n:",
  nrow(kmdat),
  "\n"
)

cat(
  "Deaths:",
  sum(kmdat$event == 1),
  "\n"
)

cat(
  "Low:",
  unname(group_counts["Low"]),
  "\n"
)

cat(
  "High:",
  unname(group_counts["High"]),
  "\n"
)

cat(
  "Median immune score:",
  format(
    immune_median,
    digits = 12
  ),
  "\n"
)

cat(
  "Exact log-rank p:",
  format(
    logrank_p,
    scientific = TRUE,
    digits = 12
  ),
  "\n"
)

cat(
  "Excluded negative survival times:",
  sum(negative_survival),
  "\n"
)

cat(
  "\nFiles saved:\n",
  outfile_cohort,
  "\n",
  outfile_summary,
  "\n",
  outfile_figure,
  "\n",
  outfile_session,
  "\n"
)

cat("\n=============================================\n")
