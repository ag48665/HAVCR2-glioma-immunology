# ============================================================
# 02_FINAL_TCGA_Cox_IDH_stratified.R
#
# FINAL REVISION ANALYSIS
# TCGA multivariable Cox proportional hazards model
#
# Model:
#   Overall survival ~ immune score + age + grade
#   stratified by IDH status
#
# Immune score expressed per 100-unit increase.
#
# Negative survival times are excluded.
# ============================================================


# ------------------------------------------------------------
# 1. Required package
# ------------------------------------------------------------

required <- c("survival")

missing <- required[
  !vapply(
    required,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing)) {
  stop(
    "Install required package(s): ",
    paste(missing, collapse = ", ")
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

outfile_cox <- file.path(
  "results",
  "Table1_FINAL_IDH_stratified_Cox.csv"
)

outfile_ph <- file.path(
  "results",
  "Table1_FINAL_Schoenfeld_test.csv"
)

outfile_cohort <- file.path(
  "results",
  "TCGA_FINAL_Cox_patient_cohort.csv"
)

outfile_excluded <- file.path(
  "results",
  "TCGA_excluded_negative_survival.csv"
)

outfile_session <- file.path(
  "results",
  "sessionInfo_FINAL_TCGA_Cox.txt"
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
cat("FINAL TCGA COX ANALYSIS\n")
cat("=============================================\n")

cat(
  "\nRows in source file:",
  nrow(d),
  "\n"
)


# ------------------------------------------------------------
# 5. Required columns
# ------------------------------------------------------------

required_columns <- c(
  "patient",
  "paper_Age..years.at.diagnosis.",
  "paper_Survival..months.",
  "paper_Vital.status..1.dead.",
  "paper_ESTIMATE.immune.score",
  "paper_IDH.status",
  "paper_Grade"
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
# Source data may contain multiple samples/aliquots
# for the same TCGA patient.
#
# Prefer a row with a non-missing ESTIMATE immune score.
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

d$age <- as.numeric(
  d[["paper_Age..years.at.diagnosis."]]
)

d$os_months <- as.numeric(
  d[["paper_Survival..months."]]
)

d$event <- as.numeric(
  d[["paper_Vital.status..1.dead."]]
)

d$immune_score <- as.numeric(
  d[["paper_ESTIMATE.immune.score"]]
)

d$IDH <- factor(
  d[["paper_IDH.status"]]
)

d$grade <- factor(
  d[["paper_Grade"]],
  levels = c(
    "G2",
    "G3",
    "G4"
  )
)


# ------------------------------------------------------------
# 8. Scale immune score
# ------------------------------------------------------------
#
# HR will correspond to a 100-unit increase
# in ESTIMATE immune score.
# ------------------------------------------------------------

d$immune_score_100 <- (
  d$immune_score / 100
)


# ------------------------------------------------------------
# 9. Complete-case filtering
# ------------------------------------------------------------

cc <- complete.cases(
  d[
    ,
    c(
      "os_months",
      "event",
      "immune_score_100",
      "age",
      "IDH",
      "grade"
    )
  ]
)

coxdat <- d[
  cc,
]


cat(
  "\nComplete cases before survival-time validation:",
  nrow(coxdat),
  "\n"
)


# ------------------------------------------------------------
# 10. Validate event variable
# ------------------------------------------------------------

if (nrow(coxdat) == 0) {
  stop(
    "No complete cases available for Cox analysis."
  )
}

if (!all(
  coxdat$event %in% c(0, 1)
)) {
  stop(
    "Event variable contains values other than 0 and 1."
  )
}


# ------------------------------------------------------------
# 11. Identify negative survival times
# ------------------------------------------------------------

negative_survival <- (
  !is.na(coxdat$os_months) &
    coxdat$os_months < 0
)

cat(
  "\nPatients with negative survival time:",
  sum(negative_survival),
  "\n"
)


# ------------------------------------------------------------
# 12. Save excluded negative-survival records
# ------------------------------------------------------------

if (any(negative_survival)) {
  
  excluded_negative <- coxdat[
    negative_survival,
  ]
  
  cat(
    "\nNegative survival values:\n"
  )
  
  print(
    excluded_negative[
      ,
      c(
        "patient",
        "os_months",
        "event",
        "age",
        "IDH",
        "grade",
        "immune_score"
      )
    ]
  )
  
  write.csv(
    excluded_negative,
    outfile_excluded,
    row.names = FALSE
  )
  
} else {
  
  excluded_negative <- coxdat[
    FALSE,
  ]
  
  cat(
    "No negative survival times detected.\n"
  )
}


# ------------------------------------------------------------
# 13. Exclude negative survival times
# ------------------------------------------------------------

coxdat <- coxdat[
  !negative_survival,
]


cat(
  "\nPatients remaining after survival-time validation:",
  nrow(coxdat),
  "\n"
)


# ------------------------------------------------------------
# 14. Check survival-time range
# ------------------------------------------------------------

cat(
  "\nSurvival-time range after validation:\n"
)

print(
  range(
    coxdat$os_months,
    na.rm = TRUE
  )
)


# ------------------------------------------------------------
# 15. Final cohort information
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("FINAL COX COHORT\n")
cat("=============================================\n")

cat(
  "Patients:",
  nrow(coxdat),
  "\n"
)

cat(
  "Deaths:",
  sum(
    coxdat$event == 1
  ),
  "\n"
)

cat(
  "Censored:",
  sum(
    coxdat$event == 0
  ),
  "\n"
)

cat(
  "\nIDH distribution:\n"
)

print(
  table(
    coxdat$IDH,
    useNA = "ifany"
  )
)

cat(
  "\nGrade distribution:\n"
)

print(
  table(
    coxdat$grade,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# 16. Save exact final cohort
# ------------------------------------------------------------

write.csv(
  coxdat,
  outfile_cohort,
  row.names = FALSE
)


# ------------------------------------------------------------
# 17. FINAL COX MODEL
# ------------------------------------------------------------
#
# IDH is a stratification variable.
#
# Covariates:
#   immune score per 100 units
#   age
#   grade
#
# Model:
#
# Surv(OS, event) ~
#   immune_score_100 +
#   age +
#   grade +
#   strata(IDH)
# ------------------------------------------------------------

fit <- coxph(
  Surv(
    os_months,
    event
  ) ~
    immune_score_100 +
    age +
    grade +
    strata(IDH),
  data = coxdat,
  ties = "efron",
  x = TRUE
)


# ------------------------------------------------------------
# 18. Cox summary
# ------------------------------------------------------------

s <- summary(fit)

cat("\n=============================================\n")
cat("FINAL COX MODEL SUMMARY\n")
cat("=============================================\n\n")

print(s)


# ------------------------------------------------------------
# 19. Prepare final Table 1
# ------------------------------------------------------------

res <- data.frame(
  
  term = rownames(
    s$coefficients
  ),
  
  beta = s$coefficients[
    ,
    "coef"
  ],
  
  HR = s$coefficients[
    ,
    "exp(coef)"
  ],
  
  CI95_low = s$conf.int[
    ,
    "lower .95"
  ],
  
  CI95_high = s$conf.int[
    ,
    "upper .95"
  ],
  
  p_value = s$coefficients[
    ,
    "Pr(>|z|)"
  ],
  
  row.names = NULL,
  check.names = FALSE
)


# ------------------------------------------------------------
# 20. Print Table 1
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("FINAL TABLE 1 RESULTS\n")
cat("=============================================\n\n")

print(
  res,
  digits = 10
)


# ------------------------------------------------------------
# 21. Save Table 1
# ------------------------------------------------------------

write.csv(
  res,
  outfile_cox,
  row.names = FALSE
)


# ------------------------------------------------------------
# 22. Schoenfeld proportional hazards test
# ------------------------------------------------------------

ph_test <- cox.zph(
  fit
)

cat("\n=============================================\n")
cat("SCHOENFELD PROPORTIONAL HAZARDS TEST\n")
cat("=============================================\n\n")

print(
  ph_test
)


# ------------------------------------------------------------
# 23. Save Schoenfeld test
# ------------------------------------------------------------

ph_table <- as.data.frame(
  ph_test$table
)

ph_table$term <- rownames(
  ph_table
)

ph_table <- ph_table[
  ,
  c(
    "term",
    setdiff(
      colnames(ph_table),
      "term"
    )
  )
]

rownames(
  ph_table
) <- NULL

write.csv(
  ph_table,
  outfile_ph,
  row.names = FALSE
)


# ------------------------------------------------------------
# 24. Extract immune-score result
# ------------------------------------------------------------

immune_row <- res[
  res$term == "immune_score_100",
  ,
  drop = FALSE
]

cat("\n=============================================\n")
cat("IMMUNE SCORE RESULT\n")
cat("=============================================\n\n")

print(
  immune_row,
  digits = 10
)


# ------------------------------------------------------------
# 25. Global Schoenfeld p-value
# ------------------------------------------------------------

if ("GLOBAL" %in% rownames(ph_test$table)) {
  
  global_ph_p <- ph_test$table[
    "GLOBAL",
    "p"
  ]
  
  cat(
    "\nGLOBAL Schoenfeld p-value:",
    format(
      global_ph_p,
      digits = 10
    ),
    "\n"
  )
}


# ------------------------------------------------------------
# 26. Concordance
# ------------------------------------------------------------

cat(
  "\nConcordance:",
  unname(
    s$concordance[1]
  ),
  "\n"
)


# ------------------------------------------------------------
# 27. Save session information
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = outfile_session
)


# ------------------------------------------------------------
# 28. Final summary
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("ANALYSIS COMPLETED\n")
cat("=============================================\n")

cat(
  "\nPatients in final model:",
  nrow(coxdat),
  "\n"
)

cat(
  "Deaths:",
  sum(
    coxdat$event == 1
  ),
  "\n"
)

cat(
  "Excluded because of negative survival time:",
  nrow(excluded_negative),
  "\n"
)

cat(
  "\nFiles saved:\n"
)

cat(
  outfile_cox,
  "\n"
)

cat(
  outfile_ph,
  "\n"
)

cat(
  outfile_cohort,
  "\n"
)

if (nrow(excluded_negative) > 0) {
  cat(
    outfile_excluded,
    "\n"
  )
}

cat(
  outfile_session,
  "\n"
)

cat("\n=============================================\n")
