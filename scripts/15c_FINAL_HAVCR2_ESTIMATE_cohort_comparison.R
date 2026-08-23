# ============================================================
# 15c_FINAL_HAVCR2_ESTIMATE_cohort_comparison.R
#
# GOAL:
# Compare HAVCR2 vs ESTIMATE immune score correlation
# across several plausible TCGA patient subsets.
#
# We will compare:
#   A. All patients with HAVCR2 + immune score
#   B. Patients with HAVCR2 + immune score + survival + event
#   C. Patients with HAVCR2 + immune score + survival + event
#      + age + grade + IDH
#   D. Same as C, but excluding negative survival times
#
# This should identify which cohort reproduces
# the manuscript value:
#   rho = 0.887
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
    "HAVCR2 file not found: ",
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
    "Missing HAVCR2 columns: ",
    paste(
      missing_hav,
      collapse = ", "
    )
  )
}

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
# 3. Load clinical data
# ------------------------------------------------------------

clin <- read.csv(
  clinical_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 4. Check required clinical columns
# ------------------------------------------------------------

required_clin <- c(
  "patient",
  "paper_ESTIMATE.immune.score",
  "paper_Survival..months.",
  "paper_Vital.status..1.dead.",
  "paper_Age..years.at.diagnosis.",
  "paper_Grade",
  "paper_IDH.status"
)

missing_clin <- setdiff(
  required_clin,
  colnames(clin)
)

if (length(missing_clin) > 0) {
  
  cat(
    "\nAvailable clinical columns:\n"
  )
  
  print(
    colnames(clin)
  )
  
  stop(
    "Missing clinical column(s): ",
    paste(
      missing_clin,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 5. Collapse clinical data to one row per patient
# ------------------------------------------------------------
#
# Prefer rows with the most complete information.
# ------------------------------------------------------------

clinical_fields <- c(
  "paper_ESTIMATE.immune.score",
  "paper_Survival..months.",
  "paper_Vital.status..1.dead.",
  "paper_Age..years.at.diagnosis.",
  "paper_Grade",
  "paper_IDH.status"
)

clin$.complete_count <- rowSums(
  !is.na(
    clin[
      ,
      clinical_fields,
      drop = FALSE
    ]
  )
)

clin$.row_order <- seq_len(
  nrow(clin)
)

clin <- clin[
  order(
    clin$patient,
    -clin$.complete_count,
    clin$.row_order
  ),
]

clin <- clin[
  !duplicated(clin$patient),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 6. Prepare clinical variables
# ------------------------------------------------------------

clin$immune_score <- as.numeric(
  clin[["paper_ESTIMATE.immune.score"]]
)

clin$os_months <- as.numeric(
  clin[["paper_Survival..months."]]
)

clin$event <- as.numeric(
  clin[["paper_Vital.status..1.dead."]]
)

clin$age <- as.numeric(
  clin[["paper_Age..years.at.diagnosis."]]
)

clin$grade <- as.character(
  clin[["paper_Grade"]]
)

clin$IDH <- as.character(
  clin[["paper_IDH.status"]]
)


# ------------------------------------------------------------
# 7. Merge HAVCR2 with clinical data
# ------------------------------------------------------------

merged <- merge(
  hav2,
  clin[
    ,
    c(
      "patient",
      "immune_score",
      "os_months",
      "event",
      "age",
      "grade",
      "IDH"
    )
  ],
  by = "patient",
  all = FALSE
)


cat("\n=============================================\n")
cat("MERGED DATA SUMMARY\n")
cat("=============================================\n")

cat(
  "Matched patients:",
  nrow(merged),
  "\n"
)


# ------------------------------------------------------------
# 8. Correlation helper
# ------------------------------------------------------------

run_cor <- function(
    data,
    label
) {
  
  ok <- complete.cases(
    data$HAVCR2,
    data$immune_score
  )
  
  d <- data[
    ok,
    ,
    drop = FALSE
  ]
  
  if (nrow(d) < 3) {
    
    return(
      data.frame(
        cohort = label,
        n = nrow(d),
        rho = NA_real_,
        p_value = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }
  
  test <- suppressWarnings(
    cor.test(
      d$HAVCR2,
      d$immune_score,
      method = "spearman",
      exact = FALSE
    )
  )
  
  data.frame(
    cohort = label,
    n = nrow(d),
    rho = unname(
      test$estimate
    ),
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 9. Cohort A:
# HAVCR2 + immune score only
# ------------------------------------------------------------

cohort_A <- merged[
  complete.cases(
    merged[
      ,
      c(
        "HAVCR2",
        "immune_score"
      )
    ]
  ),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 10. Cohort B:
# HAVCR2 + immune score + OS + event
# ------------------------------------------------------------

cohort_B <- merged[
  complete.cases(
    merged[
      ,
      c(
        "HAVCR2",
        "immune_score",
        "os_months",
        "event"
      )
    ]
  ),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 11. Cohort C:
# HAVCR2 + immune score + OS + event
# + age + grade + IDH
# ------------------------------------------------------------

cohort_C <- merged[
  complete.cases(
    merged[
      ,
      c(
        "HAVCR2",
        "immune_score",
        "os_months",
        "event",
        "age",
        "grade",
        "IDH"
      )
    ]
  ),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 12. Cohort D:
# Same as C, excluding negative OS
# ------------------------------------------------------------

cohort_D <- cohort_C[
  !is.na(cohort_C$os_months) &
    cohort_C$os_months >= 0,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 13. Run all correlations
# ------------------------------------------------------------

results <- rbind(
  run_cor(
    cohort_A,
    "A: HAVCR2 + immune score"
  ),
  run_cor(
    cohort_B,
    "B: + survival + event"
  ),
  run_cor(
    cohort_C,
    "C: + age + grade + IDH"
  ),
  run_cor(
    cohort_D,
    "D: C + non-negative OS"
  )
)


# ------------------------------------------------------------
# 14. Add rounded values
# ------------------------------------------------------------

results$rho_3dp <- round(
  results$rho,
  3
)

results$p_scientific <- format(
  results$p_value,
  scientific = TRUE,
  digits = 6
)


# ------------------------------------------------------------
# 15. Print comparison
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("HAVCR2 vs ESTIMATE COHORT COMPARISON\n")
cat("=============================================\n\n")

print(
  results[
    ,
    c(
      "cohort",
      "n",
      "rho",
      "rho_3dp",
      "p_value"
    )
  ],
  digits = 12
)


# ------------------------------------------------------------
# 16. Extra check:
# Which cohort is closest to manuscript rho = 0.887
# ------------------------------------------------------------

results$distance_from_0_887 <- abs(
  results$rho - 0.887
)

results_sorted <- results[
  order(
    results$distance_from_0_887
  ),
  ,
  drop = FALSE
]


cat("\n=============================================\n")
cat("CLOSEST TO MANUSCRIPT rho = 0.887\n")
cat("=============================================\n\n")

print(
  results_sorted[
    ,
    c(
      "cohort",
      "n",
      "rho",
      "rho_3dp",
      "distance_from_0_887"
    )
  ],
  digits = 12
)


# ------------------------------------------------------------
# 17. Additional counts
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("COHORT COUNTS\n")
cat("=============================================\n")

cat(
  "A n =",
  nrow(cohort_A),
  "\n"
)

cat(
  "B n =",
  nrow(cohort_B),
  "\n"
)

cat(
  "C n =",
  nrow(cohort_C),
  "\n"
)

cat(
  "D n =",
  nrow(cohort_D),
  "\n"
)

cat(
  "\nNegative OS in cohort C:",
  sum(
    cohort_C$os_months < 0,
    na.rm = TRUE
  ),
  "\n"
)


# ------------------------------------------------------------
# 18. Save results
# ------------------------------------------------------------

dir.create(
  "results",
  showWarnings = FALSE,
  recursive = TRUE
)

write.csv(
  results_sorted,
  file =
    "results/HAVCR2_ESTIMATE_cohort_comparison_FINAL.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 19. Save each cohort
# ------------------------------------------------------------

write.csv(
  cohort_A,
  "results/HAVCR2_ESTIMATE_cohort_A_FINAL.csv",
  row.names = FALSE
)

write.csv(
  cohort_B,
  "results/HAVCR2_ESTIMATE_cohort_B_FINAL.csv",
  row.names = FALSE
)

write.csv(
  cohort_C,
  "results/HAVCR2_ESTIMATE_cohort_C_FINAL.csv",
  row.names = FALSE
)

write.csv(
  cohort_D,
  "results/HAVCR2_ESTIMATE_cohort_D_FINAL.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 20. Final summary
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("COHORT COMPARISON COMPLETED\n")
cat("=============================================\n")

cat(
  "\nResults saved to:\n",
  "results/HAVCR2_ESTIMATE_cohort_comparison_FINAL.csv\n"
)

cat("\n=============================================\n")
