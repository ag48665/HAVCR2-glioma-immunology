# Revision analysis: stratified TCGA-LGG and TCGA-GBM analyses
# Addresses Reviewer 2 concern regarding combining LGG and GBM
# Expected working directory: repository root

required <- c("survival")
missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing)) {
  stop(
    "Install required package(s): ",
    paste(missing, collapse = ", ")
  )
}

library(survival)

# ------------------------------------------------------------
# 1. Load the patient-level cohort created in revision step 02
# ------------------------------------------------------------

infile <- file.path(
  "results",
  "TCGA_patient_level_cohort_for_revision.csv"
)

d <- read.csv(
  infile,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 2. Check cohort labels
# ------------------------------------------------------------

cat("\nProject IDs present in the dataset:\n")
print(table(d$project_id, useNA = "ifany"))

# Keep only TCGA-LGG and TCGA-GBM
d <- d[d$project_id %in% c("TCGA-LGG", "TCGA-GBM"), ]

d$project_id <- factor(
  d$project_id,
  levels = c("TCGA-LGG", "TCGA-GBM")
)

cat("\nPatients by cohort:\n")
print(table(d$project_id))

# ------------------------------------------------------------
# 3. Prepare survival variables
# ------------------------------------------------------------

d$os_months <- as.numeric(d$os_months)
d$event <- as.numeric(d$event)
d$immune_score <- as.numeric(d$immune_score)
d$age <- as.numeric(d$age)

# ------------------------------------------------------------
# 4. Descriptive statistics by cohort
# ------------------------------------------------------------

descriptive <- data.frame(
  cohort = c("TCGA-LGG", "TCGA-GBM"),
  n = c(
    sum(d$project_id == "TCGA-LGG"),
    sum(d$project_id == "TCGA-GBM")
  ),
  deaths = c(
    sum(
      d$event[d$project_id == "TCGA-LGG"] == 1,
      na.rm = TRUE
    ),
    sum(
      d$event[d$project_id == "TCGA-GBM"] == 1,
      na.rm = TRUE
    )
  ),
  median_age = c(
    median(
      d$age[d$project_id == "TCGA-LGG"],
      na.rm = TRUE
    ),
    median(
      d$age[d$project_id == "TCGA-GBM"],
      na.rm = TRUE
    )
  ),
  median_immune_score = c(
    median(
      d$immune_score[d$project_id == "TCGA-LGG"],
      na.rm = TRUE
    ),
    median(
      d$immune_score[d$project_id == "TCGA-GBM"],
      na.rm = TRUE
    )
  )
)

write.csv(
  descriptive,
  file.path(
    "results",
    "Table_revision_LGG_GBM_descriptive.csv"
  ),
  row.names = FALSE
)

cat("\nDescriptive summary:\n")
print(descriptive)

# ------------------------------------------------------------
# 5. Helper function for Cox model
# ------------------------------------------------------------

extract_cox <- function(fit, cohort_name) {
  
  s <- summary(fit)
  
  data.frame(
    cohort = cohort_name,
    term = rownames(s$coefficients),
    beta = s$coefficients[, "coef"],
    HR = s$coefficients[, "exp(coef)"],
    CI95_low = s$conf.int[, "lower .95"],
    CI95_high = s$conf.int[, "upper .95"],
    p_value = s$coefficients[, "Pr(>|z|)"],
    row.names = NULL,
    check.names = FALSE
  )
}

# ------------------------------------------------------------
# 6. LGG analysis
# ------------------------------------------------------------

lgg <- d[d$project_id == "TCGA-LGG", ]

lgg_cc <- complete.cases(
  lgg[, c(
    "os_months",
    "event",
    "immune_score",
    "age",
    "IDH",
    "grade"
  )]
)

lgg_cox <- lgg[lgg_cc, ]

cat("\nLGG complete cases:", nrow(lgg_cox), "\n")
cat(
  "LGG events:",
  sum(lgg_cox$event == 1),
  "\n"
)

# In LGG there should be no grade IV,
# but grade is retained if more than one level exists.
lgg_cox$IDH <- droplevels(factor(lgg_cox$IDH))
lgg_cox$grade <- droplevels(factor(lgg_cox$grade))

if (
  nlevels(lgg_cox$IDH) > 1 &&
  nlevels(lgg_cox$grade) > 1
) {
  
  fit_lgg <- coxph(
    Surv(os_months, event) ~
      immune_score + age + IDH + grade,
    data = lgg_cox,
    ties = "efron"
  )
  
} else if (nlevels(lgg_cox$IDH) > 1) {
  
  fit_lgg <- coxph(
    Surv(os_months, event) ~
      immune_score + age + IDH,
    data = lgg_cox,
    ties = "efron"
  )
  
} else {
  
  fit_lgg <- coxph(
    Surv(os_months, event) ~
      immune_score + age,
    data = lgg_cox,
    ties = "efron"
  )
}

res_lgg <- extract_cox(
  fit_lgg,
  "TCGA-LGG"
)

cat("\nLGG Cox results:\n")
print(res_lgg)

# ------------------------------------------------------------
# 7. GBM analysis
# ------------------------------------------------------------

gbm <- d[d$project_id == "TCGA-GBM", ]

gbm_cc <- complete.cases(
  gbm[, c(
    "os_months",
    "event",
    "immune_score",
    "age",
    "IDH"
  )]
)

gbm_cox <- gbm[gbm_cc, ]

cat("\nGBM complete cases:", nrow(gbm_cox), "\n")
cat(
  "GBM events:",
  sum(gbm_cox$event == 1),
  "\n"
)

gbm_cox$IDH <- droplevels(factor(gbm_cox$IDH))

# Grade is not included in the GBM-only model because
# TCGA-GBM is essentially grade IV by definition.
if (nlevels(gbm_cox$IDH) > 1) {
  
  fit_gbm <- coxph(
    Surv(os_months, event) ~
      immune_score + age + IDH,
    data = gbm_cox,
    ties = "efron"
  )
  
} else {
  
  fit_gbm <- coxph(
    Surv(os_months, event) ~
      immune_score + age,
    data = gbm_cox,
    ties = "efron"
  )
}

res_gbm <- extract_cox(
  fit_gbm,
  "TCGA-GBM"
)

cat("\nGBM Cox results:\n")
print(res_gbm)

# ------------------------------------------------------------
# 8. Combine and save Cox results
# ------------------------------------------------------------

stratified_results <- rbind(
  res_lgg,
  res_gbm
)

write.csv(
  stratified_results,
  file.path(
    "results",
    "Table_revision_LGG_GBM_stratified_Cox.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. Immune-score distributions in LGG vs GBM
# ------------------------------------------------------------

wilcox_result <- wilcox.test(
  immune_score ~ project_id,
  data = d,
  exact = FALSE
)

immune_comparison <- data.frame(
  comparison = "TCGA-LGG vs TCGA-GBM",
  median_LGG = median(
    d$immune_score[d$project_id == "TCGA-LGG"],
    na.rm = TRUE
  ),
  median_GBM = median(
    d$immune_score[d$project_id == "TCGA-GBM"],
    na.rm = TRUE
  ),
  p_value = wilcox_result$p.value
)

write.csv(
  immune_comparison,
  file.path(
    "results",
    "Table_revision_LGG_vs_GBM_immune_score.csv"
  ),
  row.names = FALSE
)

cat("\nLGG vs GBM immune-score comparison:\n")
print(immune_comparison)

# ------------------------------------------------------------
# 10. Save exact stratified analysis cohorts
# ------------------------------------------------------------

write.csv(
  lgg_cox,
  file.path(
    "results",
    "TCGA_LGG_revision_cohort.csv"
  ),
  row.names = FALSE
)

write.csv(
  gbm_cox,
  file.path(
    "results",
    "TCGA_GBM_revision_cohort.csv"
  ),
  row.names = FALSE
)

cat("\n--------------------------------------\n")
cat("STRATIFIED ANALYSIS COMPLETED\n")
cat("--------------------------------------\n")

cat(
  "LGG n =",
  nrow(lgg_cox),
  "; events =",
  sum(lgg_cox$event == 1),
  "\n"
)

cat(
  "GBM n =",
  nrow(gbm_cox),
  "; events =",
  sum(gbm_cox$event == 1),
  "\n"
)

cat(
  "\nResults written to the results/ folder.\n"
)
