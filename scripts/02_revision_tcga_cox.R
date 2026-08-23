# Revision analysis: TCGA multivariable Cox model including age
# Addresses Reviewer 1, Comment 3 (CIX-26-0172)
# Expected working directory: repository root

required <- c("survival")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required package(s): ", paste(missing, collapse = ", "))

library(survival)

infile <- file.path("results", "tcga_clinical_clean.csv")
outfile <- file.path("results", "Table1_REVISED_Multivariable_Cox_with_Age.csv")
cohortfile <- file.path("results", "TCGA_patient_level_cohort_for_revision.csv")

d <- read.csv(infile, check.names = FALSE, stringsAsFactors = FALSE)

# The source file is sample/aliquot-level. Survival analyses must be patient-level.
# Prefer, for each patient, a row with a non-missing ESTIMATE immune score.
d$.has_immune <- !is.na(d[["paper_ESTIMATE.immune.score"]])
d$.row_order <- seq_len(nrow(d))
d <- d[order(d$patient, -as.integer(d$.has_immune), d$.row_order), ]
d <- d[!duplicated(d$patient), ]

# Harmonize variables.
d$age <- as.numeric(d[["paper_Age..years.at.diagnosis."]])
d$os_months <- as.numeric(d[["paper_Survival..months."]])
d$event <- as.numeric(d[["paper_Vital.status..1.dead."]])
d$immune_score <- as.numeric(d[["paper_ESTIMATE.immune.score"]])
d$IDH <- factor(d[["paper_IDH.status"]])
d$grade <- factor(d[["paper_Grade"]], levels = c("G2", "G3", "G4"))

# Complete-case cohort for the prespecified revised model.
cc <- complete.cases(d[, c("os_months", "event", "immune_score", "age", "IDH", "grade")])
coxdat <- d[cc, ]

write.csv(coxdat, cohortfile, row.names = FALSE)

fit <- coxph(
  Surv(os_months, event) ~ immune_score + age + IDH + grade,
  data = coxdat,
  ties = "efron"
)

s <- summary(fit)
res <- data.frame(
  term = rownames(s$coefficients),
  beta = s$coefficients[, "coef"],
  HR = s$coefficients[, "exp(coef)"],
  CI95_low = s$conf.int[, "lower .95"],
  CI95_high = s$conf.int[, "upper .95"],
  p_value = s$coefficients[, "Pr(>|z|)"],
  row.names = NULL,
  check.names = FALSE
)

write.csv(res, outfile, row.names = FALSE)

cat("Unique patients before complete-case filtering:", nrow(d), "\n")
cat("Patients in revised Cox model:", nrow(coxdat), "\n")
cat("Events:", sum(coxdat$event == 1), "\n\n")
print(res)
cat("\nConcordance:", unname(s$concordance[1]), "\n")