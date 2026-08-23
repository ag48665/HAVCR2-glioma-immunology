# ============================================================
# 19_FINAL_CGGA_HAVCR2_grade_IDH_validation.R
#
# Validation of HAVCR2 associations with:
#   1. Tumor grade
#   2. IDH mutation status
#
# CGGA mRNAseq_693
# ============================================================


# ------------------------------------------------------------
# 1. File paths
# ------------------------------------------------------------

clinical_zip <- paste0(
  "C:/Users/PC/Downloads/",
  "CGGA.mRNAseq_693_clinical.20200506.txt.zip"
)

expr_zip <- paste0(
  "C:/Users/PC/Downloads/",
  "CGGA.mRNAseq_693.RSEM-genes.20200506.txt.zip"
)

if (!file.exists(clinical_zip)) {
  stop("Clinical CGGA ZIP file not found.")
}

if (!file.exists(expr_zip)) {
  stop("Expression CGGA ZIP file not found.")
}


# ------------------------------------------------------------
# 2. Identify TXT files inside ZIP
# ------------------------------------------------------------

clinical_inside <- unzip(
  clinical_zip,
  list = TRUE
)$Name

clinical_inside <- clinical_inside[
  grepl("\\.txt$", clinical_inside) &
    !grepl("__MACOSX", clinical_inside)
]

expr_inside <- unzip(
  expr_zip,
  list = TRUE
)$Name

expr_inside <- expr_inside[
  grepl("\\.txt$", expr_inside) &
    !grepl("__MACOSX", expr_inside)
]


# ------------------------------------------------------------
# 3. Read data
# ------------------------------------------------------------

clinical <- read.delim(
  unz(
    clinical_zip,
    clinical_inside[1]
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

expr <- read.delim(
  unz(
    expr_zip,
    expr_inside[1]
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 4. Extract HAVCR2
# ------------------------------------------------------------

havcr2_row <- grep(
  "^HAVCR2$",
  expr[[1]],
  ignore.case = TRUE
)

if (length(havcr2_row) != 1) {
  stop("HAVCR2 could not be uniquely identified.")
}

havcr2 <- data.frame(
  CGGA_ID = colnames(expr)[-1],
  HAVCR2 = as.numeric(
    expr[havcr2_row, -1]
  ),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 5. Merge
# ------------------------------------------------------------

d <- merge(
  clinical,
  havcr2,
  by = "CGGA_ID",
  all = FALSE
)

cat("\n=============================================\n")
cat("CGGA HAVCR2: GRADE AND IDH VALIDATION\n")
cat("=============================================\n")

cat(
  "Matched patients:",
  nrow(d),
  "\n"
)


# ------------------------------------------------------------
# 6. Clean grade
# ------------------------------------------------------------

d$Grade_clean <- gsub(
  "^WHO\\s*",
  "",
  trimws(d$Grade)
)

d$Grade_clean <- factor(
  d$Grade_clean,
  levels = c(
    "II",
    "III",
    "IV"
  )
)


# ------------------------------------------------------------
# 7. Clean IDH status
# ------------------------------------------------------------

d$IDH_clean <- trimws(
  d$IDH_mutation_status
)

d$IDH_clean <- factor(
  d$IDH_clean,
  levels = c(
    "Mutant",
    "Wildtype"
  )
)


# ============================================================
# A. HAVCR2 vs tumor grade
# ============================================================

grade_dat <- d[
  complete.cases(
    d[, c(
      "HAVCR2",
      "Grade_clean"
    )]
  ),
]

cat("\n=============================================\n")
cat("HAVCR2 vs TUMOR GRADE\n")
cat("=============================================\n")

cat(
  "N:",
  nrow(grade_dat),
  "\n"
)

cat("\nGrade counts:\n")
print(
  table(
    grade_dat$Grade_clean
  )
)

cat("\nHAVCR2 summary by grade:\n")

print(
  tapply(
    grade_dat$HAVCR2,
    grade_dat$Grade_clean,
    summary
  )
)

cat("\nMedian HAVCR2 by grade:\n")

grade_medians <- tapply(
  grade_dat$HAVCR2,
  grade_dat$Grade_clean,
  median,
  na.rm = TRUE
)

print(
  grade_medians
)

cat("\nIQR HAVCR2 by grade:\n")

grade_iqr <- tapply(
  grade_dat$HAVCR2,
  grade_dat$Grade_clean,
  quantile,
  probs = c(
    0.25,
    0.75
  ),
  na.rm = TRUE
)

print(
  grade_iqr
)

# Three groups -> Kruskal-Wallis
grade_test <- kruskal.test(
  HAVCR2 ~ Grade_clean,
  data = grade_dat
)

cat("\nKruskal-Wallis test:\n")
print(
  grade_test
)


# ============================================================
# B. HAVCR2 vs IDH mutation status
# ============================================================

idh_dat <- d[
  complete.cases(
    d[, c(
      "HAVCR2",
      "IDH_clean"
    )]
  ),
]

cat("\n=============================================\n")
cat("HAVCR2 vs IDH STATUS\n")
cat("=============================================\n")

cat(
  "N:",
  nrow(idh_dat),
  "\n"
)

cat("\nIDH counts:\n")

print(
  table(
    idh_dat$IDH_clean
  )
)

cat("\nHAVCR2 summary by IDH:\n")

print(
  tapply(
    idh_dat$HAVCR2,
    idh_dat$IDH_clean,
    summary
  )
)

cat("\nMedian HAVCR2 by IDH:\n")

idh_medians <- tapply(
  idh_dat$HAVCR2,
  idh_dat$IDH_clean,
  median,
  na.rm = TRUE
)

print(
  idh_medians
)

cat("\nIQR HAVCR2 by IDH:\n")

idh_iqr <- tapply(
  idh_dat$HAVCR2,
  idh_dat$IDH_clean,
  quantile,
  probs = c(
    0.25,
    0.75
  ),
  na.rm = TRUE
)

print(
  idh_iqr
)

# Two groups -> Wilcoxon rank-sum test
idh_test <- wilcox.test(
  HAVCR2 ~ IDH_clean,
  data = idh_dat,
  exact = FALSE
)

cat("\nWilcoxon rank-sum test:\n")
print(
  idh_test
)


# ------------------------------------------------------------
# 8. Save validation results
# ------------------------------------------------------------

dir.create(
  "results",
  showWarnings = FALSE
)

results <- data.frame(
  Analysis = c(
    "HAVCR2 vs tumor grade",
    "HAVCR2 vs IDH status"
  ),
  Test = c(
    "Kruskal-Wallis",
    "Wilcoxon rank-sum"
  ),
  N = c(
    nrow(grade_dat),
    nrow(idh_dat)
  ),
  p_value = c(
    grade_test$p.value,
    idh_test$p.value
  )
)

write.csv(
  results,
  "results/CGGA_HAVCR2_grade_IDH_validation.csv",
  row.names = FALSE
)

cat("\n=============================================\n")
cat("FINAL RESULTS\n")
cat("=============================================\n")

print(
  results
)

cat(
  "\nSaved:\n",
  "results/CGGA_HAVCR2_grade_IDH_validation.csv\n"
)