# ============================================================
# 18_FINAL_CGGA_baseline_characteristics.R
#
# CGGA baseline clinicopathological characteristics
# Cohort: CGGA mRNAseq_693
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


# ------------------------------------------------------------
# 2. Check that files exist
# ------------------------------------------------------------

if (!file.exists(clinical_zip)) {
  stop(
    "Clinical CGGA ZIP file not found:\n",
    clinical_zip
  )
}

if (!file.exists(expr_zip)) {
  stop(
    "Expression CGGA ZIP file not found:\n",
    expr_zip
  )
}


# ------------------------------------------------------------
# 3. Identify files inside ZIP archives
# ------------------------------------------------------------

clinical_inside <- unzip(
  clinical_zip,
  list = TRUE
)$Name

expr_inside <- unzip(
  expr_zip,
  list = TRUE
)$Name

# Exclude __MACOSX metadata files
clinical_inside <- clinical_inside[
  grepl("\\.txt$", clinical_inside) &
    !grepl("__MACOSX", clinical_inside)
]

expr_inside <- expr_inside[
  grepl("\\.txt$", expr_inside) &
    !grepl("__MACOSX", expr_inside)
]

if (length(clinical_inside) != 1) {
  stop(
    "Could not uniquely identify the clinical TXT file inside ZIP."
  )
}

if (length(expr_inside) != 1) {
  stop(
    "Could not uniquely identify the expression TXT file inside ZIP."
  )
}


# ------------------------------------------------------------
# 4. Read clinical data
# ------------------------------------------------------------

clinical <- read.delim(
  unz(
    clinical_zip,
    clinical_inside[1]
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n=============================================\n")
cat("CGGA DATA LOADING\n")
cat("=============================================\n")

cat(
  "Clinical patients:",
  nrow(clinical),
  "\n"
)

cat(
  "Clinical variables:",
  ncol(clinical),
  "\n"
)


# ------------------------------------------------------------
# 5. Read expression data
# ------------------------------------------------------------

expr <- read.delim(
  unz(
    expr_zip,
    expr_inside[1]
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat(
  "Expression genes:",
  nrow(expr),
  "\n"
)

cat(
  "Expression patient columns:",
  ncol(expr) - 1,
  "\n"
)


# ------------------------------------------------------------
# 6. Identify HAVCR2
# ------------------------------------------------------------

havcr2_rows <- grep(
  "^HAVCR2$",
  expr[[1]],
  ignore.case = TRUE
)

if (length(havcr2_rows) != 1) {
  stop(
    "HAVCR2 could not be uniquely identified. ",
    "Number of matching rows: ",
    length(havcr2_rows)
  )
}

cat(
  "HAVCR2 row:",
  havcr2_rows[1],
  "\n"
)


# ------------------------------------------------------------
# 7. Extract HAVCR2 expression
# ------------------------------------------------------------

havcr2 <- data.frame(
  CGGA_ID = colnames(expr)[-1],
  HAVCR2 = as.numeric(
    expr[havcr2_rows[1], -1]
  ),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 8. Merge clinical and HAVCR2 data
# ------------------------------------------------------------

cgga_cox <- merge(
  clinical,
  havcr2,
  by = "CGGA_ID",
  all = FALSE
)

cat(
  "Matched clinical + expression patients:",
  nrow(cgga_cox),
  "\n"
)

cat(
  "HAVCR2 missing after merge:",
  sum(is.na(cgga_cox$HAVCR2)),
  "\n"
)


# ------------------------------------------------------------
# 9. Prepare numeric variables
# ------------------------------------------------------------

cgga_cox$Age <- suppressWarnings(
  as.numeric(cgga_cox$Age)
)

cgga_cox$OS <- suppressWarnings(
  as.numeric(cgga_cox$OS)
)

cgga_cox$event <- suppressWarnings(
  as.numeric(
    cgga_cox[["Censor (alive=0; dead=1)"]]
  )
)


# ============================================================
# CGGA BASELINE CHARACTERISTICS
# ============================================================

cat("\n=============================================\n")
cat("CGGA BASELINE CHARACTERISTICS\n")
cat("=============================================\n")


# ------------------------------------------------------------
# 10. Total cohort
# ------------------------------------------------------------

cat(
  "\nTotal patients:",
  nrow(cgga_cox),
  "\n"
)


# ------------------------------------------------------------
# 11. Age
# ------------------------------------------------------------

cat("\nAGE:\n")

print(
  summary(cgga_cox$Age)
)

age_median <- median(
  cgga_cox$Age,
  na.rm = TRUE
)

age_q1 <- quantile(
  cgga_cox$Age,
  0.25,
  na.rm = TRUE
)

age_q3 <- quantile(
  cgga_cox$Age,
  0.75,
  na.rm = TRUE
)

cat(
  "Median age (IQR):",
  age_median,
  "(",
  age_q1,
  "-",
  age_q3,
  ")\n"
)

cat(
  "Age missing:",
  sum(is.na(cgga_cox$Age)),
  "\n"
)


# ------------------------------------------------------------
# 12. Gender
# ------------------------------------------------------------

cat("\nGENDER:\n")

print(
  table(
    cgga_cox$Gender,
    useNA = "ifany"
  )
)

cat(
  "Gender missing:",
  sum(is.na(cgga_cox$Gender)),
  "\n"
)


# ------------------------------------------------------------
# 13. Tumor grade
# ------------------------------------------------------------

cat("\nGRADE:\n")

print(
  table(
    cgga_cox$Grade,
    useNA = "ifany"
  )
)

cat(
  "Grade missing:",
  sum(is.na(cgga_cox$Grade)),
  "\n"
)


# ------------------------------------------------------------
# 14. IDH mutation status
# ------------------------------------------------------------

cat("\nIDH STATUS:\n")

print(
  table(
    cgga_cox$IDH_mutation_status,
    useNA = "ifany"
  )
)

cat(
  "IDH status missing:",
  sum(is.na(cgga_cox$IDH_mutation_status)),
  "\n"
)


# ------------------------------------------------------------
# 15. 1p/19q codeletion status
# ------------------------------------------------------------

cat("\n1p/19q STATUS:\n")

print(
  table(
    cgga_cox$`1p19q_codeletion_status`,
    useNA = "ifany"
  )
)

cat(
  "1p/19q status missing:",
  sum(
    is.na(
      cgga_cox$`1p19q_codeletion_status`
    )
  ),
  "\n"
)


# ------------------------------------------------------------
# 16. Vital status
# ------------------------------------------------------------

cat("\nVITAL STATUS:\n")

print(
  table(
    cgga_cox[["Censor (alive=0; dead=1)"]],
    useNA = "ifany"
  )
)

cat(
  "Vital status missing:",
  sum(is.na(cgga_cox$event)),
  "\n"
)


# ------------------------------------------------------------
# 17. Overall survival
# ------------------------------------------------------------

cat("\nOVERALL SURVIVAL:\n")

print(
  summary(cgga_cox$OS)
)

cat(
  "OS missing:",
  sum(is.na(cgga_cox$OS)),
  "\n"
)

cat(
  "Negative OS:",
  sum(
    cgga_cox$OS < 0,
    na.rm = TRUE
  ),
  "\n"
)


# ------------------------------------------------------------
# 18. HAVCR2 expression
# ------------------------------------------------------------

cat("\nHAVCR2 EXPRESSION:\n")

print(
  summary(cgga_cox$HAVCR2)
)

havcr2_median <- median(
  cgga_cox$HAVCR2,
  na.rm = TRUE
)

havcr2_q1 <- quantile(
  cgga_cox$HAVCR2,
  0.25,
  na.rm = TRUE
)

havcr2_q3 <- quantile(
  cgga_cox$HAVCR2,
  0.75,
  na.rm = TRUE
)

cat(
  "Median HAVCR2 (IQR):",
  havcr2_median,
  "(",
  havcr2_q1,
  "-",
  havcr2_q3,
  ")\n"
)

cat(
  "HAVCR2 missing:",
  sum(is.na(cgga_cox$HAVCR2)),
  "\n"
)


# ------------------------------------------------------------
# 19. Complete-case cohort for multivariable Cox model
# ------------------------------------------------------------

cgga_cox$Grade_clean <- gsub(
  "^WHO\\s*",
  "",
  cgga_cox$Grade
)

cgga_cox$Grade_clean <- factor(
  cgga_cox$Grade_clean,
  levels = c(
    "II",
    "III",
    "IV"
  )
)

cgga_cox$IDH <- factor(
  cgga_cox$IDH_mutation_status
)

cgga_cox$HAVCR2_SD <- as.numeric(
  scale(cgga_cox$HAVCR2)
)

cox_dat <- cgga_cox[
  complete.cases(
    cgga_cox[, c(
      "OS",
      "event",
      "HAVCR2_SD",
      "Age",
      "Grade_clean",
      "IDH"
    )]
  ) &
    cgga_cox$OS >= 0 &
    cgga_cox$event %in% c(0, 1),
]


# ------------------------------------------------------------
# 20. Patient flow for REMARK reporting
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("CGGA PATIENT FLOW FOR SURVIVAL ANALYSIS\n")
cat("=============================================\n")

cat(
  "Initial CGGA mRNAseq_693 cohort:",
  nrow(clinical),
  "\n"
)

cat(
  "Matched clinical + HAVCR2 expression:",
  nrow(cgga_cox),
  "\n"
)

cat(
  "Final complete-case multivariable Cox cohort:",
  nrow(cox_dat),
  "\n"
)

cat(
  "Deaths in Cox cohort:",
  sum(cox_dat$event == 1),
  "\n"
)

cat(
  "Alive/censored in Cox cohort:",
  sum(cox_dat$event == 0),
  "\n"
)

cat(
  "Excluded from multivariable Cox model:",
  nrow(cgga_cox) - nrow(cox_dat),
  "\n"
)


# ------------------------------------------------------------
# 21. Final Cox cohort distributions
# ------------------------------------------------------------

cat("\nFINAL COX COHORT - GRADE:\n")

print(
  table(
    cox_dat$Grade_clean,
    useNA = "ifany"
  )
)

cat("\nFINAL COX COHORT - IDH:\n")

print(
  table(
    cox_dat$IDH,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# 22. Save output
# ------------------------------------------------------------

dir.create(
  "results",
  showWarnings = FALSE
)

capture.output(
  {
    cat("CGGA BASELINE AND PATIENT FLOW\n\n")
    
    cat(
      "Initial cohort:",
      nrow(clinical),
      "\n"
    )
    
    cat(
      "Matched clinical + expression:",
      nrow(cgga_cox),
      "\n"
    )
    
    cat(
      "Final Cox cohort:",
      nrow(cox_dat),
      "\n"
    )
    
    cat(
      "Deaths:",
      sum(cox_dat$event == 1),
      "\n"
    )
    
    cat(
      "Alive/censored:",
      sum(cox_dat$event == 0),
      "\n\n"
    )
    
    cat(
      "Median age (IQR):",
      age_median,
      "(",
      age_q1,
      "-",
      age_q3,
      ")\n"
    )
    
    cat(
      "Median HAVCR2 (IQR):",
      havcr2_median,
      "(",
      havcr2_q1,
      "-",
      havcr2_q3,
      ")\n\n"
    )
    
    cat("Gender:\n")
    print(
      table(
        cgga_cox$Gender,
        useNA = "ifany"
      )
    )
    
    cat("\nGrade:\n")
    print(
      table(
        cgga_cox$Grade,
        useNA = "ifany"
      )
    )
    
    cat("\nIDH:\n")
    print(
      table(
        cgga_cox$IDH_mutation_status,
        useNA = "ifany"
      )
    )
    
    cat("\n1p/19q:\n")
    print(
      table(
        cgga_cox$`1p19q_codeletion_status`,
        useNA = "ifany"
      )
    )
    
    cat("\nVital status:\n")
    print(
      table(
        cgga_cox[["Censor (alive=0; dead=1)"]],
        useNA = "ifany"
      )
    )
  },
  file = "results/CGGA_baseline_characteristics_FINAL.txt"
)


# ------------------------------------------------------------
# 23. Final message
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("CGGA BASELINE ANALYSIS COMPLETED\n")
cat("=============================================\n")

cat(
  "\nSaved:\n",
  "results/CGGA_baseline_characteristics_FINAL.txt\n"
)

cat("\n")
