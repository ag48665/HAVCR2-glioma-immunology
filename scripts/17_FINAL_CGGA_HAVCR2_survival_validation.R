# ============================================================
# 17_FINAL_CGGA_HAVCR2_survival_validation.R
#
# FINAL VALIDATION
# CGGA HAVCR2 survival analysis
#
# Goals:
#   1. Read existing CGGA validation outputs
#   2. Confirm reported univariate Cox HR / CI / p
#   3. Confirm available KM validation information
#   4. Search project for raw CGGA mRNAseq_693 files
#
# Expected manuscript values:
#   KM: log-rank p < 0.0001
#   Cox HR = 1.035
#   95% CI = 1.023–1.047
#   p = 2.61e-08
# ============================================================


# ------------------------------------------------------------
# 1. Required package
# ------------------------------------------------------------

if (!requireNamespace(
  "survival",
  quietly = TRUE
)) {
  stop(
    "Install survival package first."
  )
}

library(survival)


# ------------------------------------------------------------
# 2. Existing result files
# ------------------------------------------------------------

validation_file <- file.path(
  "results",
  "Table3_CGGA_Validation.csv"
)

cox_file <- file.path(
  "results",
  "Table5_HAVCR2_Univariate_Cox_CGGA.csv"
)


cat("\n=============================================\n")
cat("FINAL CGGA HAVCR2 SURVIVAL VALIDATION\n")
cat("=============================================\n")


# ------------------------------------------------------------
# 3. Check files
# ------------------------------------------------------------

cat("\nExisting CGGA result files:\n")

cat(
  "Table3_CGGA_Validation.csv:",
  file.exists(validation_file),
  "\n"
)

cat(
  "Table5_HAVCR2_Univariate_Cox_CGGA.csv:",
  file.exists(cox_file),
  "\n"
)


# ------------------------------------------------------------
# 4. Read general CGGA validation table
# ------------------------------------------------------------

if (file.exists(validation_file)) {
  
  cgga_validation <- read.csv(
    validation_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  cat("\n=============================================\n")
  cat("EXISTING CGGA VALIDATION TABLE\n")
  cat("=============================================\n\n")
  
  print(
    cgga_validation
  )
  
  cat("\nColumns:\n")
  
  print(
    colnames(cgga_validation)
  )
  
} else {
  
  cat(
    "\nTable3_CGGA_Validation.csv not found.\n"
  )
}


# ------------------------------------------------------------
# 5. Read univariate Cox result
# ------------------------------------------------------------

if (file.exists(cox_file)) {
  
  cgga_cox <- read.csv(
    cox_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  cat("\n=============================================\n")
  cat("EXISTING CGGA UNIVARIATE COX TABLE\n")
  cat("=============================================\n\n")
  
  print(
    cgga_cox,
    digits = 12
  )
  
  cat("\nColumns:\n")
  
  print(
    colnames(cgga_cox)
  )
  
} else {
  
  cat(
    "\nTable5_HAVCR2_Univariate_Cox_CGGA.csv not found.\n"
  )
}


# ------------------------------------------------------------
# 6. Try to identify Cox columns automatically
# ------------------------------------------------------------

if (exists("cgga_cox")) {
  
  cn <- colnames(
    cgga_cox
  )
  
  hr_col <- grep(
    "^HR$|hazard",
    cn,
    ignore.case = TRUE,
    value = TRUE
  )
  
  lower_col <- grep(
    "lower|2.5|lcl",
    cn,
    ignore.case = TRUE,
    value = TRUE
  )
  
  upper_col <- grep(
    "upper|97.5|ucl",
    cn,
    ignore.case = TRUE,
    value = TRUE
  )
  
  p_col <- grep(
    "^p$|p.value|p_value|pvalue",
    cn,
    ignore.case = TRUE,
    value = TRUE
  )
  
  
  cat("\n=============================================\n")
  cat("AUTO-DETECTED COX COLUMNS\n")
  cat("=============================================\n")
  
  cat(
    "HR column:",
    paste(
      hr_col,
      collapse = " | "
    ),
    "\n"
  )
  
  cat(
    "Lower CI column:",
    paste(
      lower_col,
      collapse = " | "
    ),
    "\n"
  )
  
  cat(
    "Upper CI column:",
    paste(
      upper_col,
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
  
  
  # ----------------------------------------------------------
  # 7. Print exact Cox values when identifiable
  # ----------------------------------------------------------
  
  if (
    length(hr_col) >= 1 &&
    length(lower_col) >= 1 &&
    length(upper_col) >= 1 &&
    length(p_col) >= 1
  ) {
    
    hr_value <- as.numeric(
      cgga_cox[[hr_col[1]]][1]
    )
    
    lower_value <- as.numeric(
      cgga_cox[[lower_col[1]]][1]
    )
    
    upper_value <- as.numeric(
      cgga_cox[[upper_col[1]]][1]
    )
    
    p_value <- as.numeric(
      cgga_cox[[p_col[1]]][1]
    )
    
    
    cat("\n=============================================\n")
    cat("CGGA COX EXACT VALUES\n")
    cat("=============================================\n")
    
    cat(
      "HR =",
      format(
        hr_value,
        digits = 12
      ),
      "\n"
    )
    
    cat(
      "95% CI =",
      format(
        lower_value,
        digits = 12
      ),
      "to",
      format(
        upper_value,
        digits = 12
      ),
      "\n"
    )
    
    cat(
      "p =",
      format(
        p_value,
        scientific = TRUE,
        digits = 12
      ),
      "\n"
    )
    
    
    # --------------------------------------------------------
    # 8. Compare against manuscript values
    # --------------------------------------------------------
    
    cat("\n=============================================\n")
    cat("MANUSCRIPT COX CHECK\n")
    cat("=============================================\n")
    
    cat(
      "HR rounds to 1.035 ?",
      round(
        hr_value,
        3
      ) == 1.035,
      "\n"
    )
    
    cat(
      "Lower CI rounds to 1.023 ?",
      round(
        lower_value,
        3
      ) == 1.023,
      "\n"
    )
    
    cat(
      "Upper CI rounds to 1.047 ?",
      round(
        upper_value,
        3
      ) == 1.047,
      "\n"
    )
    
    cat(
      "p approximately 2.61e-08 ?",
      isTRUE(
        all.equal(
          p_value,
          2.61e-08,
          tolerance = 0.02
        )
      ),
      "\n"
    )
  }
}


# ------------------------------------------------------------
# 9. Search project for possible raw CGGA files
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("SEARCHING FOR RAW CGGA FILES\n")
cat("=============================================\n")

all_files <- list.files(
  path = ".",
  recursive = TRUE,
  full.names = TRUE
)

cgga_candidates <- all_files[
  grepl(
    "CGGA|cgga|mRNAseq_693|mRNAseq693|693",
    basename(all_files)
  )
]


# Remove already-known figures/scripts/results when possible

raw_candidates <- cgga_candidates[
  !grepl(
    "\\.(png|jpg|jpeg|tiff|pdf)$",
    cgga_candidates,
    ignore.case = TRUE
  )
]


cat(
  "\nPossible CGGA files found:",
  length(raw_candidates),
  "\n\n"
)

print(
  raw_candidates
)


# ------------------------------------------------------------
# 10. Search likely data locations outside project
# ------------------------------------------------------------

external_locations <- c(
  "C:/GDC",
  "C:/CGGA",
  "C:/Users/PC/Downloads",
  "C:/Users/PC/Desktop"
)

external_found <- character(0)


for (loc in external_locations) {
  
  if (dir.exists(loc)) {
    
    tmp <- list.files(
      path = loc,
      pattern =
        "CGGA|cgga|mRNAseq_693|mRNAseq693",
      recursive = TRUE,
      full.names = TRUE
    )
    
    external_found <- c(
      external_found,
      tmp
    )
  }
}


external_found <- unique(
  external_found
)


cat("\n=============================================\n")
cat("EXTERNAL CGGA FILE SEARCH\n")
cat("=============================================\n")

cat(
  "Possible external files:",
  length(external_found),
  "\n\n"
)

print(
  external_found
)


# ------------------------------------------------------------
# 11. Inspect existing validation files for KM information
# ------------------------------------------------------------

if (exists("cgga_validation")) {
  
  cat("\n=============================================\n")
  cat("SEARCHING TABLE FOR KM / LOG-RANK INFORMATION\n")
  cat("=============================================\n")
  
  text_values <- apply(
    cgga_validation,
    1,
    paste,
    collapse = " | "
  )
  
  km_rows <- grep(
    "Kaplan|log.rank|logrank|survival|HAVCR2",
    text_values,
    ignore.case = TRUE,
    value = TRUE
  )
  
  if (length(km_rows) > 0) {
    
    print(
      km_rows
    )
    
  } else {
    
    cat(
      "No explicit KM/log-rank text detected in Table3_CGGA_Validation.csv\n"
    )
  }
}


# ------------------------------------------------------------
# 12. Save validation snapshot
# ------------------------------------------------------------

dir.create(
  "results",
  showWarnings = FALSE,
  recursive = TRUE
)

capture.output(
  {
    cat(
      "CGGA HAVCR2 FINAL SURVIVAL VALIDATION\n\n"
    )
    
    if (exists("cgga_validation")) {
      cat(
        "Table3_CGGA_Validation.csv\n"
      )
      print(cgga_validation)
    }
    
    cat("\n")
    
    if (exists("cgga_cox")) {
      cat(
        "Table5_HAVCR2_Univariate_Cox_CGGA.csv\n"
      )
      print(cgga_cox)
    }
    
    cat(
      "\nPossible raw/project CGGA files:\n"
    )
    
    print(
      raw_candidates
    )
    
    cat(
      "\nPossible external CGGA files:\n"
    )
    
    print(
      external_found
    )
  },
  file =
    "results/CGGA_HAVCR2_survival_validation_FINAL.txt"
)


# ------------------------------------------------------------
# 13. Final message
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("CGGA SURVIVAL VALIDATION CHECK COMPLETED\n")
cat("=============================================\n")

cat(
  "\nSaved:\n",
  "results/CGGA_HAVCR2_survival_validation_FINAL.txt\n"
)

cat(
  "\nIMPORTANT:\n",
  "If a raw CGGA mRNAseq_693 data file is listed above,\n",
  "send me its path/name. We can then recompute KM and Cox\n",
  "directly from patient-level raw data.\n"
)

cat("\n=============================================\n")
