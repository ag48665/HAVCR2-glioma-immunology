# ============================================================
# 14_FINAL_Figure1_TCGA_KM_pretty.R
#
# FINAL PUBLICATION-STYLE FIGURE 1
# Uses already validated KM cohort:
# results/TCGA_FINAL_KM_immune_score_cohort.csv
# ============================================================

required <- c(
  "survival",
  "survminer",
  "ggplot2"
)

missing <- required[
  !vapply(
    required,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing) > 0) {
  stop(
    "Install required package(s): ",
    paste(missing, collapse = ", ")
  )
}

library(survival)
library(survminer)
library(ggplot2)


# ------------------------------------------------------------
# 1. Load final validated KM cohort
# ------------------------------------------------------------

infile <- file.path(
  "results",
  "TCGA_FINAL_KM_immune_score_cohort.csv"
)

if (!file.exists(infile)) {
  stop(
    "File not found: ",
    infile
  )
}

d <- read.csv(
  infile,
  stringsAsFactors = FALSE
)

d$immune_group <- factor(
  d$immune_group,
  levels = c(
    "Low",
    "High"
  )
)


# ------------------------------------------------------------
# 2. Confirm cohort
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("FINAL FIGURE 1 DATA CHECK\n")
cat("=============================================\n")

cat(
  "Total n:",
  nrow(d),
  "\n"
)

print(
  table(
    d$immune_group
  )
)

cat(
  "Deaths:",
  sum(d$event == 1),
  "\n"
)


# ------------------------------------------------------------
# 3. Kaplan-Meier fit
# ------------------------------------------------------------

fit <- survfit(
  Surv(
    os_months,
    event
  ) ~ immune_group,
  data = d
)


# ------------------------------------------------------------
# 4. Exact log-rank test
# ------------------------------------------------------------

lr <- survdiff(
  Surv(
    os_months,
    event
  ) ~ immune_group,
  data = d
)

logrank_p <- pchisq(
  lr$chisq,
  df = 1,
  lower.tail = FALSE
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


# ------------------------------------------------------------
# 5. Publication-style plot
# ------------------------------------------------------------

p <- ggsurvplot(
  fit,
  data = d,
  
  pval = "Log-rank p < 0.0001",
  pval.coord = c(
    8,
    0.18
  ),
  
  censor = TRUE,
  censor.shape = 3,
  censor.size = 2.5,
  
  conf.int = FALSE,
  risk.table = FALSE,
  
  legend.title = NULL,
  legend.labs = c(
    "Low immune score (n = 302)",
    "High immune score (n = 302)"
  ),
  
  xlab = "Overall survival (months)",
  ylab = "Survival probability",
  
  xlim = c(
    0,
    215
  ),
  
  break.time.by = 50,
  
  ggtheme = theme_classic(
    base_size = 13
  ),
  
  palette = c(
    "#00BFC4",
    "#F8766D"
  ),
  
  size = 1.05
)


# ------------------------------------------------------------
# 6. Additional styling
# ------------------------------------------------------------

p$plot <- p$plot +
  theme(
    legend.position = "top",
    legend.justification = "center",
    
    legend.text = element_text(
      size = 11
    ),
    
    axis.title = element_text(
      size = 13
    ),
    
    axis.text = element_text(
      size = 11
    ),
    
    plot.margin = margin(
      15,
      15,
      15,
      15
    )
  )


# ------------------------------------------------------------
# 7. Save PNG
# ------------------------------------------------------------

outfile_png <- file.path(
  "figures",
  "Figure1_FINAL_TCGA_immune_score_KM_publication.png"
)

ggsave(
  filename = outfile_png,
  plot = p$plot,
  width = 7.5,
  height = 5.5,
  units = "in",
  dpi = 600
)


# ------------------------------------------------------------
# 8. Save TIFF
# ------------------------------------------------------------

outfile_tiff <- file.path(
  "figures",
  "Figure1_FINAL_TCGA_immune_score_KM_publication.tiff"
)

ggsave(
  filename = outfile_tiff,
  plot = p$plot,
  width = 7.5,
  height = 5.5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


# ------------------------------------------------------------
# 9. Final message
# ------------------------------------------------------------

cat("\n=============================================\n")
cat("FINAL FIGURE 1 COMPLETED\n")
cat("=============================================\n")

cat(
  "\nSaved:\n",
  outfile_png,
  "\n",
  outfile_tiff,
  "\n"
)
