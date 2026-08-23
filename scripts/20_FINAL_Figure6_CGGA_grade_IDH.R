# ============================================================
# 20_FINAL_Figure6_CGGA_grade_IDH.R
#
# Final Figure 6
# HAVCR2 expression by tumor grade and IDH status
# CGGA mRNAseq_693
# ============================================================

library(ggplot2)

# ------------------------------------------------------------
# 1. Extract HAVCR2 expression
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
  HAVCR2 = as.numeric(expr[havcr2_row, -1]),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 2. Merge clinical and expression data
# ------------------------------------------------------------

d <- merge(
  clinical,
  havcr2,
  by = "CGGA_ID",
  all = FALSE
)

cat("Matched patients:", nrow(d), "\n")

# ------------------------------------------------------------
# 3. Clean grade
# ------------------------------------------------------------

d$Grade_clean <- gsub(
  "^WHO\\s*",
  "",
  trimws(d$Grade)
)

d$Grade_clean <- factor(
  d$Grade_clean,
  levels = c("II", "III", "IV"),
  labels = c("WHO II", "WHO III", "WHO IV")
)

# Complete cases only
grade_dat <- d[
  complete.cases(
    d[, c("HAVCR2", "Grade_clean")]
  ),
]

grade_dat$Grade_clean <- droplevels(
  grade_dat$Grade_clean
)

cat("\nGRADE FIGURE N =", nrow(grade_dat), "\n")
print(table(grade_dat$Grade_clean))

# ------------------------------------------------------------
# 4. Clean IDH
# ------------------------------------------------------------

d$IDH_clean <- trimws(
  d$IDH_mutation_status
)

d$IDH_clean <- factor(
  d$IDH_clean,
  levels = c("Mutant", "Wildtype")
)

# Complete cases only
idh_dat <- d[
  complete.cases(
    d[, c("HAVCR2", "IDH_clean")]
  ),
]

idh_dat$IDH_clean <- droplevels(
  idh_dat$IDH_clean
)

cat("\nIDH FIGURE N =", nrow(idh_dat), "\n")
print(table(idh_dat$IDH_clean))

# ------------------------------------------------------------
# 5. Statistical tests
# ------------------------------------------------------------

grade_test <- kruskal.test(
  HAVCR2 ~ Grade_clean,
  data = grade_dat
)

idh_test <- wilcox.test(
  HAVCR2 ~ IDH_clean,
  data = idh_dat,
  exact = FALSE
)

cat("\nGRADE TEST:\n")
print(grade_test)

cat("\nIDH TEST:\n")
print(idh_test)

# ------------------------------------------------------------
# 6. Grade panel
# ------------------------------------------------------------

p_grade <- ggplot(
  grade_dat,
  aes(
    x = Grade_clean,
    y = HAVCR2
  )
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = 16,
    outlier.size = 1.8
  ) +
  labs(
    title = "HAVCR2 expression across glioma grades",
    x = "Grade",
    y = "HAVCR2"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0,
      face = "plain"
    ),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------
# 7. IDH panel
# ------------------------------------------------------------

p_idh <- ggplot(
  idh_dat,
  aes(
    x = IDH_clean,
    y = HAVCR2
  )
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = 16,
    outlier.size = 1.8
  ) +
  labs(
    title = "HAVCR2 expression by IDH status",
    x = "IDH mutation status",
    y = "HAVCR2"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0,
      face = "plain"
    ),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------
# 8. Combine panels
# ------------------------------------------------------------

if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}

library(patchwork)

figure6 <- p_grade + p_idh +
  plot_layout(
    ncol = 2
  )

print(figure6)

# ------------------------------------------------------------
# 9. Save
# ------------------------------------------------------------

dir.create(
  "results",
  showWarnings = FALSE
)

ggsave(
  filename = "results/Figure6_CGGA_HAVCR2_grade_IDH_FINAL.png",
  plot = figure6,
  width = 11,
  height = 5.5,
  units = "in",
  dpi = 600
)

ggsave(
  filename = "results/Figure6_CGGA_HAVCR2_grade_IDH_FINAL.pdf",
  plot = figure6,
  width = 11,
  height = 5.5,
  units = "in"
)

cat("\nSaved:\n")
cat("results/Figure6_CGGA_HAVCR2_grade_IDH_FINAL.png\n")
cat("results/Figure6_CGGA_HAVCR2_grade_IDH_FINAL.pdf\n")
