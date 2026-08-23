# ============================================================
# Figure: HAVCR2-associated immune Hallmark pathways
# LGG vs GBM
# ============================================================

library(ggplot2)

# ------------------------------------------------------------
# 1. Load summary
# ------------------------------------------------------------

gsea_plot <- read.csv(
  "results/GSEA_HAVCR2_immune_pathways_summary.csv",
  stringsAsFactors = FALSE
)

# Keep only stratified cohorts for main figure
gsea_plot <- gsea_plot[
  gsea_plot$cohort %in% c("LGG", "GBM"),
]

# ------------------------------------------------------------
# 2. Clean pathway labels
# ------------------------------------------------------------

gsea_plot$pathway_clean <- gsub(
  "^HALLMARK_",
  "",
  gsea_plot$pathway
)

gsea_plot$pathway_clean <- gsub(
  "_",
  " ",
  gsea_plot$pathway_clean
)

# Order pathways by average NES
path_order <- aggregate(
  NES ~ pathway_clean,
  data = gsea_plot,
  FUN = mean
)

path_order <- path_order[
  order(path_order$NES),
]

gsea_plot$pathway_clean <- factor(
  gsea_plot$pathway_clean,
  levels = path_order$pathway_clean
)

# ------------------------------------------------------------
# 3. Convert FDR to plotting scale
# ------------------------------------------------------------

gsea_plot$minus_log10_FDR <- -log10(
  gsea_plot$padj
)

# ------------------------------------------------------------
# 4. Plot
# ------------------------------------------------------------

p <- ggplot(
  gsea_plot,
  aes(
    x = NES,
    y = pathway_clean,
    size = minus_log10_FDR,
    shape = cohort
  )
) +
  geom_point() +
  labs(
    title = "HAVCR2-associated immune pathway enrichment",
    subtitle = "Hallmark GSEA stratified by TCGA-LGG and TCGA-GBM",
    x = "Normalized enrichment score (NES)",
    y = NULL,
    size = "-log10 FDR",
    shape = "Cohort"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 10)
  )

print(p)

# ------------------------------------------------------------
# 5. Save publication-quality versions
# ------------------------------------------------------------

ggsave(
  filename =
    "figures/Figure_revision_GSEA_HAVCR2_LGG_GBM.pdf",
  plot = p,
  width = 8,
  height = 5.5
)

ggsave(
  filename =
    "figures/Figure_revision_GSEA_HAVCR2_LGG_GBM.png",
  plot = p,
  width = 8,
  height = 5.5,
  dpi = 600
)

cat(
  "\nGSEA figure saved successfully.\n"
)

print(
  file.exists(
    "figures/Figure_revision_GSEA_HAVCR2_LGG_GBM.png"
  )
)
