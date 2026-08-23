# ============================================================
# Revision figure:
# HAVCR2, macrophage/microglia signatures, xCell2 and GSEA
#
# Panels:
# A - HAVCR2 vs macrophage score
# B - HAVCR2 vs microglia score
# C - xCell2 myeloid/macrophage correlations in LGG vs GBM
# D - GSEA immune pathways in LGG vs GBM
# ============================================================

# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

required <- c(
  "ggplot2",
  "patchwork"
)

missing <- required[
  !vapply(
    required,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing)) {
  install.packages(missing)
}

library(ggplot2)
library(patchwork)


# ------------------------------------------------------------
# 2. Load patient-level macrophage/microglia data
# ------------------------------------------------------------

patient_data <- read.csv(
  "results/macrophage_microglia_LGG_GBM_patient_level.csv",
  stringsAsFactors = FALSE
)

cat("\nPatient-level data:\n")
print(table(patient_data$project_id))


# ------------------------------------------------------------
# 3. Load xCell2 stratified results
# ------------------------------------------------------------

xcell <- read.csv(
  "results/xCell2_HAVCR2_LGG_GBM_myeloid.csv",
  stringsAsFactors = FALSE
)

cat("\nxCell2 myeloid results:\n")
print(xcell)


# ------------------------------------------------------------
# 4. Load GSEA summary
# ------------------------------------------------------------

gsea <- read.csv(
  "results/GSEA_HAVCR2_immune_pathways_summary.csv",
  stringsAsFactors = FALSE
)

gsea <- gsea[
  gsea$cohort %in% c("LGG", "GBM"),
  ,
  drop = FALSE
]


# ============================================================
# PANEL A
# HAVCR2 vs macrophage marker score
# ============================================================

pA <- ggplot(
  patient_data,
  aes(
    x = macrophage_score,
    y = HAVCR2,
    shape = project_id
  )
) +
  geom_point(
    alpha = 0.45,
    size = 1.4
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.7
  ) +
  labs(
    title = "A. HAVCR2 and macrophage signature",
    x = "Macrophage marker score",
    y = "HAVCR2 expression [log2(TPM + 1)]",
    shape = "Cohort"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    legend.position = "bottom"
  )


# ============================================================
# PANEL B
# HAVCR2 vs microglia marker score
# ============================================================

pB <- ggplot(
  patient_data,
  aes(
    x = microglia_score,
    y = HAVCR2,
    shape = project_id
  )
) +
  geom_point(
    alpha = 0.45,
    size = 1.4
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.7
  ) +
  labs(
    title = "B. HAVCR2 and microglia signature",
    x = "Microglia marker score",
    y = "HAVCR2 expression [log2(TPM + 1)]",
    shape = "Cohort"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    legend.position = "bottom"
  )


# ============================================================
# PANEL C
# xCell2 myeloid results
# ============================================================

xcell$cohort <- ifelse(
  xcell$cohort == "TCGA-LGG",
  "LGG",
  "GBM"
)

xcell$cell_type_clean <- gsub(
  "_",
  " ",
  xcell$cell_type
)

xcell$cell_type_clean <- factor(
  xcell$cell_type_clean,
  levels = rev(
    c(
      "myeloid cell",
      "monocyte",
      "macrophage",
      "alternatively activated macrophage",
      "inflammatory macrophage"
    )
  )
)

pC <- ggplot(
  xcell,
  aes(
    x = rho,
    y = cell_type_clean,
    shape = cohort
  )
) +
  geom_point(
    size = 3
  ) +
  labs(
    title = "C. HAVCR2–myeloid associations by xCell2",
    x = "Spearman rho",
    y = NULL,
    shape = "Cohort"
  ) +
  xlim(
    0.5,
    1.0
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    legend.position = "bottom"
  )


# ============================================================
# PANEL D
# GSEA immune pathways
# ============================================================

gsea$pathway_clean <- gsub(
  "^HALLMARK_",
  "",
  gsea$pathway
)

gsea$pathway_clean <- gsub(
  "_",
  " ",
  gsea$pathway_clean
)

gsea_order <- aggregate(
  NES ~ pathway_clean,
  data = gsea,
  FUN = mean
)

gsea_order <- gsea_order[
  order(gsea_order$NES),
]

gsea$pathway_clean <- factor(
  gsea$pathway_clean,
  levels = gsea_order$pathway_clean
)

gsea$minus_log10_FDR <- -log10(
  gsea$padj
)

pD <- ggplot(
  gsea,
  aes(
    x = NES,
    y = pathway_clean,
    size = minus_log10_FDR,
    shape = cohort
  )
) +
  geom_point() +
  labs(
    title = "D. HAVCR2-associated immune pathway enrichment",
    x = "Normalized enrichment score (NES)",
    y = NULL,
    size = "-log10 FDR",
    shape = "Cohort"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    legend.position = "bottom"
  )

# ============================================================
# 5. Improve layout for publication
# ============================================================

# Make titles slightly smaller
pA <- pA +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 11
    ),
    legend.position = "bottom"
  )

pB <- pB +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 11
    ),
    legend.position = "bottom"
  )

pC <- pC +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 11
    ),
    legend.position = "bottom",
    axis.text.y = element_text(size = 9)
  )

pD <- pD +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 11
    ),
    legend.position = "bottom",
    axis.text.y = element_text(size = 9)
  )


# ============================================================
# 6. Combine panels
# ============================================================

top_row <- pA + pB +
  plot_layout(
    ncol = 2,
    guides = "collect"
  )

bottom_row <- pC + pD +
  plot_layout(
    ncol = 2,
    widths = c(0.85, 1.35)
  )

combined_figure <- (
  top_row /
    bottom_row
) +
  plot_layout(
    heights = c(1, 1.15)
  ) +
  plot_annotation(
    title =
      "HAVCR2-associated immune features in TCGA-LGG and TCGA-GBM",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 15,
        hjust = 0.5,
        margin = margin(b = 10)
      )
    )
  )

print(combined_figure)


# ============================================================
# 7. Save publication-quality versions
# ============================================================

dir.create(
  "figures",
  showWarnings = FALSE
)

ggsave(
  filename =
    "figures/Figure_revision_HAVCR2_immune_multiplot_FINAL.png",
  plot = combined_figure,
  width = 15,
  height = 10,
  units = "in",
  dpi = 600,
  limitsize = FALSE
)

ggsave(
  filename =
    "figures/Figure_revision_HAVCR2_immune_multiplot_FINAL.pdf",
  plot = combined_figure,
  width = 15,
  height = 10,
  units = "in",
  limitsize = FALSE
)

cat(
  "\nFINAL combined revision figure saved.\n"
)

cat(
  "PNG exists:",
  file.exists(
    "figures/Figure_revision_HAVCR2_immune_multiplot_FINAL.png"
  ),
  "\n"
)

cat(
  "PDF exists:",
  file.exists(
    "figures/Figure_revision_HAVCR2_immune_multiplot_FINAL.pdf"
  ),
  "\n"
)

figure_AB <- pA + pB +
  plot_layout(
    ncol = 2,
    guides = "collect"
  ) +
  plot_annotation(
    title = "HAVCR2 expression is associated with macrophage and microglia signatures"
  ) &
  theme(
    legend.position = "bottom",
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0.5
    )
  )

print(figure_AB)

ggsave(
  "figures/Figure_revision_HAVCR2_macrophage_microglia.png",
  figure_AB,
  width = 11,
  height = 5.5,
  dpi = 600
)

ggsave(
  "figures/Figure_revision_HAVCR2_macrophage_microglia.pdf",
  figure_AB,
  width = 11,
  height = 5.5
)
figure_CD <- pC + pD +
  plot_layout(
    ncol = 2,
    widths = c(1, 1.5),
    guides = "collect"
  ) +
  plot_annotation(
    title = "HAVCR2-associated myeloid enrichment and immune pathway activation"
  ) &
  theme(
    legend.position = "bottom",
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0.5
    )
  )

print(figure_CD)

ggsave(
  "figures/Figure_revision_HAVCR2_xCell2_GSEA.png",
  figure_CD,
  width = 13,
  height = 6,
  dpi = 600
)

ggsave(
  "figures/Figure_revision_HAVCR2_xCell2_GSEA.pdf",
  figure_CD,
  width = 13,
  height = 6
)
