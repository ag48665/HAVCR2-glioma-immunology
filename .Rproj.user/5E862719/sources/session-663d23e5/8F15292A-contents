library(TCGAbiolinks)

query <- GDCquery(
  project = c("TCGA-GBM", "TCGA-LGG"),
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "HTSeq - Counts"
)

GDCdownload(query)

exp_se <- GDCprepare(query)

saveRDS(
  exp_se,
  file = "data/tcga_glioma_expression.rds"
)

cat("DONE\n")
