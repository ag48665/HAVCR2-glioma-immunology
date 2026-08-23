library(TCGAbiolinks)

query <- GDCquery(
  project = c("TCGA-GBM", "TCGA-LGG"),
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(
  query,
  method = "api",
  files.per.chunk = 20
)

exp_se <- GDCprepare(query)

saveRDS(
  exp_se,
  file = "data/tcga_glioma_expression.rds"
)

cat("DONE\n")
