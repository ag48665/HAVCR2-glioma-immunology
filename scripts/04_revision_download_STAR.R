# Revision analysis: download TCGA-LGG + TCGA-GBM STAR-Counts
# Using GDC Data Transfer Tool (gdc-client)

library(TCGAbiolinks)

# ------------------------------------------------------------
# 1. Make local gdc-client.exe visible to R
# ------------------------------------------------------------

client_dir <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

Sys.setenv(
  PATH = paste(
    client_dir,
    Sys.getenv("PATH"),
    sep = .Platform$path.sep
  )
)

cat("GDC client location:\n")
print(Sys.which("gdc-client"))

if (Sys.which("gdc-client") == "") {
  stop("gdc-client.exe was not found.")
}

# ------------------------------------------------------------
# 2. Query TCGA STAR-Counts
# ------------------------------------------------------------

query <- GDCquery(
  project = c("TCGA-LGG", "TCGA-GBM"),
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"
)

cat(
  "\nFiles found:",
  nrow(getResults(query)),
  "\n"
)

# ------------------------------------------------------------
# 3. Download using GDC client
# Use a NEW directory, separate from failed API download
# ------------------------------------------------------------

download_dir <- "data/GDC_CLIENT"

dir.create(
  download_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

GDCdownload(
  query = query,
  method = "client",
  directory = download_dir
)

cat("\nGDC client download completed.\n")

# ------------------------------------------------------------
# 4. Prepare SummarizedExperiment
# ------------------------------------------------------------

se <- GDCprepare(
  query = query,
  directory = download_dir
)

cat("\nGDCprepare completed.\n")
cat("Samples:", ncol(se), "\n")
cat("Genes:", nrow(se), "\n")

# ------------------------------------------------------------
# 5. Save prepared object
# ------------------------------------------------------------

saveRDS(
  se,
  file = "data/TCGA_LGG_GBM_STAR_counts.rds"
)

cat(
  "\nSTAR-Counts successfully saved to:\n",
  "data/TCGA_LGG_GBM_STAR_counts.rds\n"
)

