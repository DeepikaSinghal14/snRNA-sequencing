#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 01_soupx_qc.R
# Ambient RNA removal (SoupX) + initial Seurat object creation, per sample.
# Adapted from the original Analysis_snRNAseq.Rmd exploratory notebook,
# parametrized for use as a Nextflow process.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SoupX)
})

option_list <- list(
  make_option("--cellranger_dir", type = "character",
              help = "Path to a CellRanger outs/ directory for one sample"),
  make_option("--sample_id", type = "character",
              help = "Sample identifier, used to name outputs"),
  make_option("--min_cells", type = "integer", default = 3),
  make_option("--min_features", type = "integer", default = 200),
  make_option("--outdir", type = "character", default = "results")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$cellranger_dir) || is.null(opt$sample_id)) {
  stop("--cellranger_dir and --sample_id are required")
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

message("[01_soupx_qc] Loading 10X data + estimating soup profile for ", opt$sample_id)
sc <- load10X(opt$cellranger_dir)
sc <- autoEstCont(sc)
cleaned_counts <- adjustCounts(sc)

message("[01_soupx_qc] Creating Seurat object")
seurat_obj <- CreateSeuratObject(
  counts = cleaned_counts,
  project = opt$sample_id,
  min.cells = opt$min_cells,
  min.features = opt$min_features
)

out_rds <- file.path(opt$outdir, paste0(opt$sample_id, ".soupx_clean.rds"))
saveRDS(seurat_obj, out_rds)
message("[01_soupx_qc] Wrote ", out_rds)
