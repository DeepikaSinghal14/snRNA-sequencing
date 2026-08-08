#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 02_seurat_integration.R
# Reads CellBender-denoised matrices for all samples in a samplesheet,
# builds per-sample Seurat objects annotated with sex/group metadata,
# normalizes, and integrates them.
# Adapted from the original Seurat_integration.R script, parametrized for
# use as a Nextflow process (loops over the samplesheet instead of
# hardcoding four local paths).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(Seurat)
  library(scCustomize)
})

option_list <- list(
  make_option("--samplesheet", type = "character",
              help = "CSV with columns: sample_id,cellbender_h5,sex,group"),
  make_option("--min_features", type = "integer", default = 500),
  make_option("--nfeatures", type = "integer", default = 2000),
  make_option("--outdir", type = "character", default = "results")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$samplesheet)) stop("--samplesheet is required")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
sheet <- read.csv(opt$samplesheet, stringsAsFactors = FALSE)
required_cols <- c("sample_id", "cellbender_h5", "sex", "group")
missing_cols <- setdiff(required_cols, colnames(sheet))
if (length(missing_cols) > 0) {
  stop("samplesheet is missing required columns: ", paste(missing_cols, collapse = ", "))
}

message("[02_seurat_integration] Building Seurat objects for ", nrow(sheet), " samples")

seurat_list <- lapply(seq_len(nrow(sheet)), function(i) {
  row <- sheet[i, ]
  message("  - ", row$sample_id, " (", row$sex, "/", row$group, ")")
  mat <- Read_CellBender_h5_Mat(file_name = row$cellbender_h5)
  obj <- CreateSeuratObject(counts = mat, project = row$sample_id, min.cells = 5)
  obj$sex <- row$sex
  obj$group <- row$group
  obj <- subset(obj, subset = nFeature_RNA > opt$min_features)
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst",
                               nfeatures = opt$nfeatures, verbose = FALSE)
  obj
})
names(seurat_list) <- sheet$sample_id

message("[02_seurat_integration] Finding integration anchors")
anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:20)

message("[02_seurat_integration] Integrating data")
integrated <- IntegrateData(anchorset = anchors, dims = 1:20)

out_rds <- file.path(opt$outdir, "sex.combined.rds")
saveRDS(integrated, out_rds)
message("[02_seurat_integration] Wrote ", out_rds)
