#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 03_deseq2_pseudobulk.R
# Pseudobulk differential expression on the integrated Seurat object.
# Adapted from the original DESeq2.Rmd notebook, parametrized for use as a
# Nextflow process.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
  library(Seurat)
  library(DESeq2)
  library(scater)
})

option_list <- list(
  make_option("--integrated_rds", type = "character",
              help = "Path to integrated Seurat object (from step 02)"),
  make_option("--group_var", type = "character", default = "group",
              help = "Metadata column defining the contrast, e.g. group (STRESS vs CNTRL)"),
  make_option("--sample_var", type = "character", default = "orig.ident",
              help = "Metadata column identifying the biological replicate for pseudobulk aggregation"),
  make_option("--outdir", type = "character", default = "results")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$integrated_rds)) stop("--integrated_rds is required")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

message("[03_deseq2_pseudobulk] Loading integrated object")
seurat_obj <- readRDS(opt$integrated_rds)
sce <- as.SingleCellExperiment(seurat_obj)

message("[03_deseq2_pseudobulk] Aggregating counts per ", opt$sample_var)
agg <- aggregateAcrossCells(sce, ids = colData(sce)[[opt$sample_var]])

coldata <- as.data.frame(colData(agg))
coldata[[opt$group_var]] <- factor(coldata[[opt$group_var]])

dds <- DESeqDataSetFromMatrix(
  countData = counts(agg),
  colData = coldata,
  design = as.formula(paste("~", opt$group_var))
)
dds <- DESeq(dds)

res <- results(dds)
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

out_csv <- file.path(opt$outdir, "deseq2_pseudobulk_results.csv")
write.csv(res_df, out_csv, row.names = FALSE)
message("[03_deseq2_pseudobulk] Wrote ", out_csv)
