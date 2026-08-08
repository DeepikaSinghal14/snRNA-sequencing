#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 04_ensembl_annotation.R
# Annotates a DE results table with ENSEMBL gene metadata via biomaRt.
# Adapted from the original ENSEMBL_annoations.Rmd notebook (which used
# file.choose() for interactive file selection), parametrized for use as a
# Nextflow process.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(biomaRt)
  library(dplyr)
})

option_list <- list(
  make_option("--de_results_csv", type = "character",
              help = "DE results CSV with a 'gene' column (output of step 03)"),
  make_option("--species_dataset", type = "character",
              default = "rnorvegicus_gene_ensembl",
              help = "biomaRt dataset name, e.g. rnorvegicus_gene_ensembl [default]"),
  make_option("--outdir", type = "character", default = "results")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$de_results_csv)) stop("--de_results_csv is required")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

message("[04_ensembl_annotation] Reading DE results")
de_results <- read.csv(opt$de_results_csv, stringsAsFactors = FALSE)

message("[04_ensembl_annotation] Querying Ensembl (", opt$species_dataset, ")")
ensembl <- useEnsembl(biomart = "ensembl", dataset = opt$species_dataset)
gene_info <- getBM(
  attributes = c("ensembl_gene_id_version", "ensembl_gene_id", "chromosome_name",
                 "gene_biotype", "description", "transcript_length", "external_gene_name"),
  mart = ensembl
)
gene_info <- gene_info %>% distinct(ensembl_gene_id, .keep_all = TRUE)

annotated <- de_results %>%
  left_join(gene_info, by = c("gene" = "ensembl_gene_id"))

out_csv <- file.path(opt$outdir, "de_results_annotated.csv")
write.csv(annotated, out_csv, row.names = FALSE)
message("[04_ensembl_annotation] Wrote ", out_csv)
