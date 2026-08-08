#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(biomaRt)
})

option_list <- list(
  make_option("--input", type = "character"),
  make_option("--outdir", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(opt$input, pattern = "\\.csv$", full.names = TRUE)

if (length(files) == 0) {
  stop("No CSV differential-expression files found.")
}

ensembl <- useEnsembl(
  biomart = "ensembl",
  dataset = "rnorvegicus_gene_ensembl"
)

gene_info <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "chromosome_name",
    "gene_biotype",
    "description",
    "transcript_length",
    "external_gene_name"
  ),
  mart = ensembl
) |>
  distinct(ensembl_gene_id, .keep_all = TRUE)

for (file in files) {
  dat <- read.csv(file, stringsAsFactors = FALSE)

  if (!"gene" %in% names(dat)) next

  annotated <- dat |>
    left_join(gene_info, by = c("gene" = "external_gene_name"))

  output <- file.path(
    opt$outdir,
    paste0(tools::file_path_sans_ext(basename(file)), "_annotated.csv")
  )

  write.csv(annotated, output, row.names = FALSE)
}
