#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
})

option_list <- list(
  make_option("--input", type = "character"),
  make_option("--outdir", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(opt$input)

obj$sex_cell <- paste(obj$sex, obj$seurat_clusters, sep = "_")
obj$group_cell <- paste(obj$group, obj$seurat_clusters, sep = "_")

run_de <- function(object, ident_column, ident1, ident2, label) {
  Idents(object) <- ident_column

  res <- FindMarkers(
    object,
    ident.1 = ident1,
    ident.2 = ident2,
    test.use = "wilcox",
    return.thresh = 0.05,
    only.pos = FALSE
  )

  res$gene <- rownames(res)
  res$comparison <- label
  res$direction <- ifelse(res$avg_log2FC > 0, ident1, ident2)
  res
}

clusters <- sort(unique(obj$seurat_clusters))

sex_results <- do.call(
  rbind,
  lapply(clusters, function(cluster_id) {
    run_de(
      obj, "sex_cell",
      paste("Female", cluster_id, sep = "_"),
      paste("Male", cluster_id, sep = "_"),
      paste0("Female_vs_Male_cluster_", cluster_id)
    )
  })
)

group_results <- do.call(
  rbind,
  lapply(clusters, function(cluster_id) {
    run_de(
      obj, "group_cell",
      paste("CNTRL", cluster_id, sep = "_"),
      paste("STRESS", cluster_id, sep = "_"),
      paste0("Control_vs_Stress_cluster_", cluster_id)
    )
  })
)

write.csv(sex_results, file.path(opt$outdir, "sex_comparisons_by_cluster.csv"), row.names = FALSE)
write.csv(group_results, file.path(opt$outdir, "stress_comparisons_by_cluster.csv"), row.names = FALSE)
