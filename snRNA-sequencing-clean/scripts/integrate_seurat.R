#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(scCustomize)
  library(ggplot2)
})

option_list <- list(
  make_option("--ctrl-f", type = "character"),
  make_option("--stress-f", type = "character"),
  make_option("--stress-m", type = "character"),
  make_option("--ctrl-m", type = "character"),
  make_option("--out-rds", type = "character"),
  make_option("--outdir", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

read_cellbender <- function(path, project, sex, group) {
  counts <- Read_CellBender_h5_Mat(file_name = path)

  obj <- CreateSeuratObject(
    counts = counts,
    project = project,
    min.cells = 5
  )

  obj$sex <- sex
  obj$group <- group

  obj <- subset(obj, subset = nFeature_RNA > 500)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)

  obj
}

ctrl_f   <- read_cellbender(opt$ctrl_f,   "Female_ctrl",  "Female", "CNTRL")
stress_f <- read_cellbender(opt$stress_f, "Female_stress", "Female", "STRESS")
stress_m <- read_cellbender(opt$stress_m, "Male_stress",   "Male",   "STRESS")
ctrl_m   <- read_cellbender(opt$ctrl_m,   "Male_ctrl",    "Male",   "CNTRL")

anchors <- FindIntegrationAnchors(
  object.list = list(ctrl_f, stress_f, stress_m, ctrl_m),
  dims = 1:20
)

combined <- IntegrateData(
  anchorset = anchors,
  dims = 1:20,
  k.weight = 43
)

DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:5)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:5)
combined <- FindClusters(combined, resolution = 1)

combined$sex_group <- paste(combined$sex, combined$group, sep = "_")

pdf(file.path(opt$outdir, "UMAP_sex_group.pdf"), width = 8, height = 6)
print(DimPlot(combined, reduction = "umap", group.by = "sex_group"))
dev.off()

DefaultAssay(combined) <- "RNA"
combined <- JoinLayers(combined)

markers <- FindAllMarkers(
  combined,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  assay = "RNA"
)

write.csv(
  markers,
  file.path(opt$outdir, "cluster_markers.csv"),
  row.names = FALSE
)

saveRDS(combined, opt$out_rds)
