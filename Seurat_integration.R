###This script is for the Seurat integration using the cellbender output##
##########################################################################
#########################################################################


getwd()
install.packages("Seurat")
install.packages("cowplot")

library(ggplot2)
library(dplyr)
library(magrittr)
library(patchwork)
library(viridis)
library(Seurat)
library(scCustomize)
library(qs)

cell_bender_mat_1 <- Read_CellBender_h5_Mat(file_name = "/Volumes/Singhal SSD/cellbender_output/1/output_1_filtered.h5")
cell_bender_mat_4<- Read_CellBender_h5_Mat(file_name = "/Volumes/Singhal SSD/cellbender_output/4/output_4_filtered.h5")
cell_bender_mat_A<- Read_CellBender_h5_Mat(file_name = "/Volumes/Singhal SSD/cellbender_output/A/output_A_filtered.h5")
cell_bender_mat_B <- Read_CellBender_h5_Mat(file_name = "/Volumes/Singhal SSD/cellbender_output/B/output_B_filtered.h5")

cell_bender_seurat_1 <- CreateSeuratObject(counts = cell_bender_mat_1, names.field = 1, names.delim = "_")
cell_bender_seurat_1
cell_bender_seurat_4 <- CreateSeuratObject(counts = cell_bender_mat_4, names.field = 1, names.delim = "_")
cell_bender_seurat_A <- CreateSeuratObject(counts = cell_bender_mat_A, names.field = 1, names.delim = "_")
cell_bender_seurat_B <- CreateSeuratObject(counts = cell_bender_mat_B, names.field = 1, names.delim = "_")

# Set up control object Female
ctrlF <- CreateSeuratObject(counts = cell_bender_mat_1, project = "Female_ctrl", min.cells = 5)
ctrlF$group <- "CNTRL"
ctrlF$sex <- "Female"
ctrlF <- subset(ctrlF, subset = nFeature_RNA > 500)
ctrlF <- NormalizeData(ctrlF, verbose = FALSE)
ctrlF <- FindVariableFeatures(ctrlF, selection.method = "vst", nfeatures = 2000)

# Set up stressed object Female
stressF <- CreateSeuratObject(counts = cell_bender_mat_4, project = "Female_stress", min.cells = 5)
stressF$group <- "STRESS"
stressF$sex <- "Female"
stressF <- subset(stressF, subset = nFeature_RNA > 500)
stressF <- NormalizeData(stressF, verbose = FALSE)
stressF <- FindVariableFeatures(stressF, selection.method = "vst", nfeatures = 2000)
stressF[["percent.mt"]] <- PercentageFeatureSet(stressF, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(stressF, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Set up stressed object Male
stressM <- CreateSeuratObject(counts = cell_bender_mat_A, project = "Male_stress", min.cells = 5)
stressM$group <- "STRESS"
stressM$sex <- "Male"
stressM <- subset(stressM, subset = nFeature_RNA > 500)
stressM <- NormalizeData(stressM, verbose = FALSE)
stressM <- FindVariableFeatures(stressM, selection.method = "vst", nfeatures = 2000)

# Set up control object Male
ctrlM <- CreateSeuratObject(counts = cell_bender_mat_B, project = "Male_ctrl", min.cells = 5)
ctrlM$group <- "CNTRL"
ctrlM$sex <- "Male"
ctrlM <- subset(ctrlM, subset = nFeature_RNA > 500)
ctrlM <- NormalizeData(ctrlM, verbose = FALSE)
ctrlM <- FindVariableFeatures(ctrlM, selection.method = "vst", nfeatures = 2000)

##perform intergartion
sex.anchors <- FindIntegrationAnchors(object.list = list(ctrlF, stressF, stressM, ctrlM), dims = 1:20)
sex.combined <- IntegrateData(anchorset = sex.anchors, dims = 1:20, k.weight = 43)

##Perform integration analysis
DefaultAssay(sex.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
sex.combined <- ScaleData(sex.combined, verbose = FALSE)
sex.combined <- RunPCA(sex.combined, npcs = 30, verbose = FALSE)

##dimensionality of the seurat object
ElbowPlot(sex.combined)
DimHeatmap(sex.combined, cells = 300, balanced = TRUE)

# t-SNE and Clustering
sex.combined <- RunUMAP(sex.combined, reduction = "pca", dims = 1:5)
sex.combined <- FindNeighbors(sex.combined, reduction = "pca", dims = 1:5)
sex.combined <- FindClusters(sex.combined, resolution = 1)

sex.combined$sex.group <- paste0(sex.combined$sex, "_", sex.combined$group)

##Dimplot for all the clusters
DimPlot(sex.combined, reduction = "umap", group.by = "sex.group")

# find all markers of cluster 
cluster.markers <- FindAllMarkers(sex.combined, only.pos = T)
write.csv(cluster.markers, "cluster.markers_CB.csv")

#Cluster naming
sex.combined <- RenameIdents(sex.combined, `0` = "Oligodendrocytes", `1` = "Neurons", `2` = "Oligos(Schwann Cells)", 
                             '3' = "Astrocytes", '4' = "Microglia", '5' = "Endothelial Cells")
DimPlot(sex.combined, label = TRUE)

p1 <- DimPlot(sex.combined, reduction = "umap", split.by = "sex"). ## To get the cluster distribution in males and females
plot(p1)
p2 <- DimPlot(sex.combined, reduction = "umap", split.by = "group"). ## To get the cluster distribution in controls vs stress
plot(p2)


##Identify the conserved cell type
library(BiocManager)
library(remotes)
library(devtools)

DefaultAssay(sex.combined) <- "RNA"
sex.combined <- JoinLayers(sex.combined)
sex.combined
#nk.markers <- FindConservedMarkers(sex.combined, ident.1 = "Microglia", grouping.var = "group", verbose = FALSE)
#head(nk.markers)

### Find all DEGs in the clusters - Tis will find makers differentially expressed in each identity group by comparing it to all other####
sex.combined.markers <- FindAllMarkers(sex.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay='RNA')
top_markers_cluster <-sex.combined.markers %>% group_by(cluster) %>% top_n(n = 20)
write.csv(top_markers_cluster, "top_markers _in_each_cluster.csv")

