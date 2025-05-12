## NOTES ####
#Using Seurat V5.0.1: https://satijalab.org/seurat/articles/install_v5.html
#using R 4.3.0

## Libraries ####
library(svglite)
library(tidyverse)
library(data.table)
library(readxl)
library(openxlsx)
library(dplyr)
library(glmGamPoi)
library(Seurat)
#library(Signac)
library(SeuratData)
#library(Azimuth)
library(SeuratWrappers)
library(patchwork)
library(DESeq2)
library(car)
library(sm)
library(preprocessCore)
library(reticulate)
library(anndata)
library(sceasy)
library(cowplot)
library(harmony)
library(renv)
library(BUSpaRse) #knee plot
library("loupeR")
library(scCustomize)
library(DoubletFinder) #can consider running after clustering: https://github.com/chris-mcginnis-ucsf/DoubletFinder
#likely not many doublets looking @ the knee plots


###### PART 0 - DEFINE PLOTTING FXNS ####

qc_violins <- function( seuratobj, name ) {
  
  #with dots
  qcplotname <- paste0("QCViolin_",name,".svg")
  violin <- VlnPlot(seuratobj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0.3)
  svg(qcplotname, width = 16, height = 8)
  print(violin)
  dev.off()
  
  #w/out dots
  qcplotname <- paste0("QCViolinNODOTS_", name,".svg")
  violin <- VlnPlot(seuratobj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0)
  svg(qcplotname, width = 16, height = 8)
  print(violin)
  dev.off()
}

qc_feat <- function( seuratobj, reduc, name ) {
  
  prev_assay <- DefaultAssay(seuratobj)
  
  #RNA Assay
  DefaultAssay(seuratobj) <- "RNA"
  featplotname <- paste0("QCFeatRNA_",name,".svg")
  plotheight <- ceiling(3/2)*10
  svg(featplotname, width = 20, height = plotheight)
  print(FeaturePlot(seuratobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    reduction = reduc, pt.size = 0.2, ncol = 2,
                    label = T, raster = F))
  dev.off()
  
  DefaultAssay(seuratobj) <- prev_assay
}

umaps <- function( seuratobj, reduc, name ) {
  
  #UMAP plot w/ active ident
  umapplotname <- paste0("UMAP_", name,".svg")
  umap_plot <- DimPlot(
    seuratobj, label = T,
    reduction = reduc, raster = F,
    combine = FALSE, label.size = 8)
  svg(umapplotname, width = 10, height = 10)
  print(umap_plot)
  dev.off()
  
  #UMAP plot w/ sample
  umapplotname <- paste0("UMAPsamp_", name,".svg")
  umap_plot <- DimPlot(
    seuratobj, label = F, group.by = "orig.ident",
    reduction = reduc, raster = F,
    combine = FALSE)
  svg(umapplotname, width = 10, height = 10)
  print(umap_plot)
  dev.off()
  
  #UMAP plot w/ active ident but split by sample
  umapplotname <- paste0("UMAPsplit_", name,".svg")
  umap_plot <- DimPlot(
    seuratobj, label = T,
    reduction = reduc, raster = F,
    split.by = "orig.ident",
    combine = FALSE, label.size = 2)
  svg(umapplotname, width = 25, height = 5)
  print(umap_plot)
  dev.off()
}

feat_plots <- function( seuratobj, reduc, feats, name ) {
  
  prev_assay <- DefaultAssay(seuratobj)
  
  #RNA Assay
  DefaultAssay(seuratobj) <- "RNA"
  featplotname <- paste0("FeatRNA_",name,".svg")
  plotheight <- ceiling(length(feats)/2)*10
  svg(featplotname, width = 20, height = plotheight)
  print(FeaturePlot(seuratobj, features = feats,
                    reduction = reduc, pt.size = 0.2, ncol = 2, 
                    label = F, raster = F))
  dev.off()
  
  DefaultAssay(seuratobj) <- prev_assay
}

dot_plots <- function( seuratobj, feats, name ) {
  
  prev_assay <- DefaultAssay(seuratobj)
  
  #RNA Assay
  dotplotname <- paste0("DotRNA_",name,".svg")
  plotheight <- 2+length(unique(Idents(seuratobj)))/3
  plotwidth <- 3+length(feats)/3
  svg(dotplotname, width = plotwidth, height = plotheight)
  print(DotPlot(seuratobj, features = feats, assay = "RNA", scale = F) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "italic")))
  dev.off()
  
  #RNA Assay w/scale
  dotplotname <- paste0("DotRNAscale_",name,".svg")
  plotheight <- 2+length(unique(Idents(seuratobj)))/3
  plotwidth <- 3+length(feats)/3
  svg(dotplotname, width = plotwidth, height = plotheight)
  print(DotPlot(seuratobj, features = feats, assay = "RNA", scale = T) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "italic")))
  dev.off()
  
  DefaultAssay(seuratobj) <- prev_assay
}

vln_plots <- function( seuratobj, feats, name ) {
  
  prev_assay <- DefaultAssay(seuratobj)
  
  #RNA Assay
  vlnplotname <- paste0("VlnRNA_",name,".svg")
  plotheight <- ceiling(length(feats)/4)*3
  plotwidth <- 4+length(unique(Idents(seuratobj)))
  svg(vlnplotname, width = plotwidth, height = plotheight)
  print(VlnPlot(seuratobj, features = feats, assay = "RNA", pt.size = 0.5, ncol=4))
  dev.off()
  
  DefaultAssay(seuratobj) <- prev_assay
  
}



###### PART 1 - LOAD INPUT DATASET ####

sex.combined <- readRDS("sex.combined.rds")

#Add metadata
sex.combined$sex_cell <- paste0(sex.combined$sex,"_",sex.combined$seurat_clusters)
sex.combined$group_cell <- paste0(sex.combined$group,"_",sex.combined$seurat_clusters)
Idents(sex.combined) <- "seurat_clusters"

#Working Dir
mainwd <- getwd()
if(!dir.exists("whole_dataset")){dir.create("whole_dataset")}
setwd("whole_dataset")

#Standard Plots
markers.gen <- c("Snap25", "Gfap", "Mbp", "Mog", "Tek", "Cd163")
qc_feat(sex.combined,"umap","allClusters")
qc_violins(sex.combined, "allClusters")
umaps(sex.combined,"umap","allClusters")
dot_plots(sex.combined, markers.gen,"allClusters")
vln_plots(sex.combined, markers.gen,"allClusters")
feat_plots(sex.combined,"umap", markers.gen,"allClusters")

#Make adj seurat object for DESeq2
#NOTE: cannot have a bunch of 0s or DESeq2 will fail; will add 1 count to all;
#this is a common practice for pseudobulking for DESeq2
#https://github.com/satijalab/seurat/issues/1570
adj_sex.combined <- sex.combined
adj_sex.combined[["RNA"]]$counts<-as.matrix(adj_sex.combined[["RNA"]]$counts)+1

###### PART 2 - DE PER CLUSTER ####

#Working Dir
setwd(mainwd)
if(!dir.exists("PerCluster")){dir.create("PerCluster")}
setwd("PerCluster")

#DESeq#
Idents(adj_sex.combined) <- "seurat_clusters"
sc_clust.markersDESeq <- FindAllMarkers(adj_sex.combined, test.use = "DESeq2", 
                                   return.thresh = 0.05, only.pos = T) 
#DEFAULTS: logfc.threshold=0.1; min.pct=0.1;
sc_clust.markersDESeq$negLOGp <- -log10(sc_clust.markersDESeq$p_val)
top_5 <- Extract_Top_Markers(marker_dataframe = sc_clust.markersDESeq, num_genes = 5, rank_by = "negLOGp")
head(top_5, 10)
dot_plots(adj_sex.combined,top_5,"clus_DESeq")
vln_plots(adj_sex.combined,top_5,"clus_DESeq")
feat_plots(adj_sex.combined,"umap",top_5,"clus_DESeq")
write.xlsx(sc_clust.markersDESeq, "DESeq2.xlsx")

#................................................................................#

#Wilcox#
Idents(sex.combined) <- "seurat_clusters"
sc_clust.markersWilcox <- FindAllMarkers(sex.combined, test.use = "wilcox", 
                                        return.thresh = 0.05, only.pos = T) 
#DEFAULTS: logfc.threshold=0.1; min.pct=0.1; 
sc_clust.markersWilcox$negLOGp <- -log10(sc_clust.markersWilcox$p_val)
top_5 <- Extract_Top_Markers(marker_dataframe = sc_clust.markersWilcox, num_genes = 5, rank_by = "negLOGp")
head(top_5, 10)
dot_plots(sex.combined,top_5,"clus_Wilcox")
vln_plots(sex.combined,top_5,"clus_Wilcox")
feat_plots(sex.combined,"umap",top_5,"clus_Wilcox")
write.xlsx(sc_clust.markersWilcox, "wilcox.xlsx")



###### PART 4 - DE BY SEX FOR GROUP 0 ####

#Working Dir
setwd(mainwd)
if(!dir.exists("BySexFor0")){dir.create("BySexFor0")}
setwd("BySexFor0")

#Wilcox#
Idents(sex.combined) <- "sex_cell"
sc_clust.markersWilcox <- FindMarkers(sex.combined, test.use = "wilcox", 
                                         ident.1 = "Female_0", ident.2 = "Male_0",
                                         return.thresh = 0.05, only.pos = F) 
#DEFAULTS: logfc.threshold=0.1; min.pct=0.1;
sc_clust.markersWilcox$cluster <- ifelse(sc_clust.markersWilcox$avg_log2FC>0,"Female","Male")
sc_clust.markersWilcox$gene <- rownames(sc_clust.markersWilcox)
sc_clust.markersWilcox$negLOGp <- -log10(sc_clust.markersWilcox$p_val)
top_5 <- Extract_Top_Markers(marker_dataframe = sc_clust.markersWilcox, num_genes = 5, rank_by = "negLOGp")
head(top_5, 10)
dot_plots(sex.combined,top_5,"clus_Wilcox")
vln_plots(sex.combined,top_5,"clus_Wilcox")
feat_plots(sex.combined,"umap",top_5,"clus_Wilcox")
write.xlsx(sc_clust.markersWilcox, "wilcox.xlsx")



###### PART 5 - DE BY GROUP FOR GROUP 0 ####

#Working Dir
setwd(mainwd)
if(!dir.exists("ByGroupFor0")){dir.create("ByGroupFor0")}
setwd("ByGroupFor0")

#Wilcox#
Idents(sex.combined) <- "group_cell"
sc_clust.markersWilcox <- FindMarkers(sex.combined, test.use = "wilcox", 
                                      ident.1 = "CNTRL_0", ident.2 = "STRESS_0",
                                      return.thresh = 0.05, only.pos = F) 
#DEFAULTS: logfc.threshold=0.1; min.pct=0.1;
sc_clust.markersWilcox$cluster <- ifelse(sc_clust.markersWilcox$avg_log2FC>0,"CNTRL","STRESS")
sc_clust.markersWilcox$gene <- rownames(sc_clust.markersWilcox)
sc_clust.markersWilcox$negLOGp <- -log10(sc_clust.markersWilcox$p_val)
top_5 <- Extract_Top_Markers(marker_dataframe = sc_clust.markersWilcox, num_genes = 5, rank_by = "negLOGp")
head(top_5, 10)
dot_plots(sex.combined,top_5,"clus_Wilcox")
vln_plots(sex.combined,top_5,"clus_Wilcox")
feat_plots(sex.combined,"umap",top_5,"clus_Wilcox")
write.xlsx(sc_clust.markersWilcox, "wilcox.xlsx")


###### PART 6 - DE BY SEX FOR GROUP 1 ####

#Working Dir
setwd(mainwd)
if(!dir.exists("BySexFor1")){dir.create("BySexFor1")}
setwd("BySexFor1")

#Wilcox#
Idents(sex.combined) <- "sex_cell"
sc_clust.markersWilcox <- FindMarkers(sex.combined, test.use = "wilcox", 
                                      ident.1 = "Female_1", ident.2 = "Male_1",
                                      return.thresh = 0.05, only.pos = F) 
#DEFAULTS: logfc.threshold=0.1; min.pct=0.1;
sc_clust.markersWilcox$cluster <- ifelse(sc_clust.markersWilcox$avg_log2FC>0,"Female","Male")
sc_clust.markersWilcox$gene <- rownames(sc_clust.markersWilcox)
sc_clust.markersWilcox$negLOGp <- -log10(sc_clust.markersWilcox$p_val)
top_5 <- Extract_Top_Markers(marker_dataframe = sc_clust.markersWilcox, num_genes = 5, rank_by = "negLOGp")
head(top_5, 10)
dot_plots(sex.combined,top_5,"clus_Wilcox")
vln_plots(sex.combined,top_5,"clus_Wilcox")
feat_plots(sex.combined,"umap",top_5,"clus_Wilcox")
write.xlsx(sc_clust.markersWilcox, "wilcox.xlsx")



###### PART 7 - DE BY GROUP FOR GROUP 1 ####

#Working Dir
setwd(mainwd)
if(!dir.exists("ByGroupFor1")){dir.create("ByGroupFor1")}
setwd("ByGroupFor1")

#Wilcox#
Idents(sex.combined) <- "group_cell"
sc_clust.markersWilcox <- FindMarkers(sex.combined, test.use = "wilcox", 
                                      ident.1 = "CNTRL_1", ident.2 = "STRESS_1",
                                      return.thresh = 0.05, only.pos = F) 
#DEFAULTS: logfc.threshold=0.1; min.pct=0.1;
sc_clust.markersWilcox$cluster <- ifelse(sc_clust.markersWilcox$avg_log2FC>0,"CNTRL","STRESS")
sc_clust.markersWilcox$gene <- rownames(sc_clust.markersWilcox)
sc_clust.markersWilcox$negLOGp <- -log10(sc_clust.markersWilcox$p_val)
top_5 <- Extract_Top_Markers(marker_dataframe = sc_clust.markersWilcox, num_genes = 5, rank_by = "negLOGp")
head(top_5, 10)
dot_plots(sex.combined,top_5,"clus_Wilcox")
vln_plots(sex.combined,top_5,"clus_Wilcox")
feat_plots(sex.combined,"umap",top_5,"clus_Wilcox")
write.xlsx(sc_clust.markersWilcox, "wilcox.xlsx")


###### PART 8 - DE BY SEX FOR GROUP 2 ####

#Working Dir
setwd(mainwd)
if(!dir.exists("BySexFor2")){dir.create("BySexFor2")}
setwd("BySexFor2")

#Wilcox#
Idents(sex.combined) <- "sex_cell"
sc_clust.markersWilcox <- FindMarkers(sex.combined, test.use = "wilcox", 
                                      ident.1 = "Female_2", ident.2 = "Male_2",
                                      return.thresh = 0.05, only.pos = F) 
#DEFAULTS: logfc.threshold=0.1; min.pct=0.1;
sc_clust.markersWilcox$cluster <- ifelse(sc_clust.markersWilcox$avg_log2FC>0,"Female","Male")
sc_clust.markersWilcox$gene <- rownames(sc_clust.markersWilcox)
sc_clust.markersWilcox$negLOGp <- -log10(sc_clust.markersWilcox$p_val)
top_5 <- Extract_Top_Markers(marker_dataframe = sc_clust.markersWilcox, num_genes = 5, rank_by = "negLOGp")
head(top_5, 10)
dot_plots(sex.combined,top_5,"clus_Wilcox")
vln_plots(sex.combined,top_5,"clus_Wilcox")
feat_plots(sex.combined,"umap",top_5,"clus_Wilcox")
write.xlsx(sc_clust.markersWilcox, "wilcox.xlsx")



###### PART 9 - DE BY GROUP FOR GROUP 0 ####

#Working Dir
setwd(mainwd)
if(!dir.exists("ByGroupFor2")){dir.create("ByGroupFor2")}
setwd("ByGroupFor2")

#Wilcox#
Idents(sex.combined) <- "group_cell"
sc_clust.markersWilcox <- FindMarkers(sex.combined, test.use = "wilcox", 
                                      ident.1 = "CNTRL_2", ident.2 = "STRESS_2",
                                      return.thresh = 0.05, only.pos = F) 
#DEFAULTS: logfc.threshold=0.1; min.pct=0.1;
sc_clust.markersWilcox$cluster <- ifelse(sc_clust.markersWilcox$avg_log2FC>0,"CNTRL","STRESS")
sc_clust.markersWilcox$gene <- rownames(sc_clust.markersWilcox)
sc_clust.markersWilcox$negLOGp <- -log10(sc_clust.markersWilcox$p_val)
top_5 <- Extract_Top_Markers(marker_dataframe = sc_clust.markersWilcox, num_genes = 5, rank_by = "negLOGp")
head(top_5, 10)
dot_plots(sex.combined,top_5,"clus_Wilcox")
vln_plots(sex.combined,top_5,"clus_Wilcox")
feat_plots(sex.combined,"umap",top_5,"clus_Wilcox")
write.xlsx(sc_clust.markersWilcox, "wilcox.xlsx")


###### PART 10 - DE BY SEX FOR GROUP 3 ####

#Working Dir
setwd(mainwd)
if(!dir.exists("BySexFor3")){dir.create("BySexFor3")}
setwd("BySexFor3")

#Wilcox#
Idents(sex.combined) <- "sex_cell"
sc_clust.markersWilcox <- FindMarkers(sex.combined, test.use = "wilcox", 
                                      ident.1 = "Female_3", ident.2 = "Male_3",
                                      return.thresh = 0.05, only.pos = F) 
#DEFAULTS: logfc.threshold=0.1; min.pct=0.1;
sc_clust.markersWilcox$cluster <- ifelse(sc_clust.markersWilcox$avg_log2FC>0,"Female","Male")
sc_clust.markersWilcox$gene <- rownames(sc_clust.markersWilcox)
sc_clust.markersWilcox$negLOGp <- -log10(sc_clust.markersWilcox$p_val)
top_5 <- Extract_Top_Markers(marker_dataframe = sc_clust.markersWilcox, num_genes = 5, rank_by = "negLOGp")
head(top_5, 10)
dot_plots(sex.combined,top_5,"clus_Wilcox")
vln_plots(sex.combined,top_5,"clus_Wilcox")
feat_plots(sex.combined,"umap",top_5,"clus_Wilcox")
write.xlsx(sc_clust.markersWilcox, "wilcox.xlsx")



###### PART 11 - DE BY GROUP FOR GROUP 3 ####

#Working Dir
setwd(mainwd)
if(!dir.exists("ByGroupFor3")){dir.create("ByGroupFor3")}
setwd("ByGroupFor3")

#Wilcox#
Idents(sex.combined) <- "group_cell"
sc_clust.markersWilcox <- FindMarkers(sex.combined, test.use = "wilcox", 
                                      ident.1 = "CNTRL_3", ident.2 = "STRESS_3",
                                      return.thresh = 0.05, only.pos = F) 
#DEFAULTS: logfc.threshold=0.1; min.pct=0.1;
sc_clust.markersWilcox$cluster <- ifelse(sc_clust.markersWilcox$avg_log2FC>0,"CNTRL","STRESS")
sc_clust.markersWilcox$gene <- rownames(sc_clust.markersWilcox)
sc_clust.markersWilcox$negLOGp <- -log10(sc_clust.markersWilcox$p_val)
top_5 <- Extract_Top_Markers(marker_dataframe = sc_clust.markersWilcox, num_genes = 5, rank_by = "negLOGp")
head(top_5, 10)
dot_plots(sex.combined,top_5,"clus_Wilcox")
vln_plots(sex.combined,top_5,"clus_Wilcox")
feat_plots(sex.combined,"umap",top_5,"clus_Wilcox")
write.xlsx(sc_clust.markersWilcox, "wilcox.xlsx")


###### PART 12 - DE BY SEX FOR GROUP 4 ####

#Working Dir
setwd(mainwd)
if(!dir.exists("BySexFor4")){dir.create("BySexFor4")}
setwd("BySexFor4")

#Wilcox#
Idents(sex.combined) <- "sex_cell"
sc_clust.markersWilcox <- FindMarkers(sex.combined, test.use = "wilcox", 
                                      ident.1 = "Female_4", ident.2 = "Male_4",
                                      return.thresh = 0.05, only.pos = F) 
#DEFAULTS: logfc.threshold=0.1; min.pct=0.1;
sc_clust.markersWilcox$cluster <- ifelse(sc_clust.markersWilcox$avg_log2FC>0,"Female","Male")
sc_clust.markersWilcox$gene <- rownames(sc_clust.markersWilcox)
sc_clust.markersWilcox$negLOGp <- -log10(sc_clust.markersWilcox$p_val)
top_5 <- Extract_Top_Markers(marker_dataframe = sc_clust.markersWilcox, num_genes = 5, rank_by = "negLOGp")
head(top_5, 10)
dot_plots(sex.combined,top_5,"clus_Wilcox")
vln_plots(sex.combined,top_5,"clus_Wilcox")
feat_plots(sex.combined,"umap",top_5,"clus_Wilcox")
write.xlsx(sc_clust.markersWilcox, "wilcox.xlsx")



###### PART 13 - DE BY GROUP FOR GROUP 4 ####

#Working Dir
setwd(mainwd)
if(!dir.exists("ByGroupFor4")){dir.create("ByGroupFor4")}
setwd("ByGroupFor4")

#Wilcox#
Idents(sex.combined) <- "group_cell"
sc_clust.markersWilcox <- FindMarkers(sex.combined, test.use = "wilcox", 
                                      ident.1 = "CNTRL_4", ident.2 = "STRESS_4",
                                      return.thresh = 0.05, only.pos = F) 
#DEFAULTS: logfc.threshold=0.1; min.pct=0.1;
sc_clust.markersWilcox$cluster <- ifelse(sc_clust.markersWilcox$avg_log2FC>0,"CNTRL","STRESS")
sc_clust.markersWilcox$gene <- rownames(sc_clust.markersWilcox)
sc_clust.markersWilcox$negLOGp <- -log10(sc_clust.markersWilcox$p_val)
top_5 <- Extract_Top_Markers(marker_dataframe = sc_clust.markersWilcox, num_genes = 5, rank_by = "negLOGp")
head(top_5, 10)
dot_plots(sex.combined,top_5,"clus_Wilcox")
vln_plots(sex.combined,top_5,"clus_Wilcox")
feat_plots(sex.combined,"umap",top_5,"clus_Wilcox")
write.xlsx(sc_clust.markersWilcox, "wilcox.xlsx")


###### PART 14 - DE BY SEX FOR GROUP 5 ####

#Working Dir
setwd(mainwd)
if(!dir.exists("BySexFor5")){dir.create("BySexFor5")}
setwd("BySexFor5")

#Wilcox#
Idents(sex.combined) <- "sex_cell"
sc_clust.markersWilcox <- FindMarkers(sex.combined, test.use = "wilcox", 
                                      ident.1 = "Female_5", ident.2 = "Male_5",
                                      return.thresh = 0.05, only.pos = F) 
#DEFAULTS: logfc.threshold=0.1; min.pct=0.1;
sc_clust.markersWilcox$cluster <- ifelse(sc_clust.markersWilcox$avg_log2FC>0,"Female","Male")
sc_clust.markersWilcox$gene <- rownames(sc_clust.markersWilcox)
sc_clust.markersWilcox$negLOGp <- -log10(sc_clust.markersWilcox$p_val)
top_5 <- Extract_Top_Markers(marker_dataframe = sc_clust.markersWilcox, num_genes = 5, rank_by = "negLOGp")
head(top_5, 10)
dot_plots(sex.combined,top_5,"clus_Wilcox")
vln_plots(sex.combined,top_5,"clus_Wilcox")
feat_plots(sex.combined,"umap",top_5,"clus_Wilcox")
write.xlsx(sc_clust.markersWilcox, "wilcox.xlsx")



###### PART 13 - DE BY GROUP FOR GROUP 5 ####

#Working Dir
setwd(mainwd)
if(!dir.exists("ByGroupFor5")){dir.create("ByGroupFor5")}
setwd("ByGroupFor5")

#Wilcox#
Idents(sex.combined) <- "group_cell"
sc_clust.markersWilcox <- FindMarkers(sex.combined, test.use = "wilcox", 
                                      ident.1 = "CNTRL_5", ident.2 = "STRESS_5",
                                      return.thresh = 0.05, only.pos = F) 
#DEFAULTS: logfc.threshold=0.1; min.pct=0.1;
sc_clust.markersWilcox$cluster <- ifelse(sc_clust.markersWilcox$avg_log2FC>0,"CNTRL","STRESS")
sc_clust.markersWilcox$gene <- rownames(sc_clust.markersWilcox)
sc_clust.markersWilcox$negLOGp <- -log10(sc_clust.markersWilcox$p_val)
top_5 <- Extract_Top_Markers(marker_dataframe = sc_clust.markersWilcox, num_genes = 5, rank_by = "negLOGp")
head(top_5, 10)
dot_plots(sex.combined,top_5,"clus_Wilcox")
vln_plots(sex.combined,top_5,"clus_Wilcox")
feat_plots(sex.combined,"umap",top_5,"clus_Wilcox")
write.xlsx(sc_clust.markersWilcox, "wilcox.xlsx")


