# Nasal tau vaccine elicits a distinct astrocyte population characterized by an immune-quiescent signature

# This script has been confirmed to run with Seurat version 5.1.0 and R version 4.4.0 

library(Seurat)  
library(BPCells)
library(presto)
library(glmGamPoi) 
library(tidyverse)
library(cowplot)
library(patchwork)
library(multtest)
library(metap)
library(Nebulosa)
library(ggfortify)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggupset)
library(ComplexHeatmap)
library(dendextend)


##### Preparation for integrated data generation  #####

sample <- c("nonTg_control", "nonTg_vac", "Tg_control", "Tg_vac")
condition <- c("nonTg_c", "nonTg_v", "Tg_c", "Tg_v")

for (i in seq_along(sample)) {
  
  dir_name <- paste0("./data/", sample[i], "_filtered_feature_bc_matrix/")
  project_name <- sample[i]
  
  Seurat_object <- Read10X(data.dir = dir_name)
  print(dim(Seurat_object))
  Seurat_object <- CreateSeuratObject(counts = Seurat_object, min.cells = 3, min.features = 200, project = project_name)
  print(dim(Seurat_object))
  
  Seurat_object$percent.mito <- PercentageFeatureSet(Seurat_object, pattern = '^mt-')
  Seurat_object$percent.ribo <- PercentageFeatureSet(Seurat_object, pattern = '^Rp[sl]')
  Seurat_object$condition <- condition[i]
  print(head(Seurat_object@meta.data))
  
  Seurat_object <- subset(Seurat_object, 
                          subset = nFeature_RNA > 1300 & nFeature_RNA < 8000 & nCount_RNA > 2500 & nCount_RNA < 30000 & percent.mito < 15 & percent.ribo < 20)
  print(dim(GetAssayData(Seurat_object)))

  assign(paste0(sample[i]), Seurat_object)
}

data_list <- list(nonTg_control, nonTg_vac, Tg_control, Tg_vac)

#####　Generation of integrated data  #####

data_merge <- Reduce(function(x, y) merge(x, y), data_list)

data_merge <- NormalizeData(data_merge)
data_merge <- FindVariableFeatures(data_merge)
data_merge <- ScaleData(data_merge)
data_merge <- RunPCA(data_merge, npcs = 30)

integ <- IntegrateLayers(data_merge, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
integ[["RNA"]] <- JoinLayers(integ[["RNA"]])
table(integ[['orig.ident']])
integ@assays

integ <- FindNeighbors(integ, dims = 1:30, reduction = "integrated.cca")
integ <- FindClusters(integ, resolution = 0.8, cluster.name = "integrated_clusters")
integ <- RunUMAP(integ, dims = 1:30, reduction = "integrated.cca", reduction.name = "umap")

VlnPlot(integ, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"), group.by = "seurat_clusters", ncol = 2, pt.size = 0.1)

table <- table(integ@meta.data$seurat_clusters, integ@meta.data$orig.ident)
table

DimPlot(integ, reduction = "umap", group.by = c("condition", "seurat_clusters"))
DimPlot(integ, reduction = "umap", group.by = "condition")
DimPlot(integ, reduction = "umap", label = TRUE)
DimPlot(integ, reduction = "umap", split.by = "condition", label = FALSE)

gene_names <- rownames(integ)

#####　Identification of cluster marker genes  #####

for (i in 0:23) {
  
  cluster_name <- paste0("cluster_", i)
  file_name <- paste0("./output/integ_markers_pca30_res0-8_", cluster_name, ".csv")
  
  cluster_name <- FindConservedMarkers(integ, ident.1 = as.character(i), grouping.var = "orig.ident", verbose = TRUE,
                                       min.pct = 0.1)
  write.csv(cluster_name, file_name)
}

#####　Identification of cell type  #####

FeaturePlot(integ, features = "Cx3cr1", label = FALSE, reduction = "umap")
VlnPlot(integ, features = "Cx3cr1", pt.size = 0) + 
  NoLegend() +
  theme(axis.title.x = element_blank()) +
  geom_boxplot(width = 0.3, outlier.size = 0, position = position_dodge(width = 0.9))

#####　Add cell type annotation to the UMAP  #####

# cluster13 is "undefined"
cluster.ids <- c("Microglia", "VE_cells", "Microglia", "Microglia", "Astrocytes", "CE_cells",
                 "Microglia", "VE_cells", "CE_cells", "CE_cells", "Pericytes", "VE_cells",
                 "Microglia", "undefined", "Ependymal_cells", "Astrocytes", "mature_Neurons", "immature_Neurons",
                 "Macrophages", "VL_cells", "Microglia", "VSM_cells", "VE_cells", "immature_Neurons")
names(cluster.ids) <- 0:23

integ <- RenameIdents(integ, cluster.ids)
integ$celltype <- Idents(integ)

#####　cell type extraction  #####

celltype.ids
celltype_list <- list()

for (i in seq_along(celltype.ids)) {
  
  Seurat_object <- subset(integ, idents = celltype.ids[i])
  assign(paste0(celltype.ids[i]), Seurat_object)
  
  celltype_list[[paste0(celltype.ids[i])]] <- Seurat_object
}

#####　Microglial subclustering  #####

Microglia <- FindNeighbors(Microglia, dims = 1:30, reduction = "integrated.cca")
Microglia <- FindClusters(Microglia, resolution = 0.2, cluster.name = "integrated_clusters")
Microglia <- RunUMAP(Microglia, dims = 1:30, reduction = "integrated.cca", reduction.name = "umap")

table <- table(Microglia@meta.data$seurat_clusters, Microglia@meta.data$orig.ident)
table

#####　Identification of the microglial subcluster marker genes  #####

for (i in 0:4) {
  
  cluster_name <- paste0("cluster_", i)
  file_name <- paste0("./output/Microglia_markers_pca30_res0-2_", cluster_name, ".csv") # "_posi.csv"
  
  cluster_name <- FindConservedMarkers(Microglia, ident.1 = as.character(i), grouping.var = "orig.ident", min.pct = 0.1, verbose = TRUE)
  head(cluster_name)
  write.csv(cluster_name, file_name)
}

#####　Astrocyte subclustering  #####

Astrocytes <- FindNeighbors(Astrocytes, dims = 1:30, reduction = "integrated.cca")
Astrocytes <- FindClusters(Astrocytes, resolution = 0.9, cluster.name = "integrated_clusters")
Astrocytes <- RunUMAP(Astrocytes, dims = 1:30, reduction = "integrated.cca", reduction.name = "umap")

table <- table(Astrocytes@meta.data$seurat_clusters, Astrocytes@meta.data$orig.ident)
table

#####　Identification of the astrocyte subcluster marker genes  ######

for (i in 0:8) {
  
  cluster_name <- paste0("cluster_", i)
  file_name <- paste0("./output/Astrocytes_markers_pca30_res0-9_", cluster_name, ".csv")
  
  cluster_name <- FindConservedMarkers(Astrocytes, ident.1 = as.character(i), grouping.var = "orig.ident", min.pct = 0.1, verbose = TRUE)
  head(cluster_name)
  write.csv(cluster_name, file_name)
}

# end
