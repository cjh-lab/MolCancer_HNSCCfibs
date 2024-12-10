#### PCFA construction ####
library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)

# loading in datasets, ensuring Dataset set in metadata
# making list of datasets, finding common gene space, sketch based integration
Dataset_list <- c(PDAC, Breast, HNSCC, Colon, Lung, esophageal, Gastric)

# finding common gene space between datasets and filtering to only include common genes
gene_sets <- lapply(Dataset_list, function(x) rownames(x@assays$RNA))
common_genes <- Reduce(intersect, gene_sets)

Dataset_list <- merge(Dataset_list[[1]], Dataset_list[-1])
PCFA <- subset(Dataset_list, features = common_genes)

# Ensure v5 assay for integration using sketching
PCFA[["RNA"]] <- as(object = PCFA[["RNA"]], Class = "Assay5")
PCFA <- UpdateSeuratObject(PCFA)

# split assay into n layers (n = number of datasets)
PCFA <- NormalizeData(PCFA)
PCFA[["RNA"]] <- split(PCFA[["RNA"]], f = PCFA$Dataset)
PCFA <- FindVariableFeatures(PCFA, verbose = FALSE)

# 5000 mesenchymal cells from each dataset sketched
PCFA <- SketchData(object = PCFA, ncells = 5000, method = "LeverageScore", sketched.assay = "sketch", seed = 222, over.write = T)

DefaultAssay(PCFA) <- "sketch"
PCFA <- FindVariableFeatures(PCFA, verbose = F, nfeatures = 3000)
PCFA <- ScaleData(PCFA, verbose = F)
PCFA <- RunPCA(PCFA, verbose = F)

# integrate the datasets using harmony
PCFA <- IntegrateLayers(
  object = PCFA, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE, group.by = c("Dataset", "Tech"), theta = c(1)
)

# cluster the integrated data
PCFA <- FindNeighbors(PCFA, reduction = "harmony", dims = 1:40)
PCFA <- FindClusters(PCFA, resolution = 0.5)
PCFA <- RunUMAP(PCFA, reduction = "harmony", dims = 1:40, return.model = T, verbose = F, seed.use = 22)

PCFA <- ProjectIntegration(object = PCFA, sketched.assay = "sketch", assay = "RNA", reduction = "harmony")

PCFA <- ProjectData(object = PCFA, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "harmony.full",
                    full.reduction = "harmony.full", dims = 1:40)
PCFA <- RunUMAP(PCFA, reduction = "harmony.full", dims = 1:40, reduction.name = "umap.full",
                return.model = T)

# annotation of clusters and removal of non-fibroblastic cells (Endothelial, LEC, Mural, doublets etc.)
# Repeat process using n=2500 cells
DefaultAssay(PCFA)="RNA"
PCFA[["RNA"]] <- JoinLayers(PCFA[["RNA"]])
PCFA <- NormalizeData(PCFA)

PCFA[["RNA"]] <- split(PCFA[["RNA"]], f = PCFA$Dataset)
PCFA <- FindVariableFeatures(PCFA, verbose = FALSE)
PCFA <- SketchData(object = PCFA, ncells = 2500, method = "LeverageScore", sketched.assay = "sketch", seed = 2324, over.write = T)
DefaultAssay(PCFA) <- "sketch"
PCFA <- FindVariableFeatures(PCFA, verbose = F, nfeatures = 3000)
PCFA <- ScaleData(PCFA, verbose = F)
PCFA <- RunPCA(PCFA, verbose = F)

# integrate the datasets using harmony
PCFA <- IntegrateLayers(
  object = PCFA, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE, group.by = c("Dataset", "Tech"), theta = c(1)
)

# cluster the integrated data
PCFA <- FindNeighbors(PCFA, reduction = "harmony", dims = 1:40)
PCFA <- FindClusters(PCFA, resolution = 0.5)
PCFA <- RunUMAP(PCFA, reduction = "harmony", dims = 1:40, return.model = T, verbose = F, seed.use = 22)

PCFA <- ProjectIntegration(object = PCFA, sketched.assay = "sketch", assay = "RNA", reduction = "harmony")

PCFA <- ProjectData(object = PCFA, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "harmony.full",
                    full.reduction = "harmony.full", dims = 1:40)
PCFA <- RunUMAP(PCFA, reduction = "harmony.full", dims = 1:40, reduction.name = "umap.full",
                return.model = T)

# re-cluster
DefaultAssay(PCFA)="RNA" 
PCFA <- FindNeighbors(PCFA, reduction = "harmony.full", dims = 1:40)
PCFA <- FindClusters(PCFA, resolution = 1.2)

