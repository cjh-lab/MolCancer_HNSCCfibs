#### Create seurat object ####
library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(sctransform)

seurat.obj <- CreateSeuratObject(counts, min.cells = 3, min.features = 200)
# % Mitochondrial genes
seurat.obj <- PercentageFeatureSet(seurat.obj, pattern = "^MT-", col.name = "pct.mito")
median(seurat.obj@meta.data$pct.mito)
mad(seurat.obj@meta.data$pct.mito)
thresh = median(seurat.obj@meta.data$pct.mito) +
  3*mad(seurat.obj@meta.data$pct.mito)
VlnPlot(seurat.obj, "pct.mito", group.by = "case",
        pt.size = 0.1) & geom_hline(yintercept = thresh, linetype = "dashed")
VlnPlot(seurat.obj, c("nCount_RNA", "nFeature_RNA"))
# Filtering
seurat.obj <- subset(seurat.obj, subset = nFeature_RNA <= 6000 & pct.mito < thresh)

#### Seurat RPCA integration - HPC slurm script ####
#!/usr/bin/env Rscript

# Load libraries
print("Loading Packages")
library(Seurat)
library(dplyr)
library(sctransform)

#options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

print("Loading seurat object")
# Load seurat object 
seurat.obj <- readRDS(args[1])

HNSCC.list <- SplitObject(seurat.obj, split.by = "PatientID")
HNSCC.list <- lapply(X = HNSCC.list, FUN = SCTransform)

features <- SelectIntegrationFeatures(object.list = HNSCC.list, nfeatures = 3000)
HNSCC.list <- PrepSCTIntegration(object.list = HNSCC.list, anchor.features = features)
HNSCC.list <- lapply(X = HNSCC.list, FUN = RunPCA, features = features)

HNSCC.anchors <- FindIntegrationAnchors(object.list = HNSCC.list, normalization.method = "SCT",
                                        anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 5)
HNSCC.integrated <- IntegrateData(anchorset = HNSCC.anchors, normalization.method = "SCT", dims = 1:30)
HNSCC.integrated <- RunPCA(HNSCC.integrated, verbose = FALSE)
HNSCC.integrated <- RunUMAP(HNSCC.integrated, reduction = "pca", dims = 1:20)
HNSCC.integrated <- FindNeighbors(object = HNSCC.integrated, dims = 1:20)
HNSCC.integrated <- FindClusters(object = HNSCC.integrated)

Markers <- FindAllMarkers(seurat.obj, min.pct = 0.25, logfc.threshold = 0.5, only.pos = T)

saveRDS(HNSCC.integrated, args[2])
#### Integration with GSE164690 (Kurten et al., ) ####
#!/usr/bin/env Rscript

# Iridis RPCA integration script for HNSCC EPG and Kurten datasets

# Load libraries
print("Loading Seurat")
library(Seurat)

#options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

print("Loading seurat objects (post-QC")
# Seurat_object
EPG_seurat <- readRDS(args[1])
Kurten_seurat <- readRDS(args[2])

print("Change orig.ident metadata column in EPG data to PatientID")
# Change orig.ident metadata column in EPG data to PatientID
EPG_seurat$orig.ident <- EPG_seurat$PatientID

print("remove HPV16 genes from Kurten object")
# remove HPV16 genes from Kurten object
counts <- GetAssayData(Kurten_seurat, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c("E1", "E2", "E5", "E6", "E7", "L1","L2"))),]
Kurten_seurat <- subset(Kurten_seurat, features = rownames(counts))
rm(counts)
gc()

print("merging EPG and Kurten")
# merging EPG and seurat
Merged <- merge(EPG_seurat, Kurten_seurat)

print("split into list by orig.ident")
# split into list by orig.ident
HNSCC.list <- SplitObject(Merged, split.by = "orig.ident")
rm(Merged)
gc()

print("SCTransform")
HNSCC.list <- lapply(X = HNSCC.list, FUN = SCTransform)

features <- SelectIntegrationFeatures(object.list = HNSCC.list, nfeatures = 3000)
HNSCC.list <- PrepSCTIntegration(object.list = HNSCC.list, anchor.features = features)
HNSCC.list <- lapply(X = HNSCC.list, FUN = RunPCA, features = features)

print("find integration anchors")
HNSCC.anchors <- FindIntegrationAnchors(object.list = HNSCC.list, normalization.method = "SCT",
                                        anchor.features = features, reduction = "rpca", dims = 1:30)
print("IntegrateData")
HNSCC.integrated <- IntegrateData(anchorset = HNSCC.anchors, normalization.method = "SCT", dims = 1:30)
print("RunPCA")
HNSCC.integrated <- RunPCA(HNSCC.integrated, verbose = FALSE)
HNSCC.integrated <- RunUMAP(HNSCC.integrated, reduction = "pca", dims = 1:20)
HNSCC.integrated <- FindNeighbors(object = HNSCC.integrated, dims = 1:20)
HNSCC.integrated <- FindClusters(object = HNSCC.integrated)
saveRDS(HNSCC.integrated, args[3])


#### Creation of HNSCC non-immune and immune atlases ####
# as above but for immune and non-immune cell objects separately (rather than all cells). 

