# ST Ligand Receptor Analysis
library(Seurat)
library(dplyr)
library(ggplot2)
library(nichenetr)

# read in merged visium seurat object with deconvolution metadata
Merged_Visium_srt <- readRDS("merged_visium_seurat.RDS")

Merged_Visium_srt$iCAF_niche <- 'Other'
# Identify rows where all values in the specified columns are greater than or equal to 0.05
iCAF_niche_rows <- rownames(Merged_Visium_srt@meta.data)[
  Merged_Visium_srt@meta.data$iCAF_RCTD >= 0.05 &
    (Merged_Visium_srt@meta.data$Neutrophils_RCTD >= 0.05 |
       Merged_Visium_srt@meta.data$Monocyte_RCTD >= 0.05)
]

Merged_Visium_srt@meta.data[iCAF_niche_rows, "iCAF_niche"] <- "iCAF_niche"
de_markers_iCAF_niche <- FindMarkers(Merged_Visium_srt, ident.1 = 'iCAF_niche', group.by = 'iCAF_niche', logfc.threshold = 1, min.pct = 0.25, only.pos = T)
de_markers_iCAF_niche$gene <- rownames(de_markers_iCAF_niche)

Merged_Visium_srt$FRC_like_niche <- 'Other'
FRC_like_niche_rows <- rownames(Merged_Visium_srt@meta.data)[
  Merged_Visium_srt@meta.data$FRC_like_RCTD >= 0.05 &
    (Merged_Visium_srt@meta.data$B_cell_RCTD >= 0.05 |
       Merged_Visium_srt@meta.data$CD4_T_cell_RCTD >= 0.05)
]

Merged_Visium_srt@meta.data[FRC_like_niche_rows, "FRC_like_niche"] <- "FRC_like_niche"

de_markers_FRC_likeniche <- FindMarkers(Merged_Visium_srt, ident.1 = 'FRC_like_niche', group.by = 'FRC_like_niche', logfc.threshold = 1, min.pct = 0.25, only.pos = T)
de_markers_FRC_likeniche$gene <- rownames(de_markers_FRC_likeniche)

####myCAF ####
Merged_Visium_srt$myCAF_niche <- 'Other'
myCAF_niche_rows <- rownames(Merged_Visium_srt@meta.data)[
  Merged_Visium_srt@meta.data$myCAF_RCTD >= 0.05 &
    (Merged_Visium_srt@meta.data$Macrophage_RCTD >= 0.05 |
       Merged_Visium_srt@meta.data$Epithelial_RCTD >= 0.05)
]


Merged_Visium_srt@meta.data[myCAF_niche_rows, "myCAF_niche"] <- "myCAF_niche"

de_markers_myCAFniche <- FindMarkers(Merged_Visium_srt, ident.1 = 'myCAF_niche', group.by = 'myCAF_niche', logfc.threshold = 1, min.pct = 0.25, only.pos = T)
de_markers_myCAFniche$gene <- rownames(de_markers_myCAFniche)

####PI16 ####
Merged_Visium_srt$PI16_niche <- 'Other'
PI16_niche_rows <- rownames(Merged_Visium_srt@meta.data)[
  Merged_Visium_srt@meta.data$PI16_RCTD >= 0.05 &
    (Merged_Visium_srt@meta.data$Endothelial_RCTD >= 0.05 |
       Merged_Visium_srt@meta.data$Mural_RCTD >= 0.05)
]

Merged_Visium_srt@meta.data[PI16_niche_rows, "PI16_niche"] <- "PI16_niche"

de_markers_PI16niche <- FindMarkers(Merged_Visium_srt, ident.1 = 'PI16_niche', group.by = 'PI16_niche', logfc.threshold = 1, min.pct = 0.25, only.pos = T)
de_markers_PI16niche$gene <- rownames(de_markers_PI16niche)

de_markers_FRC_like$cluster <- "FRC_like"
de_markers_iCAF$cluster <- "iCAF"
de_markers_myCAF$cluster <- "myCAF"
de_markers_PI16$cluster <- "PI16"
de_markers_all <- rbind(de_markers_FRC_like,de_markers_iCAF,de_markers_myCAF,de_markers_PI16 )

lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor)

lr_network$ligand -> ligands
lr_network$receptor -> receptors

LRs <- c(lr_network$ligand, lr_network$receptor)
# all
de_markers_all_LigRec <- de_markers_all %>%filter(gene %in% LRs)

HNSCC <- readRDS("scRNASeq_HNSCC_annotated_ref.RDS")
HNSCC <- SetIdent(HNSCC, value = 'Labels')#fine annotated cell types

FRC_like_niche_ligands <- lr_network %>%filter(ligand %in% de_markers_FRC_like$gene)

#### FRC_like ####
receiver = "FRC_like"
expressed_LRs_receiver = get_expressed_genes(receiver, HNSCC, pct = 0.25, assay_oi = 'RNA')
FRC_like_expressed_receptors <- lr_network %>%filter(receptor %in% expressed_LRs_receiver)

FRC_like_ligands_withReceptors <- FRC_like_niche_ligands%>%filter(receptor %in% FRC_like_expressed_receptors$receptor)

# finally filter for unique ligands in each niche, i.e., FRC_like_ligands_withReceptors$ligands not found in myCAF/Pi16/iCAF ST niche
iCAF_ST_ligands <- de_markers_iCAF %>%filter(avg_log2FC > 1) %>%pull(gene)
myCAF_ST_ligands <- de_markers_myCAF %>%filter(avg_log2FC > 1) %>%pull(gene)
PI16_ST_ligands <- de_markers_PI16 %>%filter(avg_log2FC > 1) %>%pull(gene)

FRC_like_ligands_withReceptors_unique <- FRC_like_ligands_withReceptors %>%filter(!ligand %in% c(iCAF_ST_ligands,myCAF_ST_ligands,PI16_ST_ligands ))

# senders=highly correlated scRNA-Seq immune cell types
sender = c("Cycling_B_cell","B_cell_FCRL4", "B_cell_YBX3", "IgM_plasma_cell", "B_cell_TNF", "NK_KIT", "CD4_TFH", "IgA_plasma_cell", "IgG_plasma_cell", "Memory_B_cell_IGKC_high", "GC_B_cell", "FRC_like")
list_expressed_LRs_sender = sender %>% unique() %>% lapply(get_expressed_genes, HNSCC, 0.1) # lapply to get the expressed LRs of every sender cell type separately here
expressed_LRs_sender = list_expressed_LRs_sender %>% unlist() %>% unique()

FRC_like_ligands_withReceptors_final <- FRC_like_ligands_withReceptors_unique %>%filter(ligand %in% expressed_LRs_sender)

#### iCAF ####
iCAF_niche_ligands <- lr_network %>%filter(ligand %in% de_markers_iCAF$gene)

receiver = "iCAF"
expressed_LRs_receiver = get_expressed_genes(receiver, HNSCC, pct = 0.25, assay_oi = 'RNA')
iCAF_expressed_receptors <- lr_network %>%filter(receptor %in% expressed_LRs_receiver)

# now filter iCAF_niche_ligands for only those that have iCAF expressed receptors
iCAF_ligands_withReceptors <- iCAF_niche_ligands%>%filter(receptor %in% iCAF_expressed_receptors$receptor)

# finally filter for unique ligands in each niche, i.e., iCAF_ligands_withReceptors$ligands not found in myCAF/Pi16/FRC-like ST niche
iCAF_ST_ligands <- de_markers_iCAF %>%filter(avg_log2FC > 1) %>%pull(gene)
myCAF_ST_ligands <- de_markers_myCAF %>%filter(avg_log2FC > 1) %>%pull(gene)
PI16_ST_ligands <- de_markers_PI16 %>%filter(avg_log2FC > 1) %>%pull(gene)
FRC_like_ST_ligands <- de_markers_FRC_like %>%filter(avg_log2FC > 1) %>%pull(gene)

iCAF_ligands_withReceptors_unique <- iCAF_ligands_withReceptors %>%filter(!ligand %in% c(FRC_like_ST_ligands,myCAF_ST_ligands,PI16_ST_ligands ))

# senders=highly correlated scRNA-Seq immune cell types
sender = c("Monocyte_IL1B", "CD4_IFN", "Monocyte_CXCL10", "Monocyte_CD14", "CD8_IFN", "cDC_LAMP3+", "CD8_GZMB", "Cycling_T_cells", "iCAF")
list_expressed_LRs_sender = sender %>% unique() %>% lapply(get_expressed_genes, HNSCC, 0.1) # lapply to get the expressed LRs of every sender cell type separately here
expressed_LRs_sender = list_expressed_LRs_sender %>% unlist() %>% unique()

iCAF_ligands_withReceptors_final <- iCAF_ligands_withReceptors_unique %>%filter(ligand %in% expressed_LRs_sender)

