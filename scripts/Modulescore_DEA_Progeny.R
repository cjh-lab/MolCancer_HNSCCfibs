#### Module scores ####
seurat.obj <- AddModuleScore(seurat.obj, assay = "RNA", features = gene_sig, name = "gene_sig_name", nbin = 24, ctrl = 100)

#### Differential expression analysis - scRNA-Seq/Pseudobulk ####
library(Seurat)
Markers <- FindAllMarkers(seurat.obj, assay = "RNA", min.pct = 0.25, logfc.threshold = 0.5, only.pos = T, test.use = 'wilcox')

#Pseudobulk
# aggregrate populations of interest (e.g., "IL11+ CAF", "IGF1+ CAF")
pseudobulk <- AggregateExpression(seurat.obj %>%subset(Labeled %in% c("IL11+ CAF", "IGF1+ CAF")), assays = "RNA", return.seurat = T, group.by = c("Labeled", "Dataset", "orig.ident"))
pseudobulk@assays[["RNA"]]@layers[["counts"]] <- as.matrix(pseudobulk@assays[["RNA"]]@layers[["counts"]])+1
UpdateSeuratObject(pseudobulk)
tail(Cells(pseudobulk))

pseudobulk$Labeled_Dataset <- rownames(pseudobulk@meta.data)

Idents(pseudobulk) <- "orig.ident"
library(DESeq2)
bulk.IL11.de <- FindMarkers(object = pseudobulk, 
                            ident.1 = "IL11+ CAF", ident.2 = "IGF1+ CAF",group.by = "Labeled",
                            test.use = "DESeq2", assay = "RNA")

#### Progeny ####
library(progeny)
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(tibble)
# read in fibroblast object
fibro <- readRDS("fibroblast_seurat_object.RDS")

DefaultAssay(fibro)="RNA"
table(fibro$Labels)
fibro <- SetIdent(fibro, value = "Labels")

CellsClusters <- data.frame(Cell = names(Idents(fibro)), 
                            CellType = as.character(Idents(fibro)),
                            stringsAsFactors = FALSE)

fibro <- progeny(fibro, scale=FALSE, organism="Human", top=500, perm=1, 
                 return_assay = TRUE)
fibro <- Seurat::ScaleData(fibro, assay = "progeny") 

progeny_scores_df <- 
  as.data.frame(t(GetAssayData(fibro, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 

progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

paletteLength = 100
myColor = colorRampPalette(c("blue", "white","red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))
progeny_hmap = pheatmap(t(summarized_progeny_scores_df[,-1]),fontsize=14, 
                        fontsize_row = 10, 
                        color=myColor, breaks = progenyBreaks, 
                        main = "PROGENy (500)", 
                        treeheight_col = 0,  border_color = NA)

