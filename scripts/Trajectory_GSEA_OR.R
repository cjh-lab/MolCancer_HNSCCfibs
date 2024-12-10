####  Slingshot ####
library(SingleCellExperiment)
library(ggplot2)
library(slingshot)
library(ggpubr)

fibro <- readRDS("fibroblast_seurat_object.RDS")

table(fibro$Labels)
sce <- as.SingleCellExperiment(fibro, assay = "integrated")

sce <- slingshot(sce, clusterLabels = sce$Labels, reducedDim = "UMAP",
                 allow.breaks = T, start.clus="PI16")
lnes <- getLineages(reducedDim(sce,"UMAP"),
                    sce$Labels, start.clus = "PI16")
lnes@metadata$lineages
crvs <- getCurves(lnes)
par(mar=c(5.1, 4.1, 4.1, 6), xpd=TRUE)

umap_coords <- as.data.frame(reducedDims(sce)$UMAP)
umap_coords$cell_type <- sce$Labels

p<- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
  geom_point(shape = 16, size = 0.8) +
  scale_color_manual(values = c("ADH1B+" = "#F8766D", "Proto-CAF" = "#619CFF", "myCAF" = "#00BA38",
                                "iCAF" = "#B79F00", "PI16" = "#00BFC4", "FRC-like" = "#F564E3")) +
  coord_fixed() +theme_pubr()

sds <- as.SlingshotDataSet(crvs)
sling_curves <- slingCurves(sds)
curve_df <- do.call(rbind, lapply(names(sling_curves), function(curve) {
  data.frame(sling_curves[[curve]]$s, curve_id = curve)
}))

p + geom_path(data = curve_df, aes(x = UMAP_1, y = UMAP_2, group = curve_id),
              color = 'black', size = 1)&NoAxes()&NoLegend()


####  Monocle3 ####
library(Signac) 
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
library(dplyr)
set.seed(22)

DefaultAssay(fibro) = "integrated"
fibro.cds <- as.cell_data_set(fibro)
fibro.cds <- cluster_cells(cds = fibro.cds, reduction_method = "UMAP")

p1 <- plot_cells(fibro.cds, color_cells_by = "cluster", reduction_method = "UMAP")
p2 <- plot_cells(fibro.cds, color_cells_by = "partition", reduction_method = "UMAP")
wrap_plots(p1, p2)

fibro.cds <- learn_graph(fibro.cds, use_partition = T, close_loop = F)
plot_cells(fibro.cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

# SELECTED PI16+ fibroblast as starting nodes
fibro.cds <- order_cells(fibro.cds, reduction_method = "UMAP")

# plot trajectories colored by pseudotime
plot_cells(
  cds = fibro.cds,
  color_cells_by = "Labels",
  show_trajectory_graph = TRUE, label_roots = F,label_branch_points = F, label_leaves = F, cell_size = 0.3, label_cell_groups = F
)

# Extract the pseudotime values and add to the Seurat object
fibro <- AddMetaData(
  object = fibro,
  metadata = fibro.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Fibro.pseudotime"
)

# now calculated trajectory with integrated assay, need to now make a new cds object with RNA assay 
# then add the trajectory constructed from integrated object
DefaultAssay(fibro) = "RNA"
RNA_from_seurat <- as.cell_data_set(fibro)
RNA_from_seurat <- estimate_size_factors(RNA_from_seurat)

RNA_from_seurat <- cluster_cells(cds = RNA_from_seurat, reduction_method = "UMAP")

fibro_combined.cds <- fibro.cds

fibro_combined.cds@assays <- RNA_from_seurat@assays
fibro_combined.cds@elementMetadata <- RNA_from_seurat@elementMetadata
fibro_combined.cds@int_elementMetadata <- RNA_from_seurat@int_elementMetadata
fibro_combined.cds@rowRanges <- RNA_from_seurat@rowRanges

#fibro_combined.cds now contains the trajectory constructed using the integrated assay
#but contains expression data from the RNA assay.
# plot trajectories colored by pseudotime

rowData(fibro.cds)$gene_name <- rownames(fibro.cds)
rowData(fibro.cds)$gene_short_name <- rowData(fibro.cds)$gene_name

fibro.cds_backup <- fibro.cds
fibro.cds <- estimate_size_factors(fibro.cds)
fibro.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(fibro[["RNA"]])

# Subset each lineage individually (e.g., FRC-like, myCAF or iCAF)
cds_subset_iCAF <- choose_cells(fibro.cds)#cells in iCAF lineage selected
## Calculate size factors using built-in function in monocle3
cds_subset_iCAF <- estimate_size_factors(cds_subset_iCAF)

# iCAF lineage
subset_pr_test_res_iCAF <- graph_test(cds_subset_iCAF, neighbor_graph="principal_graph", cores=4)
pr_deg_ids_iCAF <- row.names(subset(subset_pr_test_res_iCAF, q_value < 0.05))

#plotting heatmap of genes test
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(monocle3)

subset_pr_test_res_iCAF
subset_pr_test_res_myCAF#repeat above for myCAF
subset_pr_test_res_FRC#repeat above for FRC-like

genes <- row.names(subset(subset_pr_test_res_iCAF, q_value < 0.05 & morans_I > 0.25))
genes <- row.names(subset(subset_pr_test_res_myCAF, q_value < 0.05 & morans_I > 0.25))
genes <- row.names(subset(subset_pr_test_res_FRC, q_value < 0.05 & morans_I > 0.25))
genes

#### OR analysis of Monocle lineage genes ####
library(enrichR)
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") 
}
if (websiteLive) dbs <- listEnrichrDbs()

dbs <- c("KEGG_2021_Human", "MSigDB_Hallmark_2020")
if (websiteLive) {
  enriched <- enrichr(c(genes), dbs)
}
enriched$KEGG_2021_Human$Term <- paste(enriched$KEGG_2021_Human$Term, ": KEGG",sep = " ")
enriched$MSigDB_Hallmark_2020$Term <- paste(enriched$MSigDB_Hallmark_2020$Term, ": HM",sep = " ")

combined<-rbind(enriched$KEGG_2021_Human, enriched$MSigDB_Hallmark_2020)

#### GSEA ####
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
# repeat following for FRC-like
iCAF_markers_for_GSEA2 <- FindMarkers(fibro, ident.1 = "iCAF", min.pct = 0.1)

iCAF_markers_for_GSEA2$SYMBOL <- rownames(iCAF_markers_for_GSEA2)
iCAF_gsea <- bitr(rownames(iCAF_markers_for_GSEA2), fromType = "SYMBOL",
                  toType = c("SYMBOL", "ENTREZID"),
                  OrgDb = org.Hs.eg.db)

iCAF_GSEA <- merge(iCAF_markers_for_GSEA2, iCAF_gsea, by = "SYMBOL")
iCAF_GSEA <-iCAF_GSEA %>%filter(p_val_adj < 0.05)

d <- iCAF_GSEA
column_names <- names(d)
new_column_order <- c(column_names[7], column_names[-7])

# Rearrange the dataframe with the new column order
d <- d[, new_column_order]
d <- d[,c(1, 4)]

## 1st column is ID
## 2nd column is FC

geneList = d[,2]
names(geneList) = as.character(d[,1])
geneList = sort(geneList, decreasing = TRUE)

# hallmarks 
library(msigdbr)
msigdbr_show_species()
m_df <- msigdbr(species = "Homo sapiens")
head(m_df, 2) %>% as.data.frame

Hallmarks <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(Hallmarks)

em2 <- GSEA(geneList, TERM2GENE = Hallmarks)
head(em2)
gseaplot2(em2, geneSetID = 8, pvalue_table = F,
          color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "line", subplots = 1:5)

