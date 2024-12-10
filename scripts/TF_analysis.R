# TF analysis
## We load the required packages
library(Seurat)
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)

fibro <- readRDS("fibroblast_seurat_object.RDS")

Dorothea_decoupleR <- decoupleR::get_dorothea(organism = "human", levels = c("A", "B", "C"), )

# Extract the normalized log-transformed counts
mat <- as.matrix(fibro@assays$RNA@data)
# Run wmean
acts <- run_wmean(mat=mat, net=Dorothea_decoupleR, .source='source', .target='target',
                  .mor='mor', times = 100, minsize = 5)
acts

fibro[['tfswmean']] <- acts %>%
  filter(statistic == 'norm_wmean') %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = fibro) <- "tfswmean"
# Scale the data
fibro <- ScaleData(fibro)
fibro@assays$tfswmean@data <- fibro@assays$tfswmean@scale.data

p1 <- DimPlot(fibro, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend() + ggtitle('Cell types')
p2 <- (FeaturePlot(fibro, features = c("RELA")) & 
         scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
  ggtitle('RELA activity')
DefaultAssay(object = fibro) <- "RNA"
p3 <- FeaturePlot(fibro, features = c("RELA")) + ggtitle('RELA expression')
DefaultAssay(object = fibro) <- "tfswmean"
p1 | p2 | p3

DefaultAssay(object = fibro) <- "tfswmean"
p1 <- DimPlot(fibro, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend() + ggtitle('Cell types')
p2 <- (FeaturePlot(fibro, features = c("RELB")) & 
         scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
  ggtitle('RELB activity')
DefaultAssay(object = fibro) <- "RNA"
p3 <- FeaturePlot(fibro, features = c("RELB")) + ggtitle('RELB expression')
DefaultAssay(object = fibro) <- "tfswmean"
p1 | p2 | p3

n_tfs <- 25
# Extract activities from object as a long dataframe
df <- t(as.matrix(fibro@assays$tfswmean@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(fibro)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Get top tfs with more variable means across clusters
tfs <- df %>%
  group_by(source) %>%
  summarise(std = sd(mean)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)

# Subset long data frame to top tfs and transform to wide matrix
top_acts_mat <- df %>%
  filter(source %in% tfs) %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()


