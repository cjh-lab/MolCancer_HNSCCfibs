#### MiloR Differential abundance ####
library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(patchwork)
# read in seurat object containing all cells
HNSCC <- readRDS("scRNASeq_HNSCC_annotated_ref.RDS")
table(HNSCC$Labels2) #broad labels

lymphocyte_cells <- HNSCC$Labels2 %in% c("B_cell", "CD4_T_cell", "CD8_T_cell", "NK", "Treg", "Plasma_cell")
HNSCC$IsLymphocyte <- as.numeric(lymphocyte_cells)

HNSCC$PatientSource <- paste(HNSCC$orig.ident, HNSCC$source, sep = "_")

sample_proportions <- HNSCC@meta.data %>% filter(source == "tumour") %>% filter(Labels2 %in% c("B_cell", "CD4_T_cell", "CD8_T_cell", "NK", "Treg", "Plasma_cell", "cDC", "Macrophage", "Mast_cell", "Monocyte", "Neutrophils", "pDC")) %>%
  group_by(PatientSource, Sample.type2) %>%
  summarise(
    total_cells = n(),
    lymphocyte_count = sum(IsLymphocyte),
    lymphocyte_proportion = lymphocyte_count / total_cells
  )

expected_proportion <- mean(sample_proportions$lymphocyte_proportion)

# Compute Pearson residuals for each sample
sample_proportions <- sample_proportions %>%
  mutate(
    pearson_residual = (lymphocyte_proportion - expected_proportion) / 
      sqrt(expected_proportion * (1 - expected_proportion) / total_cells)
  )

# Assign Tertiles Based on Pearson Residuals
tertiles <- quantile(sample_proportions$pearson_residual, probs = c(0, 1/3, 2/3, 1))

sample_proportions <- sample_proportions %>%
  mutate(
    tertile = cut(pearson_residual, breaks = tertiles, labels = c("Low", "Moderate", "High"), include.lowest = TRUE)
  )

# read in fibroblasts
fibro <- readRDS("fibroblast_seurat_object.RDS")
fibro$PatientSource <- paste(fibro$orig.ident, fibro$source, sep = "_")

fibro@meta.data -> meta
meta <- meta %>%
  left_join(sample_proportions, by = "PatientSource")
table(meta$PatientSource, meta$tertile)

fibro@meta.data <- meta

fibro$Sample.type <- factor(
  fibro$hpv_stat2,
  levels = c("Positive_normal","Positive_tumour", "Negative_tumour", "Negative_normal"),
  labels = c("Normal", "HPV+ HNSCC", "HPV- HNSCC", "Normal")
)

DefaultAssay(fibro) = "integrated"
fibro_sce <- as.SingleCellExperiment(fibro, assay = "integrated")
fibro_milo <- Milo(fibro_sce)
traj_milo <- buildGraph(fibro_milo, k = 70, d = 20)
traj_milo <- makeNhoods(traj_milo, prop = 0.1, k =70, d=20, refined = TRUE)
plotNhoodSizeHist(traj_milo)

traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), sample="PatientSource")
head(nhoodCounts(traj_milo))

traj_design <- data.frame(colData(traj_milo))[,c("PatientSource", "source", "Sample.type", 'tertile')]
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$PatientSource
traj_design

traj_milo <- calcNhoodDistance(traj_milo, d=20)
# differential abundance HNSCC vs normal mucosa
da_results <- testNhoods(traj_milo, design = ~ source, design.df = traj_design)

da_results %>%
  arrange(- SpatialFDR) %>%
  head()

table(da_results$SpatialFDR < 0.05)

traj_milo <- buildNhoodGraph(traj_milo)

da_results <- annotateNhoods(traj_milo, da_results, coldata_col = "Labels")
head(da_results)
plotDAbeeswarm(da_results, group.by = "Labels", alpha = 0.05)

# differential abundance immune hot vs cold tumours
contrast.1 <- c("tertileHigh - tertileLow")
da_results <- testNhoods(traj_milo, design = ~ 0 + tertile, design.df = traj_design, model.contrasts = contrast.1, fdr.weighting="graph-overlap")

# or # fibroblast differential abundance HPV+ vs HPV- HNSCC
contrast.1 <- c("Sample.typeHPV+HNSCC - Sample.typeHPV- HNSCC")
da_results <- testNhoods(traj_milo, design = ~ 0 + Sample.type, design.df = traj_design, model.contrasts = contrast.1, fdr.weighting="graph-overlap")

table(da_results$SpatialFDR < 0.05)

traj_milo <- buildNhoodGraph(traj_milo)

da_results <- annotateNhoods(traj_milo, da_results, coldata_col = "Labels")

