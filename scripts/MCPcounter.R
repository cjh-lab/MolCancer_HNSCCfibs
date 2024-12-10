#### MCP-counter for Visium data ####
library(Seurat)
library(dpylr)
library(ggplot2)
library(ggpubr)
library(immunedeconv)
# extract Spatial_data from srt object
HNSCC_Visium_srt@assays$Spatial@data -> Spatial_data

mcp_counter_output <- immunedeconv::deconvolute(Spatial_data, 'mcp_counter')
HNSCC_Visium_srt@meta.data -> meta
meta <- cbind(meta, mcp_counter_output)

meta$T_cell <- meta$`T.cell`
meta$B_cell <- meta$`B.cell`
meta$Lymphocyte <- meta$T_cell + meta$B_cell + meta$NK.cell + meta$cytotoxicity.score

result <- meta %>%
  group_by(Source.SampleID, HPVstatus) %>%
  summarize(
    B_cell_sum = sum(B_cell),
    T_cell_sum = sum(T_cell),
    CAF_sum = sum(Cancer.associated.fibroblast),
    Lymphocyte_sum = sum(Lymphocyte),
    count = n()
  )

#### MCP-counter for Bulk RNA-Seq data ####
# read in TPM counts
tumour_only_HNSCC <- read.csv("HNSCC_TPM.csv")

mcp_counter_output <- immunedeconv::deconvolute(tumour_only_HNSCC, 'mcp_counter')
mcp_counter_output$cell_type <- paste(mcp_counter_output$cell_type, "mcp_counter", sep = "_")
names(mcp_counter_output) <- gsub("\\.", "-", names(mcp_counter_output))

