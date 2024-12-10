#### label transfer ####
# example showing CRC dataset:GSE178341 
library(Seurat)
library(dplyr)
library(ggplot2)

ColonCancer <- Read10X_h5("GSE178341_crc10x_full_c295v4_submit.h5")#.h5 data for GSE178341
ColonCancer_meta <- read.csv(file = "GSE178341_crc10x_full_c295v4_submit_metatables.csv", header = T, row.names = 1)#meta
ColonCancer_clusters <- read.csv(file = "GSE178341_crc10x_full_c295v4_submit_cluster.csv", header = T, row.names = 1)#author clusters

ColonCancerSRT <- CreateSeuratObject(counts = ColonCancer, meta.data = ColonCancer_meta, min.cells = 3, min.features = 200)
ColonCancerSRT <- AddMetaData(ColonCancerSRT, metadata = ColonCancer_clusters)
ColonCancerSRT <- AddMetaData(ColonCancerSRT, metadata = ColonCancer_meta)

ColonCancerSRT[["pct.mt"]] <- PercentageFeatureSet(ColonCancerSRT, pattern = "^MT-")
VlnPlot(ColonCancerSRT, features = c("nFeature_RNA", "nCount_RNA", "pct.mt"), ncol = 3, pt.size = 0)
ColonCancerSRT <- subset(ColonCancerSRT, subset = nFeature_RNA < 6000 & pct.mt < 20)
# subset out myeloid
ColonCancerSRT <- subset(ColonCancerSRT, subset = clTopLevel == "Myeloid")
ColonCancerSRT <- NormalizeData(ColonCancerSRT)
CRC_myeloid[["RNA"]] <- as(object = CRC_myeloid[["RNA"]], Class = "Assay")
CRC_myeloid <- UpdateSeuratObject(CRC_myeloid)

#### Label Transfer ####
# read in annotated HNSCC myeloid cells
HNSCC_myeloid <- readRDS("HNSCC_myeloid.RDS")
CRC_myeloid <- SCTransform(CRC_myeloid)
Myeloid.anchors <- FindTransferAnchors(reference = HNSCC_myeloid, query = CRC_myeloid, 
                                       dims = 1:30, normalization.method = 'SCT',
                                       query.assay = "SCT")
predictions <- TransferData(anchorset = Myeloid.anchors, refdata = HNSCC_myeloid$Labels_final, 
                            dims = 1:30)
CRC_myeloid <- AddMetaData(CRC_myeloid, metadata = predictions)
predictions <- predictions %>%filter(prediction.score.max >=0.5)

CRC_meta <- read.csv("ColonCancerSRT_QC_meta.csv", header = T, row.names = 1)

merged_df <- merge(CRC_meta, predictions, by = "row.names", all = TRUE)
# include only tumour
merged_df <- merged_df %>%filter(SPECIMEN_TYPE == "T")

PCFA <- readRDS("PCFA.RDS")
CRC_fibro <- subset(PCFA, subset = Cancer == "Colon")
CRC_fibro_meta <- CRC_fibro@meta.data
CRC_fibro_meta <- CRC_fibro_meta %>%filter(source == "tumour")

##### calculate number of cells per sample per fibroblast subset #
table(CRC_fibro_meta$Labeled)
table(CRC_fibro_meta$orig.ident)

CRC_fibro_meta.df <- as.table(as.matrix(table(CRC_fibro_meta$Labeled,
                                              CRC_fibro_meta$orig.ident)))
CRC_fibro_meta.df <- as.data.frame(CRC_fibro_meta.df)
names(CRC_fibro_meta.df) <- c("SubPop", "SampleID", "Freq")
CRC_fibro_meta.df <- reshape2::dcast(CRC_fibro_meta.df[,],
                                     formula = SampleID ~ SubPop,
                                     value.var = "Freq") 
CRC_fibro_meta.df -> df1
df1_row_sums <- rowSums(df1[, -1])  
df1_proportions <- df1[, -1] / df1_row_sums 
df1_proportions <- cbind(df1[, 1, drop = FALSE], df1_proportions)
CRC_fibro_propitions <-df1_proportions

##calculate number of cells per sample per ann.fine in all cells #
merged_df
merged_df$Labels <- paste(merged_df$predicted.id, merged_df$clTopLevel, sep = "_")
table(merged_df$Labels)

CRC_all_meta.df <- as.table(as.matrix(table(merged_df$Labels,
                                            merged_df$orig.ident)))

CRC_all_meta.df <- as.data.frame(CRC_all_meta.df)
names(CRC_all_meta.df) <- c("SubPop", "SampleID", "Freq")
CRC_all_meta.df <- reshape2::dcast(CRC_all_meta.df[,],
                                   formula = SampleID ~ SubPop,
                                   value.var = "Freq") 
CRC_all_meta.df <- CRC_all_meta.df %>%filter(SampleID %in% CRC_fibro_propitions$SampleID)

CRC_all_meta.df -> df2
df2_row_sums <- rowSums(df2[, -1])
df2_proportions <- df2[, -1] / df2_row_sums
df2_proportions <- cbind(df2[, 1, drop = FALSE], df2_proportions)
CRC_allcells_propitions <-df2_proportions
CRC_allcells_propitions <- CRC_allcells_propitions %>%filter(SampleID %in% CRC_fibro_propitions$SampleID)

Combined_CRC <- cbind(CRC_allcells_propitions, CRC_fibro_propitions)
rownames(Combined_CRC) <- Combined_CRC$SampleID
colnames(Combined_CRC)

Combined_CRC<- Combined_CRC[,-21]
colnames(Combined_CRC) <- gsub("\\+", "", colnames(Combined_CRC))
colnames(Combined_CRC) <- gsub("\\)", "", colnames(Combined_CRC))
colnames(Combined_CRC) <- gsub("\\(", "", colnames(Combined_CRC))
colnames(Combined_CRC) <- gsub(" ", "_", colnames(Combined_CRC))
colnames(Combined_CRC) <- gsub("-", "_", colnames(Combined_CRC))

# IL11_CAF cors
Combined_CRC[,c(-1)] -> df
df <- df[,c(1:19, 26)]#select columns for cor
correlation_results <- list()

# Loop through each column
for (col in colnames(df)) {
  if (col == 'IL11_CAF') next
  correlation_result <- cor.PCFA(df$IL11_CAF, df[[col]], method = "spearman")
    correlation_df <- data.frame(
    Variable1 = 'IL11_CAF',
    Variable2 = col,
    Correlation = correlation_result$estimate,
    PValue = correlation_result$p.value
  )
    correlation_results[[col]] <- correlation_df
}

result_df <- bind_rows(correlation_results)
# significant results
result_df <- result_df %>%
  filter(PValue <= 0.05 | PValue == 0)
