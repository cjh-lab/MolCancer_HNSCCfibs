#### Pan cancer TCGA analysis  ####
library(UCSCXenaTools)
data(XenaData)
library(GSVA)
FRC_like_genes <- readRDS("FRC_like_genes.RDS")
TLS_sig <- c("CD79B", "CD1D", "CCR6", "LAT","SKAP1", "CETP", "EIF1AY", "RBP5", "PTGDS")#cabrita 9 gene

# clinical
XenaGenerate(subset = XenaHostNames=="pancanAtlasHub") %>% 
  XenaFilter(filterDatasets = "TCGA_phenotype_denseDataOnlyDownload.tsv") -> df_todo

XenaQuery(df_todo) %>%
  XenaDownload() -> xe_download

cli_PanCancer = XenaPrepare(xe_download)
table(cli_PanCancer$`_primary_disease`)
table(cli_PanCancer$sample_type)
cli_PanCancer$sample -> cli_PanCancer$sampleID

# gene expression
df_todo
df_todo2 <- df_todo
df_todo2@hosts <- "https://pancanatlas.xenahubs.net"
df_todo2@cohorts <- "TCGA Pan-Cancer (PANCAN)"
df_todo2@datasets <- "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"

XenaQuery(df_todo2) %>%
  XenaDownload() -> xe_download
Sys.setenv(VROOM_CONNECTION_SIZE="252000000")
PANCAN = XenaPrepare(xe_download)

PANCAN_tumour_only <- as.data.frame(t(PANCAN))
colnames(PANCAN_tumour_only) <- PANCAN_tumour_only[1,]
PANCAN_tumour_only <- PANCAN_tumour_only[-1,]
PANCAN_tumour_only$sampleID <- rownames(PANCAN_tumour_only)
PANCAN_tumour_only <- merge(PANCAN_tumour_only, cli_PanCancer, by = "sampleID", all.x = TRUE)

PANCAN_tumour_only <- PANCAN_tumour_only %>%filter(sample_type == "Primary Tumor")
PANCAN_meta_complete <- PANCAN_tumour_only[,1:20531]
PANCAN_meta_complete <- subset(PANCAN_meta_complete, select = -c(which(colSums(is.na(PANCAN_meta_complete)) > 0)))

PANCAN_tumour_only2 <- PANCAN_tumour_only[,c(1, 20532:20535)]

PANCAN_MATRIX <- PANCAN_tumour_only[,1:20531]
PANCAN_MATRIX <- as.data.frame(t(PANCAN_MATRIX))
PANCAN_MATRIX$sample <- rownames(PANCAN_MATRIX)
PANCAN_MATRIX <- PANCAN_MATRIX[-1,]
PANCAN_MATRIX <- PANCAN_MATRIX[,-9699]
PANCAN_MATRIX <-as.matrix(PANCAN_MATRIX)
class(PANCAN_MATRIX)="numeric"
PANCAN_MATRIX <- PANCAN_MATRIX[-c(1:28),]
PANCAN_MATRIX <- PANCAN_MATRIX[-1,]

mat_pancan <- PANCAN_MATRIX[complete.cases(PANCAN_MATRIX), ]

ssgsea_results_PanCan_FRC_like <- gsva(expr = mat_pancan, gset.idx.list = list(FRC_like_genes),method = "ssgsea")
ssgsea_results_PanCan_TLS <- gsva(expr = mat_pancan, gset.idx.list = list(TLS_sig),method = "ssgsea")

PANCAN_with_ssGSEA_meta<- cbind(PANCAN_tumour_only2, t(ssgsea_results_PanCan_FRC_like),t(ssgsea_results_PanCan_TLS))

PANCAN_with_ssGSEA_meta$Cancer <- PANCAN_with_ssGSEA_meta$`_primary_disease`
PANCAN_with_ssGSEA_meta$`t(ssgsea_results_PanCan_FRC_like)` ->PANCAN_with_ssGSEA_meta$FRC_like_score
PANCAN_with_ssGSEA_meta$`t(ssgsea_results_PanCan_TLS)` ->PANCAN_with_ssGSEA_meta$TLS_score

PANCAN_with_ssGSEA_meta$'1' ->PANCAN_with_ssGSEA_meta$FRC_like_score

PANCAN_with_ssGSEA_meta$Cancer <- with(PANCAN_with_ssGSEA_meta, reorder(Cancer, FRC_like_score, mean, na.rm = TRUE))

# Adding HNSCC_HPV status
HPV_meta <- read.csv("TCGA_HNSC_HPV_metadata.csv", header = T)

merged_df <- merge(PANCAN_with_ssGSEA_meta, HPV_meta[, c("sampleID", "HPV.status")], by = "sampleID", all.x = TRUE)
merged_df$Cancer2 <- paste(merged_df$Cancer, merged_df$HPV.status, sep = "_")

merged_df$Cancer2 <- gsub("_NA", "", merged_df$Cancer2)
merged_df$Cancer2 <- with(merged_df, reorder(Cancer2, FRC_like_score, mean, na.rm = TRUE))

merged_df$TLS <- merged_df$`t(ssgsea_results_PanCan_TLS)`
