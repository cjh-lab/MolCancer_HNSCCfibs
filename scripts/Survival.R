#### survival analysis ####
# Example for GSE159067 
library(dplyr)

HNSCC_bulk <- read.delim("counts", header = T, row.names = 1)#GSE159067 counts
transposed_HNSCC_bulk <- t(HNSCC_bulk)

# GETTING FULL METADATA
library(GEOquery)

geo_accession <- "GSE159067"
gse <- getGEO(geo_accession, GSEMatrix = TRUE, AnnotGPL = FALSE)
metadata <- pData(gse[[1]])

rownames(transposed_HNSCC_bulk) <- rownames(metadata)
transposed_HNSCC_bulk_df <- as.data.frame(transposed_HNSCC_bulk)
transposed_HNSCC_bulk_df$geo_accession <- rownames(transposed_HNSCC_bulk_df)
counts_and_meta <- full_join(transposed_HNSCC_bulk_df, metadata, by = "geo_accession")

library(GSVA)
FRC_like_genes <- list(c(FRC_like_genes))
myCAF_genes <- list(c(myCAF_genes))
iCAF_genes <- list(c(iCAF_genes))
PI16_genes <- list(c(PI16_genes))

as.matrix(HNSCC_bulk) -> HNSCC_bulk

ssgsea_results <- gsva(expr = HNSCC_bulk, gset.idx.list = list(FRC_like_genes, myCAF_genes, iCAF_genes, PI16_genes),method = "ssgsea")
counts_and_meta2 <- cbind(counts_and_meta,t(ssgsea_results))
counts_and_meta2$'1' -> counts_and_meta2$FRC_like
counts_and_meta2$'2' -> counts_and_meta2$myCAF
counts_and_meta2$'3' -> counts_and_meta2$iCAF
counts_and_meta2$'4' -> counts_and_meta2$PI16

ssgsea_results <- t(ssgsea_results)
counts_and_meta2 <- cbind(counts_and_meta,ssgsea_results)

counts_and_meta2 <- counts_and_meta2 %>% dplyr::rename(os_event = 'os (event):ch1')
counts_and_meta2 <- counts_and_meta2 %>% dplyr::rename(os_time = 'os (month):ch1')
class(counts_and_meta2$os_time) <-'numeric' 
class(counts_and_meta2$os_event) <-'numeric' 

cutpoint_result <- surv_cutpoint(counts_and_meta2, time = 'os_time', event = 'os_event', variables = "FRC_like")

cut_point <- cutpoint_result$cutpoint$cutpoint
cut_point
counts_and_meta2$FRC_like <- counts_and_meta2$FRC_like
df <- counts_and_meta2

df <- df %>%
  mutate(FRC_like = ifelse(FRC_like >= cut_point, "high", "low"))
table(df$FRC_like)

cutpoint_result <- surv_cutpoint(df, time = 'os_time', event = 'os_event', variables = "myCAF")
cut_point <- cutpoint_result$cutpoint$cutpoint
cut_point

df$myCAF <- df$myCAF

df <- df %>%
  mutate(myCAF = ifelse(myCAF >= cut_point, "high", "low"))
table(df$myCAF)

cutpoint_result <- surv_cutpoint(df, time = 'os_time', event = 'os_event', variables = "iCAF")

cut_point <- cutpoint_result$cutpoint$cutpoint
cut_point
df$iCAF <- df$iCAF

df <- df %>%
  mutate(iCAF = ifelse(iCAF >= cut_point, "high", "low"))
table(df$iCAF)

cutpoint_result <- surv_cutpoint(df, time = 'os_time', event = 'os_event', variables = "PI16")

cut_point <- cutpoint_result$cutpoint$cutpoint
cut_point
df$PI16 <- df$PI16

df <- df %>%
  mutate(PI16 = ifelse(PI16 >= cut_point, "high", "low"))
table(df$PI16)

fit <- survfit(Surv(os_time, os_event) ~ FRC_like, data = df)
fit <- survfit(Surv(os_time, os_event) ~ myCAF, data = df)
fit <- survfit(Surv(os_time, os_event) ~ iCAF, data = df)
fit <- survfit(Surv(os_time, os_event) ~ FRC_like, data = df)


ggsurvplot(fit,
           pval = TRUE, conf.int = F,tables.theme = theme_pubr(),
           risk.table = T, 
           risk.table.col = "strata", 
           ggtheme = theme_pubr(), palette = c("coral2", "cadetblue", "orchid2", "pink2"),)+labs(x="Time (months)", y="Probability of OS") 

##### Cox regression survival ####
library(survival)
library(forestplot)

df$FRC_like <- factor(df$FRC_like)
df$FRC_like <- relevel(df$FRC_like, ref = "low")

df$myCAF <- factor(df$myCAF)
df$myCAF <- relevel(df$myCAF, ref = "low")

df$iCAF <- factor(df$iCAF)
df$iCAF <- relevel(df$iCAF, ref = "low")

df$PI16 <- factor(df$PI16)
df$PI16 <- relevel(df$PI16, ref = "low")

df$`age:ch1` <- as.numeric(df$`age:ch1`)
df$`Sex:ch1` <- factor(df$`Sex:ch1`)

cox_model1 <- coxph(Surv(os_time, os_event) ~ FRC_like + `Sex:ch1` +`age:ch1`, data = df, )
cox_model2 <- coxph(Surv(os_time, os_event) ~ myCAF + `Sex:ch1` +`age:ch1`, data = df, )
cox_model3 <- coxph(Surv(os_time, os_event) ~ iCAF + `Sex:ch1` +`age:ch1`, data = df, )
cox_model4 <- coxph(Surv(os_time, os_event) ~ PI16 + `Sex:ch1` +`age:ch1`, data = df, )

pdf("Fibro_MultiV_Cox.pdf", width = 6, height = 3, )
ggforest(cox_model1)
ggforest(cox_model2)
ggforest(cox_model3)
ggforest(cox_model4)
dev.off()
