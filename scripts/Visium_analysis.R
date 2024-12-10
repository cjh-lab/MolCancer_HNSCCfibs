#### Processing of Visium data ####
# Seurat integration of all ST sections

suppressPackageStartupMessages(require(Matrix))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(SeuratData))
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(patchwork))
suppressPackageStartupMessages(require(dplyr))

base_dir <- "/Visium_Data"
sample_dirs <- list.dirs(base_dir, recursive = FALSE)
seurat_objects <- list()
for (sample_dir in sample_dirs) {
  sample_name <- basename(sample_dir)
    seurat_object <- Load10X_Spatial(
    data.dir = sample_dir,
    filename = "filtered_feature_bc_matrix.h5",
    slice = sample_name,
    assay = "Spatial"
  )
    seurat_objects[[sample_name]] <- seurat_object
}

st.list = lapply(seurat_objects, SCTransform, assay = "Spatial")
st.features = SelectIntegrationFeatures(st.list, nfeatures = 3000, verbose = T)
st.list <- PrepSCTIntegration(object.list = st.list, anchor.features = st.features,
                              verbose = FALSE)
int.anchors <- FindIntegrationAnchors(object.list = st.list, normalization.method = "SCT",
                                      verbose = FALSE, anchor.features = st.features)
HNSCC.integrated <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT",
                                  verbose = FALSE)
DefaultAssay(HNSCC.integrated)="integrated"
HNSCC.integrated <- RunPCA(HNSCC.integrated, verbose = FALSE)
HNSCC.integrated <- FindNeighbors(HNSCC.integrated, dims = 1:20)
HNSCC.integrated <- FindClusters(HNSCC.integrated, verbose = FALSE, resolution = 0.5)
HNSCC.integrated <- RunUMAP(HNSCC.integrated, dims = 1:20)

# removal of poor quality clusters based on nFeature_Spatial, nCount_Spatial, % Mt, % Hb

#### Robust cell type deconvolution ####
# Spatial deconvolution using RCTD
# load in scRNA-Seq ref
HNSCC <- readRDS("scRNASeq_HNSCC_annotated_ref.RDS")

# set up reference
HNSCC <- UpdateSeuratObject(HNSCC)
HNSCC$celltype <- HNSCC$Labels

Idents(HNSCC) <- "celltype"

# extract information to pass to the RCTD Reference function
counts <- HNSCC[["RNA"]]$counts
cluster <- as.factor(HNSCC$celltype)
names(cluster) <- colnames(HNSCC)
nUMI <- HNSCC$nCount_RNA
names(nUMI) <- colnames(HNSCC)
reference <- Reference(counts, cluster, nUMI)

# for each patient Visium_srt object 
counts <- Visium_srt[["Spatial"]]$counts
coords <- GetTissueCoordinates(Visium_srt, image = "slice1")
colnames(coords) <- c("x", "y")
coords[is.na(colnames(coords))] <- NULL
query <- SpatialRNA(coords, counts, colSums(counts))

RCTD <- create.RCTD(query, reference, max_cores = 8)
RCTD <- run.RCTD(RCTD, doublet_mode = "full")

weights <- RCTD@results$weights 
norm_weights <- normalize_weights(weights) #Normalize per spot weights so cell type probabilities sum to 1 for each spot


#### SpatialDWLS ####
library(Giotto)
library(Seurat)
library(dplyr)

Visium_giotto <-seuratToGiotto(
  Visium_srt,
  spatial_assay = "Spatial",
  dim_reduction = c("pca", "umap")
)
Visium_giotto <- normalizeGiotto(gobject = Visium_giotto, scalefactor = 6000)
Visium_giotto <- addStatistics(gobject = Visium_giotto, expression_values = 'raw')
Visium_giotto <- calculateHVF(gobject = Visium_giotto)
Visium_giotto <- runPCA(gobject = Visium_giotto, center = TRUE, scale_unit = TRUE)
screePlot(Visium_giotto, ncp = 30)
#cluster and run UMAP #
Visium_giotto <- createNearestNetwork(gobject = Visium_giotto,
                                     dim_reduction_to_use = 'pca', dim_reduction_name = 'pca',
                                     dimensions_to_use = 1:10, k = 15)
Visium_giotto <- doLeidenCluster(gobject = Visium_giotto, resolution = 0.4, n_iterations = 1000)
Visium_giotto = runUMAP(Visium_giotto)

plotUMAP(gobject = Visium_giotto,
         cell_color = 'leiden_clus', show_NN_network = T, point_size = 1.5
)
spatDimPlot(gobject = Visium_giotto, cell_color = 'leiden_clus',
            dim_point_size = 2, spat_point_size = 2.5)

# reading in HNSCC sc ref and converting to giotto object and normalizing data
HNSCC_ref <- readRDS("scRNASeq_HNSCC_annotated_ref.RDS")
HNSCC_ref <- SetIdent(HNSCC_ref, value = "Labels")
set.seed(22)
HNSCC_ref_subset <- subset(x = HNSCC_ref, downsample = 500)
# remove mitochondrial and ribosomal genes
ribosomal_pattern <- "^(RPL|RPS)"
mitochondrial_pattern <- "^(MT-|MT\\.)"

# Identify ribosomal and mitochondrial genes
ribosomal_genes <- grep(ribosomal_pattern, rownames(HNSCC_ref_subset@assays$RNA), value = TRUE)
mitochondrial_genes <- grep(mitochondrial_pattern, rownames(HNSCC_ref_subset@assays$RNA), value = TRUE)
mat = HNSCC_ref_subset@assays$RNA@counts
mat <- mat[! rownames(mat) %in% mitochondrial_genes,]
mat <- mat[! rownames(mat) %in% ribosomal_genes,]

giotto_SC <- createGiottoObject(
  expression = mat)
giotto_SC <- addCellMetadata(giotto_SC,
                             new_metadata = HNSCC_ref_subset@meta.data)
giotto_SC<- normalizeGiotto(giotto_SC)
giotto_SC@cell_metadata$cell$rna@metaDT$Labels

markers_scran = findMarkers_one_vs_all(gobject=giotto_SC, method="scran",
                                       expression_values="normalized", cluster_column = "Labels", min_feats=3)
top_markers <- markers_scran[, head(.SD, 100), by="cluster"]
sc_expression_norm = getExpression(giotto_SC,
                                   values = "normalized",
                                   output = "matrix")
sc_expression_norm <- sc_expression_norm[! rownames(sc_expression_norm) %in% mitochondrial_genes,]
sc_expression_norm <- sc_expression_norm[! rownames(sc_expression_norm) %in% ribosomal_genes,]
# Create DWLS matrix
DWLS_matrix<-makeSignMatrixDWLSfromMatrix(matrix = sc_expression_norm,
                                          cell_type = pDataDT(giotto_SC)$Labels,
                                          sign_gene = top_markers$feats)

Visium_giotto = runDWLSDeconv(gobject = Visium_giotto, sign_matrix = DWLS_matrix, n_cell = 20)

DWLS <- as.data.frame(Visium_giotto@spatial_enrichment[["cell"]][["rna"]][["DWLS"]]@enrichDT)
rownames(DWLS) <- DWLS$cell_ID
DWLS$cell_ID = NULL


#### Fibroblast containing spots ####
Deconv_df## dataframe containing Deconv_output (RCTD/SpatialWDLS) for all Visium samples
Fibroblast_containing_spots <- Deconv_df %>%
  filter(iCAF >= 0.05 | myCAF >= 0.05 | FRC_like >= 0.05 | PI16 >= 0.05)

filtered_Deconv_df <- Deconv_df %>%
  mutate_all(~ ifelse(. < 0.05, 0, .))#implementing 0.05 threshold for spot cell type presence (assumed maximum of 20 cells per spot: 1/20)

#### Patient-cell-type Cor ####
compute_cor_p <- function(df, x_col, y_col) {
  cor_result <- cor.test(df[[x_col]], df[[y_col]], method = "spearman")
  return(data.frame(
    x_col = x_col,
    y_col = y_col,
    cor_coef = cor_result$estimate,
    p_value = cor_result$p.value
  ))
}

# List of column names - corresponding to cell types
column_names <- colnames(filtered_Deconv_df)[13:32]
cor_df <- data.frame()
sample_ids <- c("HN481", "HN482", "HN483", "HN485", "HN487", "HN488", "HN489", "HN490", "HN492", "HN494")
Fibro_OI <- "myCAF" # fibroblast column name of interest e.g., myCAF, iCAF, FRC_like, PI16

# Loop through each sample ID
for (sample_id in sample_ids) {
  filtered_data <- filtered_Deconv_df %>% filter(Source.SampleID == sample_id)
    for (col_name in column_names) {
    cor_result <- compute_cor_p(filtered_data, Fibro_OI, col_name)
    cor_result$sample_id <- sample_id  # Add sample_id as a new column
    cor_df <- rbind(cor_df, cor_result)
  }
}
cor_df$y_col <- with(cor_df, reorder(y_col, cor_coef, mean, na.rm = TRUE))

#### Weighted fishers test ####
cor_df <- cor_df %>%filter(!p_value == "")

cell_types <- levels(cor_df$y_col)

# Function to calculate weighted z-score Fisher's p-value
calculate_meta_p_value <- function(cor_df, cell_type) {
  # Filter dataframe to remove rows with missing p-values and for the specific cell type
  cor_df_filtered <- cor_df %>%
    filter(!is.na(p_value) & p_value != "", y_col == cell_type)
  
  if (nrow(cor_df_filtered) == 0) {
    return(NA)  
  }
    cor_df_filtered$z_score <- qnorm(1 - cor_df_filtered$p_value / 2)
    weights <- abs(cor_df_filtered$cor_coef)
    fisher_combined <- -2 * sum(weights * log(cor_df_filtered$p_value))
   total_weight <- sum(weights)
    df <- 2 * total_weight
  meta_p_value <- 1 - pchisq(fisher_combined, df)
  return(meta_p_value)
}

# Calculate meta p-value for each cell type and store the output in a dataframe
meta_p_values <- data.frame(
  cell_type = cell_types,
  meta_p_value = unlist(lapply(cell_types, function(cell_type) calculate_meta_p_value(cor_df, cell_type)))
)

