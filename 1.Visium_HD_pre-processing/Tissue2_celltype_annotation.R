# -------------------------------
# Libraries
# -------------------------------
library(SummarizedExperiment)
library(MatrixGenerics)
library(SingleCellExperiment)
library(Matrix)
library(Seurat)
library(ggplot2)
library(dplyr)
library(arrow)
library(hdf5r)
library(spacexr)
library(SeuratWrappers)
library(mclust)
library(zellkonverter)
library(SeuratDisk)

# -------------------------------
# Load & Preprocess Visium Counts
# -------------------------------

### PLACE TO CHANGE ###
counts <- readMM("/data/salomonis2/LabFiles/Shunya_Asanuma/Spatial/LungChat/Input/GSM8509590/highres_BANOSSM_SSM0003_1_LU_Whole_C1_VISHD_F07528_22MM3MLT3/outs/binned_outputs/square_008um/filtered_feature_bc_matrix/matrix.mtx.gz")
barcodes <- readLines("/data/salomonis2/LabFiles/Shunya_Asanuma/Spatial/LungChat/Input/GSM8509590/highres_BANOSSM_SSM0003_1_LU_Whole_C1_VISHD_F07528_22MM3MLT3/outs/binned_outputs/square_008um/filtered_feature_bc_matrix/barcodes.tsv.gz")
genes <- read.delim(
  "/data/salomonis2/LabFiles/Shunya_Asanuma/Spatial/LungChat/Input/GSM8509590/highres_BANOSSM_SSM0003_1_LU_Whole_C1_VISHD_F07528_22MM3MLT3/outs/binned_outputs/square_008um/filtered_feature_bc_matrix/features.tsv.gz",
  header = FALSE,
  sep = "\t")


colnames(counts) <- barcodes
rownames(counts) <- genes$V1
# Map ENSEMBL IDs to symbols
ensembl_to_symbol <- setNames(as.character(genes$V2), genes$V1)
gene_symbols <- ensembl_to_symbol[rownames(counts)]
keep <- !is.na(gene_symbols) & !duplicated(gene_symbols)
counts <- counts[keep, ]
rownames(counts) <- gene_symbols[keep]

# -------------------------------
# Load Coordinates & Subset Region
# -------------------------------
### PLACE TO CHANGE ###
coords <- read_parquet("/data/salomonis2/LabFiles/Shunya_Asanuma/Spatial/LungChat/Input/GSM8509590/highres_BANOSSM_SSM0003_1_LU_Whole_C1_VISHD_F07528_22MM3MLT3/outs/binned_outputs/square_008um/spatial/tissue_positions.parquet")
coords <- as.data.frame(coords)[, c("barcode", "pxl_col_in_fullres", "pxl_row_in_fullres")]

rownames(coords) <- coords$barcode
coords$barcode <- NULL
coords <- coords[rownames(coords) %in% colnames(counts), , drop = FALSE]

### PLACE TO CHANGE ###
tissue_2 <- coords[
  coords$pxl_col_in_fullres >= 18700 & coords$pxl_col_in_fullres <= 25000 &
    coords$pxl_row_in_fullres >= 16300 & coords$pxl_row_in_fullres <= 23500,
]
counts_2 <- counts[, rownames(tissue_2)]

# -------------------------------
# Build SpatialRNA & Reference
# -------------------------------
query_2 <- SpatialRNA(tissue_2, counts_2, colSums(counts_2))

### PLACE TO CHANGE ###
ref <- readRDS("/data/salomonis2/CCHMC-Collaborations/Minzhe/seurat4_refmap_cellrefseed/LungMAP_HumanLung_CellRef_Seed.v1.1.rds")
Idents(ref) <- "celltype_level1"

ref_filtered <- subset(ref, idents = names(table(Idents(ref)))[table(Idents(ref)) >= 25])
ref_filtered$celltype_level1 <- factor(gsub("/", "_", ref_filtered$celltype_level1))
reference <- Reference(
  counts = ref_filtered[["RNA"]]@counts,
  cell_types = ref_filtered$celltype_level1,
  nUMI = ref_filtered$nCount_RNA
)

# -------------------------------
# Run RCTD
# -------------------------------
gc()
RCTD <- create.RCTD(query_2, reference, UMI_min = 40, counts_MIN = 10)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

# -------------------------------
# Unsupervised Clustering
# -------------------------------
# load(".../Tissue2.RData")
### PLACE TO CHANGE ###
Visium2 <- RCTD

counts <- Visium2@spatialRNA@counts
coords <- Visium2@spatialRNA@coords
seurat_obj <- CreateSeuratObject(counts = counts)
seurat_obj@meta.data$x <- coords[rownames(seurat_obj@meta.data), "x"]
seurat_obj@meta.data$y <- coords[rownames(seurat_obj@meta.data), "y"]

### PLACE TO CHANGE ###
seurat_obj <- NormalizeData(seurat_obj) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters(resolution = 0.17) %>%
  RunUMAP(dims = 1:10)

# -------------------------------
# Normalize & Rotate Coordinates
# -------------------------------
df_plot <- seurat_obj@meta.data %>%
  mutate(
    X_rotated = y,          
    Y_rotated = -x,         # rotate 90° clockwise by negating x
    X_rotated_norm = (X_rotated - min(X_rotated)) / (max(X_rotated) - min(X_rotated)),
    Y_rotated_norm = (Y_rotated - min(Y_rotated)) / (max(Y_rotated) - min(Y_rotated))
  )


# Optional: Visualization
ggplot(df_plot, aes(x = X_rotated_norm, y = Y_rotated_norm, color = seurat_clusters)) +
  geom_point(size = 0.4) +
  coord_fixed() +
  theme_void() +
  labs(title = "Visium HD (Leiden Clusters, Corrected Orientation)")

# -------------------------------
# Update RCTD Object
# -------------------------------
Visium2@spatialRNA@coords <- df_plot[, c("X_rotated_norm", "Y_rotated_norm")]
colnames(Visium2@spatialRNA@coords) <- c("x", "y")

stopifnot(identical(rownames(seurat_obj@meta.data), rownames(Visium2@results$results_df)))
Visium2@results$results_df$seurat_clusters <- seurat_obj$seurat_clusters

umap_coords <- Embeddings(seurat_obj, "umap")
Visium2@results$results_df$UMAP_1 <- umap_coords[rownames(Visium2@results$results_df), 1]
Visium2@results$results_df$UMAP_2 <- umap_coords[rownames(Visium2@results$results_df), 2]

# -------------------------------
# Save
# -------------------------------
### PLACE TO CHANGE ###
save(Visium2, file = "/data/salomonis2/LabFiles/Shunya_Asanuma/Spatial/LungChat/Output/Visium_HD/GSM8509590/02.RCTD/UMI_40/Visium2.RData")
# load("/data/salomonis2/LabFiles/Shunya_Asanuma/Spatial/LungChat/Prompt_engineering/GSM8509590_GSE250346/Datasets/Tissue2/RCTD/Visium2.RData")

# -------------------------------
# Creating JSON file 
# -------------------------------
# --- Extract metadata ---
coords <- Visium2@spatialRNA@coords
n_spots <- nrow(coords)
spot_class <- Visium2@results$results_df$spot_class
first_type <- Visium2@results$results_df$first_type
seurat_clusters <- Visium2@results$results_df$seurat_clusters

# --- Unique values and counts ---
spot_class_table <- as.list(table(spot_class))
first_type_list <- sort(unique(na.omit(first_type)))
cluster_list <- sort(unique(na.omit(seurat_clusters)))

# --- Build JSON structure ---
llm_metadata <- list(
  dataset_name = "Visium2",
  description = paste0("Visium HD crop of IPF lung. Contains ", n_spots, " spatial spots."),
  species = "Homo sapiens",
  platform = "Visium HD",
  annotations = list(
    spot_class = spot_class_table,
    cell_types = first_type_list,
    seurat_clusters = cluster_list
  ),
  coordinates = list(
    x_range = range(coords[, "x"]),
    y_range = range(coords[, "y"]),
    system = "Visium raw pixel coordinates"
  ),
  available_annotations = c("first_type", "spot_class", "seurat_clusters"),
  plotting_suggestions = list(
    default_annotation = "first_type",
    color_by = c("first_type", "spot_class", "seurat_clusters"),
    spatial_axes = c("x", "y")
  )
)

# --- Export to JSON ---
jsonlite::write_json(llm_metadata, "/data/salomonis2/LabFiles/Shunya_Asanuma/Spatial/LungChat/Output/Visium_HD/GSM8509590/02.RCTD/UMI_40/JSON/Visium2.json", pretty = TRUE, auto_unbox = TRUE)

