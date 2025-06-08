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
# Load & Preprocess Xenium Counts
# -------------------------------
### PLACE TO CHANGE ###
counts <- readMM("/data/salomonis2/LabFiles/Shunya_Asanuma/Spatial/LungChat/Input/GSE250346/IPFTMA5/cell_feature_matrix/matrix.mtx.gz")
barcodes <- readLines("/data/salomonis2/LabFiles/Shunya_Asanuma/Spatial/LungChat/Input/GSE250346/IPFTMA5/cell_feature_matrix/barcodes.tsv.gz")
genes <- read.delim(
  "/data/salomonis2/LabFiles/Shunya_Asanuma/Spatial/LungChat/Input/GSE250346/IPFTMA5/cell_feature_matrix/features.tsv.gz",
  header = FALSE, sep = "\t"
)

colnames(counts) <- barcodes
rownames(counts) <- genes$V1
ensembl_to_symbol <- setNames(as.character(genes$V2), genes$V1)
gene_symbols <- ensembl_to_symbol[rownames(counts)]
keep <- !is.na(gene_symbols) & !duplicated(gene_symbols)
counts <- counts[keep, ]
rownames(counts) <- gene_symbols[keep]

# -------------------------------
# Load & Crop Cell Coordinates
# -------------------------------
### PLACE TO CHANGE ###
coords <- read_parquet("/data/salomonis2/LabFiles/Shunya_Asanuma/Spatial/LungChat/Input/GSE250346/IPFTMA5/cell_boundaries.parquet")
coords <- as.data.frame(coords[, c("cell_id", "vertex_x", "vertex_y")])
x1_pix <- 37500; x2_pix <- 54500
y1_pix <- 35900; y2_pix <- 51000
microns_per_pixel <- 0.2125


xmin <- x1_pix * microns_per_pixel; xmax <- x2_pix * microns_per_pixel
ymin <- y1_pix * microns_per_pixel; ymax <- y2_pix * microns_per_pixel
tissue_region <- coords[
  coords$vertex_x >= xmin & coords$vertex_x <= xmax &
    coords$vertex_y >= ymin & coords$vertex_y <= ymax, ]
crop_cell_ids <- unique(tissue_region$cell_id)

# -------------------------------
# Load Centroids & Subset Counts
# -------------------------------
### PLACE TO CHANGE ###
coords_centroid_all <- read_parquet("/data/salomonis2/LabFiles/Shunya_Asanuma/Spatial/LungChat/Input/GSE250346/IPFTMA5/cells.parquet")

coords_centroid_all <- as.data.frame(coords_centroid_all)
rownames(coords_centroid_all) <- coords_centroid_all$cell_id
coords_centroid_crop <- coords_centroid_all[crop_cell_ids, c("x_centroid", "y_centroid")]
colnames(coords_centroid_crop) <- c("x", "y")
common_ids <- intersect(rownames(coords_centroid_crop), colnames(counts))
coords_centroid_crop <- coords_centroid_crop[common_ids, , drop = FALSE]
counts_crop <- counts[, common_ids]
query <- SpatialRNA(coords_centroid_crop, counts_crop, colSums(counts_crop))

# -------------------------------
# Reference
# -------------------------------
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
RCTD <- create.RCTD(query, reference, UMI_min = 40, counts_MIN = 10)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
Xenium2 <- RCTD
# save(Xenium2, file = "/data/salomonis2/LabFiles/Shunya_Asanuma/Spatial/LungChat/Output/Xenium/GSE250346/02.RCTD/UMI_40/Tissue2_UMI40.RData")

# -------------------------------
# Unsupervised Clustering (Seurat)
# -------------------------------
# load("/data/salomonis2/LabFiles/Shunya_Asanuma/Spatial/LungChat/Output/Xenium/GSE250346/02.RCTD/UMI_40/Tissue2_UMI40.RData")
counts <- Xenium2@spatialRNA@counts
coords <- Xenium2@spatialRNA@coords

seurat_obj <- CreateSeuratObject(counts = counts)
seurat_obj@meta.data$x <- coords[rownames(seurat_obj@meta.data), "x"]
seurat_obj@meta.data$y <- coords[rownames(seurat_obj@meta.data), "y"]

seurat_obj <- NormalizeData(seurat_obj) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters(resolution = 0.782) %>%
  RunUMAP(dims = 1:10)

# -------------------------------
# Normalize & Rotate Coordinates
# -------------------------------
df_plot <- seurat_obj@meta.data %>%
  mutate(
    UMAP_1 = Embeddings(seurat_obj, "umap")[, 1],
    UMAP_2 = Embeddings(seurat_obj, "umap")[, 2],
    x_rotated = x,
    y_rotated = y,
    x_rotated_norm = (x_rotated - min(x_rotated)) / (max(x_rotated) - min(x_rotated)),
    y_rotated_norm = (y_rotated - min(y_rotated)) / (max(y_rotated) - min(y_rotated))
  )

# -------------------------------
# Visualization (Optional)
# -------------------------------
ggplot(df_plot, aes(x = x_rotated_norm, y = y_rotated_norm, color = seurat_clusters)) +
  geom_point(size = 0.4) +
  coord_fixed() +
  theme_void() +
  labs(title = "Xenium (Leiden Clusters")

# -------------------------------
# Update RCTD Object with Clusters & UMAP
# -------------------------------
Xenium2@spatialRNA@coords <- df_plot[, c("x_rotated_norm", "y_rotated_norm")]
colnames(Xenium2@spatialRNA@coords) <- c("x", "y")

stopifnot(identical(rownames(seurat_obj@meta.data), rownames(Xenium2@results$results_df)))
Xenium2@results$results_df$seurat_clusters <- seurat_obj$seurat_clusters
Xenium2@results$results_df$UMAP_1 <- df_plot$UMAP_1
Xenium2@results$results_df$UMAP_2 <- df_plot$UMAP_2

# -------------------------------
# Add Author Annotations from Xenium_author
# -------------------------------
### PLACE TO CHANGE ###
# Assume Xenium_author is loaded and has meta.data$cell_id
Xenium_author <- readRDS("/data/salomonis2/LabFiles/Shunya_Asanuma/Spatial/LungChat/Input/GSE250346/GSE250346_Seurat_GSE250346_CORRECTED_SEE_RDS_README_082024.rds")

seurat_meta <- Xenium_author@meta.data
cell_id_to_row <- setNames(rownames(seurat_meta), seurat_meta$cell_id)
cells_in_rctd <- rownames(Xenium2@results$results_df)
matched_cells <- intersect(cells_in_rctd, names(cell_id_to_row))
author_rows <- cell_id_to_row[matched_cells]
annotations <- seurat_meta[author_rows, c("final_lineage", "CNiche", "TNiche", "final_CT")]
rownames(annotations) <- matched_cells
Xenium2@results$results_df[matched_cells, c("final_lineage", "CNiche", "TNiche", "final_CT")] <- annotations

# -------------------------------
# Save Final Xenium2
# -------------------------------
save(Xenium2, file = "/data/salomonis2/LabFiles/Shunya_Asanuma/Spatial/LungChat/Output/Xenium/GSE250346/02.RCTD/UMI_40/Xenium2.RData")

# -------------------------------
# Creating JSON file 
# -------------------------------
# --- Extract metadata ---
coords <- Xenium2@spatialRNA@coords
n_spots <- nrow(coords)
spot_class <- Xenium2@results$results_df$spot_class
first_type <- Xenium2@results$results_df$first_type
seurat_clusters <- Xenium2@results$results_df$seurat_clusters

# --- Unique values and counts ---
spot_class_table <- as.list(table(spot_class))
first_type_list <- sort(unique(na.omit(first_type)))
cluster_list <- sort(unique(na.omit(seurat_clusters)))

# --- Build JSON structure ---
llm_metadata <- list(
  dataset_name = "Xenium2",
  description = paste0("Xenium crop of IPF lung. Contains ", n_spots, " spatial spots."),
  species = "Homo sapiens",
  platform = "Xenium",
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
jsonlite::write_json(llm_metadata, "/data/salomonis2/LabFiles/Shunya_Asanuma/Spatial/LungChat/Output/Xenium/GSE250346/02.RCTD/UMI_40/JSON/Xenium2.json", pretty = TRUE, auto_unbox = TRUE)


