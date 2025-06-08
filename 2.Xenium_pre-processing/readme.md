## Xenium pre-processing
Here we describe a series of analyses of Xenium cell type annotation.

### Outline
[2.1. Subset](#1-subset-python)

[2.2. Robust Cell Type Decomposition](#2-robust-cell-type-decomposition-cell-type-annotation-r)


---
#### 1. Subset (Python)
This Jupyter Notebook, [`Crop_high_resolution_DAPI_image.ipynb`](./Xenium_Crop_high_resolution_DAPI.ipynb), provides a pipeline to visualize, overlay grid lines, and crop high-resolution DAPI images from Xenium spatial transcriptomics output. It is especially useful for inspecting and extracting specific tissue regions for downstream segmentation or registration workflows.
- Load and visualize high-resolution DAPI `.tif` images.
- Crop regions of interest (ROIs) based on coordinate inputs.
- Save the cropped output for downstream spatial analysis (e.g., Cell type annotation, Visium HD - Xenium alignment).

#### Required input files
- High resolution DAPI tiff. 
- DAPI image from Xenium is typically located under: `morphology_focus/morphology_focus_000*.ome.tif`
---

#### 2. Robust Cell Type Decomposition, Cell type annotation (R)
This R script [`Spatial_Deconvolution_Clustering_Pipeline.R`](./Tissue2_celltype_annotation.R) performs cell type deconvolution, spatial subsetting, unsupervised clustering, and metadata annotation for Xenium spatial transcriptomics data using spacexr (RCTD) and Seurat.

##### Workflow Steps
1. Load & Preprocess Xenium Counts
- Load sparse matrix (`matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`)
- Map ENSEMBL IDs to gene symbols
- Remove duplicated or missing gene entries

2. Load & Crop Cell Boundaries
- Read `cell_boundaries.parquet`
- Subset a region of interest (ROI) using pixel-based crop bounds
- Convert pixel to micron scale using `microns_per_pixel = 0.2125`

3. Build RCTD Input
- Load cell centroids from `cells.parquet`
- Subset centroids to match cropped ROI cells
- Create SpatialRNA object for RCTD
- Load LungMAP CellRef v1.1 (or your own scRNA-seq reference)
- Remove rare cell types (<25 cells)

4. Run RCTD
- Run in `doublet` mode
- Thresholds: `UMI_min = 40`, `counts_MIN = 10`

5. Unsupervised Clustering (Seurat)
- Convert `RCTD` output to Seurat
- Normalize → PCA → Leiden clustering (`res = 0.782`) → UMAP

6. Normalize Coordinates
- No rotation applied (already aligned)
- Normalize x and y to [0, 1] range for downstream use

7. Update RCTD Object
- Replace spatial coordinates with normalized ones
- Add Seurat clusters and UMAP dimensions to `@results$results_df`

8. Merge with Author Annotations
- Add `final_lineage`, `CNiche`, `TNiche`, `final_CT` from `Xenium_author.rds` (from Vannan et al Nature Genetics manuscript)

9. Save Processed Object
- Save as `Xenium2.RData` for downstream use

10. Export Metadata (JSON)
Outputs:
- Spot class counts (singlet, doublet, reject)
- Cell type labels
- Seurat clusters
- Normalized coordinate ranges
- Plotting suggestions

Required Inputs
| Input                        | Description                                                                  |
| ---------------------------- | ---------------------------------------------------------------------------- |
| `matrix.mtx.gz`              | Sparse gene expression matrix from Xenium                                    |
| `barcodes.tsv.gz`            | Cell barcodes                                                                |
| `features.tsv.gz`            | Gene metadata (ENSEMBL ID and symbols)                                       |
| `cell_boundaries.parquet`    | Polygon-based boundaries for each cell                                       |
| `cells.parquet`              | Cell centroid coordinates                                                    |
| `scRNA-seq reference (.rds)` | Preprocessed reference (e.g., LungMAP) with celltype\_level1 annotations     |
| `ROI pixel boundaries`       | Region of interest in Xenium pixel space (converted to microns for cropping) |
| `Xenium_author.rds`          | Author-derived Seurat object with curated annotations                        |

Output Files
| Output                                | Format                          | Description                                                                |
| ------------------------------------- | ------------------------------- | -------------------------------------------------------------------------- |
| `Xenium2.RData`                       | RData                           | Final RCTD object with integrated Seurat clustering and author annotations |
| `Xenium2.json`                        | JSON                            | Summary metadata for LLM visualization or downstream tools                 |
| `seurat_clusters`, `UMAP_1`, `UMAP_2` | In `Xenium2@results$results_df` | Clustering and UMAP coordinates                                            |
| `SpatialRNA object`                   | In `Xenium2@spatialRNA`         | Coordinates normalized to \[0, 1] scale                                    |





























