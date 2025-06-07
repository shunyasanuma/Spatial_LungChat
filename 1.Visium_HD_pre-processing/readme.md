## Visium_HD_pre-processing
Here we describe a series of analyses of Visium HD cell type annotation.

### Outline
#### 1. Subset (Python)
This Jupyter Notebook, [`Crop_high_resolution_microsope_image.ipynb`](./Visium_HD_Crop_high_resolution_microsope_image.ipynb), demonstrates how to:
- Load and visualize high-resolution microscope `.tif` images.
- Crop regions of interest (ROIs) based on coordinate inputs.
- Save the cropped output for downstream spatial analysis (e.g., Cell type annotation, Visium HD - Xenium alignment).

#### Required input files
- High resolution H&E tiff


#### 2. Robust Cell Type Decomposition, Cell type annotation (R)
This Rscript, [`Spatial_Deconvolution_Clustering_Pipeline.R`](./Tissue2_celltype_annotation.R) processes Visium HD binned spatial transcriptomics data from human lung tissue to perform spatial deconvolution, clustering, coordinate normalization, and JSON export for downstream applications like visualization or LLM prompting.

##### Workflow Steps
1. Load & Preprocess Visium HD Counts
- Load filtered matrix (Space Ranger 8µm output)
- Map ENSEMBL IDs to gene symbols
- Remove duplicated or missing gene entries

2. Load Spatial Coordinates
-  Read `tissue_positions.parquet`
-  Subset a region of interest (ROI) using pixel bounds

3. Build RCTD Input
- Create `SpatialRNA` object
- Load LungMAP scRNA-seq reference (You can use your reference)
- Filter out cell types with low abundance (<25 cells)

4. Run RCTD
- Run in `doublet` mode
- Thresholds: `UMI_min = 40`, `counts_MIN = 10`


