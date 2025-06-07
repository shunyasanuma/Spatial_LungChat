## Visium_HD_pre-processing
Here we describe a series of analyses of Visium HD cell type annotation.

### Outline
#### 1. Subset (Python)
This Jupyter Notebook, [`Crop_high_resolution_microsope_image.ipynb`](./Visium_HD_Crop_high_resolution_microsope_image.ipynb), demonstrates how to:
- Load and visualize high-resolution microscope `.tif` images.
- Crop regions of interest (ROIs) based on coordinate inputs.
- Save the cropped output for downstream spatial analysis (e.g., Cell type annotation, Visium HD - Xenium alignment).

##### Required input files
- High resolution H&E tiff


#### 2. Robust Cell Type Decomposition, Cell type annotation (R)
This script processes Visium HD binned spatial transcriptomics data from human lung tissue to perform spatial deconvolution, clustering, coordinate normalization, and JSON export for downstream applications like visualization or LLM prompting.
