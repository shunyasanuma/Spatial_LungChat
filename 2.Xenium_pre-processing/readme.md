## Xenium_pre-processing
Here we describe a series of analyses of Xenium cell type annotation.

### Outline
[2.1. Subset](#1-subset-python)


---
#### 1. Subset (Python)
This Jupyter Notebook, [`Crop_high_resolution_DAPI_image.ipynb`](./Xenium_Crop_high_resolution_DAPI.ipynb), provides a pipeline to visualize, overlay grid lines, and crop high-resolution DAPI images from Xenium spatial transcriptomics output. It is especially useful for inspecting and extracting specific tissue regions for downstream segmentation or registration workflows.
- Load and visualize high-resolution DAPI `.tif` images.
- Crop regions of interest (ROIs) based on coordinate inputs.
- Save the cropped output for downstream spatial analysis (e.g., Cell type annotation, Visium HD - Xenium alignment).

#### Required input files
- High resolution DAPI tiff. 
- DAPI image from Xenium is typically located under: `morphology_focus/morphology_focus_000*.ome.tif`
