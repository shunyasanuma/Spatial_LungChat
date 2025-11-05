## 4. Spatial scTriangulate 
### Cross-Platform Cell Type Refinement
This notebook [`Spatial_scTriangulate.ipynb`](./scTriangulate_v9.ipynb) performs cell type annotation refinement across Visium HD and Xenium spatial transcriptomics platforms using `scTriangulate`. It follows up on a prior alignment step and evaluates consistency across multiple annotation sources.

---
#### Workflow Steps
#### 1. Load Aligned Expression Data
- Loads the combined expression matrix generated from the `Alignment.Rmd` step (matched Visium–Xenium spots, with `_visium` and `_xenium` gene suffixes).
- Parses cell annotations from:
- RCTD (Visium + Xenium)
  * `Visium_RCTD_LungMap_ref`
  * `Visium_RCTD_GSE250346_based`
  * `Xenium_RCTD_LungMap_ref`
- Manuscript labels
  * `Final_CT`
  * `CNiche`
  * `TNiche`

#### 2.  Build `scTriangulate` Object
- Initializes a triangulate object using:
  * Combined expression matrix
  * Multi-source annotations (RCTD, manuscript labels)
 
#### 3. Visualize Results
<img src="./figures/Spatial_scTriangulate_pruned.png" alt="Overview" width="1000"/>


