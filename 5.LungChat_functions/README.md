## 5. LungChat (LLM) functions

### Outline
[5.1. Visualization/Analysis](#1-Visualization-&-Analysis)

[5.2. CellChat](#2-CellChat)

[5.3. Neighbors enrichment_analysis](#3-neighbors-enrichment_analysis)

[5.4. Ingestion](#4-Ingestion)

---
#### 1. Visualization & Analysis
Script: [`Spatial_functions.R`](./Spatial_functions_v4.R)

#### `plot_spatial()`
- Visualize spatial annotation (categorical or numeric) from a SummarizedExperiment object.

| **Parameter**       | **Type**                                        | **Default**  | **Options**                                                              | **Description**                                                                 |
| ------------------- | ----------------------------------------------- | ------------ | ------------------------------------------------------------------------ | ------------------------------------------------------------------------------- |
| `obj`               | `SingleCellExperiment` / `SummarizedExperiment` | *(required)* | —                                                                        | H5AD-derived object loaded using `zellkonverter::readH5AD()`                    |
| `platform`          | character                                       | `"Xenium"`   | `"Xenium"`, `"Visium"`                                                   | Used to select the appropriate UMAP coordinate columns                          |
| `annotation`        | character                                       | `"Final_CT"` | Any column in `colData(obj)`, e.g., `"Final_CT"`, `"pruned"`, `"TNiche"` | Column in `colData(obj)` used to color the UMAP plot                            |
| `title`             | character                                       | `NULL`       | —                                                                        | Custom plot title; if `NULL`, uses `"<platform> UMAP: <annotation>"`            |
| `point_size`        | numeric                                         | `0.5`        | Any positive number                                                      | Size of points in the UMAP scatterplot                                          |
| `legend_point_size` | numeric                                         | `4`          | Any positive number                                                      | Size of points in the legend                                                    |
| `highlight`         | character vector                                | `NULL`       | Subset of annotation levels, e.g., `c("AT1", "AT2")`                     | If set and `annotation` is a factor, highlights selected levels, others in grey |

[`The list of available annotations is here.`](./annotation_labels.md)

#### Example usages
Generates a spatial scatterplot of scTriangulate2 cells or spots, colored by the "pruned" annotation to visualize cell type or cluster identity across tissue space.
```
plot_spatial(scTriangulate2, annotation = "pruned")
```
<img src="./figures/plot_spatial_1.png" alt="Example" width="800"/>

Generate a spatial plot highlighting `TNiche@T6` and `TNiche@T7` clusters
```
plot_spatial(scTriangulate2, 
             annotation = "pruned", 
             highlight = c("TNiche@T6", "TNiche@T7"))
```
<img src="./figures/plot_spatial_2.png" alt="Example" width="400"/>

Generates a spatial map of the IPF microenvironment where each spot is colored by confidence score
```
plot_spatial(
    obj = scTriangulate2,
    annotation = "confidence",
    title = "Spatial Confidence Map of the IPF Microenvironment",
    point_size = 0.4
)
```
<img src="./figures/plot_spatial_3.png" alt="Example" width="500"/>

---
#### `plot_umap()`
- Visualize UMAP embedding for Visium HD or Xenium platforms using metadata or highlights.

| **Parameter**       | **Type**                                        | **Default**  | **Options**                                                              | **Description**                                                                 |
| ------------------- | ----------------------------------------------- | ------------ | ------------------------------------------------------------------------ | ------------------------------------------------------------------------------- |
| `obj`               | `SingleCellExperiment` / `SummarizedExperiment` | *(required)* | —                                                                        | H5AD-derived object loaded using `zellkonverter::readH5AD()`                    |
| `platform`          | character                                       | `"Xenium"`   | `"Xenium"`, `"Visium"`                                                   | Used to select the appropriate UMAP coordinate columns                          |
| `annotation`        | character                                       | `"Final_CT"` | Any column in `colData(obj)`, e.g., `"Final_CT"`, `"pruned"`, `"TNiche"` | Column in `colData(obj)` used to color the UMAP plot                            |
| `title`             | character                                       | `NULL`       | —                                                                        | Custom plot title; if `NULL`, uses `"<platform> UMAP: <annotation>"`            |
| `point_size`        | numeric                                         | `0.5`        | Any positive number                                                      | Size of points in the UMAP scatterplot                                          |
| `legend_point_size` | numeric                                         | `4`          | Any positive number                                                      | Size of points in the legend                                                    |
| `highlight`         | character vector                                | `NULL`       | Subset of annotation levels, e.g., `c("AT1", "AT2")`                     | If set and `annotation` is a factor, highlights selected levels, others in grey |

[`The list of available annotations is here.`](./annotation_labels.md)

#### Example usages
Xenium UMAP
```
plot_umap(scTriangulate2, 
          platform = "Xenium", 
          annotation = "Final_CT")
```
Visium HD UMAP
```
plot_umap(scTriangulate2, 
          platform = "Visium", 
          annotation = "Visium_RCTD_LungMap_ref")
```
<img src="./figures/plot_umap_1.png" alt="Example" width="1000"/>

Xenium (AT1 & AT2)
```
plot_umap(scTriangulate2, 
          platform = "Visium", 
          annotation = "Visium_RCTD_LungMap_ref", 
          title = "Visium - AT1 & AT2",
          highlight = c("AT1", "AT2"))
```
Visium HD (AT1 & AT2)
```
plot_umap(scTriangulate2, 
          platform = "Xenium", 
          annotation = "Final_CT", 
          title = "Xenium - AT1 & AT2",
          highlight = c("AT1", "AT2"))
```
<img src="./figures/plot_umap_2.png" alt="Example" width="800"/>

---
#### `plot_gene_spatial()`
- Plot spatial gene expression for Visium or Xenium using platform-specific coordinates.

| **Parameter** | **Type**                                        | **Default**                   | **Options**                                                | **Description**                                                           |
| ------------- | ----------------------------------------------- | ----------------------------- | ---------------------------------------------------------- | ------------------------------------------------------------------------- |
| `obj`         | `SingleCellExperiment` / `SummarizedExperiment` | *(required)*                  | —                                                          | H5AD-derived object loaded using `zellkonverter::readH5AD()`              |
| `gene`        | character                                       | *(required)*                  | Must match a row in `rownames(assay(obj))`, e.g., `"AGER"` | Gene to visualize (suffix added based on platform)                        |
| `platform`    | character                                       | `"Xenium"`                    | `"Xenium"`, `"Visium"`                                     | Platform used to determine spatial columns and gene suffix                |
| `x_col`       | character                                       | `"X_xenium_norm"` (if Xenium) | e.g., `"X_visium_norm"`, `"X_xenium_norm"`                 | Column in `colData(obj)` for X coordinates (auto-set by platform if NULL) |
| `y_col`       | character                                       | `"Y_xenium_norm"` (if Xenium) | e.g., `"Y_visium_norm"`, `"Y_xenium_norm"`                 | Column in `colData(obj)` for Y coordinates (auto-set by platform if NULL) |
| `title`       | character                                       | `NULL`                        | —                                                          | Custom plot title; if NULL, defaults to `"Spatial Expression of <gene>"`  |
| `point_size`  | numeric                                         | `0.5`                         | Any positive number                                        | Size of points on the scatterplot                                         |
| `low_color`   | character                                       | `"grey90"`                    | Any valid R color                                          | Color used for low expression                                             |
| `high_color`  | character                                       | `"red"`                       | Any valid R color                                          | Color used for high expression                                            |

[`The list of available genes is here.`](./available_genes_by_platform.md)

#### Example usages
```
plot_gene_spatial(scTriangulate2, 
                  gene = "COL1A1", 
                  platform = "Visium", 
                  point_size = 0.4,
                  title = "Visium - Expression of COL1A1")
                  
plot_gene_spatial(scTriangulate2, 
                  gene = "COL1A1", 
                  platform = "Xenium", 
                  point_size = 0.4,
                  title = "Xenium - Expression of COL1A1")
```
<img src="./figures/plot_gene_spatial_1.png" alt="Example" width="800"/>

---
#### `rank_cross_platform_gene_correlation()`
- Rank genes by Pearson correlation across Xenium-Visium matched expression.

| **Parameter** | **Type**                                        | **Default**  | **Options**              | **Description**                                                               |
| ------------- | ----------------------------------------------- | ------------ | ------------------------ | ----------------------------------------------------------------------------- |
| `obj`         | `SingleCellExperiment` / `SummarizedExperiment` | *(required)* | —                        | Object containing gene expression from both Xenium and Visium platforms       |
| `top_genes`   | numeric                                         | `5`          | Any positive integer     | Number of top correlated gene pairs to return and visualize                   |
| `layer`       | character                                       | `"X"`        | Must match an assay name | Assay layer (e.g., `"X"`, `"counts"`, `"logcounts"`) used for expression data |

#### Example usages
Top 5 highly correlated genes (Visium HD - Xenium)
```
rank_cross_platform_gene_correlation(
    obj = scTriangulate2,
    top_genes = 5
)
```
Output
```
     Gene Correlation
1  COL1A1   0.8968483
2  COL3A1   0.6852229
3  COL1A2   0.6121173
4 SCGB3A2   0.5324030
5   MARCO   0.4863397
```
<img src="./figures/rank_cross_platform_gene_correlation_1.png" alt="Example" width="500"/>

---
#### `plot_cross_platform_gene_correlation()`
- Plot cumulative Pearson correlation across best-matched Visium-Xenium spots for a gene.

| **Parameter** | **Type**                                        | **Default**   | **Options**                                   | **Description**                                                                |
| ------------- | ----------------------------------------------- | ------------- | --------------------------------------------- | ------------------------------------------------------------------------------ |
| `obj`         | `SingleCellExperiment` / `SummarizedExperiment` | *(required)*  | —                                             | Object containing cross-platform expression matrix and metadata                |
| `gene_name`   | character                                       | *(required)*  | Any gene base name present with both suffixes | Gene to compare between platforms (must exist with appropriate suffixes)       |
| `ID1_col`     | character                                       | `"Visium_ID"` | Column name in `colData(obj)`                 | Column storing Visium-matched cell/spot IDs                                    |
| `ID2_col`     | character                                       | `"Xenium_ID"` | Column name in `colData(obj)`                 | Column storing Xenium-matched cell/spot IDs                                    |
| `hline`       | numeric                                         | `0.9`         | `NULL` or any numeric between -1 and 1        | Horizontal line for reference on the correlation plot                          |
| `y_range`     | numeric vector (length 2)                       | `c(0, 1)`     | Any valid range of Pearson correlation values | Y-axis limits for the plot                                                     |
| `title`       | character                                       | `NULL`        | —                                             | Custom plot title; if `NULL`, defaults to `"Gene: Correlation vs Top Matched"` |


[`The list of available genes is here.`](./available_genes_by_platform.md)

#### Example usages
Cumulative Pearson correlation (COL1A1)
```
plot_cross_platform_gene_correlation(
    obj = scTriangulate2,
    gene_name = "COL1A1",
    ID1_col = "Visium_ID",
    ID2_col = "Xenium_ID",
    y_range = c(0.8,1)
)
```
<img src="./figures/plot_cross_platform_gene_correlation_1.png" alt="Example" width="500"/>

---
#### `generate_stability_report()`
- Summarize per-cluster metrics (confidence, reassign, SCCAF, TF-IDF 10, TF-IDF 5m Shapley)

| **Parameter**     | **Type**                                        | **Default**  | **Options**                                          | **Description**                                                      |
| ----------------- | ----------------------------------------------- | ------------ | ---------------------------------------------------- | -------------------------------------------------------------------- |
| `sce`             | `SingleCellExperiment` / `SummarizedExperiment` | *(required)* | —                                                    | Object containing stability metric columns and cell type annotations |
| `celltype_labels` | `list` (named character vectors)                | *(required)* | Format: `list("Annotation" = c("Label1", "Label2"))` | Specifies which cell type labels (per annotation) to evaluate        |

[`The list of available annotations is here.`](./annotation_labels.md)

#### Example usages
Stability of T cells 
```
# Define Inputs and Generate Outputs
# Define the single list of clusters you want to analyze and plot
clusters_to_process <- list(
    Final_CT = c("CD4+_T-cells", "CD8+_T-cells"),
    Visium_RCTD_LungMap_ref = c("T"),
    Visium_RCTD_GSE250346_based = c("CD4+_T-cells", "CD8+_T-cells")
)

# Call the function once to get both the table and the plot
output <- generate_stability_report(
    sce = scTriangulate2, 
    celltype_labels = clusters_to_process
)


# Print Both the Table and the Plot
if (!is.null(output)) {
    
    # Print the summary table
    print(kable(output$summary_table, digits = 3, caption = "Summary of Cell Type Stability Metrics"))
    
    # Print the final plot
    print(output$stability_plot)
    
}
```
Output
```
Table: Summary of Cell Type Stability Metrics

|Cluster                                  | Reassign| SCCAF| TF-IDF 5| TF-IDF 10| Shapley| Confidence|
|:----------------------------------------|--------:|-----:|--------:|---------:|-------:|----------:|
|Final_CT@CD4+_T-cells                    |    0.562| 0.767|    0.727|     0.540|   1.783|      0.599|
|Final_CT@CD8+_T-cells                    |    0.456| 0.722|    0.837|     0.673|   6.917|      0.850|
|Visium_RCTD_LungMap_ref@T                |    0.706| 0.789|    0.623|     0.521|   8.083|      0.156|
|Visium_RCTD_GSE250346_based@CD4+_T-cells |    0.476| 0.290|    0.520|     0.469|   0.783|      0.850|
|Visium_RCTD_GSE250346_based@CD8+_T-cells |    0.000| 0.304|    0.332|     0.305|   0.067|      0.128|
```
<img src="./figures/generate_stability_report_1.png" alt="Example" width="800"/>

---
#### `rank_marker_specificity()`
- Rank how specific a gene is to a cell type using one-sided t-tests and FDR correction.

| **Parameter**  | **Type**                                        | **Default**             | **Options**                                           | **Description**                                                               |
| -------------- | ----------------------------------------------- | ----------------------- | ----------------------------------------------------- | ----------------------------------------------------------------------------- |
| `sce`          | `SingleCellExperiment` / `SummarizedExperiment` | *(required)*            | —                                                     | Object containing gene expression matrix and cluster annotations              |
| `gene_name`    | character                                       | `"COL1A1"`              | Any gene (with or without `_xenium`/`_visium` suffix) | Base name of the gene to assess for cluster-specific enrichment               |
| `cluster_cols` | character vector                                | `c("pruned")`           | One or more column names in `colData(sce)`            | Metadata columns used to define clusters or cell types                        |
| `assay_name`   | character                                       | `"X"`                   | Must exist in `assayNames(sce)`                       | Expression matrix to use                                                      |
| `platforms`    | character vector                                | `c("Xenium", "Visium")` | Subset of `c("Xenium", "Visium")`                     | Platforms to include when filtering by `meta$Platform`, if such column exists |
| `top_n`        | integer or `NULL`                               | `NULL`                  | `NULL` (all clusters) or a positive integer           | If set, returns only the top N cluster-gene combinations by adjusted p-value  |

[`The list of available genes is here.`](./available_genes_by_platform.md)

#### Example usages
COL1A1 - pruned, both Visium and Xenium
```
rank_marker_specificity(scTriangulate2, 
                        gene_name = "COL1A1",
                        top_n = 3)
```
Output
```
                    Cluster Mean_Expr Pct_Expr      P_Value Annotation_Source        P_Adj
1                 TNiche@T6 15.087205 99.83700 0.000000e+00            pruned 0.000000e+00
2 Final_CT@Inflammatory_FBs  8.892857 99.10714 2.637340e-07            pruned 3.824143e-06
3  Final_CT@Adventitial_FBs  6.216667 93.33333 3.877084e-01            pruned 1.000000e+00
```
Top 5 clusters across selected annotations
```
rank_marker_specificity(scTriangulate2,
                        gene_name = "COL1A1",
                        cluster_cols = c("Final_CT", "TNiche", "Xenium_RCTD_LungMap_ref"),
                        top_n = 5)
```
Output
```
                 Cluster Mean_Expr  Pct_Expr       P_Value       Annotation_Source         P_Adj
1                     T6 15.398747 100.00000  0.000000e+00                  TNiche  0.000000e+00
2                   SCMF 13.409579  99.00383  0.000000e+00 Xenium_RCTD_LungMap_ref  0.000000e+00
3 Activated_Fibrotic_FBs 17.800000  99.76744 6.613296e-302                Final_CT 1.851723e-300
4           Alveolar_FBs  9.208145  96.15385  2.075145e-59                Final_CT  4.357804e-58
5                    AF1 10.039419  93.77593  2.415718e-20 Xenium_RCTD_LungMap_ref  4.058406e-19
```
Visium-only, limited annotations
```
rank_marker_specificity(scTriangulate2,
                        gene_name = "COL1A1",
                        cluster_cols = c("Visium_RCTD_GSE250346_based"),
                        platforms = "Visium",
                        top_n = 3)
```
Output
```
                 Cluster Mean_Expr Pct_Expr       P_Value           Annotation_Source         P_Adj
1 Activated_Fibrotic_FBs 13.004885 97.06927 5.559134e-288 Visium_RCTD_GSE250346_based 2.501610e-286
2        Adventitial_FBs  7.719212 88.17734  1.283858e-07 Visium_RCTD_GSE250346_based  2.888679e-06
3         Subpleural_FBs  6.034483 93.10345  3.533602e-01 Visium_RCTD_GSE250346_based  1.000000e+00
```
---
#### `find_celltype_markers()`
- Identify top marker genes for a specific cell type by calculating their differential expression and statistical significance, returning both a ranked table and a volcano plot to visualize the results.

| **Parameter**     | **Type**                                        | **Default**  | **Options**                                   | **Description**                                                                |
| ----------------- | ----------------------------------------------- | ------------ | --------------------------------------------- | ------------------------------------------------------------------------------ |
| `sce`             | `SingleCellExperiment` / `SummarizedExperiment` | *(required)* | —                                             | Object containing normalized gene expression and metadata                      |
| `celltype`        | character                                       | *(required)* | Any value in `colData(sce)[[annotation_col]]` | Cell type label to compare against all others                                  |
| `annotation_col`  | character                                       | *(required)* | A column name in `colData(sce)`               | Metadata column that defines cell types or clusters                            |
| `platform`        | character                                       | `"Xenium"`   | `"Xenium"`, `"Visium"`                        | Platform used to select gene suffix (e.g., `_xenium`) for expression filtering |
| `assay_name`      | character                                       | `"X"`        | Any value in `assayNames(sce)`                | Expression assay to use                                                        |
| `top_n`           | integer                                         | `10`         | Any positive integer                          | Number of top-ranked marker genes to return                                    |
| `logfc_threshold` | numeric                                         | `0.25`       | ≥ 0                                           | Minimum log2 fold change threshold to consider as significant                  |
| `padj_threshold`  | numeric                                         | `0.05`       | (0, 1)                                        | Adjusted p-value (FDR) threshold for statistical significance                  |


#### Example usages
Identify top marker genes for T6
```
# Define Inputs and Find Markers
marker_report <- find_celltype_markers(
    sce = scTriangulate2,
    celltype = "T6",
    annotation_col = "TNiche",
    platform = "xenium"
)
# Access and Print Both the Table and the Plot
if (!is.null(marker_report)) {
    # Print the ranked table of top marker genes
    print(kable(marker_report$table, digits = 3, caption = "Top Marker Genes"))
    
    # Print the volcano plot
    print(marker_report$plot)
}
```
Output
```
Table: Top Marker Genes

|gene   | logFC| P_Value| P_Adj| log10_P_Adj|significant |
|:------|-----:|-------:|-----:|-----------:|:-----------|
|COL1A1 | 3.647|       0|     0|     335.145|TRUE        |
|MEG3   | 3.508|       0|     0|     335.145|TRUE        |
|COL3A1 | 3.403|       0|     0|     335.145|TRUE        |
|SFRP4  | 3.148|       0|     0|     335.145|TRUE        |
|LUM    | 3.067|       0|     0|     335.145|TRUE        |
|COL1A2 | 3.058|       0|     0|     335.145|TRUE        |
|CTHRC1 | 3.047|       0|     0|     335.145|TRUE        |
|DCN    | 2.268|       0|     0|     335.145|TRUE        |
|POSTN  | 1.964|       0|     0|     335.145|TRUE        |
|FN1    | 1.251|       0|     0|     335.145|TRUE        |
```
<img src="./figures/find_celltype_markers_1.png" alt="Example" width="700"/>

---
#### `generate_marker_report()`
- Analyzes the expression specificity of given genes across various cell type annotations, returning both a ranked statistical table and a dot plot summarizing the results.

| **Parameter**     | **Type**                                        | **Default**  | **Options**                                                           | **Description**                                                              |
| ----------------- | ----------------------------------------------- | ------------ | --------------------------------------------------------------------- | ---------------------------------------------------------------------------- |
| `obj`             | `SingleCellExperiment` / `SummarizedExperiment` | *(required)* | —                                                                     | Object containing expression data and cell metadata                          |
| `genes`           | `character vector`                              | *(required)* | Any set of gene names (e.g., `"COL1A1"`, `"AGER"`, `"CD3D"`)          | Gene names to test for enrichment                                            |
| `annotation_cols` | `character vector`                              | *(required)* | Must match columns in `colData(obj)` (e.g., `"Final_CT"`, `"pruned"`) | Annotation columns used to group cells for differential expression testing   |
| `assay_name`      | `character`                                     | `"X"`        | Any valid assay in `assayNames(obj)`                                  | Assay layer to pull expression data from                                     |
| `min_pct`         | `numeric`                                       | `5`          | 0–100                                                                 | Minimum percentage of expressing cells in a cluster to include the gene      |
| `max_genes_plot`  | `integer`                                       | `50`         | Positive integer or `NULL`                                            | Max number of genes to display in dot plot; `NULL` shows all passing filters |

[`The list of available genes is here.`](./available_genes_by_platform.md)

#### Example usages
List cell types that "COL1A1", "AGER", "CD3D", "SFTPC" are differentially expressed.
```
# Define Inputs and Generate Outputs
genes_of_interest <- c("COL1A1", "AGER", "CD3D", "SFTPC")
annotations_to_process <- c("pruned")

# Call the function
marker_report <- generate_marker_report(
    obj = scTriangulate2,
    genes = genes_of_interest,
    annotation_cols = annotations_to_process
)


# Access and Print Both the Table and the Plot
if (!is.null(marker_report)) {
    # Print the full, ranked data table
    print(kable(marker_report$table, digits = 3, caption = "Marker Gene Specificity Ranking"))
    
    # Print the dot plot
    print(marker_report$plot)
}
```
Output (Portion)
```
Table: Marker Gene Specificity Ranking

|Gene   |Platform |Cluster                                             | Mean_Expr| Pct_Expr| P_Value| P_Adj|
|:------|:--------|:---------------------------------------------------|---------:|--------:|-------:|-----:|
|COL1A1 |Visium   |pruned@TNiche@T6                                    |    13.681|   99.307|   0.000| 0.000|
|COL1A1 |Xenium   |pruned@TNiche@T6                                    |    15.087|   99.837|   0.000| 0.000|
|SFTPC  |Xenium   |pruned@TNiche@T11                                   |    18.464|   93.724|   0.000| 0.000|
|CD3D   |Xenium   |pruned@Xenium_RCTD_LungMap_ref@T                    |     1.008|   60.615|   0.000| 0.000|
|CD3D   |Xenium   |pruned@Final_CT@CD8+_T-cells                        |     1.052|   56.897|   0.000| 0.000|
|SFTPC  |Visium   |pruned@TNiche@T11                                   |     0.732|   31.799|   0.000| 0.000|
|AGER   |Xenium   |pruned@Final_CT@AT1                                 |     3.209|   79.104|   0.000| 0.000|
|AGER   |Xenium   |pruned@Xenium_RCTD_LungMap_ref@AT1                  |     2.037|   57.009|   0.000| 0.000|
|SFTPC  |Visium   |pruned@Visium_RCTD_LungMap_ref@AT2                  |     0.990|   46.392|   0.000| 0.000|
|AGER   |Xenium   |pruned@TNiche@T11                                   |     0.427|   27.615|   0.000| 0.000|
|COL1A1 |Xenium   |pruned@Final_CT@Inflammatory_FBs                    |     8.893|   99.107|   0.000| 0.000|
|CD3D   |Visium   |pruned@Xenium_RCTD_LungMap_ref@T                    |     0.120|    7.821|   0.000| 0.000|
|CD3D   |Xenium   |pruned@Visium_RCTD_LungMap_ref@B                    |     0.241|   14.815|   0.004| 0.065|
|COL1A1 |Visium   |pruned@Final_CT@Inflammatory_FBs                    |     7.312|   99.107|   0.010| 0.167|
|AGER   |Visium   |pruned@Xenium_RCTD_LungMap_ref@AT1                  |     0.047|    4.673|   0.017| 0.265|
|SFTPC  |Xenium   |pruned@Xenium_RCTD_LungMap_ref@AT1                  |     1.336|   35.514|   0.019| 0.282|
|CD3D   |Xenium   |pruned@Visium_RCTD_LungMap_ref@AEC                  |     0.323|   14.516|   0.021| 0.283|
|SFTPC  |Visium   |pruned@Xenium_RCTD_LungMap_ref@AT1                  |     0.178|    8.411|   0.051| 0.653|
|AGER   |Visium   |pruned@Visium_RCTD_LungMap_ref@AT2                  |     0.031|    3.093|   0.059| 0.719|
|CD3D   |Xenium   |pruned@Visium_RCTD_LungMap_ref@Pericyte             |     0.273|   10.606|   0.065| 0.759|
```
<img src="./figures/generate_marker_report_1.png" alt="Example" width="900"/>

---
#### `plot_ARI_dotplot()`
- Plot pairwise cluster overlap using dotplot and display Adjusted Rand Index (ARI).

| **Parameter** | **Type**                                         | **Default**  | **Options**                                                                                                                                    | **Description**                                                            |
| ------------- | ------------------------------------------------ | ------------ | ---------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------- |
| `obj`         | `SingleCellExperiment` or `SummarizedExperiment` | *(required)* | —                                                                                                                                              | Object containing cell metadata with clustering annotations.               |
| `annotation1` | `character`                                      | *(required)* | `pruned`, `Final_CT`, `Final_lineage`, `CNiche`, `TNiche`, `Visium_RCTD_LungMap_ref`, `Visium_RCTD_GSE250346_based`, `Xenium_RCTD_LungMap_ref` | First clustering annotation to compare.                                    |
| `annotation2` | `character`                                      | *(required)* | `pruned`, `Final_CT`, `Final_lineage`, `CNiche`, `TNiche`, `Visium_RCTD_LungMap_ref`, `Visium_RCTD_GSE250346_based`, `Xenium_RCTD_LungMap_ref` | Second clustering annotation to compare.                                   |
| `title`       | `character` or `NULL`                            | `NULL`       | Custom string or `NULL`                                                                                                                        | Optional plot title; if `NULL`, auto-generates a title with the ARI value. |
| `size_range`  | `numeric vector (length = 2)`                    | `c(1, 10)`   | Any numeric range                                                                                                                              | Range of point sizes in the dot plot based on cross-tabulated frequency.   |

#### Example usages
ARI between Visium HD and Xenium RCTD based on LungMap ref
```
plot_ARI_dotplot(
  obj = scTriangulate2,
  annotation1 = "Xenium_RCTD_LungMap_ref",
  annotation2 = "Visium_RCTD_LungMap_ref",
  title = "RCTD Matching Dot Plot"
)
```
<img src="./figures/plot_ARI_dotplot_1.png" alt="Example" width="900"/>

---
#### `plot_cluster_composition_groupedbar()`
- Compare cluster or cell-type compositions between one or two annotations using barplots.

| **Parameter** | **Type**                                        | **Default**  | **Options**                                                                                                                                    | **Description**                                                             |
| ------------- | ----------------------------------------------- | ------------ | ---------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------- |
| `obj`         | `SingleCellExperiment` / `SummarizedExperiment` | *(required)* | —                                                                                                                                              | Object containing clustering annotations and metadata in `colData`.         |
| `annotation`  | `character vector` (length 1 or 2)              | *(required)* | `pruned`, `Final_CT`, `Final_lineage`, `CNiche`, `TNiche`, `Visium_RCTD_LungMap_ref`, `Visium_RCTD_GSE250346_based`, `Xenium_RCTD_LungMap_ref` | One or two annotation columns to compare cluster composition.               |
| `palette`     | named `character` vector or `NULL`              | `NULL`       | e.g., `c("Visium" = "blue", "Xenium" = "red")`                                                                                                 | Custom fill colors for annotation sources. If `NULL`, uses default palette. |
| `title`       | `character` or `NULL`                           | `NULL`       | —                                                                                                                                              | Custom plot title. If `NULL`, uses `"Cluster Composition Comparison"`.      |
| `rotate_x`    | `logical`                                       | `TRUE`       | `TRUE`, `FALSE`                                                                                                                                | Whether to rotate x-axis labels by 45° for readability.                     |

#### Example usages
`pruned` cell type composition
```
plot_cluster_composition_groupedbar(
    obj = scTriangulate2,
    annotation = c("pruned")
)
```
<img src="./figures/plot_cluster_composition_groupedbar_1.png" alt="Example" width="900"/>

Comparison between `Visium_RCTD_LungMap_ref` and `Xenium_RCTD_LungMap_ref`
```
plot_cluster_composition_groupedbar(
    obj = scTriangulate2,
    annotation = c("Visium_RCTD_LungMap_ref", "Xenium_RCTD_LungMap_ref")
)
```
<img src="./figures/plot_cluster_composition_groupedbar_2.png" alt="Example" width="900"/>

---
#### `compare_annotations_via_heatmap()`
- Visualize overlap between two annotations using a normalized confusion matrix heatmap.

| **Parameter** | **Type**                                        | **Default**  | **Options**                                                                                                                                    | **Description**                                                       |
| ------------- | ----------------------------------------------- | ------------ | ---------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------- |
| `obj`         | `SingleCellExperiment` / `SummarizedExperiment` | *(required)* | —                                                                                                                                              | Object containing annotation metadata in `colData`.                   |
| `annotation1` | `character`                                     | *(required)* | `pruned`, `Final_CT`, `Final_lineage`, `CNiche`, `TNiche`, `Visium_RCTD_LungMap_ref`, `Visium_RCTD_GSE250346_based`, `Xenium_RCTD_LungMap_ref` | First annotation column to use as heatmap rows.                       |
| `annotation2` | `character`                                     | *(required)* | Same as `annotation1`                                                                                                                          | Second annotation column to use as heatmap columns.                   |
| `normalize`   | `character`                                     | `"none"`     | `"none"`, `"row"`, `"column"`, `"both"`                                                                                                        | Normalize values by row, column, both, or not at all (raw counts).    |
| `palette`     | `character`                                     | `"Blues"`    | Any valid palette from `RColorBrewer::brewer.pal.info` (e.g., `"Reds"`, `"Purples"`)                                                           | Color palette for the heatmap fill.                                   |
| `title`       | `character` or `NULL`                           | `NULL`       | —                                                                                                                                              | Custom plot title. If `NULL`, uses `"Annotation Comparison Heatmap"`. |

#### Example usages
Comparison between `pruned` and `Final_CT`
```
compare_annotations_via_heatmap(
    obj = scTriangulate2,
    annotation1 = "pruned",
    annotation2 = "Final_CT",
    normalize = "row"
)
```
<img src="./figures/compare_annotations_via_heatmap_1.png" alt="Example" width="900"/>

---
#### `visualize_matched_pairs_scatter()`
- Scatterplot showing gene expression correlation across Visium and Xenium matched spots.

| **Parameter**      | **Type**                                    | **Default**  | **Options / Example Values** | **Description**                                                       |
| ------------------ | ------------------------------------------- | ------------ | ---------------------------- | --------------------------------------------------------------------- |
| `obj`              | `SingleCellExperiment` or compatible object | *(required)* | —                            | Object with gene expression matrix and matched Visium/Xenium entries. |
| `gene`             | `character`                                 | *(required)* | e.g., `"COL1A1"`             | Base gene name (assumes suffixes `_visium` and `_xenium` exist).      |
| `layer`            | `character`                                 | `"X"`        | Any valid assay layer        | Layer from which to extract expression data.                          |
| `log1p`            | `logical`                                   | `TRUE`       | `TRUE`, `FALSE`              | Apply log(1 + x) transformation to expression values.                 |
| `point_size`       | `numeric`                                   | `0.4`        | e.g., `0.2`, `0.5`, `1.0`    | Size of scatter plot points.                                          |
| `alpha`            | `numeric`                                   | `0.6`        | Range: 0–1                   | Opacity of scatter plot points.                                       |
| `show_correlation` | `logical`                                   | `FALSE`      | `TRUE`, `FALSE`              | If `TRUE`, displays Pearson correlation in the title.                 |
| `title`            | `character` o                               |              |                              |                                                                       |

[`The list of available genes is here.`](./available_genes_by_platform.md)
#### Example usages
COL1A1 correlation between Visium HD and Xenium
```
visualize_matched_pairs_scatter(
  obj = scTriangulate2,
  gene = "COL1A1",
  show_correlation = TRUE
)
```
<img src="./figures/visualize_matched_pairs_scatter_1.png" alt="Example" width="700"/>

---
#### `generate_de_volcano_plot()`
- Generate table and volcano plot of DE genes
#### Example usages
Identify DEGs in 'More_Affected' T6 in IPF samples comapred to 'Unaffected' T6 control samples
```
generate_de_volcano_plot(
    obj_list = CellChat_Squidpy[["CellChat"]],
    annotation = "TNiche",
    group_1 = "More_Affected",
    group_2 = "Unaffected",
    celltype_filter = "T6",
    label_top_n = 10
)
```
Output (Top 10 genes)
```
        Gene         p_val avg_log2FC pct.1 pct.2     p_val_adj log10_padj             Sig
1      POSTN 1.782835e-270  2.1740071 0.576 0.152 6.115124e-268 267.213595   Up in Group 1
2     S100A8 2.645121e-161 -1.8177308 0.020 0.113 9.072766e-159 158.042260 Down in Group 1
3     PDGFRA 2.454663e-145 -1.3288754 0.287 0.517 8.419495e-143 142.074714 Down in Group 1
4      PLIN2 2.455126e-135 -1.6178956 0.106 0.278 8.421081e-133 132.074632 Down in Group 1
5   HIST1H1C 4.041545e-129 -1.4691020 0.140 0.328 1.386250e-126 125.858158 Down in Group 1
6      SFRP2  1.517121e-89  1.1355709 0.486 0.254  5.203725e-87  86.283686   Up in Group 1
7       SOD2  1.167657e-82 -0.9161804 0.355 0.533  4.005063e-80  79.397391 Down in Group 1
8       EGFR  1.174035e-72 -1.1470345 0.105 0.231  4.026941e-70  69.395025 Down in Group 1
9     COL1A1  5.299139e-65  0.3064564 0.959 0.933  1.817605e-62  61.740501   Up in Group 1
10     MFAP5  4.301939e-64 -2.0013742 0.016 0.065  1.475565e-61  60.831042 Down in Group 1
```
<img src="./figures/generate_de_volcano_plot_1.png" alt="Example" width="700"/>

Indentify DEGs in 'More_Affected' in IPF samples comapred to 'Unaffected' control samples
```
generate_de_volcano_plot(
    obj_list = CellChat_Squidpy[["CellChat"]],
    annotation = "TNiche",
    group_1 = "More_Affected",
    group_2 = "Unaffected"
)
```
Output (Top 10 genes)
```
$table
         Gene p_val avg_log2FC pct.1 pct.2 p_val_adj log10_padj             Sig
1       MS4A1     0   3.928060 0.051 0.008         0        300   Up in Group 1
2       CD79A     0   3.757472 0.066 0.006         0        300   Up in Group 1
3         CA4     0  -3.645654 0.011 0.111         0        300 Down in Group 1
4        MMP7     0   3.191984 0.076 0.016         0        300   Up in Group 1
5      JCHAIN     0   3.188527 0.252 0.045         0        300   Up in Group 1
6       SFRP2     0   3.173698 0.203 0.030         0        300   Up in Group 1
7   TNFRSF13C     0   3.107436 0.057 0.012         0        300   Up in Group 1
8       FOXP3     0   2.890871 0.016 0.003         0        300   Up in Group 1
9        FCN3     0  -2.884330 0.035 0.195         0        300 Down in Group 1
10       AGER     0  -2.841862 0.022 0.116         0        300 Down in Group 1
```
<img src="./figures/generate_de_volcano_plot_2.png" alt="Example" width="700"/>

---
#### 2. CellChat
Cell-Cell communication analysis pakage developed by Jin et al.

#### Pre-processing
```
# Pre-processing
library(CellChat)
# Initialize output list to store group sizes
groupSizes <- list()
# Loop over annotations and affect groups
for (annot in names(CellChat_Squidpy$CellChat)) {
  groupSizes[[annot]] <- list()
  for (affect in names(CellChat_Squidpy$CellChat[[annot]])) {
    obj <- CellChat_Squidpy$CellChat[[annot]][[affect]]
    # Store group size table
    groupSizes[[annot]][[affect]] <- as.numeric(table(obj@idents))
    # Update with centrality computation
    CellChat_Squidpy$CellChat[[annot]][[affect]] <- netAnalysis_computeCentrality(obj)
  }
}
```

#### `netVisual_circle()`
- Circle plot of cell-cell communication network
```
netVisual_circle(CellChat_Squidpy[["CellChat"]][["Final_CT"]][["Unaffected"]]@net$count, 
                 vertex.weight = groupSizes$CellChat_Squidpy[["CellChat"]][["Final_CT"]][["Unaffected"]], 
                 vertex.label.cex = 0.7, 
                 margin = 0.2, 
                 remove.isolate = TRUE)
```
<img src="./figures/netVisual_circle_1.png" alt="Example" width="500"/>

---
#### List all inferred signaling pathways
Control vs IPF (Annotations based on those annotated by Vannan et al.)
```
> # Unaffected 
> CellChat_Squidpy[["CellChat"]][["Final_CT"]][["Unaffected"]]@netP$pathways
[1] "VEGF"  "SPP1"  "UGRP1" "CCL"   "FASLG" "EGF"  
> # IPF
> CellChat_Squidpy[["CellChat"]][["Final_CT"]][["More_Affected"]]@netP$pathways
[1] "SPP1"  "UGRP1" "VEGF"  "CCL"   "EGF"   "CXCL"  "FASLG"
```

---
#### `netVisual_heatmap()`
- Generates a heatmap of intercellular communication strength
#### Example usages
Xenium, IPF
```
netVisual_heatmap(CellChat_Squidpy[["CellChat"]][["Final_CT"]][["More_Affected"]], 
                  measure = "weight", 
                  color.heatmap = "Reds", 
                  title.name = "Xenium IPF")
```
<img src="./figures/netVisual_heatmap_1.png" alt="Example" width="1000"/>

Xenium, Control
```
netVisual_heatmap(CellChat_Squidpy[["CellChat"]][["Final_CT"]][["Unaffected"]], 
                  measure = "weight", 
                  color.heatmap = "Reds", 
                  title.name = "Xenium Control")
```
<img src="./figures/netVisual_heatmap_2.png" alt="Example" width="1000"/>

---
#### `netVisual_aggregate()`
- Visualize the inferred signaling network of signaling pathways by aggregating all L-R pairs
#### Example usages
Circle plot that visualizes the intercellular communication network for the UGRP1 signaling pathway in the IPFs
```
netVisual_aggregate(CellChat_Squidpy[["CellChat"]][["Final_CT"]][["More_Affected"]], 
                    signaling = "UGRP1", 
                    layout = "circle", 
                    remove.isolate = TRUE)
```
<img src="./figures/netVisual_aggregate_1.png" alt="Example" width="500"/>

---
#### `subsetCommunication()`
- Subset the inferred cell-cell communications of interest
#### Example usages
Returns a data frame summarizing the UGRP1 pathway-level communication network from the 'Xenium_IPF_Author', author annotation CellChat object.
```
subsetCommunication(CellChat_Squidpy[["CellChat"]][["Final_CT"]][["More_Affected"]],
                    signaling = "UGRP1", 
                    slot.name = "netP")
```
Output
```
                  source                      target pathway_name         prob pval
1                    AT1        Alveolar Macrophages        UGRP1 3.572308e-04    0
2                    AT1 Macrophages - IFN-activated        UGRP1 6.142392e-05    0
3                    AT1              Monocytes/MDMs        UGRP1 1.353383e-05    0
4                    AT1       Proliferating Myeloid        UGRP1 8.546132e-05    0
5                    AT1           SPP1+ Macrophages        UGRP1 2.808011e-04    0
6                    AT2        Alveolar Macrophages        UGRP1 9.379655e-04    0
7                    AT2 Macrophages - IFN-activated        UGRP1 1.613557e-04    0
8                    AT2              Monocytes/MDMs        UGRP1 3.555505e-05    0
9                    AT2       Proliferating Myeloid        UGRP1 2.244913e-04    0
10                   AT2           SPP1+ Macrophages        UGRP1 7.373790e-04    0
11                 Basal        Alveolar Macrophages        UGRP1 1.205943e-03    0
12                 Basal Macrophages - IFN-activated        UGRP1 2.075013e-04    0
13                 Basal              Monocytes/MDMs        UGRP1 4.572496e-05    0
14                 Basal       Proliferating Myeloid        UGRP1 2.886875e-04    0
15                 Basal           SPP1+ Macrophages        UGRP1 9.481030e-04    0
16             Basophils        Alveolar Macrophages        UGRP1 1.350744e-03    0
17             Basophils Macrophages - IFN-activated        UGRP1 2.324444e-04    0
18             Basophils              Monocytes/MDMs        UGRP1 5.122242e-05    0
19             Basophils       Proliferating Myeloid        UGRP1 3.233867e-04    0
20             Basophils           SPP1+ Macrophages        UGRP1 1.061977e-03    0
21                Goblet        Alveolar Macrophages        UGRP1 1.521895e-03    0
22                Goblet Macrophages - IFN-activated        UGRP1 2.619344e-04    0
23                Goblet              Monocytes/MDMs        UGRP1 5.772228e-05    0
24                Goblet       Proliferating Myeloid        UGRP1 3.644102e-04    0
25                Goblet           SPP1+ Macrophages        UGRP1 1.196583e-03    0
26          KRT5-/KRT17+        Alveolar Macrophages        UGRP1 1.161897e-03    0
27          KRT5-/KRT17+ Macrophages - IFN-activated        UGRP1 1.999152e-04    0
28          KRT5-/KRT17+              Monocytes/MDMs        UGRP1 4.405302e-05    0
29          KRT5-/KRT17+       Proliferating Myeloid        UGRP1 2.781341e-04    0
30          KRT5-/KRT17+           SPP1+ Macrophages        UGRP1 9.134657e-04    0
31      Langerhans cells        Alveolar Macrophages        UGRP1 1.696718e-03    0
32      Langerhans cells              Monocytes/MDMs        UGRP1 6.436380e-05    0
33      Langerhans cells       Proliferating Myeloid        UGRP1 4.063248e-04    0
34      Langerhans cells           SPP1+ Macrophages        UGRP1 1.334087e-03    0
35         Multiciliated        Alveolar Macrophages        UGRP1 1.020553e-03    0
36         Multiciliated Macrophages - IFN-activated        UGRP1 1.755750e-04    0
37         Multiciliated              Monocytes/MDMs        UGRP1 3.868873e-05    0
38         Multiciliated       Proliferating Myeloid        UGRP1 2.442730e-04    0
39         Multiciliated           SPP1+ Macrophages        UGRP1 8.023189e-04    0
40        Myofibroblasts        Alveolar Macrophages        UGRP1 8.553471e-04    0
41        Myofibroblasts Macrophages - IFN-activated        UGRP1 1.471330e-04    0
42        Myofibroblasts              Monocytes/MDMs        UGRP1 3.242069e-05    0
43        Myofibroblasts       Proliferating Myeloid        UGRP1 2.047047e-04    0
44        Myofibroblasts           SPP1+ Macrophages        UGRP1 6.724169e-04    0
45                NK/NKT        Alveolar Macrophages        UGRP1 1.971445e-06    0
46                NK/NKT Macrophages - IFN-activated        UGRP1 3.388796e-07    0
47                NK/NKT              Monocytes/MDMs        UGRP1 7.466341e-08    0
48                NK/NKT       Proliferating Myeloid        UGRP1 4.715067e-07    0
49                NK/NKT           SPP1+ Macrophages        UGRP1 1.549535e-06    0
50                  PNEC        Alveolar Macrophages        UGRP1 2.951154e-03    0
51                  PNEC              Monocytes/MDMs        UGRP1 1.120854e-04    0
52                  PNEC       Proliferating Myeloid        UGRP1 7.074093e-04    0
53                  PNEC           SPP1+ Macrophages        UGRP1 2.321042e-03    0
54     Proliferating AT2        Alveolar Macrophages        UGRP1 1.649265e-03    0
55     Proliferating AT2 Macrophages - IFN-activated        UGRP1 2.838860e-04    0
56     Proliferating AT2              Monocytes/MDMs        UGRP1 6.256082e-05    0
57     Proliferating AT2       Proliferating Myeloid        UGRP1 3.949465e-04    0
58     Proliferating AT2           SPP1+ Macrophages        UGRP1 1.296762e-03    0
59  Proliferating Airway        Alveolar Macrophages        UGRP1 2.360339e-03    0
60  Proliferating Airway Macrophages - IFN-activated        UGRP1 4.065221e-04    0
61  Proliferating Airway              Monocytes/MDMs        UGRP1 8.959507e-05    0
62  Proliferating Airway       Proliferating Myeloid        UGRP1 5.655325e-04    0
63  Proliferating Airway           SPP1+ Macrophages        UGRP1 1.856139e-03    0
64 Proliferating Myeloid        Alveolar Macrophages        UGRP1 1.341633e-05    0
65 Proliferating Myeloid Macrophages - IFN-activated        UGRP1 2.306209e-06    0
66 Proliferating Myeloid              Monocytes/MDMs        UGRP1 5.081145e-07    0
67 Proliferating Myeloid       Proliferating Myeloid        UGRP1 3.208785e-06    0
68 Proliferating Myeloid           SPP1+ Macrophages        UGRP1 1.054512e-05    0
69                  RASC        Alveolar Macrophages        UGRP1 4.398891e-03    0
70                  RASC Macrophages - IFN-activated        UGRP1 7.589065e-04    0
71                  RASC              Monocytes/MDMs        UGRP1 1.673045e-04    0
72                  RASC       Proliferating Myeloid        UGRP1 1.055606e-03    0
73                  RASC           SPP1+ Macrophages        UGRP1 3.460740e-03    0
74     SPP1+ Macrophages        Alveolar Macrophages        UGRP1 1.735101e-05    0
75     SPP1+ Macrophages Macrophages - IFN-activated        UGRP1 2.982573e-06    0
76     SPP1+ Macrophages              Monocytes/MDMs        UGRP1 6.571346e-07    0
77     SPP1+ Macrophages       Proliferating Myeloid        UGRP1 4.149856e-06    0
78     SPP1+ Macrophages           SPP1+ Macrophages        UGRP1 1.363776e-05    0
79             Secretory        Alveolar Macrophages        UGRP1 4.721735e-03    0
80             Secretory Macrophages - IFN-activated        UGRP1 8.148229e-04    0
81             Secretory              Monocytes/MDMs        UGRP1 1.796394e-04    0
82             Secretory       Proliferating Myeloid        UGRP1 1.133358e-03    0
83             Secretory           SPP1+ Macrophages        UGRP1 3.714988e-03    0
84      Transitional AT2        Alveolar Macrophages        UGRP1 3.905875e-03    0
85      Transitional AT2 Macrophages - IFN-activated        UGRP1 6.735742e-04    0
86      Transitional AT2              Monocytes/MDMs        UGRP1 1.484827e-04    0
87      Transitional AT2       Proliferating Myeloid        UGRP1 9.369436e-04    0
88      Transitional AT2           SPP1+ Macrophages        UGRP1 3.072545e-03    0
89                  cDCs        Alveolar Macrophages        UGRP1 6.149080e-05    0
90                  cDCs Macrophages - IFN-activated        UGRP1 1.057042e-05    0
91                  cDCs              Monocytes/MDMs        UGRP1 2.328940e-06    0
92                  cDCs       Proliferating Myeloid        UGRP1 1.470731e-05    0
93                  cDCs           SPP1+ Macrophages        UGRP1 4.833176e-05    0
```
---
#### 3. Neighbors enrichment_analysis
- Compute enrichment score based on proximity on the connectivity graph of cell clusters by Palla, Spitzer et al. (2022) in Nature Methods.

#### `create_enrichment_heatmap()`

#### Example usages
Neighbors enrichment Unaffected (Control) vs More_Affected (IPF)
```
create_enrichment_heatmap(
    obj = CellChat_Squidpy,
    group = "Unaffected",
    annotation = "Final_CT",
    plot_style = "viridis",
    clustering = FALSE,
    base_size = 9
    )

create_enrichment_heatmap(
    obj = CellChat_Squidpy,
    group = "More_Affected",
    annotation = "Final_CT",
    plot_style = "viridis",
    clustering = FALSE,
    base_size = 9)
```
<img src="./figures/create_enrichment_heatmap_1.png" alt="Example" width="700"/>
<img src="./figures/create_enrichment_heatmap_2.png" alt="Example" width="700"/>

Clustering (IPF)
```
create_enrichment_heatmap(
    obj = CellChat_Squidpy,
    group = "More_Affected",
    annotation = "Final_CT",
    plot_style = "viridis",
    clustering = TRUE,
    base_size = 9)
```
<img src="./figures/create_enrichment_heatmap_3.png" alt="Example" width="700"/>

---
#### `analyze_and_plot_neighbors()`

#### Example usages
Cell Types Adjacent to AT1 in IPF samples
```
analyze_and_plot_neighbors(
    obj = CellChat_Squidpy,
    target_cell_type = "AT1",
    group = "More_Affected",
    annotation = "Final_CT"
)
```
Output
```
--- Cell Types Adjacent to AT1 in group: More_Affected ---
                     cell_type   log_score
1                          AT2  0.48835763
2                          AT1  0.47988747
3                    Capillary  0.46062549
4             Transitional AT2  0.45113207
5                 KRT5-/KRT17+  0.39260819
6         Alveolar Macrophages  0.37102875
7            Proliferating AT2  0.33454271
8       Activated Fibrotic FBs  0.32471795
9                       NK/NKT  0.32403627
10                        RASC  0.30969937
11              Myofibroblasts  0.30850956
12              Monocytes/MDMs  0.20140170
13                 Neutrophils  0.18579216
14           SPP1+ Macrophages  0.14433617
15       Proliferating Myeloid  0.12802797
16            Langerhans cells  0.11498948
17        Proliferating NK/NKT  0.11194960
18                        cDCs  0.03191334
19                 Mesothelial -0.02875905
20           Proliferating FBs -0.08338331
21                        PNEC -0.14449641
22       Proliferating B cells -0.16533357
23                   Basophils -0.16803249
24                   Arteriole -0.20064254
25 Macrophages - IFN-activated -0.20239647
26        Proliferating Airway -0.20580391
27       Proliferating T-cells -0.21411781
28    Interstitial Macrophages -0.23468022
29                        pDCs -0.23533174
30            Inflammatory FBs -0.26345540
31               Migratory DCs -0.26698992
32                   Secretory -0.27205697
33                       Tregs -0.30738885
34              Subpleural FBs -0.31924609
35             Adventitial FBs -0.32296172
36                   Lymphatic -0.32844154
37                        Mast -0.33705944
38                Alveolar FBs -0.34034091
39                       Basal -0.35337630
40                      Goblet -0.36002020
41              SMCs/Pericytes -0.37305499
42                      Venous -0.37539589
43                CD8+ T-cells -0.37654612
44               Multiciliated -0.39365303
45                CD4+ T-cells -0.39592181
46                     B cells -0.39938939
47                      Plasma -0.41039963
```
<img src="./figures/analyze_and_plot_neighbors_1.png" alt="Example" width="700"/>

---
#### 4. Ingestion
TO BE UPDATED 
