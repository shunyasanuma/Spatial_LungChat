## 3. Visium HD - Xenium alignment
This repository contains an R-based pipeline to align 10x Genomics Xenium and Visium HD spatial transcriptomics data from the same tissue section, with the goal of harmonizing platforms and refining cell-type annotations using spatial scTriangulate.

#### Goal
To identify the most reliable cell-type annotations by:
- Aligning matched Xenium and Visium HD spots using gene expression similarity
- Evaluating cross-platform transcript concordance
- Integrating metadata and RCTD results
- Prepareing input files for spatial scTriangulate to optimize annotations across platforms

#### Script
- [`Rmarkdown`](./Alignment.md)
- [`Rmarkdown with figures (html format)`](https://shunyasanuma.github.io/Spatial_LungChat/Alignment.html)
---

##### Workflow Steps
1. Extracting Singlet Coordinates
To ensure high-confidence cell mapping, only RCTD-classified singlets are used from each platform:
- From Visium2@results$results_df and Xenium2@results$results_df, spots with spot_class == "singlet" are extracted.
- Raw spatial coordinates (`x`, `y`) for those singlet spots are retained from:
  * `Visium2@spatialRNA@coords`
  *  `Xenium2@spatialRNA@coords`











