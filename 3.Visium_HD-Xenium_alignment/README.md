## 3. Visium HD - Xenium alignment
This repository contains an R-based pipeline to align 10x Genomics Xenium and Visium HD spatial transcriptomics data from the same tissue section, with the goal of harmonizing platforms and refining cell-type annotations using spatial scTriangulate.

#### Goal
To identify the most reliable cell-type annotations by:
- Aligning matched Xenium and Visium HD spots using gene expression similarity
- Evaluating cross-platform transcript concordance
- Integrating metadata and RCTD results
- Prepareing input files for spatial scTriangulate to optimize annotations across platforms
![Side-by-side figure](./Alignment_files/figure-gfm/unnamed-chunk-3-1.png)
![Alignment figure](.//Alignment_files/figure-gfm/unnamed-chunk-4-1.png)



#### Script
- [`Rmarkdown`](./Alignment.md)
- [`Rmarkdown with figures (html format)`](https://shunyasanuma.github.io/Spatial_LungChat/Alignment.html)
---

##### Workflow Steps
#### 1. Extracting Singlet Coordinates
To ensure high-confidence cell mapping, only RCTD-classified singlets are used from each platform:
- From Visium2@results$results_df and Xenium2@results$results_df, spots with spot_class == "singlet" are extracted.
- Raw spatial coordinates (`x`, `y`) for those singlet spots are retained from:
  * `Visium2@spatialRNA@coords`
  *  `Xenium2@spatialRNA@coords`

#### 2. Manual Alignment of Xenium to Visium Space
Due to imaging differences between Xenium and Visium HD (e.g., rotation, mirroring), direct spatial overlay is not possible. Orientation was corrected in a previous step.

Manual alignment is performed by:
- Rotating, scaling, and mirroring Xenium spatial coordinates to approximate Visium's spatial coordinate system.
- Resulting aligned coordinates (`X`, `Y`) are stored in `xenium_df` and `visium_df`.
- These aligned coordinates enable more accurate 1:1 spatial comparisons across platforms.

#### 3. Combining Raw + Aligned Spatial Coordinates
Both raw (`x`, `y`) and normalized (`X`, `Y`) spatial positions are included in the final metadata table (`combined_df`) for each matched spot:
- `X_xenium_raw`, `Y_xenium_raw`, `X_visium_raw`, `Y_visium_raw`
- `X_xenium_norm`, `Y_xenium_norm`, `X_visium_norm`, `Y_visium_norm`
- This allows downstream tools (`scTriangulate`) to use either raw or aligned coordinates for spatial integration and visualization.




