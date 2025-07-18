# This script loads and pre-processes all data, saving it into a single
# .RData file for fast loading at container runtime.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(anndata)
  library(CellChat)
})

# --- 1. Load Raw Data ---
# These paths are inside the temporary Docker build environment
adata <- read_h5ad("/Users/rajlq7/Desktop/test/data/sctriangulate.h5ad")
sce_obj <- SingleCellExperiment(
  assays = list(X = t(adata$X)),
  colData = as.data.frame(adata$obs),
  rowData = adata$var
)
CellChat_Squidpy <- readRDS("/Users/rajlq7/Desktop/test/data/CellChat_Squidpy.rds")

# --- 2. Perform CellChat Pre-processing ---
groupSizes <- list()
for (annot in names(CellChat_Squidpy$CellChat)) {
  groupSizes[[annot]] <- list()
  for (affect in names(CellChat_Squidpy$CellChat[[annot]])) {
    obj <- CellChat_Squidpy$CellChat[[annot]][[affect]]
    idents_order <- rownames(obj@net$count)
    groupSizes[[annot]][[affect]] <- as.numeric(table(factor(obj@idents, levels = idents_order)))
    withCallingHandlers({
        CellChat_Squidpy$CellChat[[annot]][[affect]] <- netAnalysis_computeCentrality(obj, slot.name = "netP")
      }, warning = function(w) {
        if (grepl("The edge weights used for computing centrality measures are from netP", w$message)) {
          invokeRestart("muffleWarning")
        }
      })
  }
}

# --- 3. Save Processed Objects ---
# Save everything into one file in the /app directory for the API to load
save(sce_obj, CellChat_Squidpy, groupSizes, file = "/Users/rajlq7/Desktop/test/processed_data.RData")