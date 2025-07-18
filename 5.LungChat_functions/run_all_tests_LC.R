# ===================================================================
# run_all_tests.R (v3 - for JSON output)
#
# Tests functions from Spatial_functions.R, captures the JSON
# output, and verifies that files were created.
# ===================================================================


# --- 1. SETUP: Load Libraries and Functions ---
cat("--- 1. SETUP: Loading Libraries and Sourcing Functions ---\n")
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scales);
  library(SingleCellExperiment); library(SummarizedExperiment); library(mclust);
  library(rlang); library(viridis); library(ggrepel); library(knitr);
  library(RColorBrewer); library(tibble); library(CellChat);
  library(ComplexHeatmap); library(circlize); library(Seurat);
  library(httr); library(jsonlite); library(anndata); library(uuid); library(tools);
})

# IMPORTANT: Update this path to the correct location of your function file
source("/Users/rajlq7/Downloads/LungChat/LungMAP_scExplore_dev/spatial_dev/Spatial_functions_LC.R")

# Create the output directory and set it as the working directory
output_dir <- "/Users/rajlq7/Downloads/LungChat/LungMAP_scExplore_dev/spatial_dev/test_scripts/test_outputs"
if (dir.exists(output_dir)) unlink(output_dir, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE)
setwd(output_dir) # Set working directory so all files are saved here
cat("✅ Setup complete. Outputs will be saved to:", getwd(), "\n")


# --- 2. DATA LOADING & CORRECTION ---
cat("\n--- 2. DATA LOADING & CORRECTION ---\n")
h5ad_path <- "/Users/rajlq7/Downloads/LungChat/LungMAP_scExplore_dev/spatial_dev/data/sctriangulate.h5ad"
rds_path <- "/Users/rajlq7/Downloads/LungChat/LungMAP_scExplore_dev/spatial_dev/data/CellChat_Squidpy.rds"

adata <- read_h5ad(h5ad_path)
sce_obj <- SingleCellExperiment(
  assays = list(X = t(adata$X)),
  colData = as.data.frame(adata$obs),
  rowData = adata$var
)
cat("✅ Loaded sctriangulate.h5ad and created sce_obj object.\n")

CellChat_Squidpy <- readRDS(rds_path)
cat("✅ Loaded CellChat_Squidpy.rds object.\n")


# --- 3. PRE-PROCESSING (FOR CELLCHAT) ---
cat("\n--- 3. PRE-PROCESSING FOR CELLCHAT ---\n")
groupSizes <- list()
for (annot in names(CellChat_Squidpy$CellChat)) {
  groupSizes[[annot]] <- list()
  for (affect in names(CellChat_Squidpy$CellChat[[annot]])) {
    obj <- CellChat_Squidpy$CellChat[[annot]][[affect]]
    idents_order <- rownames(obj@net$count)
    groupSizes[[annot]][[affect]] <- as.numeric(table(factor(obj@idents, levels = idents_order)))
    CellChat_Squidpy$CellChat[[annot]][[affect]] <- netAnalysis_computeCentrality(obj, slot.name = "netP")
  }
}
cat("✅ CellChat pre-processing completed.\n")


# --- 4. RUN ALL TESTS ---
cat("\n--- 4. EXECUTING ALL TESTS ---\n")
test_results <- list()
# Updated test runner to capture and print the JSON output
run_test <- function(test_name, test_code) {
  cat(paste0("\n▸ Testing: ", test_name, "...\n"))
  tryCatch({
    json_output <- eval(parse(text=test_code))
    cat("  Function returned JSON:\n")
    cat(json_output, "\n")
    test_results[[test_name]] <<- "✅ PASSED"
  }, error = function(e) {
    test_results[[test_name]] <<- paste("❌ FAILED:", gsub("\n", " ", e$message))
    cat(paste0("  ERROR in ", test_name, ": ", e$message, "\n"))
  })
}

# --- Visualization & Analysis Functions ---

run_test("plot_spatial", '
  plot_spatial(sce_obj, annotation = "pruned", highlight = c("TNiche@T6", "TNiche@T7"))
')
run_test("plot_umap", '
  plot_umap(sce_obj, platform = "Xenium", annotation = "Final_CT", highlight = c("AT1", "AT2"))
')
run_test("plot_gene_spatial", '
  plot_gene_spatial(sce_obj, gene = "COL1A1", platform = "Xenium", point_size = 0.4)
')
run_test("rank_cross_platform_gene_correlation", '
  rank_cross_platform_gene_correlation(sce_obj, top_genes = 5)
')
run_test("plot_cross_platform_gene_correlation", '
  plot_cross_platform_gene_correlation(sce_obj, gene_name = "COL1A1", y_range = c(0.8, 1))
')
run_test("generate_stability_report", '
  clusters_to_process <- list(Final_CT = c("CD4+_T-cells", "CD8+_T-cells"))
  generate_stability_report(sce = sce_obj, celltype_labels = clusters_to_process)
')
run_test("rank_marker_specificity", '
  rank_marker_specificity(sce_obj, gene_name = "COL1A1", cluster_cols = c("Final_CT", "TNiche"), top_n = 5)
')
run_test("find_celltype_markers", '
  find_celltype_markers(sce = sce_obj, celltype = "T6", annotation_col = "TNiche", platform = "Xenium")
')
run_test("generate_marker_report", '
  generate_marker_report(sce_obj, genes = c("COL1A1", "AGER", "SFTPC"), annotation_cols = c("pruned", "Final_CT"))
')
run_test("plot_ARI_dotplot", '
  plot_ARI_dotplot(sce_obj, annotation1 = "Xenium_RCTD_LungMap_ref", annotation2 = "Visium_RCTD_LungMap_ref")
')
run_test("plot_cluster_composition_groupedbar", '
  plot_cluster_composition_groupedbar(sce_obj, annotation = c("Visium_RCTD_LungMap_ref", "Xenium_RCTD_LungMap_ref"))
')
run_test("visualize_matched_pairs_scatter", '
  visualize_matched_pairs_scatter(sce_obj, gene = "COL1A1", show_correlation = TRUE)
')
run_test("generate_de_volcano_plot", '
  generate_de_volcano_plot(obj_list = CellChat_Squidpy[["CellChat"]], 
                           annotation = "TNiche", 
                           group_1 = "More_Affected", 
                           group_2 = "Unaffected", 
                           celltype_filter = "T6", 
                           label_top_n = 10)
')

# --- CellChat & Squidpy Functions ---

run_test("netVisual_circle", '
  weights <- groupSizes[["Final_CT"]][["Unaffected"]]
  net <- CellChat_Squidpy[["CellChat"]][["Final_CT"]][["Unaffected"]]@net$count
  netVisual_circle(net = net, vertex.weight = weights, vertex.label.cex = 0.7, remove.isolate = TRUE)
')
run_test("netVisual_heatmap", '
  netVisual_heatmap(
    object = CellChat_Squidpy[["CellChat"]][["Final_CT"]][["More_Affected"]],
    measure = "weight", color.heatmap = "Reds", title.name = "Xenium IPF"
  )
')
run_test("netVisual_aggregate", '
  netVisual_aggregate(
    object = CellChat_Squidpy$CellChat$Final_CT$More_Affected,
    signaling = "UGRP1", layout = "circle", remove.isolate = TRUE
  )
')
run_test("subsetCommunication", '
  subsetCommunication(object = CellChat_Squidpy$CellChat$Final_CT$More_Affected, signaling = "UGRP1")
')
run_test("create_enrichment_heatmap", '
  create_enrichment_heatmap(
    obj = CellChat_Squidpy, group = "More_Affected", annotation = "Final_CT",
    clustering = TRUE, base_size = 9
  )
')
run_test("analyze_and_plot_neighbors", '
  analyze_and_plot_neighbors(
    obj = CellChat_Squidpy, target_cell_type = "AT1", group = "More_Affected",
    annotation = "Final_CT"
  )
')
run_test("get_toppcell_enrichment", '
  genes <- c("SFTPC", "SFTPB", "SFTPA1", "LAMP3")
  get_toppcell_enrichment(genes)
')


# --- 5. FINAL SUMMARY ---
cat("\n\n--- 5. FINAL TEST SUMMARY ---\n")
cat(paste(rep("=", 60), collapse=""), "\n")
sorted_names <- sort(names(test_results))
for (test_name in sorted_names) {
  cat(sprintf("%-40s %s\n", test_name, test_results[[test_name]]))
}
cat(paste(rep("=", 60), collapse=""), "\n")

# --- 6. VERIFYING OUTPUT FILES ---
cat("\n--- 6. VERIFYING OUTPUT FILES ---\n")
generated_files <- list.files(".", recursive = TRUE)
if (length(generated_files) > 0) {
  cat(paste0("✅ Found ", length(generated_files), " files in the output directory:\n"))
  for(f in generated_files) {cat(paste0("  - ", f, "\n"))}
} else {
  cat("❌ No output files were generated.\n")
}

cat("\n✅✅✅ All tests completed! ✅✅✅\n")