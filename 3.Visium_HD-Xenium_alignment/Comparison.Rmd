```{r}
# Bioconductor packages (load these first)
library(SummarizedExperiment)
library(MatrixGenerics)
library(SingleCellExperiment)
# Now load Seurat and friends
library(Seurat)
library(ggplot2)
library(dplyr)
library(arrow)
library(hdf5r)
library(spacexr)
library(SeuratWrappers)
library(mclust)
library(zellkonverter)
library(SeuratDisk)
library(Matrix)
```

#### Load RCTD Objects

```{r}
load("/data/salomonis2/LabFiles/Shunya_Asanuma/Spatial/LungChat/Output/Visium_HD/GSM8509590/02.RCTD/UMI_40/Visium2.RData")
load("/data/salomonis2/LabFiles/Shunya_Asanuma/Spatial/LungChat/Output/Xenium/GSE250346/02.RCTD/UMI_40/Xenium2.RData")
load("/data/salomonis2/LabFiles/Shunya_Asanuma/Spatial/LungChat/Output/Visium_HD/GSM8509590/02.RCTD/UMI_40/Visium2_author_based.RData")
Manuscript <- readRDS("/data/salomonis2/LabFiles/Shunya_Asanuma/Spatial/LungChat/Input/GSE250346/GSE250346_Seurat_GSE250346_CORRECTED_SEE_RDS_README_082024.rds")

```

#### Extract Singlet Coordinates

```{r}
df_visium <- Visium2@results$results_df %>% filter(spot_class == "singlet")
df_xenium <- Xenium2@results$results_df %>% filter(spot_class == "singlet")

df_visium$X <- Visium2@spatialRNA@coords[rownames(df_visium), "x"]
df_visium$Y <- Visium2@spatialRNA@coords[rownames(df_visium), "y"]
df_xenium$X <- Xenium2@spatialRNA@coords[rownames(df_xenium), "x"]
df_xenium$Y <- Xenium2@spatialRNA@coords[rownames(df_xenium), "y"]

```

### Side-by-Side Plot: Visium vs. Xenium Singlets

```{r}
# Add platform labels
df_visium$Platform <- "Visium"
df_xenium$Platform <- "Xenium"

# Shift Xenium to the right for side-by-side display
df_xenium_side_by_side <- df_xenium %>%
  mutate(X = X + 1.2)

# Combine
df_combined_side_by_side <- bind_rows(df_visium, df_xenium_side_by_side)

# Plot
ggplot(df_combined_side_by_side, aes(x = X, y = Y, color = Platform)) +
  geom_point(size = 0.4) +
  scale_color_manual(values = c("Visium" = "blue", "Xenium" = "red")) +
  coord_fixed() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Side-by-Side View of Visium2 and Xenium2 Singlets",
    x = "Normalized X", y = "Normalized Y", color = "Platform"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "gray80"),
    legend.position = "right"
  ) +
  annotate("text", x = 0.5, y = 1.05, label = "Visium", size = 5, fontface = "bold") +
  annotate("text", x = 1.7, y = 1.05, label = "Xenium", size = 5, fontface = "bold")

```

#### Manual Alignment of Xenium to Visium & Combine and Plot

```{r}
# Add platform labels
df_visium$Platform <- "Visium"
df_xenium$Platform <- "Xenium"

# Copy and transform Xenium for alignment
df_xenium_aligned <- df_xenium
df_xenium_aligned$X <- df_xenium_aligned$X * 1.17 - 0.024
df_xenium_aligned$Y <- df_xenium_aligned$Y * 1.02 - 0.02

# Combine
df_combined_aligned <- bind_rows(df_visium, df_xenium_aligned)

# Plot
p_alignment <- ggplot(df_combined_aligned, aes(x = X, y = Y, color = Platform)) +
  geom_point(size = 0.4) +
  scale_color_manual(values = c("Visium" = "blue", "Xenium" = "red")) +
  coord_fixed() +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(color = "gray80", size = 0.4),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Aligned Visium2 (blue) and Xenium2 (red) - UMI > 40",
    color = "Platform"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))

print(p_alignment)
```

#### Nearest Neighbor Mapping: Xenium → Visium

```{r}
# -------------------------------
# Libraries
# -------------------------------
library(FNN)
library(dplyr)

# -------------------------------
# Step 1: Prepare aligned coordinates
# -------------------------------
# Use manually aligned Xenium coordinates for accurate 1-to-1 spatial matching
xenium_df <- df_xenium_aligned
visium_df <- df_visium  # Visium already normalized

# Extract 2D coordinates
xenium_coords <- xenium_df[, c("X", "Y")]
visium_coords <- visium_df[, c("X", "Y")]

# -------------------------------
# Step 2: Find nearest Visium neighbor for each Xenium spot
# -------------------------------
nn_assignment <- get.knnx(visium_coords, xenium_coords, k = 1)
assigned_indices <- nn_assignment$nn.index[, 1]

# Map Xenium barcodes to matched Visium barcodes
assigned_visium_spot_id <- rownames(visium_coords)[assigned_indices]
names(assigned_visium_spot_id) <- rownames(xenium_coords)

# -------------------------------
# Step 3: Save alignment map as data.frame
# -------------------------------
alignment_df <- data.frame(
  Xenium_ID = names(assigned_visium_spot_id),
  Visium_ID = unname(assigned_visium_spot_id),
  stringsAsFactors = FALSE
)


# -----------------------------------------
# Step 4: Save to CSV
# -----------------------------------------
# write.csv(
#   alignment_df,
#   "/data/salomonis2/LabFiles/Shunya_Asanuma/Spatial/LungChat/Output/Comparison/GSM8509590_GSE250346/Alignment/All/Tissue2_xenium_visium_alignment.csv",
#   row.names = FALSE
# )

```

#### Alignment Quality Evaluation

##### Plots scatterplots (log1p-transformed) for these top 3 genes.

```{r}
# -------------------------------
# Step 1: Extract raw count matrices & Get intersecting genes
# -------------------------------
expr_xenium <- Xenium2@spatialRNA@counts
expr_visium <- Visium2@spatialRNA@counts

# -------------------------------
# Step 1: Get intersecting genes
# -------------------------------
common_genes <- intersect(rownames(expr_xenium), rownames(expr_visium))

# -------------------------------
# Step 2: Compute Pearson correlation for all shared genes
# -------------------------------
cor_values <- sapply(common_genes, function(gene) {
  xen_vec <- expr_xenium[gene, names(assigned_visium_spot_id)]
  vis_vec <- expr_visium[gene, assigned_visium_spot_id]
  cor(log1p(as.numeric(xen_vec)), log1p(as.numeric(vis_vec)), method = "pearson")
})

# -------------------------------
# Step 3: Select top 3 genes by correlation
# -------------------------------
top_genes <- names(sort(cor_values, decreasing = TRUE))[1:3]

# -------------------------------
# Step 4: Generate scatterplots
# -------------------------------
plot_list <- list()

for (gene in top_genes) {
  xenium_expr_vec <- expr_xenium[gene, names(assigned_visium_spot_id)]
  visium_expr_vec <- expr_visium[gene, assigned_visium_spot_id]

  df_gene <- data.frame(
    Xenium_log = log1p(as.numeric(xenium_expr_vec)),
    Visium_log = log1p(as.numeric(visium_expr_vec)),
    Gene = gene
  )

  pearson_cor <- cor(df_gene$Xenium_log, df_gene$Visium_log, method = "pearson")

  p <- ggplot(df_gene, aes(x = Visium_log, y = Xenium_log)) +
    geom_point(alpha = 0.4, color = "blue", size = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste0(gene, " (log1p): r = ", round(pearson_cor, 2)),
      x = "Visium Expression (log1p)",
      y = "Xenium Expression (log1p)"
    )

  plot_list[[gene]] <- p
}

# -------------------------------
# Step 5: Display plots
# -------------------------------
library(cowplot)
plot_grid(plotlist = plot_list, ncol = 2)

```

##### Line plot showing the Pearson correlation (y-axis) versus the number of top aligned spot pairs (x-axis) ranked by absolute difference in expression for the most correlated gene (top_gene).

```{r}
library(ggplot2)

# Select top gene
top_gene <- top_genes[1]

# Extract aligned expression
xen_vec <- expr_xenium[top_gene, names(assigned_visium_spot_id)]
vis_vec <- expr_visium[top_gene, assigned_visium_spot_id]

# Log1p transform
xen_log <- log1p(as.numeric(xen_vec))
vis_log <- log1p(as.numeric(vis_vec))

# Absolute difference
abs_diff <- abs(xen_log - vis_log)

# Sort by similarity (smallest error first)
sorted_idx <- order(abs_diff)
xen_log_sorted <- xen_log[sorted_idx]
vis_log_sorted <- vis_log[sorted_idx]

# Evaluate correlation at increasing N
n_spots <- length(xen_log_sorted)
step_size <- 50  # you can tune this
top_n_values <- seq(step_size, n_spots, by = step_size)
cor_vals <- sapply(top_n_values, function(n) {
  cor(xen_log_sorted[1:n], vis_log_sorted[1:n], method = "pearson")
})

# Build dataframe for plotting
cor_df <- data.frame(
  Top_N_Spots = top_n_values,
  Pearson_Correlation = cor_vals
)

# Plot
ggplot(cor_df, aes(x = Top_N_Spots, y = Pearson_Correlation)) +
  geom_line(size = 1, color = "darkred") +
  geom_point(size = 2, color = "black") +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(
    title = paste("Correlation vs. Top N (", top_gene, ")", sep = ""),
    x = "Top N Paired Spots (smallest abs diff)",
    y = "Pearson Correlation"
  ) +
  theme_minimal(base_size = 14)

```

### **The top 7,500 spots achieved a Pearson correlation of 0.9 and were therefore used for downstream analysis (e.g., Spatial scTriangulate).**


##### Build Combined Data Frame for Aligned Xenium–Visium Singlets
```{r}
# -------------------------------
# Step 0: Select top N most similar spot pairs
# -------------------------------
top_n <- 7500
top_idx <- order(abs_diff)[1:top_n]  # abs_diff should be precomputed (e.g., 1 - correlation or Euclidean)

xenium_ids_top <- xenium_ids[top_idx]
visium_ids_top <- visium_ids[top_idx]

# -------------------------------
# Step 1: Map manuscript short IDs to full rownames
# -------------------------------
manuscript_map <- Manuscript@meta.data
cell_id_map <- setNames(rownames(manuscript_map), manuscript_map$cell_id)

# Only keep Xenium IDs present in the Manuscript metadata
matched_ids <- intersect(xenium_ids_top, names(cell_id_map))
mapped_row_ids <- cell_id_map[matched_ids]  # Full Xenium rownames

# -------------------------------
# Step 2: Assemble combined metadata table
# -------------------------------
combined_df <- data.frame(
  Xenium_ID = matched_ids,
  Visium_ID = visium_ids_top[match(matched_ids, xenium_ids_top)],

  # Manuscript metadata
  CNiche = manuscript_map[mapped_row_ids, "CNiche"],
  TNiche = manuscript_map[mapped_row_ids, "TNiche"],
  Final_lineage = manuscript_map[mapped_row_ids, "final_lineage"],
  Final_CT = manuscript_map[mapped_row_ids, "final_CT"],
  Disease_status = manuscript_map[mapped_row_ids, "disease_status"],
  Sample_type = manuscript_map[mapped_row_ids, "sample_type"],
  Sample_affect = manuscript_map[mapped_row_ids, "sample_affect"],

  # RCTD annotations
  Visium_RCTD_LungMap_ref = Visium2@results$results_df[visium_ids_top[match(matched_ids, xenium_ids_top)], "first_type"],
  Visium_RCTD_GSE250346_based = Visium2_author_based@results$results_df[visium_ids_top[match(matched_ids, xenium_ids_top)], "first_type"],
  Xenium_RCTD_LungMap_ref = Xenium2@results$results_df[xenium_ids_top[match(matched_ids, xenium_ids_top)], "first_type"]
)

# -------------------------------
# Step 3: Add spatial coordinates (raw)
# -------------------------------
combined_df$X_xenium_raw <- Xenium2@spatialRNA@coords[combined_df$Xenium_ID, "x"]
combined_df$Y_xenium_raw <- Xenium2@spatialRNA@coords[combined_df$Xenium_ID, "y"]
combined_df$X_visium_raw <- Visium2@spatialRNA@coords[combined_df$Visium_ID, "x"]
combined_df$Y_visium_raw <- Visium2@spatialRNA@coords[combined_df$Visium_ID, "y"]

# -------------------------------
# Step 4: Add normalized coordinates (from aligned output)
# -------------------------------
combined_df$X_xenium_norm <- xenium_df[combined_df$Xenium_ID, "X"]
combined_df$Y_xenium_norm <- xenium_df[combined_df$Xenium_ID, "Y"]
combined_df$X_visium_norm <- visium_df[combined_df$Visium_ID, "X"]
combined_df$Y_visium_norm <- visium_df[combined_df$Visium_ID, "Y"]

# -------------------------------
# Step 5: Final formatting
# -------------------------------
rownames(combined_df) <- NULL
print(dim(combined_df))
head(combined_df)
```






```{r}
# ------------ Step 1: Define genes and validate presence ------------
genes_to_plot <- "COL1A1"
stopifnot(all(genes_to_plot %in% rownames(expr_xenium)))
stopifnot(all(genes_to_plot %in% rownames(expr_visium)))

# ------------ Step 2: Select top 7500 most similar pairs based on expression difference ------------
top_n <- 7500
top_idx <- order(abs_diff)[1:top_n]
top_xenium_ids <- xenium_ids[top_idx]
top_visium_ids <- visium_ids[top_idx]

# ------------ Step 3: Filter valid IDs present in expression matrices ------------
valid_mask <- top_xenium_ids %in% colnames(expr_xenium) &
    top_visium_ids %in% colnames(expr_visium)
xenium_ids_subset <- top_xenium_ids[valid_mask]
visium_ids_subset <- top_visium_ids[valid_mask]

# ------------ Step 4: Plot expression correlations ------------
plot_list <- list()

for (gene in genes_to_plot) {
    # Extract expression vectors
    xenium_expr_vec <- expr_xenium[gene, xenium_ids_subset]
    visium_expr_vec <- expr_visium[gene, visium_ids_subset]
    
    # Log-transform and combine into a dataframe
    df_gene <- data.frame(
        Xenium_log = log1p(as.numeric(xenium_expr_vec)),
        Visium_log = log1p(as.numeric(visium_expr_vec)),
        Gene = gene
    )
    
    # Compute Pearson correlation
    pearson_cor <- cor(df_gene$Xenium_log, df_gene$Visium_log, method = "pearson")
    
    # Create scatter plot
    p <- ggplot(df_gene, aes(x = Visium_log, y = Xenium_log)) +
        geom_point(alpha = 0.4, color = "blue", size = 1) +
        geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
        theme_minimal(base_size = 14) +
        labs(
            title = paste0(gene, " (log1p): r = ", round(pearson_cor, 2)),
            x = "Visium Expression (log1p)",
            y = "Xenium Expression (log1p)"
        )
    
    plot_list[[gene]] <- p
}

# ------------ Step 5: Arrange side-by-side ------------
plot_grid(plotlist = plot_list)
```

```{r}
# --- Step 0: Get top 7500 high-correlation matched pairs ---
top_n <- 7500
top_idx <- order(abs_diff)[1:top_n]
xenium_ids_top <- xenium_ids[top_idx]
visium_ids_top <- visium_ids[top_idx]

# --- Step 1: Convert Manuscript short cell IDs to full rownames ---
manuscript_map <- Manuscript@meta.data
cell_id_map <- setNames(rownames(manuscript_map), manuscript_map$cell_id)

# --- Step 2: Only keep matched Xenium cells that exist in Manuscript meta.data ---
matched_ids <- intersect(xenium_ids_top, names(cell_id_map))
mapped_row_ids <- cell_id_map[matched_ids]

# --- Step 3: Build combined dataframe ---
combined_df <- data.frame(
  Xenium_ID = matched_ids,
  Visium_ID = visium_ids_top[match(matched_ids, xenium_ids_top)],

  # From Manuscript metadata
  CNiche = manuscript_map[mapped_row_ids, "CNiche"],
  TNiche = manuscript_map[mapped_row_ids, "TNiche"],
  Final_lineage = manuscript_map[mapped_row_ids, "final_lineage"],
  Final_CT = manuscript_map[mapped_row_ids, "final_CT"],
  Disease_status = manuscript_map[mapped_row_ids, "disease_status"],
  Sample_type = manuscript_map[mapped_row_ids, "sample_type"],
  Sample_affect = manuscript_map[mapped_row_ids, "sample_affect"],

  # RCTD annotations
  Visium_RCTD_LungMap_ref = Visium2@results$results_df[visium_ids_top[match(matched_ids, xenium_ids_top)], "first_type"],
  Visium_RCTD_GSE250346_based = Visium2_author_based@results$results_df[visium_ids_top[match(matched_ids, xenium_ids_top)], "first_type"],
  Xenium_RCTD_LungMap_ref = Xenium2@results$results_df[xenium_ids_top[match(matched_ids, xenium_ids_top)], "first_type"]
)

# --- Step 4: Add raw spatial coordinates ---
combined_df$X_xenium_raw <- Xenium2@spatialRNA@coords[combined_df$Xenium_ID, "x"]
combined_df$Y_xenium_raw <- Xenium2@spatialRNA@coords[combined_df$Xenium_ID, "y"]
combined_df$X_visium_raw <- Visium2@spatialRNA@coords[combined_df$Visium_ID, "x"]
combined_df$Y_visium_raw <- Visium2@spatialRNA@coords[combined_df$Visium_ID, "y"]

# --- Step 5: Add normalized coordinates (from aligned output) ---
combined_df$X_xenium_norm <- xenium_df[combined_df$Xenium_ID, "X"]
combined_df$Y_xenium_norm <- xenium_df[combined_df$Xenium_ID, "Y"]
combined_df$X_visium_norm <- visium_df[combined_df$Visium_ID, "X"]
combined_df$Y_visium_norm <- visium_df[combined_df$Visium_ID, "Y"]

# --- Step 6: Tidy up ---
rownames(combined_df) <- NULL
head(combined_df)
dim(combined_df)

```

```{r}
# write.csv(combined_df, "/data/salomonis2/LabFiles/Shunya_Asanuma/Spatial/LungChat/Output/Comparison/GSM8509590_GSE250346/Alignment/combined_metadata_top7500.csv", row.names = FALSE)

```

## Filter & Rename Gene Expression Matrices

```{r}
# --- Step 0: Extract barcodes from combined_df ---
top_xenium_ids <- combined_df$Xenium_ID
top_visium_ids <- combined_df$Visium_ID

# --- Step 1: Subset and rename expression matrices ---
# Visium
expr_visium_subset <- expr_visium[, top_visium_ids, drop = FALSE]
rownames(expr_visium_subset) <- paste0(rownames(expr_visium_subset), "_visium")

# Xenium
expr_xenium_subset <- expr_xenium[, top_xenium_ids, drop = FALSE]
rownames(expr_xenium_subset) <- paste0(rownames(expr_xenium_subset), "_xenium")

# --- Step 2: Transpose to spots × genes ---
expr_visium_t <- t(expr_visium_subset)  # dim: spots × genes_visium
expr_xenium_t <- t(expr_xenium_subset)  # dim: spots × genes_xenium

# --- Step 3: Combine expression matrices by column (gene-wise) ---
combined_expr <- cbind(expr_xenium_t, expr_visium_t)  # dim: 7061 × (genes_xenium + genes_visium)

# --- Step 4: Set rownames to Xenium IDs (as unique spot IDs) ---
rownames(combined_expr) <- top_xenium_ids  # or: rownames(combined_df)

# --- Step 5: Optional sanity check ---
dim(combined_expr)  # Expect: 7061 × 4827 (or actual count based on gene intersection)


```

```{r}
# Save as CSV
# write.csv(combined_expr, "/data/salomonis2/LabFiles/Shunya_Asanuma/Spatial/LungChat/Output/Comparison/GSM8509590_GSE250346/Alignment/combined_expression_matrix.csv", quote = FALSE)
```

########################### 

```{r}
plot_side_by_side <- function(Visium, Xenium, tissue_name, color_map) {
  # Step 1: Extract singlets
  df_visium <- Visium@results$results_df %>% filter(spot_class == "singlet")
  df_xenium <- Xenium@results$results_df %>% filter(spot_class == "singlet")
  
  # Step 2: Coordinates
  df_visium$X <- Visium@spatialRNA@coords[rownames(df_visium), "x"]
  df_visium$Y <- Visium@spatialRNA@coords[rownames(df_visium), "y"]
  df_xenium$X <- Xenium@spatialRNA@coords[rownames(df_xenium), "x"]
  df_xenium$Y <- Xenium@spatialRNA@coords[rownames(df_xenium), "y"]
  
  # Step 3: Rotate Visium (90° clockwise)
  df_visium <- df_visium %>%
    mutate(tempX = X, X = Y, Y = -tempX) %>%
    select(-tempX)
  
  # Normalize both
  normalize_coords <- function(df) {
    df$X <- (df$X - min(df$X)) / (max(df$X) - min(df$X))
    df$Y <- (df$Y - min(df$Y)) / (max(df$Y) - min(df$Y))
    return(df)
  }
  df_visium <- normalize_coords(df_visium)
  df_xenium <- normalize_coords(df_xenium)
  
  # Shift Xenium to the right
  df_xenium$X <- df_xenium$X + 1.2
  
  # Step 4: Label platform
  df_visium$Platform <- "Visium"
  df_xenium$Platform <- "Xenium"
  
  # Step 5: Merge
  df_all <- bind_rows(df_visium, df_xenium)
  
  # Step 6: Plot
  p <- ggplot(df_all, aes(x = X, y = Y, color = first_type)) +
    geom_point(size = 0.3) +
    coord_fixed() +
    theme_void(base_size = 14) +
    scale_color_manual(values = color_map) +
    labs(title = paste("RCTD Predicted Cell Types (", tissue_name, ")", sep = ""), color = "Cell Type") +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    ) +
    annotate("text", x = 0.5, y = 1.05, label = "Visium", size = 5, fontface = "bold") +
    annotate("text", x = 1.7, y = 1.05, label = "Xenium", size = 5, fontface = "bold") +
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  return(p)
}

```

```{r}
# Reuse the color palette from Tissue2
all_types <- sort(unique(c(
  Visium1@results$results_df$first_type,
  Xenium1@results$results_df$first_type,
  Visium2@results$results_df$first_type,
  Xenium2@results$results_df$first_type,
  Visium3@results$results_df$first_type,
  Xenium3@results$results_df$first_type,
  Visium4@results$results_df$first_type,
  Xenium4@results$results_df$first_type
)))
color_map <- Seurat:::DiscretePalette(length(all_types))
names(color_map) <- all_types

# Plot all tissues side-by-side
plot_side_by_side(Visium1, Xenium1, "Tissue1", color_map)
plot_side_by_side(Visium2, Xenium2, "Tissue2", color_map)
plot_side_by_side(Visium3, Xenium3, "Tissue3", color_map)
plot_side_by_side(Visium4, Xenium4, "Tissue4", color_map)

```

```{r}
plot_bar_chart <- function(Visium, Xenium, tissue_name, color_map) {
  # Step 1: Extract singlets
  df_visium <- Visium@results$results_df %>% filter(spot_class == "singlet")
  df_xenium <- Xenium@results$results_df %>% filter(spot_class == "singlet")
  
  # Step 2: Add platform labels
  df_visium$Platform <- "Visium"
  df_xenium$Platform <- "Xenium"
  
  # Step 3: Merge
  df_all <- bind_rows(df_visium, df_xenium)
  
  # Step 4: Summarize proportions
  df_summary <- df_all %>%
    group_by(Platform, first_type) %>%
    summarize(n = n()) %>%
    group_by(Platform) %>%
    mutate(Fraction = n / sum(n))
  
  # Step 5: Bar plot
  p <- ggplot(df_summary, aes(x = Platform, y = Fraction, fill = first_type)) +
    geom_bar(stat = "identity", position = "fill", width = 0.6) +
    scale_fill_manual(values = color_map) +
    labs(title = paste("Cell Type Proportions -", tissue_name),
         y = "Fraction", x = NULL, fill = "Cell Type") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 12),
      legend.position = "right"
    )
  
  return(p)
}

```

```{r}
# Again reusing the color_map created earlier

# Bar plots for all tissues
plot_bar_chart(Visium1, Xenium1, "Tissue1", color_map)
plot_bar_chart(Visium2, Xenium2, "Tissue2", color_map)
plot_bar_chart(Visium3, Xenium3, "Tissue3", color_map)
plot_bar_chart(Visium4, Xenium4, "Tissue4", color_map)

```

# Tissue 2 (TEST)

```{r}
# ------------ Step 1: Extract singlets ------------
df_visium <- Visium2@results$results_df %>% filter(spot_class == "singlet")
df_xenium <- Xenium2@results$results_df %>% filter(spot_class == "singlet")

# ------------ Step 2: Pull coordinates ------------
df_visium$X <- Visium2@spatialRNA@coords[rownames(df_visium), "x"]
df_visium$Y <- Visium2@spatialRNA@coords[rownames(df_visium), "y"]
df_xenium$X <- Xenium2@spatialRNA@coords[rownames(df_xenium), "x"]
df_xenium$Y <- Xenium2@spatialRNA@coords[rownames(df_xenium), "y"]

# ------------ Step 3: Rotate Visium by 90° clockwise ------------
rotate90_clockwise <- function(df) {
  df_rot <- df
  df_rot$X <- df$Y
  df_rot$Y <- -df$X
  return(df_rot)
}
df_visium <- rotate90_clockwise(df_visium)

# Normalize both to [0, 1]
normalize_coords <- function(df) {
  df$X <- (df$X - min(df$X)) / (max(df$X) - min(df$X))
  df$Y <- (df$Y - min(df$Y)) / (max(df$Y) - min(df$Y))
  return(df)
}
df_visium <- normalize_coords(df_visium)
df_xenium <- normalize_coords(df_xenium)

# Shift Xenium for side-by-side view
df_xenium$X <- df_xenium$X + 1.2

# ------------ Step 4: Add platform labels ------------
df_visium$Platform <- "Visium"
df_xenium$Platform <- "Xenium"

# ------------ Step 5: Merge and assign consistent colors ------------
df_all <- bind_rows(df_visium, df_xenium)
cell_types <- sort(unique(df_all$first_type))
colors <- Seurat:::DiscretePalette(length(cell_types))
names(colors) <- cell_types

# ------------ Step 6: Plot ------------
ggplot(df_all, aes(x = X, y = Y, color = first_type)) +
  geom_point(size = 0.5) +
  coord_fixed() +
  theme_void(base_size = 14) +
  scale_color_manual(values = colors) +
  labs(title = "RCTD Predicted Cell Types (Tissue2)", color = "Cell Type") +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  ) +
  annotate("text", x = 0.5, y = 1.05, label = "Visium", size = 5, fontface = "bold") +
  annotate("text", x = 1.7, y = 1.05, label = "Xenium", size = 5, fontface = "bold") +
  guides(color = guide_legend(override.aes = list(size = 4)))
```

```{r}
plot_bar_chart_dodge <- function(Visium, Xenium, tissue_name, color_map) {
  # Step 1: Extract singlets
  df_visium <- Visium@results$results_df %>% filter(spot_class == "singlet")
  df_xenium <- Xenium@results$results_df %>% filter(spot_class == "singlet")
  
  # Step 2: Add platform labels
  df_visium$Platform <- "Visium"
  df_xenium$Platform <- "Xenium"
  
  # Step 3: Merge
  df_all <- bind_rows(df_visium, df_xenium)
  
  # Step 4: Summarize counts and fractions
  df_summary <- df_all %>%
    group_by(Platform, first_type) %>%
    summarize(n = n(), .groups = "drop") %>%
    group_by(Platform) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()
  
  # Step 5: Bar plot
  p <- ggplot(df_summary, aes(x = first_type, y = proportion, fill = Platform)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
    scale_fill_manual(values = c("Visium" = "blue", "Xenium" = "red")) +
    labs(
      title = paste("Cell Type Proportions -", tissue_name),
      x = "Cell Type", y = "Proportion of Singlet Spots", fill = "Platform"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  
  return(p)
}

```

```{r}
plot_bar_chart_dodge(Visium1, Xenium1, "Tissue1", color_map)
plot_bar_chart_dodge(Visium2, Xenium2, "Tissue2", color_map)
plot_bar_chart_dodge(Visium3, Xenium3, "Tissue3", color_map)
plot_bar_chart_dodge(Visium4, Xenium4, "Tissue4", color_map)
```

```{r}
# Expression matrices
expr_visium <- Visium2@spatialRNA@counts
expr_xenium <- Xenium2@spatialRNA@counts

genes_of_interest <- c(
  "MFAP5", "SCARA5", "FGF18", "DACH2", "ACTA2",
  "ITGBL1", "PDGFRA", "POSTN", "COL1A1", "IGFBP5",
  "MYH11", "FBNL2", "WNT2", "COL13A1", "TAGLN"
)

# Subset to genes that actually exist
genes_visium <- genes_of_interest[genes_of_interest %in% rownames(expr_visium)]
genes_xenium <- genes_of_interest[genes_of_interest %in% rownames(expr_xenium)]
common_genes <- intersect(genes_visium, genes_xenium)

cat("Common genes to plot:", common_genes, "\n")

library(cowplot)  # for nice multi-panel layout
library(viridis)  # for nice expression colors

plot_gene_expression <- function(counts, coords, genes, title_prefix) {
  plots <- list()
  
  for (g in genes) {
    expr <- counts[g, ]
    df <- data.frame(
      X = coords[colnames(counts), "x"],
      Y = coords[colnames(counts), "y"],
      expr = as.numeric(expr)
    )
    
    # Scale expression for better color dynamic
    df$expr <- pmin(df$expr, quantile(df$expr, 0.99, na.rm = TRUE))  # clip top 1%
    
    p <- ggplot(df, aes(x = X, y = Y, color = expr)) +
      geom_point(size = 0.5) +
      coord_fixed() +
      theme_void() +
      scale_color_viridis(option = "magma", direction = -1) +
      labs(title = paste0(title_prefix, " - ", g), color = "Expr") +
      theme(plot.title = element_text(size = 10, face = "bold"))
    
    plots[[g]] <- p
  }
  
  return(plots)
}

# Visium2
plots_visium2 <- plot_gene_expression(
  counts = expr_visium,
  coords = Visium2@spatialRNA@coords,
  genes = common_genes,
  title_prefix = "Visium2"
)

# Xenium2
plots_xenium2 <- plot_gene_expression(
  counts = expr_xenium,
  coords = Xenium2@spatialRNA@coords,
  genes = common_genes,
  title_prefix = "Xenium2"
)

# Visium2 grid
plot_grid(plotlist = plots_visium2, ncol = 5)

# Xenium2 grid
plot_grid(plotlist = plots_xenium2, ncol = 5)

```

#################### 

```{r}
# Each RCTD object: @results$results_df
df1 <- Visium1@results$results_df %>% mutate(platform = "Visium", tissue = "Tissue1")
df2 <- Xenium1@results$results_df %>% mutate(platform = "Xenium", tissue = "Tissue1")
df3 <- Visium2@results$results_df %>% mutate(platform = "Visium", tissue = "Tissue2")
df4 <- Xenium2@results$results_df %>% mutate(platform = "Xenium", tissue = "Tissue2")

# Combine all into one
combined_df <- bind_rows(df1, df2, df3, df4)

celltype_props <- combined_df %>%
  filter(spot_class == "singlet") %>%
  count(tissue, platform, first_type) %>%
  group_by(tissue, platform) %>%
  mutate(proportion = n / sum(n))

# Make a separate plot for each tissue
unique_tissues <- unique(celltype_props$tissue)

for (t in unique_tissues) {
  df_tissue <- celltype_props %>% filter(tissue == t)
  
  p <- ggplot(df_tissue, aes(x = first_type, y = proportion, fill = platform)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
    labs(
      title = paste("Cell Type Proportions -", t),
      x = "Cell Type", y = "Proportion of Singlet Spots", fill = "Platform"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  print(p)  # Display the plot

  # Optionally: Save each plot as a file
  ggsave(paste0("CellTypeProportions_", t, ".pdf"), plot = p, width = 8, height = 5)
}
```

### Visium HD Spatial

```{r}
# Extract singlet RCTD results
df_visium <- Visium1@results$results_df %>%
  filter(spot_class == "singlet")

# Map spatial coordinates from SpatialRNA
df_visium$X <- Visium1@spatialRNA@coords[rownames(df_visium), "x"]
df_visium$Y <- Visium1@spatialRNA@coords[rownames(df_visium), "y"]

# Plot
ggplot(df_visium, aes(x = X, y = Y, color = first_type)) +
  geom_point(size = 1) +
  coord_fixed() +
  theme_void() +
  labs(title = "Visium HD Tissue1") +
  theme(legend.position = "right") +
  scale_color_manual(values = Seurat:::DiscretePalette(length(unique(df_visium$first_type))))
```

```{r}
df_visium <- Visium2@results$results_df %>%
  filter(spot_class == "singlet")

# Assign X and Y coordinates using rownames
df_visium$X <- Visium2@spatialRNA@coords[rownames(df_visium), "x"]
df_visium$Y <- Visium2@spatialRNA@coords[rownames(df_visium), "y"]

ggplot(df_visium, aes(x = X, y = Y, color = first_type)) +
  geom_point(size = 1) +
  coord_fixed() +
  theme_void() +
  labs(title = "Visium HD Tissue2") +
  theme(legend.position = "right") +
  scale_color_manual(values = Seurat:::DiscretePalette(length(unique(df_visium$first_type))))
```

### Xenium Spatial

```{r}
# Extract results and coordinates for singlets
df <- Xenium1@results$results_df %>%
  filter(spot_class == "singlet")

# Add X/Y spatial coordinates from SpatialRNA
df$X <- Xenium1@spatialRNA@coords[rownames(df), "x"]
df$Y <- Xenium1@spatialRNA@coords[rownames(df), "y"]

# Plot with ggplot
library(ggplot2)
library(Seurat)

ggplot(df, aes(x = X, y = Y, color = first_type)) +
  geom_point(size = 1) +
  coord_fixed() +
  theme_void() +
  labs(title = "Xenium Tissue1") +
  theme(legend.position = "right") +
  scale_color_manual(values = Seurat:::DiscretePalette(length(unique(df$first_type))))


```

```{r}
df$X <- Xenium2@spatialRNA@coords[rownames(df), "x"]
df$Y <- Xenium2@spatialRNA@coords[rownames(df), "y"]

ggplot(df, aes(x = X, y = Y, color = first_type)) +
  geom_point(size = 1) +
  coord_fixed() +
  theme_void() +
  labs(title = "Xenium Tissue2") +
  theme(legend.position = "right") +
  scale_color_manual(values = Seurat:::DiscretePalette(length(unique(df$first_type))))

```

```{r}
library(ggplot2)
library(dplyr)
library(Seurat)

# ------------ Step 1: Extract singlets ------------
df_visium <- Visium2@results$results_df %>% filter(spot_class == "singlet")
df_xenium <- Xenium2@results$results_df %>% filter(spot_class == "singlet")

# ------------ Step 2: Pull coordinates ------------
df_visium$X <- Visium2@spatialRNA@coords[rownames(df_visium), "x"]
df_visium$Y <- Visium2@spatialRNA@coords[rownames(df_visium), "y"]
df_xenium$X <- Xenium2@spatialRNA@coords[rownames(df_xenium), "x"]
df_xenium$Y <- Xenium2@spatialRNA@coords[rownames(df_xenium), "y"]

# ------------ Step 3: Rotate and normalize ------------
rotate_visium <- function(df) {
  df_rot <- df
  df_rot$X <- -df$X
  df_rot$Y <- -df$Y
  return(df_rot)
}

df_visium <- rotate90(df_visium)

# Normalize to [0, 1] and rescale to same range
normalize_coords <- function(df) {
  df$X <- (df$X - min(df$X)) / (max(df$X) - min(df$X))
  df$Y <- (df$Y - min(df$Y)) / (max(df$Y) - min(df$Y))
  return(df)
}
df_visium <- normalize_coords(df_visium)
df_xenium <- normalize_coords(df_xenium)

# Shift side by side
df_visium$X <- df_visium$X
df_xenium$X <- df_xenium$X + 1.2  # spacing between panels

# ------------ Step 4: Add platform labels ------------
df_visium$Platform <- "Visium"
df_xenium$Platform <- "Xenium"

# ------------ Step 5: Merge and map consistent colors ------------
df_all <- bind_rows(df_visium, df_xenium)

cell_types <- sort(unique(df_all$first_type))
colors <- Seurat:::DiscretePalette(length(cell_types))
names(colors) <- cell_types

# ------------ Step 6: Plot ------------
ggplot(df_all, aes(x = X, y = Y, color = first_type)) +
  geom_point(size = 1.2) +
  coord_fixed() +
  theme_void(base_size = 14) +
  scale_color_manual(values = colors) +
  labs(title = "RCTD Predicted Cell Types (Tissue2)", color = "Cell Type") +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  ) +
  annotate("text", x = 0.5, y = 1.05, label = "Visium", size = 5, fontface = "bold") +
  annotate("text", x = 1.7, y = 1.05, label = "Xenium", size = 5, fontface = "bold")

```

```{r}
plot_bar_chart_dodge <- function(author_labels, manuscript_labels, color_map, title = "Cell Type Proportions") {
  library(ggplot2)
  library(dplyr)
  
  # Step 1: Combine with labels
  df_author <- data.frame(first_type = author_labels, Platform = "Visium (RCTD based on Xenium annotation)")
  df_manuscript <- data.frame(first_type = manuscript_labels, Platform = "Xenium (Author's annotation)")
  df_all <- bind_rows(df_author, df_manuscript)
  
  # Step 2: Summarize proportions
  df_summary <- df_all %>%
    group_by(Platform, first_type) %>%
    summarize(n = n(), .groups = "drop") %>%
    group_by(Platform) %>%
    mutate(Proportion = n / sum(n)) %>%
    ungroup()
  
  # Step 3: Consistent cell type order
  df_summary$first_type <- factor(df_summary$first_type, levels = names(color_map))
  
  # Step 4: Plot
  p <- ggplot(df_summary, aes(x = first_type, y = Proportion, fill = Platform)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
    scale_fill_manual(values = c(
      "Visium (RCTD based on Xenium annotation)" = "blue",
      "Xenium (Author's annotation)" = "red"
    )) +
    labs(
      title = title,
      x = "Cell Type", y = "Proportion of Singlet Spots",
      fill = "Source"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  
  return(p)
}

# Call with correct color map
author_labels <- Visium2_author_based@results$results_df$first_type
manuscript_labels <- plot_df$final_CT

# Unified color map (used for cell type x-axis consistency)
all_types <- sort(unique(c(author_labels, manuscript_labels)))
color_map <- Seurat:::DiscretePalette(length(all_types))
names(color_map) <- all_types

# Generate plot
p <- plot_bar_chart_dodge(author_labels, manuscript_labels, color_map)
print(p)


```
