# ------------------------------------------------------------
# Combined Spatial LLM Functions for Visualization and Analysis
# ------------------------------------------------------------

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(mclust)
library(rlang)
library(viridis)
library(ggrepel)
library(knitr)
library(RColorBrewer)
library(tibble)
library(CellChat)

# ------------------------------------------------------------
# 1. plot_spatial: General spatial plot from colData(obj)
# ------------------------------------------------------------
plot_spatial <- function(obj,
                         annotation,
                         title = NULL,
                         point_size = 0.5,
                         base_size = 12,
                         coord_flip = FALSE,
                         highlight = NULL) {
  # Validate that the annotation column exists
  stopifnot(annotation %in% colnames(colData(obj)))
  df <- as.data.frame(colData(obj))
  
  # Try to find known spatial coordinate sets
  coord_sets <- list(
    c("X", "Y"),
    c("X_visium", "Y_visium"),
    c("X_xenium", "Y_xenium"),
    c("X_visium_norm", "Y_visium_norm"),
    c("X_xenium_norm", "Y_xenium_norm")
  )
  
  coord_found <- FALSE
  for (coords in coord_sets) {
    if (all(coords %in% colnames(df))) {
      df$X <- as.numeric(df[[coords[1]]])
      df$Y <- as.numeric(df[[coords[2]]])
      coord_found <- TRUE
      break
    }
  }
  
  if (!coord_found) {
    stop("Spatial coordinates not found in colData(). Tried: X/Y, *_visium, *_xenium, *_norm")
  }
  
  # --- Process annotation column ---
  raw_vec <- df[[annotation]]
  
  if (is.numeric(raw_vec)) {
    df$Annotation <- raw_vec
    annotation_type <- "numeric"
  } else {
    # Attempt to convert to numeric; if it fails, treat as categorical
    try_num <- suppressWarnings(as.numeric(as.character(raw_vec)))
    if (all(!is.na(try_num))) {
      df$Annotation <- try_num
      annotation_type <- "numeric"
    } else {
      df$Annotation <- as.character(raw_vec)
      annotation_type <- "categorical"
    }
  }
  
  # Determine if highlighting should be used
  use_highlight <- annotation_type == "categorical" && !is.null(highlight)
  if (use_highlight) {
    df$HighlightGroup <- ifelse(df$Annotation %in% highlight, df$Annotation, "Other")
    df$HighlightGroup[is.na(df$Annotation)] <- NA
    df$HighlightGroup <- factor(df$HighlightGroup, levels = c(highlight, "Other"))
  }
  
  # --- Build plot ---
  p <- ggplot(df, aes(x = X, y = Y)) +
    coord_fixed() +
    theme_void(base_size = base_size) +
    labs(title = title, color = annotation) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size * 0.8)
    )
  
  # Add layers based on annotation type
  if (annotation_type == "numeric") {
    p <- p +
      geom_point(aes(color = Annotation), size = point_size, na.rm = TRUE) +
      scale_color_gradient(low = "blue", high = "red")
  } else if (use_highlight) {
    pal <- c(
      setNames(scales::hue_pal()(length(highlight)), highlight),
      "Other" = "grey80"
    )
    p <- p +
      geom_point(aes(color = HighlightGroup), size = point_size, na.rm = TRUE) +
      scale_color_manual(values = pal, drop = TRUE, na.translate = FALSE) +
      guides(color = guide_legend(override.aes = list(size = 4)))
  } else {
    df$Annotation <- factor(df$Annotation)
    p <- p +
      geom_point(aes(color = Annotation), size = point_size, na.rm = TRUE) +
      scale_color_manual(
        values = scales::hue_pal()(length(unique(na.omit(df$Annotation)))),
        na.translate = FALSE
      ) +
      guides(color = guide_legend(override.aes = list(size = 4)))
  }
  
  if (coord_flip) {
    p <- p + scale_y_reverse()
  }
  
  return(p)
}


# ------------------------------------------------------------
# 2. plot_umap: Platform-specific UMAP plot from colData(obj)
# ------------------------------------------------------------
plot_umap <- function(obj,
                      platform = c("Xenium", "Visium"),
                      annotation = "Final_CT",
                      title = NULL,
                      point_size = 0.5,
                      legend_point_size = 4,
                      highlight = NULL) {
  platform <- match.arg(platform)
  df <- as.data.frame(colData(obj))
  
  # Validate annotation and UMAP columns
  if (!(annotation %in% colnames(df))) {
    stop(paste("Annotation column", annotation, "not found in colData."))
  }
  
  umap_cols <- switch(platform,
                      Visium = c("UMAP_1_Visium", "UMAP_2_Visium"),
                      Xenium = c("UMAP_1_Xenium", "UMAP_2_Xenium")
  )
  
  if (!all(umap_cols %in% colnames(df))) {
    stop(paste("UMAP columns", paste(umap_cols, collapse = ", "), "not found in colData."))
  }
  
  x_col <- umap_cols[1]
  y_col <- umap_cols[2]
  
  plot_title <- if (!is.null(title)) title else paste(platform, "UMAP:", annotation)
  
  # Handle highlighting
  if (!is.null(highlight) && is.factor(df[[annotation]])) {
    df$highlight_flag <- ifelse(df[[annotation]] %in% highlight, as.character(df[[annotation]]), "Other")
    df$highlight_flag <- factor(df$highlight_flag, levels = c(highlight, "Other"))
    color_aes <- "highlight_flag"
    
    # Generate colors for highlighted groups
    if (length(highlight) <= 9) {
      base_colors <- RColorBrewer::brewer.pal(n = max(length(highlight), 3), "Set1")[seq_along(highlight)]
    } else {
      base_colors <- grDevices::rainbow(length(highlight))
    }
    names(base_colors) <- highlight
    colors <- c(base_colors, Other = "grey80")
  } else {
    color_aes <- annotation
    colors <- NULL
  }
  
  # Create plot
  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[color_aes]])) +
    geom_point(size = point_size) +
    coord_fixed() +
    theme_void() +
    labs(title = plot_title, color = annotation) +
    guides(color = guide_legend(override.aes = list(size = legend_point_size))) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  if (!is.null(colors)) {
    p <- p + scale_color_manual(values = colors)
  }
  
  return(p)
}


# ------------------------------------------------------------
# 3. plot_gene_spatial: Spatial gene expression visualization
# ------------------------------------------------------------
plot_gene_spatial <- function(obj,
                              gene,
                              platform = c("Xenium", "Visium"),
                              x_col = NULL,
                              y_col = NULL,
                              title = NULL,
                              point_size = 0.5,
                              low_color = "grey90",
                              high_color = "red") {
  platform <- match.arg(platform)
  
  # Define default spatial coordinates based on platform
  if (is.null(x_col)) x_col <- ifelse(platform == "Visium", "X_visium_norm", "X_xenium_norm")
  if (is.null(y_col)) y_col <- ifelse(platform == "Visium", "Y_visium_norm", "Y_xenium_norm")
  
  # Construct the full gene name with platform suffix
  gene_full <- paste0(gene, "_", tolower(platform))
  
  # Extract expression data
  expr_matrix <- assays(obj)$X
  if (!(gene_full %in% rownames(expr_matrix))) {
    stop(paste("Gene", gene_full, "not found in expression matrix."))
  }
  expr_vec <- expr_matrix[gene_full, ]
  
  df <- as.data.frame(colData(obj))
  df$Expr <- log1p(expr_vec) # Use log(1+x) transform
  
  # Validate coordinate columns
  if (!all(c(x_col, y_col) %in% colnames(df))) {
    stop(paste("Spatial coordinates", x_col, "and/or", y_col, "not found."))
  }
  
  plot_title <- if (!is.null(title)) title else paste("Spatial Expression of", gene)
  
  # Create plot
  p <- ggplot(df, aes_string(x = x_col, y = y_col, color = "Expr")) +
    geom_point(size = point_size) +
    scale_color_gradient(low = low_color, high = high_color) +
    coord_fixed() +
    theme_void() +
    labs(title = plot_title, color = "log1p(expr)") +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  return(p)
}


# ------------------------------------------------------------
# 4. rank_cross_platform_gene_correlation: Find top correlated genes
# ------------------------------------------------------------
rank_cross_platform_gene_correlation <- function(obj,
                                                 top_genes = 5,
                                                 layer = "X") {
  # Ensure the specified assay layer exists
  if (!layer %in% assayNames(obj)) {
    stop(paste("Assay '", layer, "' not found in the object."))
  }
  expr <- assay(obj, layer)
  
  # Identify genes present in both Xenium and Visium platforms
  all_genes <- rownames(expr)
  xenium_genes <- all_genes[grepl("_xenium$", all_genes)]
  visium_genes <- all_genes[grepl("_visium$", all_genes)]
  
  base_genes <- intersect(
    gsub("_xenium$", "", xenium_genes),
    gsub("_visium$", "", visium_genes)
  )
  
  if (length(base_genes) == 0) {
    stop("No matched gene pairs found between xenium and visium platforms.")
  }
  
  # Correlate expression of each gene pair across all cells
  gene_corrs <- sapply(base_genes, function(gene) {
    xenium_row <- paste0(gene, "_xenium")
    visium_row <- paste0(gene, "_visium")
    
    if (xenium_row %in% all_genes && visium_row %in% all_genes) {
      x <- expr[xenium_row, ]
      y <- expr[visium_row, ]
      # Use log1p transformation for sparse count data
      cor(log1p(x), log1p(y), method = "pearson", use = "pairwise.complete.obs")
    } else {
      NA # Return NA if a pair is missing
    }
  })
  
  # Remove NAs and create a clean results data frame
  gene_corrs <- gene_corrs[!is.na(gene_corrs)]
  results_df <- data.frame(Correlation = gene_corrs) %>%
    tibble::rownames_to_column("Gene") %>%
    arrange(desc(Correlation))
  
  top_results_df <- head(results_df, n = top_genes)
  
  # --- Generate Bar Plot ---
  correlation_plot <- ggplot(top_results_df, aes(x = reorder(Gene, -Correlation), y = Correlation)) +
    geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
    coord_cartesian(ylim = c(0, 1)) + # Fixed y-axis for consistency
    labs(
      title = "Top Cross-Platform Gene Correlations (Xenium vs. Visium)",
      x = "Gene",
      y = "Pearson Correlation"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank()
    )
  
  # Return both the table and the plot
  return(list(
    results_table = top_results_df,
    correlation_plot = correlation_plot
  ))
}


# ------------------------------------------------------------
# 5. plot_cross_platform_gene_correlation: Cumulative correlation plot
# ------------------------------------------------------------
plot_cross_platform_gene_correlation <- function(obj,
                                                 gene_name,
                                                 ID1_col = "Visium_ID",
                                                 ID2_col = "Xenium_ID",
                                                 hline = 0.9,
                                                 y_range = c(0, 1),
                                                 title = NULL) {
  coldata <- as.data.frame(colData(obj))
  if (!all(c(ID1_col, ID2_col) %in% colnames(coldata))) {
    stop("colData must include both ID1_col and ID2_col.")
  }
  
  # Infer platform suffixes from ID column names
  suffix1 <- paste0("_", tolower(gsub("_ID$", "", ID1_col)))
  suffix2 <- paste0("_", tolower(gsub("_ID$", "", ID2_col)))
  
  gene1 <- paste0(gene_name, suffix1)
  gene2 <- paste0(gene_name, suffix2)
  
  expr <- assays(obj)$X
  if (!(gene1 %in% rownames(expr)) || !(gene2 %in% rownames(expr))) {
    stop(paste("Gene", gene_name, "not found with expected suffixes in expression matrix."))
  }
  
  # Filter for matched cells
  matched_df <- coldata[!is.na(coldata[[ID1_col]]) & !is.na(coldata[[ID2_col]]), ]
  if (nrow(matched_df) < 2) stop("Fewer than 2 matched spots available.")
  matched_cells <- rownames(matched_df)
  
  # Extract and transform expression data
  log1 <- log1p(as.numeric(expr[gene1, matched_cells]))
  log2 <- log1p(as.numeric(expr[gene2, matched_cells]))
  
  # Order cells by absolute difference in expression
  abs_diff <- abs(log1 - log2)
  order_idx <- order(abs_diff)
  n_total <- length(order_idx)
  
  # Calculate cumulative correlation
  cor_values <- sapply(1:n_total, function(i) {
    if (i < 2) return(NA)
    x_sub <- log1[order_idx[1:i]]
    y_sub <- log2[order_idx[1:i]]
    cor(x_sub, y_sub, method = "pearson")
  })
  
  df_plot <- data.frame(
    Top_N = 1:n_total,
    Correlation = cor_values
  )
  
  plot_title <- if (is.null(title)) paste0(gene_name, ": Correlation vs Top Matched Spots") else title
  
  # Create plot
  p <- ggplot(df_plot, aes(x = Top_N, y = Correlation)) +
    geom_line(color = "steelblue") +
    geom_point(size = 0.6, color = "black") +
    labs(
      title = plot_title,
      x = "Number of Top Matched Spots (N)",
      y = "Pearson Correlation"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.05))) +
    scale_y_continuous(limits = y_range) +
    theme_minimal(base_size = 14)
  
  if (!is.null(hline)) {
    p <- p + geom_hline(yintercept = hline, linetype = "dashed", color = "red")
  }
  
  return(p)
}


# ------------------------------------------------------------
# 7. rank_marker_specificity: Rank gene enrichment per cluster
# ------------------------------------------------------------
rank_marker_specificity <- function(sce,
                                    gene_name = "COL1A1",
                                    cluster_cols = c("pruned"),
                                    assay_name = "X",
                                    platforms = c("Xenium", "Visium"),
                                    top_n = NULL) {
  # --- Step 1: Select appropriate gene variant based on platform ---
  gene_variants <- rownames(sce)
  full_gene <- if (paste0(gene_name, "_xenium") %in% gene_variants && "Xenium" %in% platforms) {
    paste0(gene_name, "_xenium")
  } else if (paste0(gene_name, "_visium") %in% gene_variants && "Visium" %in% platforms) {
    paste0(gene_name, "_visium")
  } else if (gene_name %in% gene_variants) {
    gene_name
  } else {
    stop(paste("Gene", gene_name, "not found for selected platform(s)."))
  }
  
  # --- Step 2: Extract expression and filter metadata ---
  expr_vector <- assay(sce, assay_name)[full_gene, ]
  meta <- as.data.frame(colData(sce))
  meta$expr <- expr_vector
  
  if ("Platform" %in% colnames(meta)) {
    meta <- meta[meta$Platform %in% platforms, ]
  }
  
  # --- Step 3: Compute specificity stats for each cluster column ---
  get_specificity_stats <- function(cluster_col) {
    if (!(cluster_col %in% colnames(meta))) return(NULL)
    
    df <- data.frame(expr = meta$expr, cluster = as.character(meta[[cluster_col]]))
    df$cluster[is.na(df$cluster)] <- "NA"
    
    result_list <- lapply(unique(df$cluster), function(cl) {
      in_group <- df$expr[df$cluster == cl]
      out_group <- df$expr[df$cluster != cl]
      if (length(in_group) < 2 || length(out_group) < 2) return(NULL)
      
      # Perform one-sided t-test
      t_result <- tryCatch(t.test(in_group, out_group, alternative = "greater"), error = function(e) NULL)
      
      data.frame(
        Cluster = cl,
        Mean_Expr = mean(in_group),
        Pct_Expr = mean(in_group > 0) * 100,
        P_Value = if (!is.null(t_result)) t_result$p.value else NA,
        Annotation_Source = cluster_col
      )
    })
    
    do.call(rbind, result_list)
  }
  
  # --- Step 4: Aggregate results and rank ---
  stats_list <- lapply(cluster_cols, get_specificity_stats)
  combined_stats <- bind_rows(stats_list)
  
  if (nrow(combined_stats) == 0) {
    warning("No valid clusters found to calculate statistics.")
    return(tibble())
  }
  
  ranked <- combined_stats %>%
    mutate(P_Adj = p.adjust(P_Value, method = "fdr")) %>%
    arrange(P_Adj, desc(Mean_Expr))
  
  if (!is.null(top_n)) {
    return(head(ranked, top_n))
  } else {
    return(ranked)
  }
}


# ------------------------------------------------------------
# 9. plot_ARI_dotplot: Dot plot for annotation concordance (ARI)
# ------------------------------------------------------------
plot_ARI_dotplot <- function(obj,
                             annotation1,
                             annotation2,
                             title = NULL,
                             size_range = c(1, 10)) {
  # --- Step 1: Extract and clean annotation data ---
  meta <- as.data.frame(colData(obj))
  if (!(annotation1 %in% colnames(meta)) || !(annotation2 %in% colnames(meta))) {
    stop("One or both annotation columns not found.")
  }
  
  labels1 <- as.character(meta[[annotation1]])
  labels2 <- as.character(meta[[annotation2]])
  
  # Filter to common, non-NA cells
  mask <- !is.na(labels1) & !is.na(labels2)
  labels1 <- labels1[mask]
  labels2 <- labels2[mask]
  
  # --- Step 2: Compute Adjusted Rand Index (ARI) ---
  ari_score <- adjustedRandIndex(labels1, labels2)
  
  # --- Step 3: Cross-tabulation ---
  cross_tab <- table(X = labels1, Y = labels2)
  df <- as.data.frame(cross_tab) %>% filter(Freq > 0)
  
  # --- Step 4: Create Dot Plot ---
  plot_title <- if (is.null(title)) {
    sprintf("%s vs %s (ARI = %.3f)", annotation1, annotation2, ari_score)
  } else {
    sprintf("%s (ARI = %.3f)", title, ari_score)
  }
  
  p <- ggplot(df, aes(x = X, y = Y, size = Freq)) +
    geom_point(alpha = 0.8, color = "steelblue") +
    scale_size_continuous(range = size_range) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major = element_line(color = "gray80"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    ) +
    labs(
      title = plot_title,
      x = annotation1,
      y = annotation2,
      size = "Count"
    )
  
  return(p)
}


# ------------------------------------------------------------
# 11. plot_cluster_composition_groupedbar: Grouped proportion bar plot
# ------------------------------------------------------------
plot_cluster_composition_groupedbar <- function(obj,
                                                annotation,
                                                palette = NULL,
                                                title = NULL,
                                                rotate_x = TRUE) {
  meta <- as.data.frame(colData(obj))
  
  # Validate number and existence of annotations
  stopifnot(length(annotation) %in% c(1, 2))
  stopifnot(all(annotation %in% colnames(meta)))
  
  # Handle one or two annotations
  if (length(annotation) == 1) {
    df <- meta %>%
      filter(!is.na(.data[[annotation[1]]])) %>%
      count(annotation_value = .data[[annotation[1]]]) %>%
      mutate(Proportion = n / sum(n), annotation_source = annotation[1])
  } else { # length(annotation) == 2
    # Calculate proportions for the first annotation
    df1 <- meta %>%
      filter(!is.na(.data[[annotation[1]]])) %>%
      count(annotation_value = .data[[annotation[1]]]) %>%
      mutate(Proportion = n / sum(n), annotation_source = annotation[1])
    
    # Calculate proportions for the second annotation
    df2 <- meta %>%
      filter(!is.na(.data[[annotation[2]]])) %>%
      count(annotation_value = .data[[annotation[2]]]) %>%
      mutate(Proportion = n / sum(n), annotation_source = annotation[2])
    
    df <- bind_rows(df1, df2)
  }
  
  # Order x-axis by total abundance across annotations
  annotation_order <- df %>%
    group_by(annotation_value) %>%
    summarize(total = sum(Proportion), .groups = "drop") %>%
    arrange(desc(total)) %>%
    pull(annotation_value)
  
  df <- df %>%
    mutate(
      annotation_value = factor(annotation_value, levels = annotation_order),
      annotation_source = factor(annotation_source)
    )
  
  # Use a default palette if comparing Visium and Xenium RCTD results
  if (is.null(palette) && all(sort(annotation) == sort(c("Visium_RCTD_LungMap_ref", "Xenium_RCTD_LungMap_ref")))) {
    palette <- c("Visium_RCTD_LungMap_ref" = "blue", "Xenium_RCTD_LungMap_ref" = "red")
  }
  
  plot_title <- if (is.null(title)) "Cluster Composition Comparison" else title
  
  # Create plot
  p <- ggplot(df, aes(x = annotation_value, y = Proportion, fill = annotation_source)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    labs(
      x = "Cell Type",
      y = "Proportion",
      fill = "Annotation Source",
      title = plot_title
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = if (rotate_x) 45 else 0, hjust = if (rotate_x) 1 else 0.5),
      panel.grid.major.x = element_blank()
    )
  
  if (!is.null(palette)) {
    p <- p + scale_fill_manual(values = palette)
  }
  
  return(p)
}


# ------------------------------------------------------------
# 12. compare_annotations_via_heatmap: Confusion matrix as heatmap
# ------------------------------------------------------------
compare_annotations_via_heatmap <- function(obj,
                                            annotation1,
                                            annotation2,
                                            normalize = c("none", "row", "column", "both"),
                                            palette = "Blues",
                                            title = NULL) {
  normalize <- match.arg(normalize)
  
  # Extract metadata and count overlaps
  meta <- as.data.frame(colData(obj))
  stopifnot(annotation1 %in% colnames(meta), annotation2 %in% colnames(meta))
  
  df <- meta %>%
    filter(!is.na(.data[[annotation1]]), !is.na(.data[[annotation2]])) %>%
    count(Annotation1 = .data[[annotation1]], Annotation2 = .data[[annotation2]])
  
  # Pivot to a wide matrix format
  mat <- df %>%
    pivot_wider(names_from = Annotation2, values_from = n, values_fill = 0) %>%
    column_to_rownames("Annotation1") %>%
    as.matrix()
  
  # Normalize the matrix based on user's choice
  if (normalize == "row") {
    mat <- sweep(mat, 1, rowSums(mat), FUN = "/")
  } else if (normalize == "column") {
    mat <- sweep(mat, 2, colSums(mat), FUN = "/")
  } else if (normalize == "both") {
    mat <- mat / sum(mat)
  }
  
  # Convert back to long format for plotting
  df_plot <- as.data.frame(as.table(mat))
  colnames(df_plot) <- c("Annotation1", "Annotation2", "Value")
  
  plot_title <- if (is.null(title)) "Annotation Comparison Heatmap" else title
  
  # Create plot
  p <- ggplot(df_plot, aes(x = Annotation2, y = Annotation1, fill = Value)) +
    geom_tile(color = "white") +
    scale_fill_distiller(palette = palette, direction = 1, na.value = "grey90") +
    geom_text(aes(label = ifelse(Value == 0, "", scales::comma(round(Value, 2)))), size = 3.2) +
    theme_minimal(base_size = 12) +
    labs(
      x = annotation2,
      y = annotation1,
      fill = ifelse(normalize == "none", "Count", "Proportion"),
      title = plot_title
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = element_blank()
    )
  
  return(p)
}


# ------------------------------------------------------------
# 13. visualize_matched_pairs_scatter: Scatter plot for matched gene expression
# ------------------------------------------------------------
visualize_matched_pairs_scatter <- function(obj,
                                            gene,
                                            layer = "X",
                                            log1p_transform = TRUE,
                                            point_size = 0.4,
                                            alpha = 0.6,
                                            show_correlation = FALSE,
                                            title = NULL) {
  # Extract expression data
  expr_mat <- assay(obj, layer)
  gene_visium <- paste0(gene, "_visium")
  gene_xenium <- paste0(gene, "_xenium")
  
  if (!(gene_visium %in% rownames(expr_mat)) || !(gene_xenium %in% rownames(expr_mat))) {
    stop(paste("Gene", gene, "not found in both Visium and Xenium formats."))
  }
  
  visium_expr <- expr_mat[gene_visium, ]
  xenium_expr <- expr_mat[gene_xenium, ]
  
  if (log1p_transform) {
    visium_expr <- log1p(visium_expr)
    xenium_expr <- log1p(xenium_expr)
  }
  
  df <- data.frame(
    Xenium = xenium_expr,
    Visium = visium_expr
  )
  
  # Compute correlation if requested
  plot_title <- if (!is.null(title)) title else paste0("Xenium vs Visium: ", gene)
  if (show_correlation) {
    corr <- cor(df$Xenium, df$Visium, method = "pearson", use = "pairwise.complete.obs")
    corr_text <- paste0("Pearson r = ", round(corr, 3))
    plot_title <- paste0(plot_title, "\n", corr_text)
  }
  
  # Create plot
  p <- ggplot(df, aes(x = Xenium, y = Visium)) +
    geom_point(size = point_size, alpha = alpha) +
    labs(
      x = if (log1p_transform) "Xenium Expression (log1p)" else "Xenium Expression",
      y = if (log1p_transform) "Visium Expression (log1p)" else "Visium Expression",
      title = plot_title
    ) +
    theme_minimal(base_size = 13)
  
  return(p)
}


# ------------------------------------------------------------
# 14. generate_stability_report: Report on cluster stability metrics
# ------------------------------------------------------------
generate_stability_report <- function(sce, celltype_labels) {
  # --- Part 1: Extract and clean data from the sce object ---
  metadata <- as.data.frame(colData(sce))
  all_cluster_scores <- list()
  
  # Iterate through the user-provided list of annotations and labels
  for (annot in names(celltype_labels)) {
    for (label in celltype_labels[[annot]]) {
      cluster_string <- paste(annot, label, sep = "@")
      
      # Find a representative cell for the cluster
      rep_cell <- metadata[metadata[[annot]] == label, ][1, ]
      if (nrow(rep_cell) == 0) {
        warning(paste("Could not find cluster:", cluster_string))
        next
      }
      
      # Define expected column names for stability scores
      col_reassign <- paste0("reassign.", annot)
      col_sccaf <- paste0("SCCAF.", annot)
      col_tfidf5 <- paste0("tfidf5.", annot)
      col_tfidf10 <- paste0("tfidf10.", annot)
      col_shapley <- paste0(annot, "_shapley")
      col_confidence <- "confidence"
      expected_cols <- c(col_reassign, col_sccaf, col_tfidf5, col_tfidf10, col_shapley, col_confidence)
      
      if (!all(expected_cols %in% names(rep_cell))) {
        missing <- expected_cols[!expected_cols %in% names(rep_cell)]
        warning(paste("For '", annot, "', columns not found: ", paste(missing, collapse = ", "), ". Skipping.", sep = ""))
        next
      }
      
      scores <- rep_cell[1, expected_cols]
      cluster_df <- tibble(
        Cluster = cluster_string,
        Metric = c("Reassign", "SCCAF", "TF-IDF 5", "TF-IDF 10", "Shapley", "Confidence"),
        Score = as.numeric(scores)
      )
      all_cluster_scores[[cluster_string]] <- cluster_df
    }
  }
  
  long_stability_data <- bind_rows(all_cluster_scores) %>% filter(!is.na(Score))
  if (nrow(long_stability_data) == 0) {
    message("No stability data could be extracted based on your input.")
    return(NULL)
  }
  
  # --- Part 2: Create a wide-format summary table ---
  summary_table <- long_stability_data %>%
    pivot_wider(
      names_from = Metric,
      values_from = Score
    )
  
  # --- Part 3: Create the stability bar plot ---
  metric_order <- c("Confidence", "Reassign", "SCCAF", "TF-IDF 10", "TF-IDF 5", "Shapley")
  long_stability_data$Metric <- factor(long_stability_data$Metric, levels = metric_order)
  
  stability_plot <- ggplot(long_stability_data, aes(x = Metric, y = Score, fill = Cluster)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.9) +
    scale_fill_viridis_d(name = "Cluster Annotation", guide = guide_legend(nrow = 2)) +
    labs(y = "Score", x = NULL, title = "Cluster Stability Metrics") +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )
  
  # --- Part 4: Return both table and plot ---
  return(list(summary_table = summary_table, stability_plot = stability_plot))
}


# ------------------------------------------------------------
# 15. generate_marker_report: Comprehensive marker analysis and dot plot
# ------------------------------------------------------------
generate_marker_report <- function(obj,
                                   genes,
                                   annotation_cols,
                                   assay_name = "X",
                                   min_pct = 5,
                                   max_genes_plot = 50) {
  # Step 1: Define platform-specific rules for annotations
  rule_map <- tibble(
    annotation = c("Final_CT", "Xenium_RCTD_LungMap_ref", "Visium_RCTD_LungMap_ref", "Visium_RCTD_GSE250346_based", "pruned", "pruned"),
    platform = c("Xenium", "Xenium", "Visium", "Visium", "Xenium", "Visium")
  )
  
  # Create tasks for platform-specific and base genes
  rules_to_run <- rule_map %>%
    filter(annotation %in% annotation_cols) %>%
    crossing(gene = genes) %>%
    mutate(full_gene = paste0(gene, "_", tolower(platform)))
  base_gene_tasks <- tibble(annotation = annotation_cols, platform = "Base", gene = genes, full_gene = genes)
  tasks <- bind_rows(rules_to_run, base_gene_tasks)
  
  # Step 2: Run analysis for each task
  all_results <- lapply(seq_len(nrow(tasks)), function(i) {
    task <- tasks[i, ]
    if (!task$full_gene %in% rownames(obj)) return(NULL)
    
    expr_vector <- assay(obj, assay_name)[task$full_gene, ]
    meta <- as.data.frame(colData(obj))
    
    df <- data.frame(expr = expr_vector, cluster = as.character(meta[[task$annotation]]))
    df$cluster[is.na(df$cluster)] <- "NA"
    
    # Calculate stats for each cluster
    cluster_stats <- lapply(unique(df$cluster), function(cl) {
      in_group_expr <- df$expr[df$cluster == cl]
      out_group_expr <- df$expr[df$cluster != cl]
      if (length(in_group_expr) < 2 || length(out_group_expr) < 2) return(NULL)
      
      t_result <- tryCatch(t.test(in_group_expr, out_group_expr, alternative = "greater"), error = function(e) NULL)
      
      data.frame(
        Gene = task$gene,
        Platform = task$platform,
        Cluster = paste(task$annotation, cl, sep = "@"),
        Mean_Expr = mean(in_group_expr, na.rm = TRUE),
        Pct_Expr = mean(in_group_expr > 0, na.rm = TRUE) * 100,
        P_Value = if (!is.null(t_result)) t_result$p.value else NA_real_
      )
    })
    bind_rows(cluster_stats)
  })
  
  # Step 3: Create final table and dot plot
  final_table <- bind_rows(all_results)
  if (nrow(final_table) == 0) {
    message("No results could be calculated.")
    return(NULL)
  }
  
  final_table <- final_table %>%
    mutate(P_Adj = p.adjust(P_Value, method = "BH")) %>%
    arrange(P_Adj, P_Value)
  
  # Prepare data for dot plot
  dot_data <- final_table %>%
    filter(Pct_Expr >= min_pct) %>%
    mutate(
      gene_label = if_else(Platform != "Base", paste0(Gene, " (", Platform, ")"), Gene),
      Mean_Expr = ifelse(is.nan(Mean_Expr), 0, Mean_Expr)
    )
  
  # Limit number of genes in the plot for readability
  if (!is.null(max_genes_plot) && n_distinct(dot_data$gene_label) > max_genes_plot) {
    top_genes <- dot_data %>%
      group_by(gene_label) %>%
      summarize(score = -log10(min(P_Adj, na.rm = TRUE) + 1e-300) * max(Pct_Expr, na.rm = TRUE)) %>%
      arrange(desc(score)) %>%
      slice_head(n = max_genes_plot) %>%
      pull(gene_label)
    dot_data <- dot_data %>% filter(gene_label %in% top_genes)
  }
  
  p <- ggplot(dot_data, aes(x = gene_label, y = Cluster)) +
    geom_point(aes(size = Pct_Expr, color = Mean_Expr)) +
    scale_size(name = "% Expressing", range = c(1, 8)) +
    scale_color_gradient(name = "Avg. Expression", low = "grey90", high = "red") +
    theme_minimal(base_size = 14) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    )
  
  return(list(plot = p, table = final_table))
}


# ------------------------------------------------------------
# 16. find_celltype_markers: Find DE genes for a cell type (Volcano Plot)
# ------------------------------------------------------------
find_celltype_markers <- function(sce,
                                  celltype,
                                  annotation_col,
                                  platform = "Xenium",
                                  assay_name = "X",
                                  top_n = 10,
                                  logfc_threshold = 0.25,
                                  padj_threshold = 0.05) {
  # Step 1: Filter genes by platform and get expression data
  expr_mat <- assay(sce, assay_name)
  suffix <- paste0("_", tolower(platform), "$")
  gene_mask <- grepl(suffix, rownames(expr_mat))
  
  # Fallback to unique base gene names if no platform-specific genes found
  if (sum(gene_mask) == 0) {
    message("No platform-specific genes found, using unique base gene names.")
    base_genes <- gsub(suffix, "", rownames(expr_mat))
    gene_mask <- !duplicated(base_genes)
    if (sum(gene_mask) == 0) stop("No suitable genes found.")
  }
  expr_mat <- expr_mat[gene_mask, , drop = FALSE]
  
  # Step 2: Define target and background cell groups
  meta <- as.data.frame(colData(sce))
  target_mask <- meta[[annotation_col]] == celltype
  other_mask <- meta[[annotation_col]] != celltype
  
  if (sum(target_mask, na.rm = TRUE) < 3 || sum(other_mask, na.rm = TRUE) < 3) {
    stop("Not enough cells in target or other groups for comparison.")
  }
  
  # Step 3: Compute logFC and p-values using Wilcoxon test
  mean_target <- rowMeans(expr_mat[, target_mask], na.rm = TRUE)
  mean_other <- rowMeans(expr_mat[, other_mask], na.rm = TRUE)
  logFC <- log2(mean_target + 1e-9) - log2(mean_other + 1e-9)
  
  p_values <- apply(expr_mat, 1, function(gene_expr) {
    tryCatch(
      wilcox.test(gene_expr[target_mask], gene_expr[other_mask], alternative = "greater")$p.value,
      error = function(e) NA_real_
    )
  })
  
  # Step 4: Create full results table with significance flags
  full_results_df <- tibble(
    gene = gsub(suffix, "", rownames(expr_mat)),
    logFC = logFC,
    P_Value = p_values
  ) %>%
    filter(!is.na(logFC) & !is.na(P_Value)) %>%
    mutate(
      P_Adj = p.adjust(P_Value, method = "BH"),
      significant = P_Adj < padj_threshold & logFC > logfc_threshold
    ) %>%
    mutate(log10_P_Adj = -log10(P_Adj))
  
  
  # Handle infinite -log10(P_Adj) values resulting from p-values of 0
  if (any(is.infinite(full_results_df$log10_P_Adj))) {
    max_finite_val <- max(full_results_df$log10_P_Adj[is.finite(full_results_df$log10_P_Adj)], na.rm = TRUE)
    full_results_df$log10_P_Adj[is.infinite(full_results_df$log10_P_Adj)] <- max_finite_val * 1.1
  }
  
  # Step 5: Create a filtered table of top markers
  top_markers_table <- full_results_df %>%
    filter(significant == TRUE) %>%
    arrange(desc(log10_P_Adj), desc(logFC)) %>%
    slice_head(n = top_n)
  
  # Step 6: Generate the Volcano Plot
  volcano_plot <- ggplot(full_results_df, aes(x = logFC, y = log10_P_Adj)) +
    geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
    geom_text_repel(
      data = top_markers_table,
      aes(label = gene),
      size = 4, box.padding = 0.5, max.overlaps = Inf
    ) +
    scale_color_manual(
      name = "",
      values = c("TRUE" = "red", "FALSE" = "grey50"),
      labels = c("TRUE" = "Significant", "FALSE" = "Not Significant")
    ) +
    geom_vline(xintercept = logfc_threshold, linetype = "dashed", color = "blue") +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "blue") +
    # Expand axis limits to give labels more space
    coord_cartesian(
      xlim = c(min(full_results_df$logFC, na.rm = TRUE) - 1, max(full_results_df$logFC, na.rm = TRUE) * 1.15),
      ylim = c(-1, max(full_results_df$log10_P_Adj, na.rm = TRUE) * 1.1)
    ) +
    labs(
      title = paste("Volcano Plot for", celltype),
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value"
    ) +
    theme_classic(base_size = 14) +
    theme(legend.position = "top")
  
  # Step 7: Return both the table and the plot
  return(list(table = top_markers_table, plot = volcano_plot))
}

# ------------------------------------------------------------
# Squidpy
# ------------------------------------------------------------
# ------------------------------------------------------------
# 17. create_enrichment_heatmap
# ------------------------------------------------------------
create_enrichment_heatmap <- function(obj,
                                      group = "Unaffected",
                                      annotation = "Final_CT",
                                      plot_style = "viridis",
                                      clustering = TRUE,
                                      title = NULL,
                                      base_size = 10) {
  
  # Load required packages
  require(ComplexHeatmap)
  require(circlize)
  require(viridis)
  require(dplyr)
  
  # Extract enrichment matrix
  enrichment_df <- obj$Squidpy[[annotation]][[group]]
  if (is.null(enrichment_df)) stop("No enrichment matrix found for annotation: ", annotation, ", group: ", group)
  
  enrichment_mat <- as.matrix(enrichment_df)
  if (!is.numeric(enrichment_mat)) stop("Enrichment matrix must be numeric.")
  
  # Replace NaN/Inf with NA
  enrichment_mat[!is.finite(enrichment_mat)] <- NA
  num_na <- sum(is.na(enrichment_mat))
  if (num_na > 0) {
    warning(paste("Matrix contains", num_na, "non-finite values. Replacing with 0."))
    enrichment_mat[is.na(enrichment_mat)] <- 0
  }
  
  # Set default title
  if (is.null(title)) {
    title <- paste0("Enrichment Heatmap: ", group, " (", annotation, ")")
  }
  
  # Determine viridis scale
  min_val <- min(enrichment_mat, na.rm = TRUE)
  max_val <- max(enrichment_mat, na.rm = TRUE)
  if (!is.finite(min_val) || !is.finite(max_val) || min_val == max_val) {
    min_val <- 0
    max_val <- 1
    warning("Invalid viridis range. Using [0, 1].")
  }
  
  # Color function (viridis only)
  col_fun <- circlize::colorRamp2(seq(min_val, max_val, length.out = 100), viridis::viridis(100))
  
  # Plot heatmap
  Heatmap(
    enrichment_mat,
    name = "Z-score",
    col = col_fun,
    cluster_rows = clustering,
    cluster_columns = clustering,
    row_names_gp = grid::gpar(fontsize = base_size),
    column_names_gp = grid::gpar(fontsize = base_size),
    heatmap_legend_param = list(
      title = "Z-score",
      title_gp = grid::gpar(fontsize = base_size + 1),
      labels_gp = grid::gpar(fontsize = base_size)
    ),
    column_title = title,
    column_title_gp = grid::gpar(fontsize = base_size + 2, fontface = "bold")
  )
}


# ------------------------------------------------------------
# 18. analyze_and_plot_neighbors: rank and visualize neighboring cell types
# ------------------------------------------------------------
analyze_and_plot_neighbors <- function(obj,
                                       target_cell_type,
                                       group = "More_Affected",
                                       annotation = "Final_CT",
                                       base_size = 11) {
  
  # Load required packages
  require(ggplot2)
  require(scales)
  require(dplyr)
  
  # Step 1: Extract the enrichment matrix
  enrichment_df <- obj$Squidpy[[annotation]][[group]]
  
  if (is.null(enrichment_df)) {
    stop(paste0("No enrichment matrix found for group = '", group, "' and annotation = '", annotation, "'"))
  }
  
  if (!(target_cell_type %in% rownames(enrichment_df))) {
    stop(paste0("Cell type '", target_cell_type, "' not found in rownames. Available types: ", 
                paste(rownames(enrichment_df), collapse = ", ")))
  }
  
  # Step 2: Rank and transform
  ranked_values <- sort(enrichment_df[target_cell_type, ], decreasing = TRUE)
  log_transformed <- sign(ranked_values) * log10(abs(ranked_values) + 1)
  log_df <- data.frame(
    cell_type = factor(names(log_transformed), levels = rev(names(log_transformed))),
    log_score = as.numeric(log_transformed)
  )
  
  # Step 3: Print transformed scores
  cat("\n--- Cell Types Adjacent to", target_cell_type, "in group:", group, "---\n")
  print(log_df)
  cat("\n")
  
  # Step 4: Plot
  p <- ggplot(log_df, aes(x = log_score, y = cell_type, fill = log_score >= 0)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("red", "royalblue"), guide = "none") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    labs(
      title = paste0("Neighborhood Enrichment: ", target_cell_type, " in ", gsub("_", " ", group)),
      x = "Signed log10(z-score + 1)",
      y = "Neighboring Cell Type"
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.y = element_text(color = "black"),
      axis.text.x = element_text(color = "black")
    )
  
  print(p)
}

# ------------------------------------------------------------
# 19. generate_de_volcano_plot: returns table + volcano plot of DE genes
# ------------------------------------------------------------
generate_de_volcano_plot <- function(obj_list,
                                     annotation = "TNiche",
                                     group_1 = "More_Affected",
                                     group_2 = "Unaffected",
                                     celltype_filter = NULL,
                                     top_n = 100,
                                     label_top_n = 20,
                                     logfc_threshold = 0.25,
                                     padj_threshold = 0.05) {
  # Load required libraries
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(tibble)
  
  # Step 1: Extract expression matrices and metadata
  mat1 <- obj_list[[annotation]][[group_1]]@data
  mat2 <- obj_list[[annotation]][[group_2]]@data
  meta1 <- obj_list[[annotation]][[group_1]]@meta
  meta2 <- obj_list[[annotation]][[group_2]]@meta
  
  if (!is.null(celltype_filter)) {
    keep1 <- which(meta1[[annotation]] == celltype_filter)
    keep2 <- which(meta2[[annotation]] == celltype_filter)
    mat1 <- mat1[, keep1, drop = FALSE]
    mat2 <- mat2[, keep2, drop = FALSE]
    message(paste0("Subsetting to cell type '", celltype_filter, "' — ",
                   ncol(mat1), " in ", group_1, ", ", ncol(mat2), " in ", group_2))
  }
  
  # Combine expression and metadata
  combined_mat <- cbind(mat1, mat2)
  labels <- c(rep(group_1, ncol(mat1)), rep(group_2, ncol(mat2)))
  meta_df <- data.frame(annotation = labels)
  colnames(meta_df) <- annotation
  rownames(meta_df) <- colnames(combined_mat)
  
  # Step 2: Create Seurat object and normalize
  seu <- CreateSeuratObject(counts = combined_mat)
  seu@meta.data[[annotation]] <- meta_df[[annotation]]
  seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Step 3: Differential expression analysis
  deg_result <- FindMarkers(
    object = seu,
    ident.1 = group_1,
    ident.2 = group_2,
    group.by = annotation,
    logfc.threshold = 0,
    min.pct = 0.0,
    min.diff.pct = -Inf,
    test.use = "wilcox",
    verbose = FALSE
  ) %>%
    rownames_to_column("Gene") %>%
    mutate(
      log10_padj = -log10(p_val_adj + 1e-300),
      log10_padj = ifelse(log10_padj > 300, 300, log10_padj),
      Sig = case_when(
        p_val_adj < padj_threshold & avg_log2FC > logfc_threshold ~ "Up in Group 1",
        p_val_adj < padj_threshold & avg_log2FC < -logfc_threshold ~ "Down in Group 1",
        TRUE ~ "Not significant"
      )
    )
  
  # Step 4: Sort by log10_padj and then abs(avg_log2FC)
  deg_result <- deg_result %>%
    arrange(desc(log10_padj), desc(abs(avg_log2FC)))
  
  top_table <- deg_result %>%
    slice(1:top_n)
  
  label_genes <- deg_result %>%
    slice(1:label_top_n) %>%
    pull(Gene)
  
  # Step 5: Volcano plot
  volcano <- ggplot(deg_result, aes(x = avg_log2FC, y = log10_padj, color = Sig)) +
    geom_point(alpha = 0.7, size = 1) +
    scale_color_manual(values = c(
      "Up in Group 1" = "red",
      "Down in Group 1" = "blue",
      "Not significant" = "gray"
    )) +
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black") +
    geom_text_repel(
      data = subset(deg_result, Gene %in% label_genes),
      aes(label = Gene),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.5,
      segment.size = 0.2
    ) +
    labs(
      title = paste(group_1, "vs", group_2, if (!is.null(celltype_filter)) paste("in", celltype_filter) else ""),
      x = "Average Log2 Fold Change",
      y = "-log10(Adjusted P-value)",
      color = "Significance"
    ) +
    theme_minimal(base_size = 12)
  
  return(list(table = top_table, plot = volcano, full_result = deg_result))
}


# ------------------------------------------------------------
# CellChat
# ------------------------------------------------------------
# ------------------------------------------------------------
# 20. netVisual_circle: Wrapper with all parameters for CellChat's netVisual_circle
# ------------------------------------------------------------
# Description:
# This function is a wrapper for `CellChat::netVisual_circle`, which visualizes the 
# cell-cell communication network as a circular graph, showing direction and strength of interactions.
#
# Main argument:
#   net: The communication network matrix (e.g., cellchat@net$count).
#
# Example Usage:
#   netVisual_circle(
#     net = CellChat_Squidpy[["CellChat"]][["Final_CT"]][["Unaffected"]]@net$count,
#     vertex.weight = groupSizes$CellChat_Squidpy[["CellChat"]][["Final_CT"]][["Unaffected"]],
#     vertex.label.cex = 0.7,
#     margin = 0.2,
#     remove.isolate = TRUE
#   )
#
netVisual_circle <- function(
    net,
    color.use = NULL,
    title.name = NULL,
    sources.use = NULL,
    targets.use = NULL,
    idents.use = NULL,
    remove.isolate = FALSE,
    top = 1,
    weight.scale = FALSE,
    vertex.weight = 20,
    vertex.weight.max = NULL,
    vertex.size.max = NULL,
    vertex.label.cex = 1,
    vertex.label.color = "black",
    edge.weight.max = NULL,
    edge.width.max = 8,
    alpha.edge = 0.6,
    label.edge = FALSE,
    edge.label.color = "black",
    edge.label.cex = 0.8,
    edge.curved = 0.2,
    shape = "circle",
    layout = igraph::in_circle(),
    margin = 0.2,
    vertex.size = NULL,
    arrow.width = 1,
    arrow.size = 0.2,
    text.x = 0,
    text.y = 1.5
) {
  # --- Input Validation ---
  if (missing(net) || (!is.data.frame(net) && !is.matrix(net))) {
    stop("The 'net' argument is required and must be a data frame or matrix.", call. = FALSE)
  }
  
  # --- Call the original CellChat function ---
  CellChat::netVisual_circle(
    net = net,
    color.use = color.use,
    title.name = title.name,
    sources.use = sources.use,
    targets.use = targets.use,
    idents.use = idents.use,
    remove.isolate = remove.isolate,
    top = top,
    weight.scale = weight.scale,
    vertex.weight = vertex.weight,
    vertex.weight.max = vertex.weight.max,
    vertex.size.max = vertex.size.max,
    vertex.label.cex = vertex.label.cex,
    vertex.label.color = vertex.label.color,
    edge.weight.max = edge.weight.max,
    edge.width.max = edge.width.max,
    alpha.edge = alpha.edge,
    label.edge = label.edge,
    edge.label.color = edge.label.color,
    edge.label.cex = edge.label.cex,
    edge.curved = edge.curved,
    shape = shape,
    layout = layout,
    margin = margin,
    vertex.size = vertex.size,
    arrow.width = arrow.width,
    arrow.size = arrow.size,
    text.x = text.x,
    text.y = text.y
  )
}


# ------------------------------------------------------------
# 21. netVisual_heatmap: Wrapper with all parameters for CellChat's netVisual_heatmap
# ------------------------------------------------------------
# Description:
# This function wraps `CellChat::netVisual_heatmap`, which generates a heatmap visualization
# of the communication network. It supports plotting interaction count or strength (weight)
# in a single CellChat object, or differential interactions across two conditions.
#
# Main argument:
#   object: A CellChat object (e.g., CellChat_Squidpy[["CellChat"]][["Final_CT"]][["More_Affected"]])
#
# Example Usage:
#   netVisual_heatmap(
#     object = CellChat_Squidpy[["CellChat"]][["Final_CT"]][["More_Affected"]],
#     measure = "weight",
#     color.heatmap = "Reds",
#     title.name = "Xenium IPF"
#   )
#
netVisual_heatmap <- function(
    object,
    comparison = c(1, 2),
    measure = c("count", "weight"),
    signaling = NULL,
    slot.name = c("netP", "net"),
    color.use = NULL,
    color.heatmap = NULL,
    title.name = NULL,
    width = NULL,
    height = NULL,
    font.size = 8,
    font.size.title = 10,
    cluster.rows = FALSE,
    cluster.cols = FALSE,
    sources.use = NULL,
    targets.use = NULL,
    remove.isolate = FALSE,
    row.show = NULL,
    col.show = NULL
) {
  # --- Input Validation ---
  if (!inherits(object, "CellChat")) {
    stop("The 'object' must be a valid CellChat object.", call. = FALSE)
  }
  
  # --- Call original CellChat function ---
  CellChat::netVisual_heatmap(
    object = object,
    comparison = comparison,
    measure = measure,
    signaling = signaling,
    slot.name = slot.name,
    color.use = color.use,
    color.heatmap = color.heatmap,
    title.name = title.name,
    width = width,
    height = height,
    font.size = font.size,
    font.size.title = font.size.title,
    cluster.rows = cluster.rows,
    cluster.cols = cluster.cols,
    sources.use = sources.use,
    targets.use = targets.use,
    remove.isolate = remove.isolate,
    row.show = row.show,
    col.show = col.show
  )
}

# ------------------------------------------------------------
# 22. netVisual_aggregate: Wrapper with all parameters for CellChat's netVisual_aggregate
# ------------------------------------------------------------
# Description:
# This function wraps `CellChat::netVisual_aggregate`, which visualizes the inferred signaling
# network by aggregating all ligand–receptor pairs within a specified signaling pathway.
# Multiple layouts (circle, hierarchy, chord, spatial) are supported.
#
# Main argument:
#   object: A CellChat object (e.g., CellChat_Squidpy[["CellChat"]][["Final_CT"]][["More_Affected"]])
#   signaling: The signaling pathway to visualize (e.g., "UGRP1")
#
# Example Usage:
#   netVisual_aggregate(
#     object = CellChat_Squidpy[["CellChat"]][["Final_CT"]][["More_Affected"]],
#     signaling = "UGRP1",
#     layout = "circle",
#     remove.isolate = TRUE
#   )
#
netVisual_aggregate <- function(
    object,
    signaling,
    signaling.name = NULL,
    color.use = NULL,
    thresh = 0.05,
    vertex.receiver = NULL,
    sources.use = NULL,
    targets.use = NULL,
    idents.use = NULL,
    top = 1,
    remove.isolate = FALSE,
    vertex.weight = 1,
    vertex.weight.max = NULL,
    vertex.size.max = NULL,
    weight.scale = TRUE,
    edge.weight.max = NULL,
    edge.width.max = 8,
    layout = c("circle", "hierarchy", "chord", "spatial"),
    pt.title = 12,
    title.space = 6,
    vertex.label.cex = 0.8,
    sample.use = NULL,
    alpha.image = 0.15,
    point.size = 1.5,
    group = NULL,
    cell.order = NULL,
    small.gap = 1,
    big.gap = 10,
    scale = FALSE,
    reduce = -1,
    show.legend = FALSE,
    legend.pos.x = 20,
    legend.pos.y = 20,
    ...
) {
  # --- Input Validation ---
  if (missing(object) || !inherits(object, "CellChat")) {
    stop("The 'object' must be a valid CellChat object.", call. = FALSE)
  }
  if (missing(signaling)) {
    stop("The 'signaling' pathway must be specified.", call. = FALSE)
  }
  
  # --- Call original CellChat function ---
  CellChat::netVisual_aggregate(
    object = object,
    signaling = signaling,
    signaling.name = signaling.name,
    color.use = color.use,
    thresh = thresh,
    vertex.receiver = vertex.receiver,
    sources.use = sources.use,
    targets.use = targets.use,
    idents.use = idents.use,
    top = top,
    remove.isolate = remove.isolate,
    vertex.weight = vertex.weight,
    vertex.weight.max = vertex.weight.max,
    vertex.size.max = vertex.size.max,
    weight.scale = weight.scale,
    edge.weight.max = edge.weight.max,
    edge.width.max = edge.width.max,
    layout = layout,
    pt.title = pt.title,
    title.space = title.space,
    vertex.label.cex = vertex.label.cex,
    sample.use = sample.use,
    alpha.image = alpha.image,
    point.size = point.size,
    group = group,
    cell.order = cell.order,
    small.gap = small.gap,
    big.gap = big.gap,
    scale = scale,
    reduce = reduce,
    show.legend = show.legend,
    legend.pos.x = legend.pos.x,
    legend.pos.y = legend.pos.y,
    ...
  )
}


# ------------------------------------------------------------
# 23. subsetCommunication: Wrapper with all parameters for CellChat's subsetCommunication
# ------------------------------------------------------------
# Description:
# This function wraps `CellChat::subsetCommunication` to extract specific subsets of inferred
# cell-cell communication events based on sender/receiver identities, signaling pathway, 
# ligand-receptor pairs, statistical thresholds, or expression-based filters.
#
# If all arguments are NULL, it returns all inferred communications in the specified slot.
#
# Main argument:
#   object: A CellChat object.
#
# Example Usage:
#   subsetCommunication(
#     object = CellChat_Squidpy[["CellChat"]][["Final_CT"]][["More_Affected"]],
#     signaling = "UGRP1",
#     slot.name = "netP"
#   )
#
subsetCommunication <- function(
    object = NULL,
    net = NULL,
    slot.name = "net",
    sources.use = NULL,
    targets.use = NULL,
    signaling = NULL,
    pairLR.use = NULL,
    thresh = 0.05,
    datasets = NULL,
    ligand.pvalues = NULL,
    ligand.logFC = NULL,
    ligand.pct.1 = NULL,
    ligand.pct.2 = NULL,
    receptor.pvalues = NULL,
    receptor.logFC = NULL,
    receptor.pct.1 = NULL,
    receptor.pct.2 = NULL
) {
  # --- Input Validation ---
  if (is.null(object) && is.null(net)) {
    stop("You must provide either a 'CellChat' object or a communication data frame in 'net'.", call. = FALSE)
  }
  
  # --- Call original CellChat function ---
  CellChat::subsetCommunication(
    object = object,
    net = net,
    slot.name = slot.name,
    sources.use = sources.use,
    targets.use = targets.use,
    signaling = signaling,
    pairLR.use = pairLR.use,
    thresh = thresh,
    datasets = datasets,
    ligand.pvalues = ligand.pvalues,
    ligand.logFC = ligand.logFC,
    ligand.pct.1 = ligand.pct.1,
    ligand.pct.2 = ligand.pct.2,
    receptor.pvalues = receptor.pvalues,
    receptor.logFC = receptor.logFC,
    receptor.pct.1 = receptor.pct.1,
    receptor.pct.2 = receptor.pct.2
  )
}


# ------------------------------------------------------------
# Accessing Active Signaling Pathways from netP
# ------------------------------------------------------------
# Description:
# Each CellChat object stores the list of inferred signaling pathways in:
#   object@netP$pathways
# This vector represents the signaling pathways detected after statistical filtering.

# Example:
# Unaffected sample pathways
# CellChat_Squidpy[["CellChat"]][["Final_CT"]][["Unaffected"]]@netP$pathways
# [1] "VEGF"  "SPP1"  "UGRP1" "CCL"   "FASLG" "EGF"

# More Affected (IPF) sample pathways
# CellChat_Squidpy[["CellChat"]][["Final_CT"]][["More_Affected"]]@netP$pathways
# [1] "SPP1"  "UGRP1" "VEGF"  "CCL"   "EGF"   "CXCL"  "FASLG"




