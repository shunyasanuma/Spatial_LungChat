
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

# ------------------------------------------------------------
# 1. plot_spatial: general spatial plot from colData(obj)
# ------------------------------------------------------------
plot_spatial <- function(obj,
                         annotation,
                         title = NULL,
                         point_size = 0.5,
                         base_size = 12,
                         coord_flip = FALSE,
                         highlight = NULL) {
  stopifnot(annotation %in% colnames(colData(obj)))
  df <- as.data.frame(colData(obj))

  # Try known spatial coordinate sets
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
    stop(" Spatial coordinates not found in colData(). Tried: X/Y, *_visium, *_xenium, *_norm")
  }

  # -------------------------------
  # Process annotation column
  # -------------------------------
  raw_vec <- df[[annotation]]

  if (is.numeric(raw_vec)) {
    df$Annotation <- raw_vec
    annotation_type <- "numeric"
  } else {
    try_num <- suppressWarnings(as.numeric(as.character(raw_vec)))
    if (all(!is.na(try_num))) {
      df$Annotation <- try_num
      annotation_type <- "numeric"
    } else {
      df$Annotation <- as.character(raw_vec)
      annotation_type <- "categorical"
    }
  }

  # Determine whether to highlight
  use_highlight <- annotation_type == "categorical" && !is.null(highlight)
  if (use_highlight) {
    df$HighlightGroup <- ifelse(df$Annotation %in% highlight, df$Annotation, "Other")
    df$HighlightGroup[is.na(df$Annotation)] <- NA
    df$HighlightGroup <- factor(df$HighlightGroup, levels = c(highlight, "Other"))
  }

  # -------------------------------
  # Build plot
  # -------------------------------
  p <- ggplot(df, aes(x = X, y = Y)) +
    coord_fixed() +
    theme_void(base_size = base_size) +
    labs(title = title, color = annotation) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size * 0.8)
    )

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
        values = scales::hue_pal()(length(unique(df$Annotation))),
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
# 2. plot_umap: platform-specific UMAP plot from colData(obj)
# ------------------------------------------------------------
plot_umap <- function(obj,
                      platform = c("Visium", "Xenium"),
                      annotation = "Final_CT",
                      title = NULL,
                      point_size = 0.5,
                      legend_point_size = 4,
                      highlight = NULL) {
  platform <- match.arg(platform)
  df <- as.data.frame(colData(obj))
  df$cell <- rownames(df)

  # Validate annotation
  if (!(annotation %in% colnames(df))) {
    stop(paste("Annotation column", annotation, "not found in colData."))
  }

  # UMAP columns based on platform
  umap_cols <- switch(platform,
                      Visium = c("UMAP_1_Visium", "UMAP_2_Visium"),
                      Xenium = c("UMAP_1_Xenium", "UMAP_2_Xenium"))

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

  # Plot
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

  if (!is.null(highlight) && is.factor(df[[annotation]])) {
    p <- p + scale_color_manual(values = colors)
  }

  return(p)
}

# ------------------------------------------------------------
# 3. plot_gene_spatial: spatial gene expression visualization
# ------------------------------------------------------------
plot_gene_spatial <- function(obj,
                              gene,
                              platform = c("Visium", "Xenium"),
                              x_col = NULL,
                              y_col = NULL,
                              title = NULL,
                              point_size = 0.5,
                              low_color = "grey90",
                              high_color = "red") {
  platform <- match.arg(platform)

  # Spatial coordinate defaults
  if (is.null(x_col)) x_col <- ifelse(platform == "Visium", "X_visium_norm", "X_xenium_norm")
  if (is.null(y_col)) y_col <- ifelse(platform == "Visium", "Y_visium_norm", "Y_xenium_norm")

  # Construct gene name with suffix
  gene_full <- paste0(gene, "_", tolower(platform))

  # Extract expression values
  expr_matrix <- assays(obj)$X  # assumes counts are in assays(obj)$X
  if (!(gene_full %in% rownames(expr_matrix))) {
    stop(paste("Gene", gene_full, "not found in expression matrix."))
  }

  expr_vec <- expr_matrix[gene_full, ]
  df <- as.data.frame(colData(obj))
  df$Expr <- log1p(expr_vec)  # log(1 + x)
  df$cell <- rownames(df)

  # Validate coordinates
  if (!all(c(x_col, y_col) %in% colnames(df))) {
    stop(paste("Spatial coordinates", x_col, "and/or", y_col, "not found."))
  }

  plot_title <- if (!is.null(title)) title else paste("Spatial Expression of", gene)

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
# 4. rank_cross_platform_gene_correlation: top correlated genes
# ------------------------------------------------------------
rank_cross_platform_gene_correlation <- function(
  obj,
  top_genes = 5,
  layer = "X"
) {
  expr <- assay(obj, layer)
  
  # Identify base gene names with both _xenium and _visium rows
  all_genes <- rownames(expr)
  base_genes <- intersect(
    gsub("_xenium$", "", all_genes[grepl("_xenium$", all_genes)]),
    gsub("_visium$", "", all_genes[grepl("_visium$", all_genes)])
  )

  if (length(base_genes) == 0) stop("No matched gene pairs found.")

  # For each gene, correlate Xenium vs Visium row across all cells (columns)
  gene_corrs <- sapply(base_genes, function(gene) {
    x <- expr[paste0(gene, "_xenium"), ]
    y <- expr[paste0(gene, "_visium"), ]
    cor(log1p(x), log1p(y), method = "pearson", use = "pairwise.complete.obs")
  })

  # Sort and return
  gene_corrs <- sort(gene_corrs, decreasing = TRUE, na.last = NA)
  return(gene_corrs[1:min(top_genes, length(gene_corrs))])
}

# ------------------------------------------------------------
# 5. plot_cross_platform_gene_correlation: cumulative correlation
# ------------------------------------------------------------
plot_cross_platform_gene_correlation <- function(
  obj,
  gene_name,
  ID1_col = "Visium_ID",
  ID2_col = "Xenium_ID",
  hline = 0.9,
  y_range = c(0, 1),
  title = NULL
) {
  coldata <- colData(obj)
  if (!all(c(ID1_col, ID2_col) %in% colnames(coldata))) {
    stop("colData must include both ID1_col and ID2_col.")
  }

  # Suffix inference
  suffix1 <- paste0("_", tolower(gsub("_ID$", "", ID1_col)))
  suffix2 <- paste0("_", tolower(gsub("_ID$", "", ID2_col)))

  gene1 <- paste0(gene_name, suffix1)
  gene2 <- paste0(gene_name, suffix2)

  expr <- assays(obj)$X
  if (!(gene1 %in% rownames(expr)) || !(gene2 %in% rownames(expr))) {
    stop(paste("Gene", gene_name, "not found with expected suffixes in expression matrix."))
  }

  # Matched cells
  matched_df <- coldata[!is.na(coldata[[ID1_col]]) & !is.na(coldata[[ID2_col]]), ]
  matched_cells <- rownames(matched_df)
  if (length(matched_cells) < 2) stop("Fewer than 2 matched spots available.")

  # Extract log expression
  log1 <- log1p(as.numeric(expr[gene1, matched_cells]))
  log2 <- log1p(as.numeric(expr[gene2, matched_cells]))

  abs_diff <- abs(log1 - log2)
  order_idx <- order(abs_diff)
  n_total <- length(order_idx)

  cor_values <- numeric(n_total)
  for (i in 1:n_total) {
    x_sub <- log1[order_idx[1:i]]
    y_sub <- log2[order_idx[1:i]]
    cor_values[i] <- if (length(x_sub) >= 2) cor(x_sub, y_sub, method = "pearson") else NA
  }

  df_plot <- data.frame(
    Top_N = 1:n_total,
    Correlation = cor_values
  )

  p <- ggplot(df_plot, aes(x = Top_N, y = Correlation)) +
    geom_line(color = "steelblue") +
    geom_point(size = 0.6, color = "black") +
    theme_minimal(base_size = 14) +
    labs(
      title = ifelse(is.null(title), paste0(gene_name, ": Correlation vs Top Matched Spots"), title),
      x = "Number of Top Matched Spots (N)",
      y = "Pearson Correlation"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.05))) +
    scale_y_continuous(limits = y_range)

  if (!is.null(hline)) {
    p <- p + geom_hline(yintercept = hline, linetype = "dashed", color = "red")
  }

  return(p)
}

# ------------------------------------------------------------
# 6. scTriangulate_summarize_celltype_stability: cell type metrics
# ------------------------------------------------------------
scTriangulate_summarize_celltype_stability <- function(sce, celltype_labels, show_all_metrics = FALSE) {
  
  stopifnot(inherits(sce, "SingleCellExperiment"))
  stopifnot(is.list(celltype_labels))

  df_list <- lapply(names(celltype_labels), function(annot_col) {
    if (!(annot_col %in% colnames(colData(sce)))) return(NULL)
    labels <- as.character(colData(sce)[[annot_col]])
    selected_cells <- labels %in% celltype_labels[[annot_col]]
    if (sum(selected_cells, na.rm = TRUE) == 0) return(NULL)
    meta <- as.data.frame(colData(sce)[selected_cells, , drop = FALSE])
    meta$Cluster <- paste0(annot_col, "@", labels[selected_cells])
    meta$Platform <- ifelse(grepl("visium", rownames(meta), ignore.case = TRUE), "Visium", "Xenium")
    return(meta)
  })

  df_combined <- bind_rows(df_list)

  if (nrow(df_combined) == 0) return(tibble())

  summary_df <- df_combined %>%
    group_by(Cluster, Platform)

  if (show_all_metrics) {
    summary_df <- summary_df %>%
      select(Cluster, Platform, confidence, starts_with("SCCAF"), starts_with("tfidf"), starts_with("doublet"), starts_with("reassign."), ends_with("_shapley"))
  } else {
    summary_df <- summary_df %>%
      select(Cluster, Platform, confidence, ends_with("_shapley"))
  }

  summary_df <- summary_df %>%
    summarize(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
    
  # --- CORRECTED: Use order() for robust sorting ---
  if (nrow(summary_df) > 0) {
    # Isolate the shapley columns
    shapley_cols <- names(summary_df)[grepl("_shapley$", names(summary_df))]
    # Isolate the other, non-shapley columns
    other_cols <- setdiff(names(summary_df), shapley_cols)
    
    # Get the shapley values from the first row as a numeric vector
    shapley_values_vector <- as.numeric(summary_df[1, shapley_cols])
    
    # Get the sorting order (indices) of the columns
    col_order <- order(shapley_values_vector, decreasing = TRUE)
    
    # Use the ordering indices to sort the original column names
    sorted_shapley_names <- shapley_cols[col_order]
    
    # Create the final, sorted column order
    final_col_order <- c(other_cols, sorted_shapley_names)
    
    # Reorder the dataframe
    summary_df <- summary_df[, final_col_order]
  }
  
  return(summary_df)
}

# ------------------------------------------------------------
# 7. rank_marker_specificity: enrichment of gene per cluster
# ------------------------------------------------------------
rank_marker_specificity <- function(sce,
                                    gene_name = "COL1A1",
                                    cluster_cols = c("pruned"),
                                    assay_name = "X",
                                    platforms = c("Xenium", "Visium"),
                                    top_n = NULL) {
  # -------------------------------
  # Step 1: Select appropriate gene variant
  # -------------------------------
  gene_variants <- rownames(sce)
  full_gene <- if (paste0(gene_name, "_xenium") %in% gene_variants && "Xenium" %in% platforms) {
    paste0(gene_name, "_xenium")
  } else if (paste0(gene_name, "_visium") %in% gene_variants && "Visium" %in% platforms) {
    paste0(gene_name, "_visium")
  } else if (gene_name %in% gene_variants) {
    gene_name
  } else {
    stop(paste("Gene", gene_name, "not found in the dataset for selected platform(s)."))
  }

  # -------------------------------
  # Step 2: Extract expression & filter by platform
  # -------------------------------
  expr_vector <- assay(sce, assay_name)[full_gene, ]
  meta <- as.data.frame(colData(sce))
  meta$expr <- expr_vector

  if ("Platform" %in% colnames(meta)) {
    meta <- meta[meta$Platform %in% platforms, ]
  }

  # -------------------------------
  # Step 3: Compute specificity per cluster column
  # -------------------------------
  get_specificity_stats <- function(cluster_col) {
    if (!(cluster_col %in% colnames(meta))) return(NULL)

    df <- data.frame(expr = meta$expr, cluster = as.character(meta[[cluster_col]]))
    df$cluster[is.na(df$cluster)] <- "NA"

    result_list <- lapply(unique(df$cluster), function(cl) {
      in_group <- df$expr[df$cluster == cl]
      out_group <- df$expr[df$cluster != cl]
      if (length(in_group) < 2 || length(out_group) < 2) return(NULL)

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

  # -------------------------------
  # Step 4: Aggregate & rank
  # -------------------------------
  stats_list <- lapply(cluster_cols, get_specificity_stats)
  combined_stats <- bind_rows(stats_list)
  combined_stats$P_Adj <- p.adjust(combined_stats$P_Value, method = "fdr")

  ranked <- combined_stats %>%
    arrange(P_Adj, desc(Mean_Expr))

  if (!is.null(top_n)) {
    return(head(ranked, top_n))
  } else {
    return(ranked)
  }
}

# ------------------------------------------------------------
# 8. find_top_markers_for_celltype: markers per celltype/platform
# ------------------------------------------------------------
find_top_markers_for_celltype <- function(obj,
                                          celltype,
                                          annotation_col,
                                          platform_suffixes,
                                          top_n = 10) {
  expr_mat <- assay(obj, "X")  # genes x cells
  meta <- as.data.frame(colData(obj))
  meta$annotation <- meta[[annotation_col]]
  
  # Step 1: Filter genes by suffix match
  pattern <- paste0("_(", paste(platform_suffixes, collapse = "|"), ")$")
  gene_mask <- grepl(pattern, rownames(expr_mat))
  if (sum(gene_mask) == 0) stop("No genes found for platform suffixes: ", paste(platform_suffixes, collapse = ", "))
  expr_mat <- expr_mat[gene_mask, , drop = FALSE]
  
  # Step 2: Group definition
  target_mask <- meta$annotation == celltype
  other_mask  <- meta$annotation != celltype
  
  if (sum(target_mask) < 5 || sum(other_mask) < 5) {
    stop("Not enough cells for comparison.")
  }
  
  # Step 3: Compute logFC and t-test
  results <- apply(expr_mat, 1, function(gene_expr) {
    expr_b <- gene_expr[target_mask]
    expr_o <- gene_expr[other_mask]
    if (mean(expr_b) > 0 && mean(expr_o) > 0) {
      logFC <- log1p(mean(expr_b)) - log1p(mean(expr_o))
      pval <- tryCatch(t.test(expr_b, expr_o, alternative = "greater")$p.value, error = function(e) NA)
      return(c(logFC, pval))
    } else {
      return(c(NA, NA))
    }
  })
  
  result_df <- data.frame(
    gene = gsub(pattern, "", rownames(expr_mat)),
    logFC = results[1, ],
    P_Value = results[2, ]
  ) %>%
    filter(!is.na(logFC) & !is.na(P_Value)) %>%
    mutate(P_Adj = p.adjust(P_Value, method = "fdr")) %>%
    arrange(P_Adj, desc(logFC)) %>%
    slice_head(n = top_n)
  
  return(result_df)
}

# ------------------------------------------------------------
# 9. plot_ARI_dotplot: dotplot for annotation concordance (ARI)
# ------------------------------------------------------------
plot_ARI_dotplot <- function(obj,
                             annotation1,
                             annotation2,
                             title = NULL,
                             size_range = c(1, 10)) {
  # -------------------------------
  # Step 1: Extract annotations
  # -------------------------------
  meta <- as.data.frame(colData(obj))
  if (!(annotation1 %in% colnames(meta)) || !(annotation2 %in% colnames(meta))) {
    stop("One or both annotation columns not found.")
  }

  labels1 <- as.character(meta[[annotation1]])
  labels2 <- as.character(meta[[annotation2]])
  mask <- !is.na(labels1) & !is.na(labels2)
  labels1 <- labels1[mask]
  labels2 <- labels2[mask]

  # -------------------------------
  # Step 2: Compute ARI
  # -------------------------------
  ari_score <- adjustedRandIndex(labels1, labels2)

  # -------------------------------
  # Step 3: Cross-tabulation
  # -------------------------------
  cross_tab <- table(X = labels1, Y = labels2)
  df <- as.data.frame(cross_tab) %>% filter(Freq > 0)

  # -------------------------------
  # Step 4: Dot plot
  # -------------------------------
  plot_title <- ifelse(
    is.null(title),
    sprintf("%s vs %s (ARI = %.3f)", annotation1, annotation2, ari_score),
    sprintf("%s (ARI = %.3f)", title, ari_score)
  )

  p <- ggplot(df, aes(x = X, y = Y, size = Freq)) +
    geom_point(alpha = 0.8, color = "steelblue") +
    scale_size(range = size_range) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
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
# 10. plot_marker_expression_dotplot: gene expression dotplot
# ------------------------------------------------------------
plot_marker_expression_dotplot <- function(obj,
                                           genes,
                                           platforms = c("visium", "xenium"),
                                           annotation = "pruned",
                                           assay_name = "X",
                                           min_pct = 0.05,
                                           max_genes = 50,
                                           rotate_x = TRUE) {
  platforms <- match.arg(platforms, choices = c("visium", "xenium"), several.ok = TRUE)
  suffixes <- paste0("_", platforms)

  stopifnot(annotation %in% colnames(colData(obj)))

  expr <- assays(obj)[[assay_name]]
  coldata <- as.data.frame(colData(obj))

  # Build platform-specific gene names
  gene_map <- expand.grid(gene = genes, suffix = suffixes, stringsAsFactors = FALSE)
  gene_map$full_gene <- paste0(gene_map$gene, gene_map$suffix)
  gene_map$platform <- gsub("_", "", gene_map$suffix)
  gene_map <- gene_map[gene_map$full_gene %in% rownames(expr), ]

  if (nrow(gene_map) == 0) stop("No matching Visium/Xenium genes found.")

  # Assemble long format expression
  expr_list <- lapply(seq_len(nrow(gene_map)), function(i) {
    g <- gene_map[i, ]
    vec <- as.numeric(expr[g$full_gene, ])
    data.frame(
      expr = vec,
      annotation = coldata[[annotation]],
      gene = paste0(g$gene, " (", tools::toTitleCase(g$platform), ")"),
      stringsAsFactors = FALSE
    )
  })

  df_long <- bind_rows(expr_list)

  # Dotplot stats
  dot_data <- df_long %>%
    group_by(annotation, gene) %>%
    summarize(
      avg_expr = mean(log1p(expr), na.rm = TRUE),
      pct_expr = mean(expr > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(pct_expr >= min_pct)

  if (!is.null(max_genes) && length(unique(dot_data$gene)) > max_genes) {
    top_genes <- dot_data %>%
      group_by(gene) %>%
      summarize(score = max(avg_expr * pct_expr)) %>%
      arrange(desc(score)) %>%
      slice_head(n = max_genes) %>%
      pull(gene)

    dot_data <- dot_data %>% filter(gene %in% top_genes)
  }

  # Plot
  p <- ggplot(dot_data, aes(x = gene, y = annotation)) +
    geom_point(aes(size = pct_expr * 100, color = avg_expr)) +
    scale_size(name = "% Expressed", range = c(1, 6)) +
    scale_color_gradient(name = "Avg log1p(expr)", low = "grey90", high = "red") +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = if (rotate_x) 45 else 0, hjust = 1),
      axis.title = element_blank()
    )

  return(p)
}

# ------------------------------------------------------------
# 11. plot_cluster_composition_groupedbar: proportion barplot
# ------------------------------------------------------------
plot_cluster_composition_groupedbar <- function(obj,
                                                annotation,
                                                palette = NULL,
                                                title = NULL,
                                                rotate_x = TRUE) {
  meta <- as.data.frame(colData(obj))
  
  # Validate annotations
  stopifnot(length(annotation) %in% c(1, 2))
  stopifnot(all(annotation %in% colnames(meta)))
  
  # Handle one annotation
  if (length(annotation) == 1) {
    df_tmp <- meta %>% filter(!is.na(.data[[annotation[1]]]))
    df <- df_tmp %>%
      dplyr::count(annotation_value = !!sym(annotation[1])) %>%
      mutate(Proportion = n / sum(n), annotation_source = annotation[1])
  }
  
  # Handle two annotations
  if (length(annotation) == 2) {
    meta1 <- meta %>% filter(!is.na(.data[[annotation[1]]]))
    meta2 <- meta %>% filter(!is.na(.data[[annotation[2]]]))
    
    labels1 <- unique(trimws(as.character(meta1[[annotation[1]]])))
    labels2 <- unique(trimws(as.character(meta2[[annotation[2]]])))
    shared <- intersect(labels1, labels2)
    
    message(length(shared), " overlapping labels: ", paste(shared, collapse = ", "))
    
    df_tmp1 <- meta1
    df_tmp2 <- meta2
    
    df1 <- df_tmp1 %>%
      dplyr::count(annotation_value = !!sym(annotation[1])) %>%
      mutate(Proportion = n / sum(n), annotation_source = annotation[1])
    
    df2 <- df_tmp2 %>%
      dplyr::count(annotation_value = !!sym(annotation[2])) %>%
      mutate(Proportion = n / sum(n), annotation_source = annotation[2])
    
    df <- bind_rows(df1, df2)
  }
  
  # Sort x-axis by total abundance
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
  
  # Default palette for Visium vs Xenium RCTD
  if (is.null(palette) &&
      all(sort(annotation) == sort(c("Visium_RCTD_LungMap_ref", "Xenium_RCTD_LungMap_ref")))) {
    palette <- c(
      "Visium_RCTD_LungMap_ref" = "blue",
      "Xenium_RCTD_LungMap_ref" = "red"
    )
  }
  
  # Plot
  p <- ggplot(df, aes(x = annotation_value, y = Proportion, fill = annotation_source)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    labs(
      x = "Cell Type",
      y = "Proportion",
      fill = "Annotation Source",
      title = ifelse(is.null(title), "Cluster Composition Comparison", title)
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = if (rotate_x) 45 else 0, hjust = 1),
      panel.grid.major.x = element_blank()
    )
  
  if (!is.null(palette)) {
    p <- p + scale_fill_manual(values = palette)
  }
  
  return(p)
}


# ------------------------------------------------------------
# 12. compare_annotations_via_heatmap: confusion matrix heatmap
# ------------------------------------------------------------
compare_annotations_via_heatmap <- function(obj,
                                            annotation1,
                                            annotation2,
                                            normalize = c("none", "row", "column", "both"),
                                            palette = "Blues",
                                            title = NULL) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)

  # Match normalization option
  normalize <- match.arg(normalize)

  # Extract metadata
  meta <- as.data.frame(SummarizedExperiment::colData(obj))
  stopifnot(annotation1 %in% colnames(meta), annotation2 %in% colnames(meta))

  df <- meta %>%
    filter(!is.na(.data[[annotation1]]), !is.na(.data[[annotation2]])) %>%
    count(Annotation1 = .data[[annotation1]], Annotation2 = .data[[annotation2]])

  # Pivot to wide matrix
  mat <- tidyr::pivot_wider(df, names_from = Annotation2, values_from = n, values_fill = 0) %>%
    tibble::column_to_rownames("Annotation1") %>%
    as.matrix()

  # Normalize
  if (normalize == "row") {
    mat <- sweep(mat, 1, rowSums(mat), FUN = "/")
  } else if (normalize == "column") {
    mat <- sweep(mat, 2, colSums(mat), FUN = "/")
  } else if (normalize == "both") {
    mat <- mat / sum(mat)
  }

  # Convert back to long format
  df_plot <- as.data.frame(as.table(mat))
  colnames(df_plot) <- c("Annotation1", "Annotation2", "Value")

  p <- ggplot(df_plot, aes(x = Annotation2, y = Annotation1, fill = Value)) +
    geom_tile(color = "white") +
    scale_fill_distiller(palette = palette, direction = 1, na.value = "grey90") +
    geom_text(aes(label = ifelse(Value == 0, "", comma(round(Value, 2)))), size = 3.2) +
    theme_minimal(base_size = 12) +
    labs(
      x = annotation2,
      y = annotation1,
      fill = ifelse(normalize == "none", "Count", "Proportion"),
      title = ifelse(is.null(title), "Annotation Comparison Heatmap", title)
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )

  return(p)
}

# ------------------------------------------------------------
# 13. visualize_matched_pairs_scatter: matched gene expression
# ------------------------------------------------------------
visualize_matched_pairs_scatter <- function(obj,
                                            gene,
                                            layer = "X",
                                            log1p = TRUE,
                                            point_size = 0.4,
                                            alpha = 0.6,
                                            show_correlation = FALSE,
                                            title = NULL) {
  # Extract expression matrix
  expr_mat <- assay(obj, layer)
  
  gene_visium <- paste0(gene, "_visium")
  gene_xenium <- paste0(gene, "_xenium")
  
  if (!(gene_visium %in% rownames(expr_mat)) || !(gene_xenium %in% rownames(expr_mat))) {
    stop(paste("Gene", gene, "not found in both Visium and Xenium formats"))
  }
  
  visium_expr <- expr_mat[gene_visium, ]
  xenium_expr <- expr_mat[gene_xenium, ]
  
  if (log1p) {
    visium_expr <- log1p(visium_expr)
    xenium_expr <- log1p(xenium_expr)
  }
  
  df <- data.frame(
    Xenium = xenium_expr,
    Visium = visium_expr
  )
  
  # Compute correlation
  corr <- cor(df$Xenium, df$Visium, method = "pearson")
  corr_text <- paste0("Pearson r = ", round(corr, 3))
  
  # Title logic
  plot_title <- if (!is.null(title)) title else paste0("Xenium vs Visium: ", gene)
  if (show_correlation) {
    plot_title <- paste0(plot_title, "\n", corr_text)
  }
  
  # Plot
  p <- ggplot(df, aes(x = Xenium, y = Visium)) +
    geom_point(size = point_size, alpha = alpha) +
    labs(
      x = "Xenium Expression",
      y = "Visium Expression",
      title = plot_title
    ) +
    theme_minimal(base_size = 13)
  
  return(p)
}
