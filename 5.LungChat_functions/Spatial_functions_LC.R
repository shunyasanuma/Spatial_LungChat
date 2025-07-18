# Spatial_functions_LC.R
# ============================================================================ #
# Combined Spatial Functions for Visualization and Analysis
# ============================================================================ #

# ---------------------------------------------------------------------------- #
# Load Required Libraries
# ---------------------------------------------------------------------------- #

suppressPackageStartupMessages({
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
  library(uuid)
  library(jsonlite)
  library(tools)
  library(ComplexHeatmap)
  library(circlize)
  library(Seurat)
  library(httr)
  library(patchwork)
  library(magick)
  library(cowplot)
})

# ---------------------------------------------------------------------------- #
# Internal Helper Functions
# ---------------------------------------------------------------------------- #

#' Generate a dynamic and unique filename base.
#' @param base_name A string for the base of the filename.
#' @param ... Key arguments from the parent function to include for descriptiveness.
#' @param max_items The maximum number of elements from a vector argument to include.
#' @return A sanitized, unique base string for filenames.
.generate_filename_base <- function(base_name, ..., max_items = 2) {
  args <- list(...)
  descriptive_parts <- sapply(names(args), function(arg_name) {
    x <- args[[arg_name]]
    if (is.null(x) || length(x) == 0 || all(is.na(x)) || !nzchar(x[[1]])) return(NULL)
    val <- if (length(x) > max_items) {
      paste(c(as.character(head(x, max_items)), "etc"), collapse = "-")
    } else {
      paste(as.character(x), collapse = "-")
    }
    gsub("[^a-zA-Z0-9_.-]", "_", val)
  })

  descriptive_string <- paste(unlist(descriptive_parts), collapse = "_")
  unique_id <- substr(uuid::UUIDgenerate(), 1, 5)

  filename_base <- if (nzchar(descriptive_string)) {
    paste(base_name, descriptive_string, unique_id, sep = "_")
  } else {
    paste(base_name, unique_id, sep = "_")
  }

  sanitized_base <- gsub("[[:punct:]]{2,}", "_", gsub(" ", "_", filename_base))
  sub("_$", "", sanitized_base)
}

#' Saves plots generated from base R or grid graphics to PNG and PDF.
#' @param plot_fun A zero-argument function that generates the plot.
#' @param base_filename The base name for the output files.
#' @param width The width of the plot.
#' @param height The height of the plot.
#' @param dpi The resolution for the PNG file.
.save_base_graphics_plot <- function(plot_fun, base_filename, width, height, dpi) {
  png_file <- paste0(base_filename, ".png")
  pdf_file <- paste0(base_filename, ".pdf")

  png(png_file, width = width, height = height, units = "in", res = dpi)
  plot_fun()
  dev.off()

  pdf(pdf_file, width = width, height = height)
  plot_fun()
  dev.off()

  return(list(
    pdf_path = tools::file_path_as_absolute(pdf_file),
    png_path = tools::file_path_as_absolute(png_file)
  ))
}

#' Create a ggplot object from an image file.
#' @param image_path Path to the image file (e.g., TIFF, PNG).
#' @param title The title for the plot.
#' @return A ggplot object displaying the image, or a placeholder on error.
.create_image_plot <- function(image_path, title) {
    # Check if file exists
    if (!file.exists(image_path)) {
        error_msg <- paste("File not found:\n", image_path)
        return(
            ggplot() +
            annotate("text", x = 0, y = 0, label = error_msg, size = 4, color = "red", lineheight = .9) +
            labs(title = title) +
            theme_void() +
            theme(
                plot.title = element_text(hjust = 0.5, face = "bold"),
                plot.background = element_rect(fill = "white", color = NA),
                panel.background = element_rect(fill = "#f0f0f0", color = "grey50", linetype = "dashed")
            )
        )
    }

    # Read image using magick
    img <- tryCatch({
        magick::image_read(image_path)
    }, error = function(e) {
        error_msg <- paste("Failed to read image:\n", basename(image_path))
        return(
            ggplot() +
            annotate("text", x = 0, y = 0, label = error_msg, size = 4, color = "red", lineheight = .9) +
            labs(title = title) +
            theme_void() +
            theme(
                plot.title = element_text(hjust = 0.5, face = "bold"),
                plot.background = element_rect(fill = "white", color = NA),
                panel.background = element_rect(fill = "#f0f0f0", color = "grey50", linetype = "dashed")
            )
        )
    })

    # Convert to ggplot object
    p <- magick::image_ggplot(img, interpolate = TRUE) +
        labs(title = title) +
        theme_void() +
        theme(
            plot.title = element_text(hjust = 0.5),
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA)
        )
    return(p)
}

# ---------------------------------------------------------------------------- #
# 1. General Spatial Plot
# ---------------------------------------------------------------------------- #
plot_spatial <- function(obj,
                         annotation,
                         title = NULL,
                         point_size = 0.5,
                         base_size = 12,
                         coord_flip = FALSE,
                         highlight = NULL,
                         plot_width = 12,
                         plot_height = 10,
                         plot_dpi = 300,
                         return_plot = FALSE) {
  # Validate annotation exists
  stopifnot(annotation %in% colnames(colData(obj)))
  df <- as.data.frame(colData(obj))

  # Find spatial coordinates
  coord_sets <- list(c("X", "Y"), c("X_visium", "Y_visium"), c("X_xenium", "Y_xenium"),
                     c("X_visium_norm", "Y_visium_norm"), c("X_xenium_norm", "Y_xenium_norm"))
  coord_found <- FALSE
  for (coords in coord_sets) {
    if (all(coords %in% colnames(df))) {
      df$X <- as.numeric(df[[coords[1]]])
      df$Y <- as.numeric(df[[coords[2]]])
      coord_found <- TRUE
      break
    }
  }
  if (!coord_found) stop("Spatial coordinates not found in colData().")

  # Process annotation column type (numeric vs categorical)
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

  # Handle categorical highlighting
  use_highlight <- annotation_type == "categorical" && !is.null(highlight)
  if (use_highlight) {
    df$HighlightGroup <- ifelse(df$Annotation %in% highlight, df$Annotation, "Other")
    df$HighlightGroup[is.na(df$Annotation)] <- NA
    df$HighlightGroup <- factor(df$HighlightGroup, levels = c(highlight, "Other"))
  }

  # Build the ggplot object
  p <- ggplot(df, aes(x = X, y = Y)) +
    coord_fixed() + theme_bw(base_size = base_size) +
    labs(title = title, color = annotation) +
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_text(size = base_size),
          legend.text = element_text(size = base_size * 0.8),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))

  # Add point layers based on annotation type
  if (annotation_type == "numeric") {
    p <- p + geom_point(aes(color = Annotation), size = point_size, na.rm = TRUE) +
      scale_color_gradient(low = "blue", high = "red")
  } else if (use_highlight) {
    pal <- c(setNames(scales::hue_pal()(length(highlight)), highlight), "Other" = "grey80")
    p <- p + geom_point(aes(color = HighlightGroup), size = point_size, na.rm = TRUE) +
      scale_color_manual(values = pal, drop = TRUE, na.translate = FALSE) +
      guides(color = guide_legend(override.aes = list(size = 4)))
  } else {
    df$Annotation <- factor(df$Annotation)
    p <- p + geom_point(aes(color = Annotation), size = point_size, na.rm = TRUE) +
      scale_color_manual(values = scales::hue_pal()(length(unique(na.omit(df$Annotation)))),
                         na.translate = FALSE) +
      guides(color = guide_legend(override.aes = list(size = 4)))
  }
  if (coord_flip) p <- p + scale_y_reverse()

  # Check if the plot object itself should be returned
  if (return_plot) {
    return(p)
  }

  # Save outputs and generate JSON
  base_filename <- .generate_filename_base("spatial", annotation = annotation, highlight = highlight)
  png_file <- paste0(base_filename, ".png")
  pdf_file <- paste0(base_filename, ".pdf")
  ggsave(png_file, plot = p, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
  ggsave(pdf_file, plot = p, width = plot_width, height = plot_height)

  output_paths <- list(plot_type = "spatial_plot",
                       pdf_path = tools::file_path_as_absolute(pdf_file),
                       png_path = tools::file_path_as_absolute(png_file))

  return(jsonlite::toJSON(output_paths, auto_unbox = TRUE, pretty = TRUE))
}

# ---------------------------------------------------------------------------- #
# 2. Plot Quad View: Spatial Plots with H&E and DAPI (REPAIRED)
# ---------------------------------------------------------------------------- #
#' Generates an aesthetically improved 2x2 plot with a shared legend panel.
#'
#' @param visium_obj A SingleCellExperiment object for Visium data.
#' @param xenium_obj A SingleCellExperiment object for Xenium data.
#' @param he_path A string path to the H&E image file (.tif).
#' @param dapi_path A string path to the DAPI image file (.png).
#' @param visium_annotation A string for the Visium metadata column.
#' @param xenium_annotation A string for the Xenium metadata column.
#' @param point_size The size of points in the spatial plots.
#' @param coord_flip_visium Logical, whether to flip the Y-axis for the Visium plot.
#' @param coord_flip_xenium Logical, whether to flip the Y-axis for the Xenium plot.
#' @param title An overall title for the combined plot.
#' @param plot_width The width of the output plot in inches.
#' @param plot_height The height of the output plot in inches.
#' @param plot_dpi The resolution of the output PNG file.
#' @return A JSON string containing the absolute paths to the generated plot files.
plot_quad_spatial_view <- function(visium_obj,
                                     xenium_obj,
                                     he_path,
                                     dapi_path,
                                     visium_annotation = "Visium_RCTD_LungMap_ref",
                                     xenium_annotation = "Final_CT",
                                     point_size = 0.5,
                                     coord_flip_visium = TRUE,
                                     coord_flip_xenium = FALSE,
                                     title = "Comparative Spatial View",
                                     plot_width = 16,
                                     plot_height = 10,
                                     plot_dpi = 300) {

    # 1. Generate Visium spatial plot (no legend, tight margin)
    p_visium <- plot_spatial(
        obj = visium_obj,
        annotation = visium_annotation,
        title = "Visium HD",
        point_size = point_size,
        coord_flip = coord_flip_visium,
        return_plot = TRUE
    ) + theme(legend.position = "none", plot.margin = margin(b = 0.1, unit = "cm"))

    # 2. Generate Xenium spatial plot (no legend, tight margin)
    p_xenium <- plot_spatial(
        obj = xenium_obj,
        annotation = xenium_annotation,
        title = "Xenium",
        point_size = point_size,
        coord_flip = coord_flip_xenium,
        return_plot = TRUE
    ) + theme(legend.position = "none", plot.margin = margin(b = 0.1, unit = "cm"))

    # 3. Create H&E and DAPI image plots (tight margins)
    p_he <- .create_image_plot(image_path = he_path, title = "H&E") +
        theme(plot.margin = margin(t = 0.1, unit = "cm"))
    p_dapi <- .create_image_plot(image_path = dapi_path, title = "DAPI") +
        theme(plot.margin = margin(t = 0.1, unit = "cm"))

    # 4. Create separate plots just to extract their legends
    p_vis_legend_dummy <- plot_spatial(visium_obj, visium_annotation, return_plot = TRUE) +
        guides(color = guide_legend(title = "Visium Annotation", ncol = 1, override.aes = list(size = 3)))
    vis_legend <- get_legend(p_vis_legend_dummy)

    p_xen_legend_dummy <- plot_spatial(xenium_obj, xenium_annotation, return_plot = TRUE) +
        guides(color = guide_legend(title = "Xenium Annotation", ncol = 1, override.aes = list(size = 3)))
    xen_legend <- get_legend(p_xen_legend_dummy)

    # 5. Combine the plots and legends using patchwork
    plot_grid <- (p_visium | p_xenium) / (p_he | p_dapi)

    # FIX: Use wrap_plots() to combine legend grobs correctly before adding to the layout
    legend_panel <- wrap_plots(vis_legend, xen_legend, ncol = 1)

    combined_plot <- (plot_grid | legend_panel) +
        plot_layout(widths = c(3.5, 1)) + # Assign relative widths (plots vs legends)
        plot_annotation(
            title = title,
            theme = theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
        )

    # 6. Save final combined plot and generate JSON output
    base_filename <- .generate_filename_base("quad_view", v_annot = visium_annotation, x_annot = xenium_annotation)
    png_file <- paste0(base_filename, ".png")
    pdf_file <- paste0(base_filename, ".pdf")

    ggsave(png_file, plot = combined_plot, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
    ggsave(pdf_file, plot = combined_plot, width = plot_width, height = plot_height)

    output_paths <- list(
        plot_type = "quad_spatial_view",
        pdf_path = tools::file_path_as_absolute(pdf_file),
        png_path = tools::file_path_as_absolute(png_file)
    )

    return(jsonlite::toJSON(output_paths, auto_unbox = TRUE, pretty = TRUE))
}

# ---------------------------------------------------------------------------- #
# 3. Platform-Specific UMAP Plot
# ---------------------------------------------------------------------------- #
plot_umap <- function(obj,
                      platform = c("Xenium", "Visium"),
                      annotation = "Final_CT",
                      title = NULL,
                      point_size = 0.5,
                      legend_point_size = 4,
                      highlight = NULL,
                      plot_width = 12,
                      plot_height = 10,
                      plot_dpi = 300) {
  platform <- match.arg(platform)
  df <- as.data.frame(colData(obj))

  # Validate columns
  if (!(annotation %in% colnames(df))) stop(paste("Annotation column", annotation, "not in colData."))
  umap_cols <- switch(platform, Visium = c("UMAP_1_Visium", "UMAP_2_Visium"),
                      Xenium = c("UMAP_1_Xenium", "UMAP_2_Xenium"))
  if (!all(umap_cols %in% colnames(df))) stop(paste("UMAP columns not found for", platform))

  # Handle highlighting
  if (!is.null(highlight) && is.factor(df[[annotation]])) {
    df$highlight_flag <- ifelse(df[[annotation]] %in% highlight, as.character(df[[annotation]]), "Other")
    df$highlight_flag <- factor(df$highlight_flag, levels = c(highlight, "Other"))
    color_aes <- "highlight_flag"
    n_colors <- length(highlight)
    base_colors <- if (n_colors <= 9) RColorBrewer::brewer.pal(n = max(n_colors, 3), "Set1")[seq_len(n_colors)] else grDevices::rainbow(n_colors)
    colors <- c(setNames(base_colors, highlight), Other = "grey80")
  } else {
    color_aes <- annotation
    colors <- NULL
  }

  # Build plot
  plot_title <- if (!is.null(title)) title else paste(platform, "UMAP:", annotation)
  p <- ggplot(df, aes(x = .data[[umap_cols[1]]], y = .data[[umap_cols[2]]], color = .data[[color_aes]])) +
    geom_point(size = point_size) + coord_fixed() + theme_bw() +
    labs(title = plot_title, color = annotation) +
    guides(color = guide_legend(override.aes = list(size = legend_point_size))) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right",
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))
  if (!is.null(colors)) p <- p + scale_color_manual(values = colors)

  # Save outputs and generate JSON
  base_filename <- .generate_filename_base("umap", platform = platform, annotation = annotation, highlight = highlight)
  png_file <- paste0(base_filename, ".png")
  pdf_file <- paste0(base_filename, ".pdf")
  ggsave(png_file, plot = p, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
  ggsave(pdf_file, plot = p, width = plot_width, height = plot_height)

  output_paths <- list(plot_type = "umap_plot",
                       pdf_path = tools::file_path_as_absolute(pdf_file),
                       png_path = tools::file_path_as_absolute(png_file))

  return(jsonlite::toJSON(output_paths, auto_unbox = TRUE, pretty = TRUE))
}

# ---------------------------------------------------------------------------- #
# 4. Spatial Gene Expression Visualization
# ---------------------------------------------------------------------------- #
plot_gene_spatial <- function(obj,
                              gene,
                              platform = c("Xenium", "Visium"),
                              x_col = NULL,
                              y_col = NULL,
                              title = NULL,
                              point_size = 0.5,
                              low_color = "grey90",
                              high_color = "red",
                              plot_width = 12,
                              plot_height = 10,
                              plot_dpi = 300) {
  platform <- match.arg(platform)

  # Define default coordinates
  if (is.null(x_col)) x_col <- ifelse(platform == "Visium", "X_visium_norm", "X_xenium_norm")
  if (is.null(y_col)) y_col <- ifelse(platform == "Visium", "Y_visium_norm", "Y_xenium_norm")

  # Construct gene name with platform suffix
  platform_suffix <- paste0("_", tolower(platform))
  gene_full <- if (endsWith(gene, platform_suffix)) gene else paste0(gene, platform_suffix)

  # Extract data
  expr_matrix <- assays(obj)$X
  if (!(gene_full %in% rownames(expr_matrix))) stop(paste("Gene", gene_full, "not in expression matrix."))
  df <- as.data.frame(colData(obj))
  df$Expr <- log1p(as.numeric(expr_matrix[gene_full, , drop = TRUE]))
  if (!all(c(x_col, y_col) %in% colnames(df))) stop(paste("Spatial coordinates not found."))

  # Build plot
  plot_title <- if (!is.null(title)) title else paste("Spatial Expression of", gene)
  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[["Expr"]])) + # Modernized aes() call
    geom_point(size = point_size) +
    scale_color_gradient(low = low_color, high = high_color) +
    coord_fixed() + theme_bw() +
    labs(title = plot_title, color = "log1p(expr)") +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right",
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))

  # Save outputs and generate JSON
  base_filename <- .generate_filename_base("gene_spatial", gene = gene, platform = platform)
  png_file <- paste0(base_filename, ".png")
  pdf_file <- paste0(base_filename, ".pdf")
  ggsave(png_file, plot = p, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
  ggsave(pdf_file, plot = p, width = plot_width, height = plot_height)

  output_paths <- list(plot_type = "gene_spatial_plot",
                       pdf_path = tools::file_path_as_absolute(pdf_file),
                       png_path = tools::file_path_as_absolute(png_file))

  return(jsonlite::toJSON(output_paths, auto_unbox = TRUE, pretty = TRUE))
}

# ---------------------------------------------------------------------------- #
# 5. Rank Cross-Platform Gene Correlation
# ---------------------------------------------------------------------------- #
rank_cross_platform_gene_correlation <- function(obj,
                                                 top_genes = 5,
                                                 layer = "X",
                                                 plot_width = 12,
                                                 plot_height = 10,
                                                 plot_dpi = 300) {
  # Get expression data and find shared genes
  if (!layer %in% assayNames(obj)) stop(paste("Assay '", layer, "' not found."))
  expr <- assay(obj, layer)
  all_genes <- rownames(expr)
  base_genes <- intersect(gsub("_xenium$", "", all_genes[grepl("_xenium$", all_genes)]),
                          gsub("_visium$", "", all_genes[grepl("_visium$", all_genes)]))
  if (length(base_genes) == 0) stop("No matched gene pairs found between platforms.")

  # Calculate correlations, suppressing warnings for zero-variance genes
  gene_corrs <- suppressWarnings({
    sapply(base_genes, function(gene) {
      x <- expr[paste0(gene, "_xenium"), ]
      y <- expr[paste0(gene, "_visium"), ]
      cor(log1p(x), log1p(y), method = "pearson", use = "pairwise.complete.obs")
    })
  })
  results_df <- data.frame(Correlation = gene_corrs[!is.na(gene_corrs)]) %>%
    tibble::rownames_to_column("Gene") %>%
    arrange(desc(Correlation))
  top_results_df <- head(results_df, n = top_genes)

  # Build plot
  p <- ggplot(top_results_df, aes(x = reorder(Gene, -Correlation), y = Correlation)) +
    geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "Top Cross-Platform Gene Correlations (Xenium vs. Visium)",
         x = "Gene", y = "Pearson Correlation") +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))

  # Save outputs and generate JSON
  base_filename <- .generate_filename_base("cross_platform_corr", top_n = top_genes)
  tsv_file <- paste0(base_filename, ".tsv")
  png_file <- paste0(base_filename, ".png")
  pdf_file <- paste0(base_filename, ".pdf")
  write.table(top_results_df, tsv_file, sep = "\t", row.names = FALSE, quote = FALSE)
  ggsave(png_file, plot = p, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
  ggsave(pdf_file, plot = p, width = plot_width, height = plot_height)

  output_paths <- list(plot_type = "cross_platform_correlation_rank",
                       tsv_path = tools::file_path_as_absolute(tsv_file),
                       pdf_path = tools::file_path_as_absolute(pdf_file),
                       png_path = tools::file_path_as_absolute(png_file))

  return(jsonlite::toJSON(output_paths, auto_unbox = TRUE, pretty = TRUE))
}

# ---------------------------------------------------------------------------- #
# 6. Plot Cumulative Cross-Platform Correlation
# ---------------------------------------------------------------------------- #
plot_cross_platform_gene_correlation <- function(obj,
                                                 gene_name,
                                                 ID1_col = "Visium_ID",
                                                 ID2_col = "Xenium_ID",
                                                 hline = 0.9,
                                                 y_range = c(0, 1),
                                                 title = NULL,
                                                 plot_width = 10,
                                                 plot_height = 8,
                                                 plot_dpi = 300) {
  # Extract data and find matched cells
  coldata <- as.data.frame(colData(obj))
  if (!all(c(ID1_col, ID2_col) %in% colnames(coldata))) stop("ID columns not found.")
  matched_df <- coldata[!is.na(coldata[[ID1_col]]) & !is.na(coldata[[ID2_col]]), ]
  if (nrow(matched_df) < 2) stop("Fewer than 2 matched spots available.")

  # Construct platform-specific gene names
  suffix1 <- paste0("_", tolower(gsub("_ID$", "", ID1_col)))
  suffix2 <- paste0("_", tolower(gsub("_ID$", "", ID2_col)))
  base_gene <- if (endsWith(gene_name, suffix1) || endsWith(gene_name, suffix2)) {
    sub(paste0(suffix1, "|", suffix2), "", gene_name)
  } else {
    gene_name
  }
  gene1 <- paste0(base_gene, suffix1)
  gene2 <- paste0(base_gene, suffix2)
  expr <- assays(obj)$X
  if (!all(c(gene1, gene2) %in% rownames(expr))) stop("Gene pair not found in expression matrix.")

  # Calculate cumulative correlation
  log1 <- log1p(as.numeric(expr[gene1, rownames(matched_df)]))
  log2 <- log1p(as.numeric(expr[gene2, rownames(matched_df)]))
  order_idx <- order(abs(log1 - log2))
  n_total <- length(order_idx)
  cor_values <- suppressWarnings({
    sapply(seq_len(n_total), function(i) {
      if (i < 2) return(NA)
      cor(log1[order_idx[1:i]], log2[order_idx[1:i]], method = "pearson")
    })
  })
  df_plot <- data.frame(Top_N = seq_len(n_total), Correlation = cor_values)

  # Build plot
  plot_title <- if (is.null(title)) paste0(gene_name, ": Correlation vs Top Matched Spots") else title
  p <- ggplot(df_plot, aes(x = Top_N, y = Correlation)) +
    geom_line(color = "steelblue") + geom_point(size = 0.6, color = "black") +
    labs(title = plot_title, x = "Number of Top Matched Spots (N)", y = "Pearson Correlation") +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.05))) +
    scale_y_continuous(limits = y_range) +
    theme_bw(base_size = 14) +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))
  if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, linetype = "dashed", color = "red")

  # Save outputs and generate JSON
  base_filename <- .generate_filename_base("cumulative_corr", gene = gene_name)
  png_file <- paste0(base_filename, ".png")
  pdf_file <- paste0(base_filename, ".pdf")
  
  # Suppress harmless warning about removing rows with NA values
  suppressWarnings({
    ggsave(png_file, plot = p, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
    ggsave(pdf_file, plot = p, width = plot_width, height = plot_height)
  })

  output_paths <- list(plot_type = "cumulative_correlation_plot",
                       pdf_path = tools::file_path_as_absolute(pdf_file),
                       png_path = tools::file_path_as_absolute(png_file))

  return(jsonlite::toJSON(output_paths, auto_unbox = TRUE, pretty = TRUE))
}

# ---------------------------------------------------------------------------- #
# 7. Rank Marker Specificity
# ---------------------------------------------------------------------------- #
rank_marker_specificity <- function(sce,
                                    gene_name = "COL1A1",
                                    cluster_cols = c("pruned"),
                                    assay_name = "X",
                                    platforms = c("Xenium", "Visium"),
                                    top_n = NULL) {
  # Select appropriate gene variant based on platform
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

  # Extract data
  meta <- as.data.frame(colData(sce))
  meta$expr <- assay(sce, assay_name)[full_gene, ]
  if ("Platform" %in% colnames(meta)) meta <- meta[meta$Platform %in% platforms, ]

  # Compute specificity stats for each cluster column
  get_specificity_stats <- function(cluster_col) {
    if (!(cluster_col %in% colnames(meta))) return(NULL)
    df <- data.frame(expr = meta$expr, cluster = as.character(meta[[cluster_col]]))
    df$cluster[is.na(df$cluster)] <- "NA"
    result_list <- lapply(unique(df$cluster), function(cl) {
      in_group <- df$expr[df$cluster == cl]
      if (length(in_group) < 2) return(NULL)
      t_result <- tryCatch(t.test(in_group, df$expr[df$cluster != cl], alternative = "greater"), error = function(e) NULL)
      data.frame(Cluster = cl, Mean_Expr = mean(in_group), Pct_Expr = mean(in_group > 0) * 100,
                 P_Value = if (!is.null(t_result)) t_result$p.value else NA,
                 Annotation_Source = cluster_col)
    })
    do.call(rbind, result_list)
  }

  # Aggregate results and rank
  combined_stats <- bind_rows(lapply(cluster_cols, get_specificity_stats))
  if (nrow(combined_stats) == 0) {
    warning("No valid clusters found to calculate statistics."); return(NULL)
  }
  ranked <- combined_stats %>%
    mutate(P_Adj = p.adjust(P_Value, method = "fdr")) %>%
    arrange(P_Adj, desc(Mean_Expr))
  final_result <- if (!is.null(top_n)) head(ranked, top_n) else ranked

  # Save output and generate JSON
  base_filename <- .generate_filename_base("marker_specificity", gene = gene_name, platform = platforms)
  tsv_file <- paste0(base_filename, ".tsv")
  write.table(final_result, tsv_file, sep = "\t", row.names = FALSE, quote = FALSE)

  output_paths <- list(data_type = "marker_specificity_table",
                       tsv_path = tools::file_path_as_absolute(tsv_file))

  return(jsonlite::toJSON(output_paths, auto_unbox = TRUE, pretty = TRUE))
}

# ---------------------------------------------------------------------------- #
# 8. Plot Annotation Concordance (ARI)
# ---------------------------------------------------------------------------- #
plot_ARI_dotplot <- function(obj,
                             annotation1,
                             annotation2,
                             title = NULL,
                             size_range = c(1, 10),
                             plot_width = 12,
                             plot_height = 10,
                             plot_dpi = 300) {
  # Extract and clean annotation data
  meta <- as.data.frame(colData(obj))
  if (!all(c(annotation1, annotation2) %in% colnames(meta))) stop("Annotation columns not found.")
  labels1 <- as.character(meta[[annotation1]])
  labels2 <- as.character(meta[[annotation2]])
  mask <- !is.na(labels1) & !is.na(labels2)
  if (sum(mask) == 0) stop("No overlapping cells between annotations.")

  # Compute ARI and cross-tabulation
  ari_score <- adjustedRandIndex(labels1[mask], labels2[mask])
  df <- as.data.frame(table(X = labels1[mask], Y = labels2[mask])) %>% filter(Freq > 0)

  # Build plot
  plot_title <- if (is.null(title)) {
    sprintf("%s vs %s (ARI = %.3f)", annotation1, annotation2, ari_score)
  } else {
    sprintf("%s (ARI = %.3f)", title, ari_score)
  }
  p <- ggplot(df, aes(x = X, y = Y, size = Freq)) +
    geom_point(alpha = 0.8, color = "steelblue") +
    scale_size_continuous(range = size_range) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.major = element_line(color = "gray80"), panel.grid.minor = element_blank(),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA)) +
    labs(title = plot_title, x = annotation1, y = annotation2, size = "Count")

  # Save outputs and generate JSON
  base_filename <- .generate_filename_base("ARI_dotplot", annot1 = annotation1, annot2 = annotation2)
  png_file <- paste0(base_filename, ".png")
  pdf_file <- paste0(base_filename, ".pdf")
  ggsave(png_file, plot = p, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
  ggsave(pdf_file, plot = p, width = plot_width, height = plot_height)

  output_paths <- list(plot_type = "ARI_dotplot",
                       pdf_path = tools::file_path_as_absolute(pdf_file),
                       png_path = tools::file_path_as_absolute(png_file))

  return(jsonlite::toJSON(output_paths, auto_unbox = TRUE, pretty = TRUE))
}

# ---------------------------------------------------------------------------- #
# 9. Plot Cluster Composition
# ---------------------------------------------------------------------------- #
plot_cluster_composition_groupedbar <- function(obj,
                                                annotation,
                                                palette = NULL,
                                                title = NULL,
                                                rotate_x = TRUE,
                                                plot_width = 12,
                                                plot_height = 8,
                                                plot_dpi = 300) {
  # Validate inputs
  meta <- as.data.frame(colData(obj))
  stopifnot(length(annotation) %in% c(1, 2), all(annotation %in% colnames(meta)))

  # Calculate proportions for one or two annotations
  if (length(annotation) == 1) {
    df <- meta %>% filter(!is.na(.data[[annotation[1]]])) %>%
      dplyr::count(annotation_value = .data[[annotation[1]]]) %>%
      mutate(Proportion = n / sum(n), annotation_source = annotation[1])
  } else {
    df1 <- meta %>% filter(!is.na(.data[[annotation[1]]])) %>%
      dplyr::count(annotation_value = .data[[annotation[1]]]) %>%
      mutate(Proportion = n / sum(n), annotation_source = annotation[1])
    df2 <- meta %>% filter(!is.na(.data[[annotation[2]]])) %>%
      dplyr::count(annotation_value = .data[[annotation[2]]]) %>%
      mutate(Proportion = n / sum(n), annotation_source = annotation[2])
    df <- bind_rows(df1, df2)
  }

  # Order categories by abundance
  annotation_order <- df %>% group_by(annotation_value) %>%
    summarize(total = sum(Proportion), .groups = "drop") %>%
    arrange(desc(total)) %>% pull(annotation_value)
  df <- df %>% mutate(annotation_value = factor(annotation_value, levels = annotation_order),
                      annotation_source = factor(annotation_source))

  # Build plot
  plot_title <- if (is.null(title)) "Cluster Composition Comparison" else title
  p <- ggplot(df, aes(x = annotation_value, y = Proportion, fill = annotation_source)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    labs(x = "Cell Type", y = "Proportion", fill = "Annotation Source", title = plot_title) +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = if (rotate_x) 45 else 0, hjust = if (rotate_x) 1 else 0.5),
          panel.grid.major.x = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))
  if (!is.null(palette)) p <- p + scale_fill_manual(values = palette)

  # Save outputs and generate JSON
  base_filename <- .generate_filename_base("composition", annotation = annotation)
  png_file <- paste0(base_filename, ".png")
  pdf_file <- paste0(base_filename, ".pdf")
  ggsave(png_file, plot = p, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
  ggsave(pdf_file, plot = p, width = plot_width, height = plot_height)

  output_paths <- list(plot_type = "composition_barplot",
                       pdf_path = tools::file_path_as_absolute(pdf_file),
                       png_path = tools::file_path_as_absolute(png_file))

  return(jsonlite::toJSON(output_paths, auto_unbox = TRUE, pretty = TRUE))
}

# ---------------------------------------------------------------------------- #
# 10. Visualize Matched Pairs Scatter Plot
# ---------------------------------------------------------------------------- #
visualize_matched_pairs_scatter <- function(obj,
                                            gene,
                                            layer = "X",
                                            log1p_transform = TRUE,
                                            point_size = 0.4,
                                            alpha = 0.6,
                                            show_correlation = FALSE,
                                            title = NULL,
                                            plot_width = 8,
                                            plot_height = 8,
                                            plot_dpi = 300) {
  # Extract expression data
  expr_mat <- assay(obj, layer)
  gene_visium <- paste0(gene, "_visium")
  gene_xenium <- paste0(gene, "_xenium")
  if (!all(c(gene_visium, gene_xenium) %in% rownames(expr_mat))) {
    stop(paste("Gene", gene, "not found in both Visium and Xenium formats."))
  }
  visium_expr <- expr_mat[gene_visium, ]
  xenium_expr <- expr_mat[gene_xenium, ]
  if (log1p_transform) {
    visium_expr <- log1p(visium_expr)
    xenium_expr <- log1p(xenium_expr)
  }
  df <- data.frame(Xenium = xenium_expr, Visium = visium_expr)

  # Build plot
  plot_title <- if (!is.null(title)) title else paste0("Xenium vs Visium: ", gene)
  if (show_correlation) {
    corr <- cor(df$Xenium, df$Visium, method = "pearson", use = "pairwise.complete.obs")
    plot_title <- paste0(plot_title, "\n", paste0("Pearson r = ", round(corr, 3)))
  }
  p <- ggplot(df, aes(x = Xenium, y = Visium)) +
    geom_point(size = point_size, alpha = alpha) +
    labs(x = if (log1p_transform) "Xenium Expression (log1p)" else "Xenium Expression",
         y = if (log1p_transform) "Visium Expression (log1p)" else "Visium Expression",
         title = plot_title) +
    theme_bw(base_size = 13) +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))

  # Save outputs and generate JSON
  base_filename <- .generate_filename_base("scatter_matched_pairs", gene = gene)
  png_file <- paste0(base_filename, ".png")
  pdf_file <- paste0(base_filename, ".pdf")
  ggsave(png_file, plot = p, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
  ggsave(pdf_file, plot = p, width = plot_width, height = plot_height)

  output_paths <- list(plot_type = "matched_pair_scatter",
                       pdf_path = tools::file_path_as_absolute(pdf_file),
                       png_path = tools::file_path_as_absolute(png_file))

  return(jsonlite::toJSON(output_paths, auto_unbox = TRUE, pretty = TRUE))
}

# ---------------------------------------------------------------------------- #
# 11. Generate Cluster Stability Report
# ---------------------------------------------------------------------------- #
generate_stability_report <- function(sce,
                                      celltype_labels,
                                      plot_width = 10,
                                      plot_height = 8,
                                      plot_dpi = 300) {
  # Extract stability scores from metadata
  metadata <- as.data.frame(colData(sce))
  all_cluster_scores <- list()
  for (annot in names(celltype_labels)) {
    for (label in celltype_labels[[annot]]) {
      cluster_string <- paste(annot, label, sep = "@")
      rep_cell <- metadata[metadata[[annot]] == label, ][1, ]
      if (nrow(rep_cell) == 0) { warning(paste("Could not find cluster:", cluster_string)); next }
      expected_cols <- c(paste0("reassign.", annot), paste0("SCCAF.", annot),
                         paste0("tfidf5.", annot), paste0("tfidf10.", annot),
                         paste0(annot, "_shapley"), "confidence")
      if (!all(expected_cols %in% names(rep_cell))) {
        warning(paste("Stability columns not found for", annot)); next
      }
      scores <- rep_cell[1, expected_cols]
      all_cluster_scores[[cluster_string]] <- tibble(
        Cluster = cluster_string,
        Metric = c("Reassign", "SCCAF", "TF-IDF 5", "TF-IDF 10", "Shapley", "Confidence"),
        Score = as.numeric(scores)
      )
    }
  }
  long_stability_data <- bind_rows(all_cluster_scores) %>% filter(!is.na(Score))
  if (nrow(long_stability_data) == 0) {
    warning("No stability data could be extracted."); return(NULL)
  }
  summary_table <- long_stability_data %>% pivot_wider(names_from = Metric, values_from = Score)

  # Build plot
  metric_order <- c("Confidence", "Reassign", "SCCAF", "TF-IDF 10", "TF-IDF 5", "Shapley")
  long_stability_data$Metric <- factor(long_stability_data$Metric, levels = metric_order)
  p <- ggplot(long_stability_data, aes(x = Metric, y = Score, fill = Cluster)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.9) +
    scale_fill_viridis_d(name = "Cluster Annotation", guide = guide_legend(nrow = 2)) +
    labs(y = "Score", x = NULL, title = "Cluster Stability Metrics") +
    theme_bw(base_size = 14) +
    theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))

  # Save outputs and generate JSON
  base_filename <- .generate_filename_base("stability_report", annotation = names(celltype_labels))
  tsv_file <- paste0(base_filename, "_table.tsv")
  png_file <- paste0(base_filename, "_plot.png")
  pdf_file <- paste0(base_filename, "_plot.pdf")
  write.table(summary_table, tsv_file, sep = "\t", row.names = FALSE, quote = FALSE)
  ggsave(png_file, plot = p, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
  ggsave(pdf_file, plot = p, width = plot_width, height = plot_height)

  output_paths <- list(report_type = "stability_report",
                       tsv_path = tools::file_path_as_absolute(tsv_file),
                       pdf_path = tools::file_path_as_absolute(pdf_file),
                       png_path = tools::file_path_as_absolute(png_file))

  return(jsonlite::toJSON(output_paths, auto_unbox = TRUE, pretty = TRUE))
}

# ---------------------------------------------------------------------------- #
# 12. Generate Comprehensive Marker Report
# ---------------------------------------------------------------------------- #
generate_marker_report <- function(obj,
                                   genes,
                                   annotation_cols,
                                   assay_name = "X",
                                   min_pct = 5,
                                   max_genes_plot = 50,
                                   plot_width = 14,
                                   plot_height = 10,
                                   plot_dpi = 300) {
  # Define rules and create analysis tasks
  rule_map <- tibble(annotation = c("Final_CT", "Xenium_RCTD_LungMap_ref", "Visium_RCTD_LungMap_ref",
                                    "Visium_RCTD_GSE250346_based", "pruned", "pruned"),
                     platform = c("Xenium", "Xenium", "Visium", "Visium", "Xenium", "Visium"))
  rules_to_run <- rule_map %>% filter(annotation %in% annotation_cols) %>%
    tidyr::crossing(gene = genes) %>% mutate(full_gene = paste0(gene, "_", tolower(platform)))
  base_gene_tasks <- tidyr::crossing(annotation = annotation_cols, platform = "Base", gene = genes) %>%
    mutate(full_gene = gene)
  tasks <- bind_rows(rules_to_run, base_gene_tasks)

  # Run analysis for each task
  all_results <- lapply(seq_len(nrow(tasks)), function(i) {
    task <- tasks[i, ]
    if (!task$full_gene %in% rownames(obj)) return(NULL)
    meta <- as.data.frame(colData(obj))
    df <- data.frame(expr = assay(obj, assay_name)[task$full_gene, ],
                     cluster = as.character(meta[[task$annotation]]))
    df$cluster[is.na(df$cluster)] <- "NA"
    cluster_stats <- lapply(unique(df$cluster), function(cl) {
      in_group_expr <- df$expr[df$cluster == cl]
      if (length(in_group_expr) < 2) return(NULL)
      t_result <- tryCatch(t.test(in_group_expr, df$expr[df$cluster != cl], alternative = "greater"), error = function(e) NULL)
      data.frame(Gene = task$gene, Platform = task$platform, Cluster = paste(task$annotation, cl, sep = "@"),
                 Mean_Expr = mean(in_group_expr, na.rm = TRUE), Pct_Expr = mean(in_group_expr > 0, na.rm = TRUE) * 100,
                 P_Value = if (!is.null(t_result)) t_result$p.value else NA_real_)
    })
    bind_rows(cluster_stats)
  })

  # Aggregate results and prepare for plotting
  final_table <- bind_rows(all_results)
  if (nrow(final_table) == 0) { warning("No marker results could be calculated."); return(NULL) }
  final_table <- final_table %>% mutate(P_Adj = p.adjust(P_Value, method = "BH")) %>% arrange(P_Adj, P_Value)
  dot_data <- final_table %>% filter(Pct_Expr >= min_pct) %>%
    mutate(gene_label = if_else(Platform != "Base", paste0(Gene, " (", Platform, ")"), Gene),
           Mean_Expr = ifelse(is.nan(Mean_Expr), 0, Mean_Expr))

  # Limit genes for plot readability
  if (!is.null(max_genes_plot) && n_distinct(dot_data$gene_label) > max_genes_plot) {
    top_genes <- dot_data %>% group_by(gene_label) %>%
      summarize(score = -log10(min(P_Adj, na.rm = TRUE) + 1e-300) * max(Pct_Expr, na.rm = TRUE)) %>%
      arrange(desc(score)) %>% slice_head(n = max_genes_plot) %>% pull(gene_label)
    dot_data <- dot_data %>% filter(gene_label %in% top_genes)
  }

  # Build plot
  p <- ggplot(dot_data, aes(x = gene_label, y = Cluster)) +
    geom_point(aes(size = Pct_Expr, color = Mean_Expr)) +
    scale_size(name = "% Expressing", range = c(1, 8)) +
    scale_color_gradient(name = "Avg. Expression", low = "grey90", high = "red") +
    theme_bw(base_size = 14) +
    theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))

  # Save outputs and generate JSON
  base_filename <- .generate_filename_base("marker_report", genes = genes, annotations = annotation_cols)
  tsv_file <- paste0(base_filename, "_table.tsv")
  png_file <- paste0(base_filename, "_dotplot.png")
  pdf_file <- paste0(base_filename, "_dotplot.pdf")
  write.table(final_table, tsv_file, sep = "\t", row.names = FALSE, quote = FALSE)
  ggsave(png_file, plot = p, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
  ggsave(pdf_file, plot = p, width = plot_width, height = plot_height)

  output_paths <- list(report_type = "marker_report",
                       tsv_path = tools::file_path_as_absolute(tsv_file),
                       pdf_path = tools::file_path_as_absolute(pdf_file),
                       png_path = tools::file_path_as_absolute(png_file))

  return(jsonlite::toJSON(output_paths, auto_unbox = TRUE, pretty = TRUE))
}

# ---------------------------------------------------------------------------- #
# 13. Find Cell Type Markers (Volcano Plot)
# ---------------------------------------------------------------------------- #
find_celltype_markers <- function(sce,
                                  celltype,
                                  annotation_col,
                                  platform = "Xenium",
                                  assay_name = "X",
                                  top_n = 10,
                                  logfc_threshold = 0.25,
                                  padj_threshold = 0.05,
                                  plot_width = 10,
                                  plot_height = 8,
                                  plot_dpi = 300) {
  # Filter genes by platform
  expr_mat <- assay(sce, assay_name)
  suffix <- paste0("_", tolower(platform), "$")
  gene_mask <- grepl(suffix, rownames(expr_mat))
  if (sum(gene_mask) == 0) { # Fallback to base gene names if no platform-specific ones
    gene_mask <- !duplicated(gsub(suffix, "", rownames(expr_mat)))
    if(sum(gene_mask) == 0) stop("No suitable genes found.")
  }
  expr_mat <- expr_mat[gene_mask, , drop = FALSE]

  # Define cell groups
  meta <- as.data.frame(colData(sce))
  target_mask <- meta[[annotation_col]] == celltype
  other_mask <- meta[[annotation_col]] != celltype
  if (sum(target_mask, na.rm = TRUE) < 3 || sum(other_mask, na.rm = TRUE) < 3) {
    stop("Not enough cells in target or other groups for comparison.")
  }

  # Compute logFC and p-values
  mean_target <- rowMeans(expr_mat[, target_mask], na.rm = TRUE)
  mean_other <- rowMeans(expr_mat[, other_mask], na.rm = TRUE)
  logFC <- log2(mean_target + 1e-9) - log2(mean_other + 1e-9)
  p_values <- apply(expr_mat, 1, function(gene_expr) {
    tryCatch(wilcox.test(gene_expr[target_mask], gene_expr[other_mask], alternative = "greater")$p.value,
             error = function(e) NA_real_)
  })

  # Create results table
  full_results_df <- tibble(gene = gsub(suffix, "", rownames(expr_mat)), logFC = logFC, P_Value = p_values) %>%
    filter(!is.na(logFC) & !is.na(P_Value)) %>%
    mutate(P_Adj = p.adjust(P_Value, method = "BH"),
           significant = P_Adj < padj_threshold & logFC > logfc_threshold,
           log10_P_Adj = -log10(P_Adj))
  if (any(is.infinite(full_results_df$log10_P_Adj))) { # Handle p-values of 0
    max_finite <- max(full_results_df$log10_P_Adj[is.finite(full_results_df$log10_P_Adj)], na.rm = TRUE)
    full_results_df$log10_P_Adj[is.infinite(full_results_df$log10_P_Adj)] <- max_finite * 1.1
  }
  top_markers_table <- full_results_df %>% filter(significant == TRUE) %>%
    arrange(desc(log10_P_Adj), desc(logFC)) %>% slice_head(n = top_n)

  # Build plot
  p <- ggplot(full_results_df, aes(x = logFC, y = log10_P_Adj)) +
    geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
    geom_text_repel(data = top_markers_table, aes(label = gene), size = 4, box.padding = 0.5, max.overlaps = Inf) +
    scale_color_manual(name = "", values = c("TRUE" = "red", "FALSE" = "grey50"),
                       labels = c("TRUE" = "Significant", "FALSE" = "Not Significant")) +
    geom_vline(xintercept = logfc_threshold, linetype = "dashed", color = "blue") +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "blue") +
    labs(title = paste("Volcano Plot for", celltype), x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
    theme_bw(base_size = 14) +
    theme(legend.position = "top", plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))

  # Save outputs and generate JSON
  base_filename <- .generate_filename_base("celltype_markers", celltype = celltype, platform = platform)
  tsv_file <- paste0(base_filename, "_table.tsv")
  png_file <- paste0(base_filename, "_plot.png")
  pdf_file <- paste0(base_filename, "_plot.pdf")
  write.table(top_markers_table, tsv_file, sep = "\t", row.names = FALSE, quote = FALSE)
  ggsave(png_file, plot = p, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
  ggsave(pdf_file, plot = p, width = plot_width, height = plot_height)

  output_paths <- list(report_type = "celltype_marker_volcano",
                       tsv_path = tools::file_path_as_absolute(tsv_file),
                       pdf_path = tools::file_path_as_absolute(pdf_file),
                       png_path = tools::file_path_as_absolute(png_file))

  return(jsonlite::toJSON(output_paths, auto_unbox = TRUE, pretty = TRUE))
}

# ---------------------------------------------------------------------------- #
# 14. Squidpy: Create Neighborhood Enrichment Heatmap
# ---------------------------------------------------------------------------- #
create_enrichment_heatmap <- function(obj,
                                      group = "Unaffected",
                                      annotation = "Final_CT",
                                      clustering = TRUE,
                                      title = NULL,
                                      base_size = 10,
                                      plot_width = 10,
                                      plot_height = 8,
                                      plot_dpi = 300) {
  # Extract and clean enrichment matrix
  enrichment_df <- obj$Squidpy[[annotation]][[group]]
  if (is.null(enrichment_df)) stop("Enrichment matrix not found for specified group/annotation.")
  enrichment_mat <- as.matrix(enrichment_df)
  if (!is.numeric(enrichment_mat)) stop("Enrichment matrix must be numeric.")
  enrichment_mat[!is.finite(enrichment_mat)] <- 0

  # Define plot aesthetics
  plot_title <- if (is.null(title)) paste0("Enrichment Heatmap: ", group, " (", annotation, ")") else title
  min_val <- min(enrichment_mat, na.rm = TRUE)
  max_val <- max(enrichment_mat, na.rm = TRUE)
  if (!is.finite(min_val) || !is.finite(max_val) || min_val == max_val) {
    min_val <- 0; max_val <- 1
  }
  col_fun <- circlize::colorRamp2(seq(min_val, max_val, length.out = 100), viridis::viridis(100))

  # Create Heatmap object
  hm <- Heatmap(enrichment_mat, name = "Z-score", col = col_fun,
                cluster_rows = clustering, cluster_columns = clustering,
                row_names_gp = grid::gpar(fontsize = base_size),
                column_names_gp = grid::gpar(fontsize = base_size),
                heatmap_legend_param = list(title = "Z-score",
                                            title_gp = grid::gpar(fontsize = base_size + 1),
                                            labels_gp = grid::gpar(fontsize = base_size)),
                column_title = plot_title,
                column_title_gp = grid::gpar(fontsize = base_size + 2, fontface = "bold"))

  # Save outputs and generate JSON
  base_filename <- .generate_filename_base("enrichment_heatmap", group = group, annotation = annotation)
  png_file <- paste0(base_filename, ".png")
  pdf_file <- paste0(base_filename, ".pdf")
  png(png_file, width = plot_width, height = plot_height, units = "in", res = plot_dpi)
  draw(hm)
  dev.off()
  pdf(pdf_file, width = plot_width, height = plot_height)
  draw(hm)
  dev.off()

  output_paths <- list(plot_type = "enrichment_heatmap",
                       pdf_path = tools::file_path_as_absolute(pdf_file),
                       png_path = tools::file_path_as_absolute(png_file))

  return(jsonlite::toJSON(output_paths, auto_unbox = TRUE, pretty = TRUE))
}

# ---------------------------------------------------------------------------- #
# 15. Squidpy: Analyze and Plot Cell Type Neighbors
# ---------------------------------------------------------------------------- #
analyze_and_plot_neighbors <- function(obj,
                                       target_cell_type,
                                       group = "More_Affected",
                                       annotation = "Final_CT",
                                       base_size = 11,
                                       plot_width = 10,
                                       plot_height = 10,
                                       plot_dpi = 300) {
  # Extract and validate data
  enrichment_df <- obj$Squidpy[[annotation]][[group]]
  if (is.null(enrichment_df)) stop("Enrichment matrix not found for specified group/annotation.")
  if (!(target_cell_type %in% rownames(enrichment_df))) {
    stop(paste0("Target cell type '", target_cell_type, "' not found in matrix."))
  }

  # Rank and transform scores
  ranked_values <- sort(enrichment_df[target_cell_type, ], decreasing = TRUE)
  log_transformed <- sign(ranked_values) * log10(abs(ranked_values) + 1)
  log_df <- data.frame(cell_type = factor(names(log_transformed), levels = rev(names(log_transformed))),
                       log_score = as.numeric(log_transformed))

  # Build plot
  p <- ggplot(log_df, aes(x = log_score, y = cell_type, fill = log_score >= 0)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("red", "royalblue"), guide = "none") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    labs(title = paste0("Neighborhood Enrichment: ", target_cell_type, " in ", gsub("_", " ", group)),
         x = "Signed log10(z-score + 1)", y = "Neighboring Cell Type") +
    theme_bw(base_size = base_size) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          axis.text = element_text(color = "black"),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))

  # Save outputs and generate JSON
  base_filename <- .generate_filename_base("neighbors", target = target_cell_type, group = group, annotation = annotation)
  png_file <- paste0(base_filename, ".png")
  pdf_file <- paste0(base_filename, ".pdf")
  ggsave(png_file, plot = p, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
  ggsave(pdf_file, plot = p, width = plot_width, height = plot_height)

  output_paths <- list(plot_type = "neighbor_enrichment_plot",
                       pdf_path = tools::file_path_as_absolute(pdf_file),
                       png_path = tools::file_path_as_absolute(png_file))

  return(jsonlite::toJSON(output_paths, auto_unbox = TRUE, pretty = TRUE))
}

# ---------------------------------------------------------------------------- #
# 16. Generate DE Volcano Plot
# ---------------------------------------------------------------------------- #
generate_de_volcano_plot <- function(obj_list,
                                     annotation = "TNiche",
                                     group_1 = "More_Affected",
                                     group_2 = "Unaffected",
                                     celltype_filter = NULL,
                                     top_n = 100,
                                     label_top_n = 20,
                                     logfc_threshold = 0.25,
                                     padj_threshold = 0.05,
                                     plot_width = 10,
                                     plot_height = 8,
                                     plot_dpi = 300) {
  # Extract data from CellChat objects
  cc_obj1 <- obj_list[[annotation]][[group_1]]
  cc_obj2 <- obj_list[[annotation]][[group_2]]
  mat1 <- cc_obj1@data
  mat2 <- cc_obj2@data
  if (!is.null(celltype_filter)) {
    mat1 <- mat1[, cc_obj1@idents == celltype_filter, drop = FALSE]
    mat2 <- mat2[, cc_obj2@idents == celltype_filter, drop = FALSE]
  }
  combined_mat <- cbind(mat1, mat2)
  meta_df <- data.frame(group = c(rep(group_1, ncol(mat1)), rep(group_2, ncol(mat2))))
  rownames(meta_df) <- colnames(combined_mat)

  # Create Seurat object and find markers, suppressing warnings
  suppressWarnings({
    seu <- CreateSeuratObject(counts = combined_mat, meta.data = meta_df)
    seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    deg_result <- FindMarkers(object = seu, ident.1 = group_1, ident.2 = group_2, group.by = "group",
                              logfc.threshold = 0, min.pct = 0.0, test.use = "wilcox", verbose = FALSE)
  })

  # Process results
  deg_result <- deg_result %>% as.data.frame() %>% rownames_to_column("Gene") %>%
    mutate(log10_padj = -log10(p_val_adj + 1e-300),
           log10_padj = ifelse(log10_padj > 300, 300, log10_padj),
           Sig = case_when(
             p_val_adj < padj_threshold & avg_log2FC > logfc_threshold ~ paste("Up in", group_1),
             p_val_adj < padj_threshold & avg_log2FC < -logfc_threshold ~ paste("Up in", group_2),
             TRUE ~ "Not significant"
           ))
  top_table <- deg_result %>% arrange(desc(log10_padj), desc(abs(avg_log2FC))) %>% head(top_n)
  label_genes <- head(top_table$Gene, label_top_n)

  # Build plot
  plot_title <- paste(group_1, "vs", group_2, if (!is.null(celltype_filter)) paste("-", celltype_filter) else "")
  p <- ggplot(deg_result, aes(x = avg_log2FC, y = log10_padj, color = Sig)) +
    geom_point(alpha = 0.7, size = 1) +
    scale_color_manual(values = c("red", "blue", "gray"), name = "Significance") +
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black") +
    geom_text_repel(data = subset(deg_result, Gene %in% label_genes), aes(label = Gene),
                    size = 3, max.overlaps = Inf, box.padding = 0.5, segment.size = 0.2) +
    labs(title = plot_title, x = "Average Log2 Fold Change", y = "-log10(Adjusted P-value)") +
    theme_bw(base_size = 12) +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))

  # Save outputs and generate JSON
  base_filename <- .generate_filename_base("de_volcano", g1 = group_1, g2 = group_2, filter = celltype_filter)
  tsv_file <- paste0(base_filename, "_table.tsv")
  png_file <- paste0(base_filename, "_plot.png")
  pdf_file <- paste0(base_filename, "_plot.pdf")
  write.table(top_table, tsv_file, sep = "\t", row.names = FALSE, quote = FALSE)
  ggsave(png_file, plot = p, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
  ggsave(pdf_file, plot = p, width = plot_width, height = plot_height)

  output_paths <- list(report_type = "de_volcano",
                       tsv_path = tools::file_path_as_absolute(tsv_file),
                       pdf_path = tools::file_path_as_absolute(pdf_file),
                       png_path = tools::file_path_as_absolute(png_file))

  return(jsonlite::toJSON(output_paths, auto_unbox = TRUE, pretty = TRUE))
}

# ---------------------------------------------------------------------------- #
# 17. CellChat: Circle Plot of Communication Network
# ---------------------------------------------------------------------------- #
netVisual_circle <- function(net, ..., plot_width = 12, plot_height = 12, plot_dpi = 300) {
  if (missing(net)) stop("The 'net' argument is required.")
  plot_args <- list(...)
  plot_fun <- function() do.call(CellChat::netVisual_circle, c(list(net = net), plot_args))

  base_filename <- .generate_filename_base("netVisual_circle", title = plot_args$title.name)
  saved_paths <- .save_base_graphics_plot(plot_fun, base_filename, plot_width, plot_height, plot_dpi)

  output_paths <- list(plot_type = "cellchat_circle_plot",
                       pdf_path = saved_paths$pdf_path,
                       png_path = saved_paths$png_path)

  return(jsonlite::toJSON(output_paths, auto_unbox = TRUE, pretty = TRUE))
}

# ---------------------------------------------------------------------------- #
# 18. CellChat: Heatmap of Communication Network
# ---------------------------------------------------------------------------- #
netVisual_heatmap <- function(object, ..., plot_width = 10, plot_height = 12, plot_dpi = 300) {
  # Validate input
  if (!inherits(object, "CellChat")) {
    stop("The 'object' must be a valid CellChat object.")
  }

  # Create the heatmap object, suppressing console messages
  hm <- suppressMessages({
    do.call(CellChat::netVisual_heatmap, c(list(object = object), list(...)))
  })

  # Save outputs and generate JSON
  plot_args <- list(...)
  base_filename <- .generate_filename_base("netVisual_heatmap", title = plot_args$title.name, measure = plot_args$measure[1], signaling = plot_args$signaling)
  png_file <- paste0(base_filename, ".png")
  pdf_file <- paste0(base_filename, ".pdf")
  png(png_file, width = plot_width, height = plot_height, units = "in", res = plot_dpi)
  draw(hm)
  dev.off()
  pdf(pdf_file, width = plot_width, height = plot_height)
  draw(hm)
  dev.off()

  output_paths <- list(plot_type = "cellchat_heatmap",
                       pdf_path = tools::file_path_as_absolute(pdf_file),
                       png_path = tools::file_path_as_absolute(png_file))

  return(jsonlite::toJSON(output_paths, auto_unbox = TRUE, pretty = TRUE))
}

# ---------------------------------------------------------------------------- #
# 19. CellChat: Aggregate Signaling Pathway Plot
# ---------------------------------------------------------------------------- #
netVisual_aggregate <- function(object, signaling, ..., plot_width = 10, plot_height = 10, plot_dpi = 300) {
  if (!inherits(object, "CellChat")) stop("The 'object' must be a valid CellChat object.")
  if (missing(signaling)) stop("The 'signaling' pathway must be specified.")
  plot_args <- list(...)
  plot_fun <- function() do.call(CellChat::netVisual_aggregate, c(list(object = object, signaling = signaling), plot_args))

  base_filename <- .generate_filename_base("netVisual_aggregate", signaling = signaling, layout = plot_args$layout[1])
  saved_paths <- .save_base_graphics_plot(plot_fun, base_filename, plot_width, plot_height, plot_dpi)

  output_paths <- list(plot_type = "cellchat_aggregate_plot",
                       pdf_path = saved_paths$pdf_path,
                       png_path = saved_paths$png_path)

  return(jsonlite::toJSON(output_paths, auto_unbox = TRUE, pretty = TRUE))
}

# ---------------------------------------------------------------------------- #
# 20. CellChat: Subset Communication Data
# ---------------------------------------------------------------------------- #
subsetCommunication <- function(object = NULL, ...) {
  # Validate input
  if (is.null(object) && is.null(list(...)$net)) {
    stop("You must provide either a 'CellChat' object or a 'net' data frame.")
  }

  # Subset the data and save it
  df_result <- do.call(CellChat::subsetCommunication, c(list(object = object), list(...)))
  base_filename <- .generate_filename_base("subset_communication", signaling = list(...)$signaling)
  tsv_file <- paste0(base_filename, ".tsv")
  write.table(df_result, tsv_file, sep = "\t", row.names = FALSE, quote = FALSE)

  output_paths <- list(data_type = "communication_subset",
                       tsv_path = tools::file_path_as_absolute(tsv_file))

  return(jsonlite::toJSON(output_paths, auto_unbox = TRUE, pretty = TRUE))
}

# ---------------------------------------------------------------------------- #
# 21. ToppCell Enrichment Analysis
# ---------------------------------------------------------------------------- #
get_toppcell_enrichment <- function(gene_symbols) {
  # Helper to return a standardized error JSON
  .return_error <- function(msg) {
    return(jsonlite::toJSON(list(data_type = "toppcell_enrichment_table", tsv_path = NA, error = msg), auto_unbox = TRUE, pretty = TRUE))
  }
  
  # Convert gene symbols to Entrez IDs via ToppGene API
  lookup_url <- "https://toppgene.cchmc.org/API/lookup"
  lookup_body <- list(Symbols = gene_symbols)
  lookup_response <- POST(url = lookup_url, body = toJSON(lookup_body, auto_unbox = TRUE), add_headers("Content-Type" = "text/json"))
  if (status_code(lookup_response) != 200) return(.return_error("ToppGene API lookup failed"))
  
  lookup_results <- fromJSON(content(lookup_response, "text", encoding = "UTF-8"))
  entrez_ids <- as.numeric(lookup_results$Genes$Entrez[!is.na(lookup_results$Genes$Entrez)])
  if (length(entrez_ids) < 5) return(.return_error("Fewer than 5 valid Entrez IDs found for enrichment"))

  # Submit Entrez IDs to ToppCell API for enrichment
  enrich_url <- "https://toppgene.cchmc.org/API/enrich"
  enrich_body <- list(Genes = entrez_ids, Categories = list(list(Type = "ToppCell", PValue = 0.05, MaxResults = 5, Correction = "FDR")))
  enrich_response <- POST(url = enrich_url, body = toJSON(enrich_body, auto_unbox = TRUE), add_headers("Content-Type" = "text/json"))
  if (status_code(enrich_response) != 200) return(.return_error("ToppGene API enrich failed"))
  
  enrich_results <- fromJSON(content(enrich_response, "text", encoding = "UTF-8"))

  # Process and save results
  results_df <- NULL
  if (!is.null(enrich_results$Annotations) && length(enrich_results$Annotations) > 0) {
    results_df <- as.data.frame(enrich_results$Annotations)
  }
  if (is.null(results_df) || nrow(results_df) == 0) return(.return_error("No significant enrichment results found"))

  base_filename <- .generate_filename_base("toppcell_enrichment", genes = gene_symbols)
  tsv_file <- paste0(base_filename, ".tsv")
  write.table(results_df, tsv_file, sep = "\t", row.names = FALSE, quote = FALSE)

  output_paths <- list(data_type = "toppcell_enrichment_table",
                       tsv_path = tools::file_path_as_absolute(tsv_file))

  return(jsonlite::toJSON(output_paths, auto_unbox = TRUE, pretty = TRUE))
}