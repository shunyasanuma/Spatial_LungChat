import os
import json
from typing import Literal, Optional, List, Union, Dict
from pydantic import BaseModel, Field

# -----------------------------
# Base Configuration Classes #Pankaj
# -----------------------------

class BaseSpatialPlotConfig(BaseModel):
    """
    Base configuration for spatial plots.
    All spatial plot configs inherit from this.
    """
    spatial_data_file: str = Field(..., description="Path to the spatial .h5ad file.")
    plot_type: str = Field(..., description="Type of spatial plot to generate.")

# -----------------------------
# Individual Spatial Plot Configurations
# -----------------------------

class PlotSpatialConfig(BaseSpatialPlotConfig):
    """
    Configuration for plot_spatial function.
    Shows spatial distribution of annotations on tissue coordinates.
    R signature: plot_spatial(obj, annotation, title=NULL, point_size=0.5, base_size=12, coord_flip=FALSE, highlight=NULL)
    """
    plot_type: Literal["plot_spatial"] = "plot_spatial"
    annotation: str = Field(..., description="Column name in colData to visualize (e.g., 'Final_CT', 'CNiche', 'TNiche').")
    point_size: Optional[float] = Field(0.5, description="Size of points in the plot.")
    title: Optional[str] = Field(None, description="Custom title for the plot.")
    base_size: Optional[int] = Field(12, description="Base font size for the plot.")
    coord_flip: Optional[bool] = Field(False, description="Whether to flip coordinates.")
    highlight: Optional[List[str]] = Field(None, description="Cell types to highlight (others will be gray).")

class PlotUmapConfig(BaseSpatialPlotConfig):
    """
    Configuration for plot_umap function.
    Shows UMAP embedding colored by annotations.
    R signature: plot_umap(obj, platform=c("Visium", "Xenium"), annotation="Final_CT", title=NULL, point_size=0.5, legend_point_size=4, highlight=NULL)
    """
    plot_type: Literal["plot_umap"] = "plot_umap"
    platform: str = Field("Visium", description="Platform to use ('Visium' or 'Xenium').")
    annotation: str = Field("Final_CT", description="Column name in colData to color by (e.g., 'Final_CT', 'platform').")
    title: Optional[str] = Field(None, description="Custom title for the plot.")
    point_size: Optional[float] = Field(0.5, description="Size of points in the plot.")
    legend_point_size: Optional[int] = Field(4, description="Size of points in legend.")
    highlight: Optional[List[str]] = Field(None, description="Cell types to highlight (others will be gray).")

class PlotGeneSpatialConfig(BaseSpatialPlotConfig):
    """
    Configuration for plot_gene_spatial function.
    Shows spatial expression patterns of specific genes.
    R signature: plot_gene_spatial(obj, gene, platform=c("Visium", "Xenium"), x_col=NULL, y_col=NULL, title=NULL, point_size=0.5, low_color="grey90", high_color="red")
    """
    plot_type: Literal["plot_gene_spatial"] = "plot_gene_spatial"
    gene: str = Field(..., description="Gene name to visualize (base name without platform suffix).")
    platform: str = Field("Visium", description="Platform to use ('Visium' or 'Xenium').")
    x_col: Optional[str] = Field(None, description="X coordinate column name (auto-detected if None).")
    y_col: Optional[str] = Field(None, description="Y coordinate column name (auto-detected if None).")
    title: Optional[str] = Field(None, description="Custom title for the plot.")
    point_size: Optional[float] = Field(0.5, description="Size of points in the plot.")
    low_color: Optional[str] = Field("grey90", description="Color for low expression.")
    high_color: Optional[str] = Field("red", description="Color for high expression.")

class RankCrossPlatformGeneCorrelationConfig(BaseSpatialPlotConfig):
    """
    Configuration for rank_cross_platform_gene_correlation function.
    Identifies genes with highest correlation between platforms.
    R signature: rank_cross_platform_gene_correlation(obj, top_genes=5, layer="X")
    """
    plot_type: Literal["rank_cross_platform_gene_correlation"] = "rank_cross_platform_gene_correlation"
    top_genes: Optional[int] = Field(5, description="Number of top correlated genes to return.")
    layer: Optional[str] = Field("X", description="Assay layer to use.")

class PlotCrossPlatformGeneCorrelationConfig(BaseSpatialPlotConfig):
    """
    Configuration for plot_cross_platform_gene_correlation function.
    Visualizes correlation curves for specific genes across platforms.
    R signature: plot_cross_platform_gene_correlation(obj, gene_name, ID1_col="Visium_ID", ID2_col="Xenium_ID", hline=0.9, y_range=c(0,1), title=NULL)
    """
    plot_type: Literal["plot_cross_platform_gene_correlation"] = "plot_cross_platform_gene_correlation"
    gene_name: str = Field(..., description="Gene name to visualize correlation for.")
    ID1_col: Optional[str] = Field("Visium_ID", description="Column name for first platform IDs.")
    ID2_col: Optional[str] = Field("Xenium_ID", description="Column name for second platform IDs.")
    hline: Optional[float] = Field(0.9, description="Horizontal line to draw on plot.")
    y_range: Optional[List[float]] = Field([0, 1], description="Y-axis range as [min, max].")
    title: Optional[str] = Field(None, description="Custom title for the plot.")

class RankMarkerSpecificityConfig(BaseSpatialPlotConfig):
    """
    Configuration for rank_marker_specificity function.
    Ranks genes by their specificity to cell types.
    """
    plot_type: Literal["rank_marker_specificity"] = "rank_marker_specificity"
    platform: str = Field(..., description="Platform to analyze ('Visium' or 'Xenium').")
    annotation: str = Field(..., description="Column name for cell type annotation (e.g., 'Final_CT').")
    n_genes: Optional[int] = Field(10, description="Number of top marker genes per cell type.")

class FindTopMarkersForCelltypeConfig(BaseSpatialPlotConfig):
    """
    Configuration for find_top_markers_for_celltype function.
    Finds top marker genes for specific cell types.
    R signature: find_top_markers_for_celltype(obj, celltype, annotation_col, platform_suffixes, top_n=10)
    """
    plot_type: Literal["find_top_markers_for_celltype"] = "find_top_markers_for_celltype"
    celltype: str = Field(..., description="Cell type to find markers for (e.g., 'AT1', 'AT2', 'C1').")
    annotation_col: str = Field(..., description="Column name for cell type annotation (e.g., 'Final_CT').")
    platform_suffixes: List[str] = Field(..., description="List of platform suffixes (e.g., ['visium'], ['xenium'], or ['visium', 'xenium']).")
    top_n: Optional[int] = Field(10, description="Number of top marker genes to return.")

class ScTriangulateSummarizeCelltypeStabilityConfig(BaseSpatialPlotConfig):
    """
    Configuration for scTriangulate_summarize_celltype_stability function.
    Analyzes stability of cell type annotations.
    """
    plot_type: Literal["scTriangulate_summarize_celltype_stability"] = "scTriangulate_summarize_celltype_stability"
    labels_list: Dict[str, List[str]] = Field(..., description="Dictionary mapping annotation names to lists of cell types.")

class PlotARIDotplotConfig(BaseSpatialPlotConfig):
    """
    Configuration for plot_ARI_dotplot function.
    Shows adjusted rand index between annotations.
    """
    plot_type: Literal["plot_ARI_dotplot"] = "plot_ARI_dotplot"
    annotation1: str = Field(..., description="First annotation to compare (e.g., 'pruned').")
    annotation2: str = Field(..., description="Second annotation to compare (e.g., 'Final_CT').")

class PlotMarkerExpressionDotplotConfig(BaseSpatialPlotConfig):
    """
    Configuration for plot_marker_expression_dotplot function.
    Dotplot showing marker gene expression across cell types.
    """
    plot_type: Literal["plot_marker_expression_dotplot"] = "plot_marker_expression_dotplot"
    genes: List[str] = Field(..., description="List of genes to visualize.")
    platform: str = Field(..., description="Platform to analyze ('Visium' or 'Xenium').")
    annotation: str = Field(..., description="Column name for cell type annotation (e.g., 'Final_CT').")

class PlotClusterCompositionGroupedbarConfig(BaseSpatialPlotConfig):
    """
    Configuration for plot_cluster_composition_groupedbar function.
    Bar plot showing composition of clusters.
    """
    plot_type: Literal["plot_cluster_composition_groupedbar"] = "plot_cluster_composition_groupedbar"
    annotation: str = Field(..., description="Column name for annotation to visualize (e.g., 'Final_CT').")

class CompareAnnotationsViaHeatmapConfig(BaseSpatialPlotConfig):
    """
    Configuration for compare_annotations_via_heatmap function.
    Heatmap comparing different annotation methods.
    """
    plot_type: Literal["compare_annotations_via_heatmap"] = "compare_annotations_via_heatmap"
    annotation1: str = Field(..., description="First annotation to compare (e.g., 'pruned').")
    annotation2: str = Field(..., description="Second annotation to compare (e.g., 'Final_CT').")
    normalize: Optional[str] = Field("row", description="Normalization method ('row', 'column', or 'none').")
    title: Optional[str] = Field(None, description="Custom title for the heatmap.")

class VisualizeMatchedPairsScatterConfig(BaseSpatialPlotConfig):
    """
    Configuration for visualize_matched_pairs_scatter function.
    Scatter plot for matched pairs analysis.
    """
    plot_type: Literal["visualize_matched_pairs_scatter"] = "visualize_matched_pairs_scatter"
    gene: str = Field(..., description="Gene name to visualize.")
    annotation: str = Field(..., description="Column name for cell type annotation.")

# -----------------------------
# Plot Class Mapping Function
# -----------------------------

def get_spatial_plot_class(plot_type: str):
    """Map the plot_type string to the correct spatial Pydantic class."""
    plot_class_map = {
        "plot_spatial": PlotSpatialConfig,
        "plot_umap": PlotUmapConfig,
        "plot_gene_spatial": PlotGeneSpatialConfig,
        "rank_cross_platform_gene_correlation": RankCrossPlatformGeneCorrelationConfig,
        "plot_cross_platform_gene_correlation": PlotCrossPlatformGeneCorrelationConfig,
        "rank_marker_specificity": RankMarkerSpecificityConfig,
        "find_top_markers_for_celltype": FindTopMarkersForCelltypeConfig,
        "scTriangulate_summarize_celltype_stability": ScTriangulateSummarizeCelltypeStabilityConfig,
        "plot_ARI_dotplot": PlotARIDotplotConfig,
        "plot_marker_expression_dotplot": PlotMarkerExpressionDotplotConfig,
        "plot_cluster_composition_groupedbar": PlotClusterCompositionGroupedbarConfig,
        "compare_annotations_via_heatmap": CompareAnnotationsViaHeatmapConfig,
        "visualize_matched_pairs_scatter": VisualizeMatchedPairsScatterConfig,
    }
    
    if plot_type not in plot_class_map:
        raise ValueError(f"Unsupported spatial plot type: {plot_type}")
    
    return plot_class_map[plot_type]

# -----------------------------
# Spatial Plot Guides
# -----------------------------

SPATIAL_PLOT_GUIDES = {
    "plot_spatial": (
        "Spatial Plot:\n"
        "Shows spatial distribution of cell types or annotations on tissue coordinates.\n"
        "Required fields:\n"
        "  - annotation: Column name in colData to visualize (e.g., 'Final_CT', 'CNiche', 'TNiche')\n"
        "Optional fields:\n"
        "  - platform: Filter by platform ('Visium' or 'Xenium')\n"
        "  - point_size: Size of points (default 0.5)\n"
        "  - alpha: Transparency (default 0.7)\n"
        "  - title: Custom plot title\n"
    ),
    
    "plot_umap": (
        "UMAP Plot:\n"
        "Shows cells in UMAP embedding colored by annotations.\n"
        "Required fields:\n"
        "  - annotation: Column name to color by (e.g., 'Final_CT', 'platform')\n"
        "Optional fields:\n"
        "  - platform: Filter by platform ('Visium' or 'Xenium')\n"
        "  - facet: Create separate facets for platforms (default False)\n"
        "  - point_size: Size of points (default 0.5)\n"
        "  - alpha: Transparency (default 0.7)\n"
    ),
    
    "plot_gene_spatial": (
        "Gene Spatial Plot:\n"
        "Shows spatial expression patterns of specific genes.\n"
        "Required fields:\n"
        "  - gene: Gene name to visualize (must exist in data)\n"
        "Optional fields:\n"
        "  - platform: Filter by platform ('Visium' or 'Xenium')\n"
        "  - point_size: Size of points (default 0.5)\n"
        "  - alpha: Transparency (default 0.7)\n"
        "  - title: Custom plot title\n"
    ),
    
    "rank_cross_platform_gene_correlation": (
        "Cross-Platform Gene Correlation Ranking:\n"
        "Identifies genes with highest correlation between Visium and Xenium platforms.\n"
        "Optional fields:\n"
        "  - n_genes: Number of top correlated genes to return (default 20)\n"
    ),
    
    "plot_cross_platform_gene_correlation": (
        "Cross-Platform Gene Correlation Plot:\n"
        "Visualizes correlation curves for specific genes across platforms.\n"
        "Required fields:\n"
        "  - genes: List of gene names to visualize correlations for\n"
    ),
    
    "rank_marker_specificity": (
        "Marker Specificity Ranking:\n"
        "Ranks genes by their specificity to cell types in spatial context.\n"
        "Required fields:\n"
        "  - platform: Platform to analyze ('Visium' or 'Xenium')\n"
        "  - annotation: Cell type annotation column (e.g., 'Final_CT')\n"
        "Optional fields:\n"
        "  - n_genes: Number of top markers per cell type (default 10)\n"
    ),
    
    "find_top_markers_for_celltype": (
        "Find Top Markers for Cell Type:\n"
        "Finds top marker genes for specific cell types.\n"
        "Required fields:\n"
        "  - celltype: Cell type to find markers for (e.g., 'AT1', 'AT2')\n"
        "  - platform: Platform to analyze ('Visium' or 'Xenium')\n"
        "  - annotation: Cell type annotation column (e.g., 'Final_CT')\n"
        "Optional fields:\n"
        "  - n_genes: Number of top markers to return (default 10)\n"
    ),
    
    "scTriangulate_summarize_celltype_stability": (
        "Cell Type Stability Summary:\n"
        "Analyzes stability of cell type annotations across different methods.\n"
        "Required fields:\n"
        "  - labels_list: Dictionary mapping annotation names to cell type lists\n"
    ),
    
    "plot_ARI_dotplot": (
        "ARI Dotplot:\n"
        "Shows adjusted rand index between different annotations.\n"
        "Required fields:\n"
        "  - annotation1: First annotation to compare (e.g., 'pruned')\n"
        "  - annotation2: Second annotation to compare (e.g., 'Final_CT')\n"
    ),
    
    "plot_marker_expression_dotplot": (
        "Marker Expression Dotplot:\n"
        "Dotplot showing marker gene expression across cell types.\n"
        "Required fields:\n"
        "  - genes: List of genes to visualize\n"
        "  - platform: Platform to analyze ('Visium' or 'Xenium')\n"
        "  - annotation: Cell type annotation column (e.g., 'Final_CT')\n"
    ),
    
    "plot_cluster_composition_groupedbar": (
        "Cluster Composition Bar Plot:\n"
        "Bar plot showing composition of clusters or annotations.\n"
        "Required fields:\n"
        "  - annotation: Annotation column to visualize (e.g., 'Final_CT')\n"
    ),
    
    "compare_annotations_via_heatmap": (
        "Annotation Comparison Heatmap:\n"
        "Heatmap comparing different annotation methods.\n"
        "Required fields:\n"
        "  - annotation1: First annotation to compare (e.g., 'pruned')\n"
        "  - annotation2: Second annotation to compare (e.g., 'Final_CT')\n"
        "Optional fields:\n"
        "  - normalize: Normalization method ('row', 'column', 'none', default 'row')\n"
        "  - title: Custom heatmap title\n"
    ),
    
    "visualize_matched_pairs_scatter": (
        "Matched Pairs Scatter Plot:\n"
        "Scatter plot for matched pairs analysis.\n"
        "Required fields:\n"
        "  - gene: Gene name to visualize\n"
        "  - annotation: Cell type annotation column\n"
    ),
}
