#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def create_enrichment_heatmap(group: str, cluster_key: str, hdf_path: str, plot_style: str = 'blue_white_red'):
    """
    Loads RAW neighborhood enrichment data, log-transforms it, and generates a heatmap.
    This final version handles potential NaN values in the data.
    """
    # 1. Load the RAW z-score data
    key_to_load = f"/{group}_{cluster_key}_raw_zscore"
    try:
        with pd.HDFStore(hdf_path, mode='r') as store:
            raw_zscore_df = store[key_to_load]
    except KeyError:
        print(f"Error: Could not find data for key '{key_to_load}'.")
        return

    # 2. Perform the signed log10 transformation
    data_df = np.sign(raw_zscore_df) * np.log10(np.abs(raw_zscore_df) + 1)
    
    # 3. Replace any NaN values with 0
    data_df = data_df.fillna(0)
    
    # 4. Create the plot based on the chosen style
    fig, ax = plt.subplots(figsize=(12, 10))
    title = f"Neighborhood Enrichment for {group.replace('_', ' ')} ({cluster_key})"

    if plot_style == 'blue_white_red':
        plot_df = data_df.copy()
        cmap = LinearSegmentedColormap.from_list("bwr", ["blue", "white", "red"])
        max_val = np.abs(plot_df.values).max()
        im = ax.imshow(plot_df.values, cmap=cmap, vmin=-max_val, vmax=max_val, aspect="auto")
        cbar_label = "Signed log10(z-score + 1)"
        
    elif plot_style == 'grey_red':
        plot_df = data_df.copy()
        plot_df[plot_df < 0] = 0
        cmap = LinearSegmentedColormap.from_list("grey_red", ["lightgrey", "red"])
        im = ax.imshow(plot_df.values, cmap=cmap, aspect="auto")
        cbar_label = "log10(z-score + 1) for Enrichment"
        
    else:
        print(f"Error: Invalid plot_style '{plot_style}'. Choose 'grey_red' or 'blue_white_red'.")
        return

    # 5. Set final plot details
    ax.set_xticks(np.arange(len(data_df.columns)))
    ax.set_yticks(np.arange(len(data_df.index)))
    ax.set_xticklabels(data_df.columns, rotation=90, fontsize=12)
    ax.set_yticklabels(data_df.index, fontsize=12)
    ax.set_title(title, fontsize=16)
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(cbar_label, fontsize=14)
    plt.tight_layout()
    plt.show()


# In[2]:


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def analyze_and_plot_neighbors(target_cell_type: str, group: str, cluster_key: str, hdf_path: str):
    """
    Ranks neighboring cell types and visualizes the result as a bar plot.
    This version prints the log-transformed scores.

    Args:
        target_cell_type (str): The cell type to analyze (e.g., 'AT1', 'T6').
        group (str): The sample group to use (e.g., 'More_Affected').
        cluster_key (str): The cell annotation set ('TNiche' or 'final_CT').
        hdf_path (str): The full path to the input HDF5 file.
    """
    # --- Step 1: Load and Rank Data ---
    key_to_load = f"/{group}_{cluster_key}_raw_zscore"
    try:
        with pd.HDFStore(hdf_path, mode='r') as store:
            zscore_df = store[key_to_load]
    except KeyError:
        print(f"Error: Could not find data for key '{key_to_load}'.")
        return

    try:
        # Ranking is still done on raw scores to get the correct order
        ranked_neighbors = zscore_df.loc[target_cell_type].sort_values(ascending=False)
    except KeyError:
        print(f"Error: Cell type '{target_cell_type}' not found for {group} using {cluster_key}.")
        print(f"Available types: {zscore_df.index.to_list()}")
        return

    # --- Step 2: Log-transform the scores BEFORE printing ---
    log_transformed_scores = np.sign(ranked_neighbors) * np.log10(np.abs(ranked_neighbors) + 1)

    # --- Step 3: Print the TRANSFORMED results ---
    print(f"--- Cell Types Adjacent to {target_cell_type} in '{group}' ---")
    print(f"(Ranked by Z-score, showing Signed Log10 Values)")
    print(log_transformed_scores) # Now prints the log-transformed series
    print("\n")

    # --- Step 4: Create the Bar Plot ---
    palette = ["royalblue" if x >= 0 else "red" for x in log_transformed_scores.values]

    plt.figure(figsize=(10, 8))
    sns.barplot(
        x=log_transformed_scores.values,
        y=log_transformed_scores.index,
        palette=palette
    )

    plt.title(f"Neighborhood Enrichment around {target_cell_type} in {group.replace('_', ' ')}", fontsize=16)
    plt.xlabel("Signed log10(z-score + 1)", fontsize=14)
    plt.ylabel("Neighboring Cell Type", fontsize=14)
    plt.axvline(x=0, color='black', linestyle='--', linewidth=1.5)

    plt.tight_layout()
    plt.show()

