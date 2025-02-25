import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import muon as mu
import os
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Process a single sample for QC analysis.")
parser.add_argument('--sample_name', type=str, required=True, help="Name of the sample to process.")
parser.add_argument('--output', type=str, required=True, help="Path to save the output .h5ad file.")

# Parse arguments
args = parser.parse_args()

sample_name = args.sample_name
output_file = args.output

# Define the input and output directories for the current sample
input_folder = f"/path/to/sample_directory/{sample_name}/post-cellbender_QC"
h5_file = f"/path/to/sample_directory/{sample_name}/post-cellbender_QC/{sample_name}-raw.h5ad"

# Check if the input file exists
if not os.path.exists(h5_file):
    print(f"File not found: {h5_file}")
    exit(1)

# Load the data
adata = sc.read_h5ad(h5_file)
print(f"Processing {sample_name}...")


# Filter data
mito_percent_thresh = 15
ribo_percent_thresh = 10
doublet_thresh = 0.2
min_genes_per_cell = 300

# Apply filtering using muon
mu.pp.filter_obs(adata, 'pct_counts_mt', lambda x: x <= mito_percent_thresh)  # Keep cells with low mito percent
mu.pp.filter_obs(adata, 'pct_counts_rb', lambda x: x <= ribo_percent_thresh)  # Keep cells with low ribo percent
mu.pp.filter_obs(adata, 'doublet_scores', lambda x: x < doublet_thresh)  # Keep cells with low doublet score
mu.pp.filter_obs(adata, 'n_genes_by_counts', lambda x: x >= min_genes_per_cell)  # Keep cells with enough genes

# Save the filtered counts (after filtering)
adata.layers['filtered_counts'] = adata.X.copy()

# Print remaining cells and genes after filtering
print(f"Remaining cells: {adata.shape[0]}")
print(f"Number of cells after filtering: {adata.n_obs}")
print(f"Number of genes after filtering: {adata.n_vars}")

# Normalize data
sc.pp.normalize_total(adata, target_sum=1e6)
adata.layers['cpm'] = adata.X.copy()  # Save normalized counts
sc.pp.log1p(adata)
adata.layers['data'] = adata.X.copy()  # Save log-transformed data

#Reduce Dimensionality
sc.pp.pca(adata, n_comps=50)

# Calculate nearest neighbors
sc.pp.neighbors(adata)

# Visualize with UMAP
sc.tl.umap(adata)

# Cell cycle analysis
try:
    # Load cell cycle gene list
    cell_cycle_genes = [x.strip() for x in open('/path/to/lab_cell_cycle_genes.txt')]  # Adjust path if necessary
    s_genes = cell_cycle_genes[:43]  # S-phase genes
    g2m_genes = cell_cycle_genes[43:]  # G2/M-phase genes
    
    # Calculate cell cycle phase scores
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    
    # Check top cells in S and G2M phases
    high_s_phase = adata.obs[adata.obs['S_score'] > 0.5]
    high_g2m_phase = adata.obs[adata.obs['G2M_score'] > 0.5]
    
    print("Cells with high S-phase score:")
    print(high_s_phase[['S_score', 'G2M_score']].head())
    
    print("\nCells with high G2/M-phase score:")
    print(high_g2m_phase[['S_score', 'G2M_score']].head())
except Exception as e:
    print(f"Error in cell cycle analysis: {e}")

#Save adata to the output folder
adata.write(output_file, compression='gzip')

# Check the contents of each layer
print("\n--- Checking the contents of each layer ---")
for layer_name in adata.layers:
    print(f"Layer: {layer_name}")
    print(f"  Shape: {adata.layers[layer_name].shape}")
    print(f"  Preview (first 5 values): {adata.layers[layer_name].toarray()[:5]}")
    print()

# Plot UMAP for the entire dataset
sc.pl.umap(
    adata,
    color=['pct_counts_mt', 'pct_counts_rb', 'doublet_scores', 'total_counts', 'n_genes_by_counts', 'S_score', 'G2M_score'],
    size=2,
    ncols=3,
    cmap='viridis',
    show=False
)
# Save UMAP plot
umap_plot_file = f"/path/to/sample_directory/{sample_name}/post-cellbender_QC/{sample_name}_umap.png"
plt.savefig(umap_plot_file)
plt.close()

# Plot cell cycle phases
fig, axes = plt.subplots(1, 2, figsize=(12, 6))
sc.pl.scatter(adata, x='S_score', y='G2M_score', color='phase', title=f"S-phase vs G2/M-phase ({sample_name})", show=False, ax=axes[0])

phase_counts = adata.obs['phase'].value_counts()
phase_counts.plot(kind='bar', color=['#1f77b4', '#2ca02c', '#ff7f0e'], ax=axes[1])
axes[1].set_title(f'Cell Distribution Across Phases ({sample_name})')
axes[1].set_xlabel('Phase')
axes[1].set_ylabel('Number of Cells')
axes[1].tick_params(axis='x', rotation=0)

plt.tight_layout()
# Save cell cycle phase plot
phase_plot_file = f"/path/to/sample_directory/{sample_name}/post-cellbender_QC/{sample_name}_cell_cycle_phases.png"
plt.savefig(phase_plot_file)
plt.close()

print(f"Finished processing {sample_name}")

