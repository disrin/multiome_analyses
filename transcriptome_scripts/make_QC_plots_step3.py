import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import numpy as np

# Parse command-line arguments
sample_name = sys.argv[1]
output_dir = sys.argv[2]

# Debugging: Confirm output directory
print(f"Saving plots to the following directory: {output_dir}")

# Check if the output directory exists
if not os.path.exists(output_dir):
    print(f"Error: The directory does not exist: {output_dir}")
    sys.exit(1)  # Exit the script if the directory does not exist

# Read the .h5ad file for the specific sample
adata_raw = sc.read(f'/path/to/sample_directory/{sample_name}/post-cellbender_QC/{sample_name}-raw.h5ad')

# Get the values for mitochondrial, ribosomal, n_genes, and doublet scores
mito_pct = adata_raw.obs['pct_counts_mt']
ribosomal_pct = adata_raw.obs['pct_counts_rb']
n_genes = adata_raw.obs['n_genes_by_counts']
doublet_scores = adata_raw.obs['doublet_scores']  # Corrected to match your adata key

# Count the total number of cell barcodes (cells) in the AnnData object
total_cells = adata_raw.obs.shape[0]  # Number of observations (cells)

# Define colors for the plots
mt_color = 'skyblue'  # Light blue for mitochondrial
rb_color = '#1f3c75'  # Darker blue for ribosomal
genes_color = 'lightgreen'  # Color for number of genes
doublets_color = 'salmon'  # Color for doublet score

# Threshold values
mito_threshold = 15  # 15% threshold for mitochondrial percentage
ribo_threshold = 10  # 10% threshold for ribosomal percentage
genes_threshold = 300  # 300 genes threshold for n_genes_by_counts
doublets_threshold = 0.2  # 0.2 threshold for doublet score

# Create the output directory for plots
# os.makedirs(output_dir, exist_ok=True)

# Plot 1: Mitochondrial Percentage
fig, ax = plt.subplots(1, 2, figsize=(12, 6))

# Mitochondrial Percentage Distribution (Histogram)
ax[0].hist(mito_pct, bins=30, color=mt_color, edgecolor='black', alpha=0.7)
ax[0].set_title('Mitochondrial Percentage Distribution')
ax[0].set_xlabel('Percent Mitochondria')
ax[0].set_ylabel('Number of Cells')
ax[0].set_xlim(0, 100)
ax[0].plot([mito_threshold, mito_threshold], [0, ax[0].get_ylim()[1]], '--r', label=f'{mito_threshold}% Threshold')
ax[0].legend(loc="upper right")

# Mitochondrial Percentage Violin Plot
sc.pl.violin(adata_raw, ['pct_counts_mt'], jitter=0.5, ax=ax[1], show=False, color=mt_color)
ax[1].plot([-.5, .5], [mito_threshold, mito_threshold], '--r', label=f'{mito_threshold}% Threshold')
ax[1].set_ylabel('Percent Mitochondria')
ax[1].set_title('Mitochondrial Percentage per Cell')
ax[1].legend()

# Add total number of cells (barcodes) on top of the plots
for axis in ax:
    axis.text(0.5, 1.05, f'Total cells: {total_cells}', ha='center', va='bottom', transform=axis.transAxes, fontsize=10)

# Save the mitochondrial distribution plot
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'mito_distribution.png'), dpi=300)

# Plot 2: Ribosomal Percentage
fig, ax = plt.subplots(1, 2, figsize=(12, 6))

# Ribosomal Percentage Distribution (Histogram)
ax[0].hist(ribosomal_pct, bins=30, color=rb_color, edgecolor='black', alpha=0.7)
ax[0].set_title('Ribosomal Percentage Distribution')
ax[0].set_xlabel('Percent Ribosomal')
ax[0].set_ylabel('Number of Cells')
ax[0].set_xlim(0, 100)
ax[0].plot([ribo_threshold, ribo_threshold], [0, ax[0].get_ylim()[1]], '--r', label=f'{ribo_threshold}% Threshold')
ax[0].legend(loc="upper right")

# Ribosomal Percentage Violin Plot
sc.pl.violin(adata_raw, ['pct_counts_rb'], jitter=0.5, ax=ax[1], show=False, color=rb_color)
ax[1].plot([-.5, .5], [ribo_threshold, ribo_threshold], '--r', label=f'{ribo_threshold}% Threshold')
ax[1].set_ylabel('Percent Ribosomal')
ax[1].set_title('Ribosomal Percentage per Cell')
ax[1].legend()

# Add total number of cells (barcodes) on top of the plots
for axis in ax:
    axis.text(0.5, 1.05, f'Total cells: {total_cells}', ha='center', va='bottom', transform=axis.transAxes, fontsize=10)

# Save the ribosomal distribution plot
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'ribo_distribution.png'), dpi=300)

# Plot 3: Number of Genes Detected
fig, ax = plt.subplots(1, 2, figsize=(12, 6))

# Number of Genes Detected Distribution (Histogram)
ax[0].hist(n_genes, bins=30, color=genes_color, edgecolor='black', alpha=0.7)
ax[0].set_title('Number of Genes Detected Distribution')
ax[0].set_xlabel('Number of Genes Detected')
ax[0].set_ylabel('Number of Cells')
ax[0].set_xlim(0, max(n_genes))
ax[0].plot([genes_threshold, genes_threshold], [0, ax[0].get_ylim()[1]], '--r', label=f'{genes_threshold} Genes Threshold')
ax[0].legend(loc="upper right")

# Number of Genes Detected Violin Plot
sc.pl.violin(adata_raw, ['n_genes_by_counts'], jitter=0.5, ax=ax[1], show=False, color=genes_color)
ax[1].plot([-.5, .5], [genes_threshold, genes_threshold], '--r', label=f'{genes_threshold} Genes Threshold')
ax[1].set_ylabel('Number of Genes Detected')
ax[1].set_title('Number of Genes Detected per Cell')
ax[1].legend()

# Add total number of cells (barcodes) on top of the plots
for axis in ax:
    axis.text(0.5, 1.05, f'Total cells: {total_cells}', ha='center', va='bottom', transform=axis.transAxes, fontsize=10)

# Save the number of genes distribution plot
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'genes_distribution.png'), dpi=300)

# Plot 4: Doublet Score
fig, ax = plt.subplots(1, 2, figsize=(12, 6))

# Doublet Score Distribution (Histogram)
ax[0].hist(doublet_scores, bins=30, color=doublets_color, edgecolor='black', alpha=0.7)
ax[0].set_title('Doublet Score Distribution')
ax[0].set_xlabel('Doublet Score')
ax[0].set_ylabel('Number of Cells')
ax[0].set_xlim(0, 1)
ax[0].plot([doublets_threshold, doublets_threshold], [0, ax[0].get_ylim()[1]], '--r', label=f'{doublets_threshold} Threshold')
ax[0].legend(loc="upper right")

# Doublet Score Violin Plot
sc.pl.violin(adata_raw, ['doublet_scores'], jitter=0.5, ax=ax[1], show=False, color=doublets_color)
ax[1].plot([-.5, .5], [doublets_threshold, doublets_threshold], '--r', label=f'{doublets_threshold} Threshold')
ax[1].set_ylabel('Doublet Score')
ax[1].set_title('Doublet Score per Cell')
ax[1].legend()

# Add total number of cells (barcodes) on top of the plots
for axis in ax:
    axis.text(0.5, 1.05, f'Total cells: {total_cells}', ha='center', va='bottom', transform=axis.transAxes, fontsize=10)

# Save the doublet score distribution plot
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'doublet_distribution.png'), dpi=300)
