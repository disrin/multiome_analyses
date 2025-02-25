import scanpy as sc
import argparse
import os

parser = argparse.ArgumentParser(description="QC for each sample")
parser.add_argument('--sample_name', type=str, help="Name of the sample (xxxx-xx format)")
parser.add_argument('--output', type=str, help="Path to save the output .h5ad file")

# Parse the arguments
args = parser.parse_args()

# Load the .h5 file
input_file = f"/path/to/input/directory/with/Samples/{args.sample_name}/cellbender_gex_counts_filtered.h5"
adata = sc.read_10x_h5(input_file)

# Ensure unique gene names
adata.var_names_make_unique()

# Add mitochondrial and ribosomal gene markers
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['rb'] = adata.var_names.str.startswith(('RPL', 'RPS'))

# Calculate QC metrics (e.g., number of genes expressed, mitochondrial percentage, ribosomal percentage)
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'rb'], percent_top=None, log1p=False, inplace=True)

# Run Scrublet for doublet detection
import scrublet as scr
doublet_detector = scr.Scrublet(adata.X, expected_doublet_rate=(adata.n_obs / 1000) * 0.008)  # Adjust rate if needed
doublet_scores, predicted_doublets = doublet_detector.scrub_doublets()

# Add doublet scores and predictions to adata.obs
adata.obs['cell_barcode'] = adata.obs_names
adata.obs['doublet_scores'] = doublet_scores
adata.obs['predicted_doublet'] = predicted_doublets

# Create the output directory if it doesn't exist
output_dir = os.path.dirname(args.output)
os.makedirs(output_dir, exist_ok=True)

# Save the processed .h5ad file
adata.write(args.output)
print(f"Processed data saved to: {args.output}")
