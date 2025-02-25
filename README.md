# Multiome Analysis for 10X Data (RNA + ATAC)
## This is under active development
This repository contains code for performing multiome analysis on single-cell RNA (scRNA-seq) and ATAC-seq (chromatin accessibility) data from 10X Genomics. The code is designed to work after the CellRanger-ARC-2.0 pipeline has been run on individual samples and is built to handle multiple samples. For now, to run the respective Python scripts as batch jobs on an HPC system, submit the corresponding .sh files

Currently, the code is still under development for automation with Snakemake, and for now, the analysis uses conda environments.

## Transcriptome

### Input
The input paths - 
"/path/to/batch#/sample_directory/"
within the sample_directory are all Samples - Sample_001, Sample_002 and so on. 

### Step1: Cellbender
Removes ambient RNA from scRNAseq data. Although CellRanger pipeline performs its own QC filtering, we use cellbender to remove any lingering RNA backgroung. CellRanger output `raw_feature_bc_matrix.h5` serves as CellBender input.

### Step2: QC
CellBender output `cellbender_gex_counts_filtered.h5` serves as input. 

Summary of Step 2:
- **Input=`cellbender_gex_counts_filtered.h5`**
- **Reads the input into an AnnData object.**
- **Ensures gene names are unique.**
- **Identifies mitochondrial and ribosomal genes for QC.**
- **Calculates quality control metrics (like mitochondrial and ribosomal gene expression).**
- **Runs Scrublet for detecting doublets and adds the results to the dataset.**
- **Creates a directory for saving output files.**
- **Saves the entire AnnData object, including expression matrix, annotations, QC metrics, and other calculations (like doublet scores), into a new .h5ad file for next steps.**
- **Output=`{sample_name}-raw.h5ad`**

### Step3: Make QC plots
Here we set QC thresholds (can be modifed as needed in the .py file) and make plots for every QC parameter.

Summary of step 3:
- **Input= `{sample_name}-raw.h5ad`.**
- **Generates histograms and violin plots for mitochondrial percentage, ribosomal percentage, gene detection, and doublet score.**
- **Adds threshold lines to the plots to highlight cells that may be outliers based on QC metrics (e.g., high mitochondrial content, low gene detection).**
- **Displays the total number of cells (barcodes) on top of each plot for reference.**
- **Customizable Output: Saves the generated QC plots to the specified output directory as PNG files.**
- **Output = .png files

### Step4: Filter and Normalize
This code processes single-cell RNA-seq data by filtering, normalizing, and performing dimensionality reduction, followed by cell cycle analysis and visualization, saving results to output files.

Summary of step 4:
- **Input = `{sample_name}-raw.h5ad`**
- **Filters out cells with muon based on thresholds for mitochondrial percentage, ribosomal percentage, doublet score, and minimum gene count.**
- **Normalizes the data, applies log transformation, and saves these processed versions in different layers of the AnnData object.**
- **Performs PCA for dimensionality reduction, calculates nearest neighbors, and visualizes the data with UMAP.**
- **Loads a predefined list of cell cycle genes `lab_cell_cycle_genes.txt`, calculates S-phase and G2/M-phase scores, and generates visualizations for the UMAP and cell cycle phase distribution, saving them as image files.**

