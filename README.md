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
- **Reads the input into an AnnData object.**
- **Ensures gene names are unique.**
- **Identifies mitochondrial and ribosomal genes for QC.**
- **Calculates quality control metrics (like mitochondrial and ribosomal gene expression).**
- **Runs Scrublet for detecting doublets and adds the results to the dataset.**
- **Creates a directory for saving output files.**
- **Saves the entire AnnData object, including expression matrix, annotations, QC metrics, and other calculations (like doublet scores), into a new .h5ad file for next steps.**
