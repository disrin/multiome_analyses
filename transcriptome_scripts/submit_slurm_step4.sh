#!/bin/bash
#SBATCH --job-name=unfilter_norm_rna
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32G
#SBATCH --time=36:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100x:1
#SBATCH --array=0-63

# Load the conda environment
source myconda
conda activate /path/to/conda_envs/PSP_RNA_batch1

# Define the sample array
SAMPLES=($(ls /path/to/sample_directory/)) ## | grep -E '^[A-Za-z0-9]{4}-[0-9]{2}-[A-Za-z]{3}$'))

# Get the sample name based on the SLURM_ARRAY_TASK_ID
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Print the sample being processed
echo "Processing sample: $SAMPLE"

# Define input file paths and output file paths
h5_file="/path/to/sample_directory/$SAMPLE/post-cellbender_QC/${SAMPLE}-raw.h5ad"
OUTPUT_DIR="/path/to/sample_directory/$SAMPLE/post-cellbender_QC"
OUTPUT_FILE="${OUTPUT_DIR}/${SAMPLE}_norm_filtered.h5ad"

# Check if the h5ad file exists
if [[ -f "$h5_file" ]]; then
    echo "Found input file: $h5_file"
    
    # Run the Python script with sample name and output file path
    python /path/to/scripts/Norm_unfiltered_rna_step4.py --sample_name "$SAMPLE" --output "$OUTPUT_FILE"
    
    echo "Finished processing sample: $SAMPLE"
else
    echo "File not found: $h5_file"
fi
