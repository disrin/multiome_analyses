#!/bin/bash
#SBATCH --job-name=QC_analysis
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --cpus-per-task 4
#SBATCH --mem-per-cpu=32G
#SBATCH --time 36:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100x:1
#SBATCH --array=0-63

#Source conda environment
source myconda

#Activate conda environment
mamba activate /path/to/conda_envs/PSP_RNA_batch1

#Define sample directory and sample name
SAMPLES=($(ls /path/to/input/directory/with/Samples))

#Get the sample index based on job array
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

#print the sample name being processed
echo "Processing sample: $SAMPLE"

# Define the path to the .h5 file for the current sample
h5_file="/path/to/input/directory/with/Samples/$SAMPLE/cellbender_gex_counts_filtered.h5"


# Check if the file exists before proceeding
if [[ -f "$h5_file" ]]; then
    echo "Processing file: $h5_file"

    OUTPUT_DIR="/path/to/input/directory/with/Samples/$SAMPLE/post-cellbender_QC"
    
    # Define the output .h5ad file path (with the '-raw' suffix)
    OUTPUT_FILE="$OUTPUT_DIR/${SAMPLE}-raw.h5ad"
    
    # Run the Python script for the current sample, passing the sample name and output file path
    python /path/to/scripts/QC_step2.py --sample_name "$SAMPLE" --output "$OUTPUT_FILE"
    
    echo "Finished processing sample: $SAMPLE"
else
    echo "File not found: $h5_file"
fi
