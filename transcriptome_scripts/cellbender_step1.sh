#!/bin/bash
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=32G
#SBATCH --time 24:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100x:1
#SBATCH --array=0-39
1:28
module load cellbender
module load CUDA/12.1
# Run cellbender within this base directory
output_file_base=$(echo /path/to/input/directory/with/Samples/batch-num_or_Controls/Multiome/*/)
# Convert directory locations into an array
out_dirs=($(echo ${output_file_base}))
# Iterate through the array of sample directories(/base directory/outs)
cd ${out_dirs[$SLURM_ARRAY_TASK_ID]}/outs 
cellbender remove-background --input raw_feature_bc_matrix.h5 --output cellbender_gex_counts.h5 --fpr 0 --cuda
