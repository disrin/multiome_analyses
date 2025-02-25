#!/bin/bash
#SBATCH --job-name=make_QC_plots
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32G
#SBATCH --time 10:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100x:1
#SBATCH --array=0-63

#Source conda environment
source myconda

#Activate conda environment
mamba activate /path/to/conda_envs/PSP_RNA_batch1

#list samples
SAMPLES=($(ls /path/to/sample_directory)) ## | grep -E '^[A-Za-z0-9]{4}-[0-9]{2}-[A-Za-z]{3}$'))

#Get sample name - ensures the script processes one sample at a time 
sample_name=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
 
echo "Processing sample: $sample_name"

#Specify output directory
output_dir="/path/to/sample_directory/${sample_name}/post-cellbender_QC/"

# Run the Python script 
python /path/to/scripts/make_QC_plots_step3.py $sample_name $output_dir

echo "Completed processing for: $sample_name"
