#!/bin/bash

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --mem=40G
#SBATCH -c 10
#SBATCH --job-name=blaze
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/blaze/raw_high_sensitivity.%a.txt
#SBATCH -e logs/blaze/raw_high_sensitivity.%a.txt
#SBATCH --array=1-4
#SBATCH --time=7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

CONFIG=/users/sparthib/retina_lrs/raw_data/single_cell.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
num_cells=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $5}' $CONFIG)
path=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $3}' $CONFIG)
seq_sum=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $4}' $CONFIG)
echo $sample

source activate flair

input_fastq=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/01_input_fastqs/${sample}.fastq.gz
output_dir=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/03_blaze_processed/raw/high_sensitivity
mkdir -p $output_dir
output_prefix=${output_dir}/${sample}_

echo "processing input fastq"
blaze --expect-cells $num_cells --output-prefix $output_prefix \
--threads $SLURM_CPUS_PER_TASK  $input_fastq  --high-sensitivity-mode

echo "**** Job ends ****"
date +"%Y-%m-%d %T"

