#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=40G
#SBATCH -c 10
#SBATCH --job-name=check_sequence
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --array=9-12
#SBATCH --output=logs/check_sequence.%a.log
#SBATCH --error=logs/check_sequence.%a.log
#SBATCH -t 7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config

sample=$(awk -v Index=${SLURM_ARRAY_TASK_ID} '$1==Index {print $2}' $CONFIG)
echo "$sample"
input_file="/dcs04/hicks/data/sparthib/retina_lrs/01_input_fastqs/$sample.fastq.gz"

source activate check_seq
python check_sequences.py $input_file

conda deactivate 
echo "**** Job ends ****"
date +"%Y-%m-%d %T"
