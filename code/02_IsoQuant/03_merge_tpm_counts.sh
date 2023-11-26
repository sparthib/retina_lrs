#!/bin/bash

#SBATCH -p shared
#SBATCH -p shared
#SBATCH --mem=25G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=counts_tpm
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/counts_tpm.%a.txt
#SBATCH -e logs/counts_tpm.%a.txt
#SBATCH --array=1-4

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.3

Rscript 03_merge_tpm_counts.R


echo "**** Job ends ****"
date +"%Y-%m-%d %T"



