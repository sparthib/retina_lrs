#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=30G
#SBATCH -c 20
#SBATCH --job-name=gene_counts
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/gene_counts.%a.txt
#SBATCH -e logs/gene_counts.%a.txt
#SBATCH --time=7-00:00:00
#SBATCH --array=1-7

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Array ID: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.3.x
# Rscript 06_gene_counts.R 
Rscript 06_gene_counts.R 

echo "**** Job ends ****"
date
