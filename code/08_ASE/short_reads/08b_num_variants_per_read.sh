#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH -c 20
#SBATCH --job-name=num_var_per_read
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/num_var_per_read.%a.txt
#SBATCH -e logs/num_var_per_read.%a.txt
#SBATCH --time=7-00:00:00
#SBATCH --array=1-11

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
Rscript 08b_num_variants_per_read.R 

echo "**** Job ends ****"
date
