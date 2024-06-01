#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=200G
#SBATCH -c 10
#SBATCH --job-name=test_bambu
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/extended_annotation.txt
#SBATCH -e logs/extended_annotation.txt
#SBATCH --time=7-00:00:00


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"


module load conda_R/4.3.x
Rscript 02a_bambu_generate_rcs_ROs.R 
Rscript 02b_bambu_generate_rcs_RGC_FT.R

echo "**** Job ends ****"
date