#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=200G
#SBATCH -c 10
#SBATCH --job-name=extended_annotation
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
Rscript 02_bambu_generate_rcs.R 
echo "**** Job ends ****"
date +"%Y-%m-%d %T"

