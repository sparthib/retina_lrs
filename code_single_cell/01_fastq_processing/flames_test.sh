#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH --job-name=flames
#SBATCH -c 20
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/flames/10x_D200-EP1-1_B2.txt
#SBATCH -e logs/flames/10x_D200-EP1-1_B2.txt
#SBATCH --time=7-00:00:00


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

ml load conda_R/4.4.x
Rscript flames_test.R

echo "**** Job ends ****"
date +"%Y-%m-%d %T"