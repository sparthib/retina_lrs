#!/bin/bash

#SBATCH -p shared
#SBATCH -p shared
#SBATCH --mem=15G
#SBATCH --cpus-per-task=5
#SBATCH --job-name=switch
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/switch_alignment_mode.text
#SBATCH -e logs/switch_alignment_mode.text

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

ml load conda_R/4.3.x
Rscript 01_IsoformSwitchAnalyzer.R 

echo "**** Job ends ****"
date +"%Y-%m-%d %T"