#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=10G
#SBATCH --job-name=change_fa_header
#SBATCH -o logs/change_fa_header.txt
#SBATCH -e logs/change_fa_header.txt


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

ml purge
ml load python/3.11.8
python change_fa_header.py

echo "**** Job ends ****"
date +"%Y-%m-%d %T"