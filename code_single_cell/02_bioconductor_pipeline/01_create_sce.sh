#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH -c 5
#SBATCH --job-name=create_sce
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/create_sce.txt
#SBATCH -e logs/create_sce.txt
#SBATCH --time=7-00:00:00


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Array ID: ${SLURM_ARRAY_TASK_ID}"


ml load conda_R/4.4.x



echo "**** Job ends ****"
date