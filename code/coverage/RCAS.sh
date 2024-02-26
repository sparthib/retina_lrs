#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH --cpus-per-task=10
#SBATCH --job-name=test_RCAS
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/test_RCAS.txt
#SBATCH -e logs/test_RCAS.txt
#SBATCH --array=1
#SBATCH -t 4-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)

module load conda_R/4.3.x
Rscript RCAS.R 

echo "**** Job ends ****"
date