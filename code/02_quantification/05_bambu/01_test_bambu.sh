#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=250G
#SBATCH --cpus-per-task=20
#SBATCH --job-name=test_bambu
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/test_bambu.txt
#SBATCH -e logs/test_bambu.txt
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
mkdir -p /dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/$sample

module load conda_R/4.3.x
Rscript 01_test_bambu.R 

echo "**** Job ends ****"
date