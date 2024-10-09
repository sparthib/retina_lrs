#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH -c 10
#SBATCH --job-name=transcriptome_coverage
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/coverage.%a.txt
#SBATCH -e logs/coverage.%a.txt
#SBATCH --array=1-12

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
echo "$sample"


ml load conda_R/4.4.x
# Rscript transcriptome_coverage.R $sample
Rscript cov_by_categories.R $sample

echo "**** Job ends ****"
date +"%Y-%m-%d %T"

