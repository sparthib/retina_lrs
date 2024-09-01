#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH -c 5
#SBATCH --job-name=plot
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/plots.%a.txt
#SBATCH -e logs/plots.%a.txt
#SBATCH --array=9
#SBATCH --time=7-00:00:00


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "****"


CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=${SLURM_ARRAY_TASK_ID} '$1==Index {print $2}' $CONFIG)
echo $sample

ml load conda_R/4.4.x
Rscript 01f_plot_attributes.R $sample

echo "**** Job ends ****"
date +"%Y-%m-%d %T"