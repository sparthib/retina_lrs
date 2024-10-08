#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=200G
#SBATCH -c 20
#SBATCH --job-name=generate_RCs
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/generate_RCs.%a.txt
#SBATCH -e logs/generate_RCs.%a.txt
#SBATCH --time=7-00:00:00
#SBATCH --array=1-15

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Array ID: ${SLURM_ARRAY_TASK_ID}"

CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
mkdir /dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/primary_assembly/$sample/

echo "**** Processing sample $sample ****"

module load conda_R/4.3.x
Rscript 01_generate_sample_wise_read_class.R $sample

echo "**** Job ends ****"
date
