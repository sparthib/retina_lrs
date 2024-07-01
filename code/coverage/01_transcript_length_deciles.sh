#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=90G
#SBATCH -c 10
#SBATCH --job-name=transcript_length_deciles
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/transcript_length_deciles.%a.txt
#SBATCH -e logs/transcript_length_deciles.%a.txt
#SBATCH --array=1-15

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
BAM_FOLDER=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only/
sample=$(awk -v Index=${SLURM_ARRAY_TASK_ID} '$1==Index {print $2}' $CONFIG)
echo "$sample"

ml load conda_R/4.3.x

Rscript 01_transcript_length_deciles.R $sample

echo "**** Job ends ****"
date +"%Y-%m-%d %T"
