#!/bin/bash

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH -c 5
#SBATCH --job-name=pomoxis
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/pomoxis.%a.log
#SBATCH -e logs/pomoxis.%a.log
#SBATCH --array=1-15


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

source activate pomoxis

CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=${SLURM_ARRAY_TASK_ID} '$1==Index {print $2}' $CONFIG)
echo $sample

BAM_FILE=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly${sample}_sorted.bam
OUT_FILE=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/$sample.stats

stats_from_bam $BAM_FILE > $OUT_FILE

BAM_FILE=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality/${sample}_primary_over_30_chr_only_sorted.bam
OUT_FILE=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality/$sample.stats

stats_from_bam $BAM_FILE > $OUT_FILE

echo "**** Job ends ****"
date +"%Y-%m-%d %T"