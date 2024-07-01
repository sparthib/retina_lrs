#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=10G
#SBATCH -c 5
#SBATCH --job-name=extract_read_ids
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/extract_read_ids.%a.txt
#SBATCH -e logs/extract_read_ids.%a.txt
#SBATCH --array=1-15


CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
BAM_FOLDER=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only/
sample=$(awk -v Index=${SLURM_ARRAY_TASK_ID} '$1==Index {print $2}' $CONFIG)
echo "$sample"

ml load samtools 
# Extract read ids from bam files
samtools view $BAM_FOLDER/${sample}_primary_over_30_chr_only_sorted.bam | awk '{if($1 ~ /^@/ || $1 ~ /NH:i:1/) print $1}' > $BAM_FOLDER/${sample}_read_ids.txt