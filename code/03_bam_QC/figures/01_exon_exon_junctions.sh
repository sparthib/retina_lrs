#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH -c 20
#SBATCH --job-name=exon_exon_junctions
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/exon_exon/sample.%a.txt
#SBATCH -e logs/exon_exon/sample.%a.txt
#SBATCH --array=9




echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

LOGS_FOLDER=/users/sparthib/retina_lrs/code/01_fastq_processing/logs/bam_stats_transcriptome_gencode_splice
CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=${SLURM_ARRAY_TASK_ID} '$1==Index {print $2}' $CONFIG)
echo "$sample"

ml load conda_R/4.4.x
Rscript 01_exon_exon_junctions.R $sample

echo "**** Job ends ****"
date +"%Y-%m-%d %T"
