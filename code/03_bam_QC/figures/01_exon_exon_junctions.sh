#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH -c 20
#SBATCH --job-name=exon_exon_junctions
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/exon_exon/sample.%a.txt
#SBATCH -e logs/exon_exon/sample.%a.txt
#SBATCH --array=1-15
#SBATCH --time=7-00:00:00



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

ml load samtools

alignment_dir='/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only'
alignment="${alignment_dir}/${sample}_primary_over_30_chr_only_sorted.bam"

# 2. Compute number of junctions in read
output_file="/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only/CIGAR_sequences/${sample}_cigar_sequences.txt"
nums=()

samtools view $alignment | awk '{print $6}' > $output_file

echo "**** Job ends ****"
date +"%Y-%m-%d %T"
