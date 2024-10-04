#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH --job-name=flames
#SBATCH -c 20
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/flames_min_count2.%a.txt
#SBATCH -e logs/flames_min_count2.%a.txt
#SBATCH --array=1-15
#SBATCH --time=7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"

CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
mkdir /dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flames/$sample/

#create symlink for FASTQ and bam file 
ln -s /dcs04/hicks/data/sparthib/retina_lrs/03a_nanofilt_fastqs/${sample}.fastq /dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flames/$sample/matched_reads.fastq
ln -s /dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality/${sample}_primary_over_30_sorted.bam /dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flames/$sample/align2genome.bam
ln -s /dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality/${sample}_primary_over_30_sorted.bam.bai /dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flames/$sample/align2genome.bam.bai


echo "**** Processing sample $sample ****"
module load conda_R/4.4.x
Rscript 01_flames.R $sample
