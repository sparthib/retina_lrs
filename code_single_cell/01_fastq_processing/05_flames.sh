#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=5G
#SBATCH --job-name=flames
#SBATCH -c 2
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/flames/flames.%a.txt
#SBATCH -e logs/flames/flames.%a.txt
#SBATCH --array=1-12
#SBATCH --time=7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

CONFIG=/users/sparthib/retina_lrs/raw_data/single_cell.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
num_cells=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $5}' $CONFIG)
path=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $3}' $CONFIG)
seq_sum=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $4}' $CONFIG)
echo $sample

#sce pipeline failed because it didn't find the matched_reads.fastq file
  #https://github.com/mritchielab/FLAMES/issues/30
  #work around is to create a symlink to the blaze processed fastq file 
#create symlink for the fastq files 
# mkdir -p /dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/$sample

#symlink for bam files 
# rm -r /dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/$sample
# mkdir -p /dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/$sample
# gunzip -c /dcs04/hicks/data/sparthib/retina_single_cell_lrs/03_blaze_processed/raw/high_sensitivity/${sample}_matched_reads.fastq.gz > /dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/$sample/matched_reads.fastq
# 
# 
unlink /dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/$sample/align2genome.bam
unlink /dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/$sample/align2genome.bam.bai

ln -s /dcs04/hicks/data/sparthib/retina_single_cell_lrs/04_minimap2_output/genome/bams/primary_over_30_chr_only/${sample}_primary_over_30_chr_only.bam /dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/$sample/align2genome.bam
ln -s /dcs04/hicks/data/sparthib/retina_single_cell_lrs/04_minimap2_output/genome/bams/primary_over_30_chr_only/${sample}_primary_over_30_chr_only.bam.bai /dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/$sample/align2genome.bam.bai


# ml load conda_R/4.4.x
# Rscript 05_flames.R "$sample"

echo "**** Job ends ****"
date +"%Y-%m-%d %T"