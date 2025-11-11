#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=150G
#SBATCH --job-name=flames_multisample
#SBATCH -c 20
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/multisample.%a.txt
#SBATCH -e logs/multisample.%a.txt
#SBATCH --time=7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
mkdir /dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flames/$sample/

#create symlink for FASTQ and bam file 
# ln -s /dcs04/hicks/data/sparthib/retina_lrs/03a_nanofilt_fastqs/${sample}.fastq /dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flames/$sample/matched_reads.fastq
# ln -s /dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality/${sample}_primary_over_30_sorted.bam /dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flames/$sample/matched_reads_align2genome.bam
# ln -s /dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality/${sample}_primary_over_30_sorted.bam.bai /dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flames/$sample/matched_reads_align2genome.bam.bai

# RO_samples=("EP1-BRN3B-RO" "EP1-WT_hRO_2" "EP1-WT_ROs_D45" "H9-BRN3B_hRO_2" "H9-BRN3B-RO" "H9-CRX_hRO_2" "H9-CRX_ROs_D45")
# 
#   
# for sample in "${RO_samples[@]}"; do
#   # Construct source file path
#   source_file="/dcs04/hicks/data/sparthib/retina_lrs/03a_nanofilt_fastqs/${sample}.fastq"
#   # Construct target symlink path
#   target_link="/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flames/ROs/fastqs/$sample.fastq"
#   # Create the symlink
#   ln -s "$source_file" "$target_link"
#   
#   # Optional: echo to confirm what's happening
#   echo "Created symlink for $sample"
# done


# 
# FT_RGC_samples=("H9-FT_1" "H9-hRGC_1" "H9-FT_2" "H9-hRGC_2" )
# 
# for sample in "${FT_RGC_samples[@]}"; do
#   # Construct source file path
#   source_file="/dcs04/hicks/data/sparthib/retina_lrs/03a_nanofilt_fastqs/${sample}.fastq"
#   # Construct target symlink path
#   target_link="/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flames/FT_vs_RGC/fastqs/$sample.fastq"
#   # Create the symlink
#   ln -s "$source_file" "$target_link"
#   
#   # Optional: echo to confirm what's happening
#   echo "Created symlink for $sample"
# done

# 
# FT_RGC_samples=("H9-FT_1" "H9-hRGC_1" "H9-FT_2" "H9-hRGC_2" )
# 
# for sample in "${FT_RGC_samples[@]}"; do
#   # Construct source file path
#   source_file="/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality/${sample}_primary_over_30_sorted.bam"
#   # Construct target symlink path
#   target_link="/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flames/FT_vs_RGC/${sample}_align2genome.bam"
#   # Create the symlink
#   ln -s "$source_file" "$target_link"
#   
#   # Optional: echo to confirm what's hppening
#   echo "Created symlink for $sample"
#   
#   source_file="/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality/${sample}_primary_over_30_sorted.bam.bai"
#   # Construct target symlink path
#   target_link="/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flames/FT_vs_RGC/${sample}_align2genome.bam.bai"
#   # Create the symlink
#   ln -s "$source_file" "$target_link"
#   
# done


echo "**** Processing sample $sample ****"
module load conda_R/4.4.x
Rscript 02_multisample.R

echo "**** Job ends ****"
date +"%Y-%m-%d %T"