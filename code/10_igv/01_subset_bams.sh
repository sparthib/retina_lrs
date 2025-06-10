#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=30G
#SBATCH -c 10
#SBATCH --job-name=igv
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1,2,5-6,9-15
#SBATCH -o logs/igv.%a.txt
#SBATCH -e logs/igv.%a.txt
#SBATCH --time=7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Array job ID: ${SLURM_ARRAY_JOB_ID}"


CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=${SLURM_ARRAY_TASK_ID} '$1==Index {print $2}' $CONFIG)
echo $sample

bam_dir=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality
output_dir=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality/igv

mkdir -p $output_dir

# Load the required module
ml load samtools 

#REEP6 chromosome 	NC_000019.10 (1491181..1497927)
#ADD1 chromosome 	NC_000004.12 (2843844..2930062)

samtools view -b $bam_dir/${sample}_primary_over_30_chr_only_sorted.bam "chr19:1491181-1497927" > $output_dir/REEP/${sample}_REEP6.bam
samtools view -b $bam_dir/${sample}_primary_over_30_chr_only_sorted.bam "chr4:2843844-2930062" > $output_dir/ADD1/${sample}_ADD1.bam

## index bam files
for dir in REEP ADD1; do
    mkdir -p $output_dir/$dir
    samtools index $output_dir/$dir/${sample}_$dir.bam > $output_dir/$dir/${sample}_$dir.bam.bai
done


