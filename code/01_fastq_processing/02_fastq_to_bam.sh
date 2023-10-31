#!/bin/bash

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --job-name=fastq2bam
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/fastq2bam.%a.txt
#SBATCH -e logs/fastq2bam.%a.txt
#SBATCH --array=1-4

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"


CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
INPUT_FOLDER=/dcs04/hicks/data/sparthib/casey/fastqs
REFERENCE_FASTA=/dcs04/hicks/data/sparthib/GENCODE_FASTA.fa
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
SAM_FOLDER=/dcs04/hicks/data/sparthib/casey/sams
BAM_FOLDER=/dcs04/hicks/data/sparthib/casey/bams

mkdir -p $SAM_FOLDER
mdir -p $BAM_FOLDER

cd ~/minimap2
./minimap2 -ax splice -uf map-ont --secondary=no $REFERENCE_FASTA ${INPUT_FOLDER}/${sample}.fastq.gz > ${SAM_FOLDER}/${sample}.sam

ml load samtools
samtools view -bS -@ $SLURM_NTASKS_PER_NODE  ${SAM_FOLDER}/${sample}.sam -o ${BAM_FOLDER}/${sample}.bam
samtools sort ${BAM_FOLDER}/${sample}.bam -o ${BAM_FOLDER}/${sample}_sorted.bam
samtools index ${BAM_FOLDER}/${sample}_sorted.bam ${BAM_FOLDER}/${sample}_sorted.bam.bai
samtools idxstats ${BAM_FOLDER}/${sample}_sorted.bam

echo "**** Job ends ****"
date +"%Y-%m-%d %T"