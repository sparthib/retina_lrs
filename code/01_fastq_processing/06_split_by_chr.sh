#!/bin/bash

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH -c 5
#SBATCH --job-name=fastq_qc
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/06_split/gencode_splice/split.%a.txt
#SBATCH -e logs/06_split/gencode_splice/split.%a.txt
#SBATCH --array=13-15

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

ml load samtools
CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
echo "$sample"

output=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/${sample}_chromosome_level/
mkdir $output
bam=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/${sample}_sorted.bam


for chr in chr{1..22} chrX chrY chrM
do
    echo "Processing $chr"
    samtools view -b $bam $chr > ${output}/${sample}_${chr}.bam
    samtools index ${output}/${sample}_${chr}.bam ${output}/${sample}_${chr}.bam.bai
done

echo "**** Job ends ****"
date +"%Y-%m-%d %T"

