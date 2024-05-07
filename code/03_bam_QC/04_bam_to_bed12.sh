#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=40G
#SBATCH -c 20
#SBATCH --job-name=bamtobed
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/bam2bed12/genome_splice/bam2bed12.%a.txt
#SBATCH -e logs/bam2bed12/genome_splice/bam2bed12.%a.txt
#SBATCH --array=1-15


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "****"


CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=${SLURM_ARRAY_TASK_ID} '$1==Index {print $2}' $CONFIG)
echo $sample

input_bam=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/${sample}_sorted.bam
bed12_folder=/dcs04/hicks/data/sparthib/retina_lrs/05b_beds/genome/GENCODE_splice/


#use bedtools to convert bam to bed files 
bedtools bamtobed -bed12 -i $input_bam  > $bed12_folder/${sample}.bed12


echo "**** Job ends ****"
date +"%Y-%m-%d %T"

