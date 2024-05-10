#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH -c 10
#SBATCH --job-name=transcriptome_cov
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/transcriptome_cov/transcriptome_cov.%a.txt
#SBATCH -e logs/transcriptome_cov/transcriptome_cov.%a.txt
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

transcriptome_input_bam=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/transcriptome/GENCODE/sorted/primary_over_30/${sample}_sorted.bam
transcriptome_coverage=/dcs04/hicks/data/sparthib/retina_lrs/05c_coverage/transcriptome

#use bedtools to convert bam to transcritpome coverage files 
bedtools genomecov -ibam $transcriptome_input_bam -d  > $transcriptome_coverage/${sample}.coverage.txt

echo "**** Job ends ****"
date +"%Y-%m-%d %T"

