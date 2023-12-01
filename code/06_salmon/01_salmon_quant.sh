#!/bin/bash

#SBATCH -p shared
#SBATCH -p shared
#SBATCH --mem=150G
#SBATCH --cpus-per-task=20
#SBATCH --job-name=salmon
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/salmon.%a.txt
#SBATCH -e logs/salmon.%a.txt
#SBATCH --array=1-4

#try running for all chromosomes

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

source activate salmon

CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
echo "${sample}"
BAM_FOLDER=/dcs04/hicks/data/sparthib/casey/bams
REFERENCE_FASTA=/dcs04/hicks/data/sparthib/GENCODE_FASTA.fa 

rm -r /dcs04/hicks/data/sparthib/casey/salmon_output
OUTPUT_FOLDER=/dcs04/hicks/data/sparthib/casey/salmon_output/$sample
mkdir -p $OUTPUT_FOLDER
salmon quant -t $REFERENCE_FASTA -a ${BAM_FOLDER}/${sample}_sorted.bam --libType A \
 -o $OUTPUT_FOLDER -p ${SLURM_CPUS_PER_TASK} --dumpEq --numBootstraps 100 

 
conda deactivate

echo "**** Job ends ****"

date +"%Y-%m-%d %T"