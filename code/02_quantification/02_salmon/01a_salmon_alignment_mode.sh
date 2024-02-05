#!/bin/bash

#SBATCH -p shared
#SBATCH -p shared
#SBATCH --mem=10G
#SBATCH --cpus-per-task=10
#SBATCH --job-name=salmon
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/salmon.%a.txt
#SBATCH -e logs/salmon.%a.txt
#SBATCH --array=1-12

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
BAM_FOLDER=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/transcriptome/GENCODE
REFERENCE_FASTA=/dcs04/hicks/data/sparthib/references/transcriptome/GENCODE/gencode.v44.transcripts_short_header.fa


OUTPUT_FOLDER=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/salmon/alignment_mode/$sample
rm -r $OUTPUT_FOLDER
mkdir -p $OUTPUT_FOLDER

salmon quant -t $REFERENCE_FASTA --libType U -a $BAM_FOLDER/$sample.bam  -o $OUTPUT_FOLDER --ont -p 10
 
conda deactivate

echo "**** Job ends ****"

date +"%Y-%m-%d %T"