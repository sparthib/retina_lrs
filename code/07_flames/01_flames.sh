#!/bin/bash

#SBATCH -p shared
#SBATCH -p shared
#SBATCH --mem=150G
#SBATCH --cpus-per-task=20
#SBATCH --job-name=flames
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/flames.%a.txt
#SBATCH -e logs/flames.%a.txt
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

source activate FLAMES
ml load samtools

CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
echo "${sample}"
BAM_FOLDER=/dcs04/hicks/data/sparthib/casey/bams
REFERENCE_FASTA=/dcs04/hicks/data/sparthib/ENSEMBLE_CDNA.fa.gz       
OUTPUT_FOLDER=/dcs04/hicks/data/sparthib/casey/flames_output/$sample
INPUT_FASTQ=/dcs04/hicks/data/sparthib/casey/fastqs_post_qc/$sample
mkdir -p $OUTPUT_FOLDER


~/flames/python/bulk_long_pipeline.py \
    -i $INPUT_FASTQ \
    --gff3 /dcs04/hicks/data/sparthib/ENSEMBL_GTF.gtf.gz \
    --genomefa $REFERENCE_FASTA \
    --config_file /users/sparthib/retina_lrs/code/07_flames/config.json \
    --outdir $OUTPUT_FOLDER
   

conda deactivate

echo "**** Job ends ****"

date +"%Y-%m-%d %T"