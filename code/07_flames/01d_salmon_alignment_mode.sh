#!/bin/bash

#SBATCH -p shared
#SBATCH -p shared
#SBATCH --mem=150G
#SBATCH --nodes=1
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
FASTQ=/dcs04/hicks/data/sparthib/casey/fastqs_post_qc
REFERENCE_FASTA=/dcs04/hicks/data/sparthib/ENSEMBLE_CDNA.fa.gz  

OUTPUT_FOLDER=/dcs04/hicks/data/sparthib/casey/salmon_outputs_transcript_level_01d/$sample
mkdir -p $OUTPUT_FOLDER

salmon quant -t transcripts.fa -l U -a aln.bam -o salmon_quant \
--dumpEq --useEM -o $OUTPUT_FOLDER 
 
conda deactivate

echo "**** Job ends ****"

date +"%Y-%m-%d %T"