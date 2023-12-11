#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=100G,h_vmem=100G,h_fsize=100G
#$ -N salmon_transcript_level
#$ -o logs/salmon_transcript_level.\$TASK_ID.txt
#$ -e logs/salmon_transcript_level.\$TASK_ID.txt
#$ -m e
#$ -M sparthi1@jhu.edu
#$ -t 1-4
#$ -tc 20

echo "**** Job starts ****"
date +"%Y-%m-%d %T"

echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${SGE_TASK_ID}"

ml load conda
source activate salmon

CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=${SGE_TASK_ID} '$1==Index {print $2}' $CONFIG)
echo "${sample}"
BAM_FOLDER=/dcs04/hicks/data/sparthib/casey/bams_3
REFERENCE_FASTA=/dcs04/hicks/data/sparthib/ENSEMBLE_CDNA.fa 


OUTPUT_FOLDER=/dcs04/hicks/data/sparthib/casey/salmon_outputs_transcript_level/ensembl_fa/$sample
mkdir -p $OUTPUT_FOLDER

salmon quant -t $REFERENCE_FASTA --libType U -a $BAM_FOLDER/$sample.bam  -o $OUTPUT_FOLDER --ont --noErrorModel -p 10
 
conda deactivate

echo "**** Job ends ****"

date +"%Y-%m-%d %T"