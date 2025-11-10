#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=80G
#SBATCH --job-name=seqkit
#SBATCH -c 10
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/seqkit/seqkit.%a.txt
#SBATCH -e logs/seqkit/seqkit.%a.txt
#SBATCH --array=9-12

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

source activate seqkit

ENV_FILE="../../.env"
if [ -f $ENV_FILE ]; then
    set -a
    source $ENV_FILE
    set +a
fi

sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
echo $sample
INPUT_FOLDER=$retina_lrs_dir/01_input_fastqs
OUTPUT_FOLDER=$retina_lrs_dir/01a_seqkit_processed


zcat $INPUT_FOLDER/${sample}.fastq.gz |  seqkit rmdup --dup-num-file $OUTPUT_FOLDER/$sample.dup_ids.txt  \
--dup-seqs-file $OUTPUT_FOLDER/$sample.duplicated.fa.gz \
-o $OUTPUT_FOLDER/$sample.cleaned.fastq.gz \

conda deactivate 

echo "**** Job ends ****"
date +"%Y-%m-%d %T"