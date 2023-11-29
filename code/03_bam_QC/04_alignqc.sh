#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=80G
#SBATCH --job-name=alignqc
#SBATCH -c 10
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/alignqc.%a.txt
#SBATCH -e logs/alignqc.%a.txt
#SBATCH --array=1-4

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

source activate alignqc

# python3 -m pip list
# Package    Version
# ---------- -------
# AlignQC    2.0.5
# pip        23.3
# seq-tools  1.0.10
# setuptools 68.0.0
# wheel      0.41.2
#Successfully installed AlignQC-2.0.5 seq-tools-1.0.10
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
BAM_FOLDER=/dcs04/hicks/data/sparthib/casey/bams

REFERENCE_FASTA=/dcs04/hicks/data/sparthib/GENCODE_FASTA.fa
REFERENCE_GTF=/dcs04/hicks/data/sparthib/GENCODE_GTF.gtf
ALIGNQC_OUTPUT=/dcs04/hicks/data/sparthib/alignqc_output


alignqc analysis ${BAM_FOLDER}/${sample}_sorted.bam -g $REFERENCE_FASTA -t $REFERENCE_GTF --output_folder $ALIGNQC_OUTPUT/${sample}

conda deactivate 

echo "**** Job ends ****"
date +"%Y-%m-%d %T"
