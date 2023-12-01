#!/bin/bash

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH -c 10
#SBATCH --job-name=liqa
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/liqa.%a.txt
#SBATCH -e logs/liqa.%a.txt
#SBATCH --array=1-4

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

source activate LIQA

CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
echo "${sample}"
BAM_FOLDER=/dcs04/hicks/data/sparthib/casey/bams
REFERENCE_GTF=/dcs04/hicks/data/sparthib/GENCODE_GTF.gtf
STEP1_REFGENE=/dcs04/hicks/data/sparthib/LIQA_STEP_1.refgene
STEP2_OUTPUT=/dcs04/hicks/data/sparthib/casey/liqa/step2_outputs
STEP_3_OUTPUT=/dcs04/hicks/data/sparthib/casey/liqa/step3_outputs
mkdir -p $STEP2_OUTPUT
mkdir -p $STEP_3_OUTPUT

##STEP 1 - preprocess GTF
liqa -task refgene -ref $REFERENCE_GTF -format gtf -out $STEP1_REFGENE

##STEP 2 - quantify isoform expression 

liqa -task quantify -refgene $STEP1_REFGENE -bam ${BAM_FOLDER}/${sample}_sorted.bam -out $STEP2_OUTPUT/$sample -max_distance 20 -f_weight 1

conda deactivate

echo "**** Job ends ****"
date +"%Y-%m-%d %T"