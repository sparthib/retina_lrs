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
#SBATCH --array=1
#SBATCH -t 4-00:00:00

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
BAM_FOLDER=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE/MAPQ_FILTERED
# REFERENCE_GTF=/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf
REFERENCE_GTF=/users/sparthib/LIQA_example/sample.gtf
STEP1_OUTPUT=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/LIQA/LIQA_STEP_1
STEP2_OUTPUT=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/LIQA/step2_outputs
STEP_3_OUTPUT=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/LIQA/step3_outputs


rm $STEP1_OUTPUT/step1.refgene
touch $STEP1_OUTPUT/step1.refgene
# mkdir -p $STEP2_OUTPUT
# mkdir -p $STEP_3_OUTPUT

##STEP 1 - preprocess GTF
liqa -task refgene -ref $REFERENCE_GTF -format gtf -out $STEP1_REFGENE/step1.refgene

##STEP 2 - quantify isoform expression 

# liqa -task quantify -refgene $STEP1_OUTPUT/step1.refgene -bam ${BAM_FOLDER}/${sample}_sorted.bam -out $STEP2_OUTPUT/$sample -max_distance 20 -f_weight 1

conda deactivate

echo "**** Job ends ****"
date +"%Y-%m-%d %T"