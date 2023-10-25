#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=25G
#SBATCH --job-name=concat_fqs
#SBATCH -c 2
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/concat_fqs.4.txt
#SBATCH -e logs/concat_fqs.4.txt


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"


CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
#INPUT_FOLDER=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $3}' $CONFIG)
#sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
sample=DG-WT-hRGC
INPUT_FOLDER=/dcs04/hicks/data/globus/casey/230920_Casey/DG-WT-hRGCs/20230920_1343_3D_PAQ42790_791d991e/fastq_pass
OUTPUT_FOLDER=/dcs04/hicks/data/sparthib/casey/fastqs
mkdir -p $OUTPUT_FOLDER

cd $INPUT_FOLDER
cat *.fastq.gz > $OUTPUT_FOLDER/${sample}.fastq.gz

echo "**** Job ends ****"
date +"%Y-%m-%d %T"