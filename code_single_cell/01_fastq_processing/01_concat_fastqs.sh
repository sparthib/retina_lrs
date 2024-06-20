#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=25G
#SBATCH --job-name=concat_fqs
#SBATCH -c 2
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/concat/concat_fqs.%a.txt
#SBATCH -e logs/concat/concat_fqs.%a.txt
#SBATCH --array=5-12

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"


CONFIG=/users/sparthib/retina_lrs/raw_data/single_cell.config
INPUT_FOLDER=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $3}' $CONFIG)
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
echo $sample
OUTPUT_FOLDER=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/01_input_fastqs

cd $INPUT_FOLDER/fastq_pass/
cat *.fastq.gz > $OUTPUT_FOLDER/${sample}.fastq.gz

echo "**** Job ends ****"
date +"%Y-%m-%d %T"