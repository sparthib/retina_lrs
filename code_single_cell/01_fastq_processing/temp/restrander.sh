#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH --job-name=restrander
#SBATCH -c 10
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/restrander/restrander.%a.txt
#SBATCH -e logs/restrander/restrander.%a.txt
#SBATCH --array=1-4


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"


CONFIG=/users/sparthib/retina_lrs/raw_data/single_cell.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
INPUT_FILE=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/02_nanofilt_processed/${sample}.fastq.gz
OUTPUT_FILE=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/restrander_output/${sample}.fastq.gz
STATS_FILE=/users/sparthib/retina_lrs/code_single_cell/01_fastq_processing/logs/restrander/${sample}_stats.json

cd /users/sparthib/restrander
./restrander \
    $INPUT_FILE \
    $OUTPUT_FILE \
    /users/sparthib/restrander/config/10X-3prime.json \
        > $STATS_FILE
        
echo "**** Job ends ****"
date +"%Y-%m-%d %T"
