#!/bin/bash

#SBATCH -p shared
#SBATCH --mem-per-cpu=50G
#SBATCH --job-name=restrander
#SBATCH -c 10
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/lib1_restrander_split.%a.txt
#SBATCH -e logs/lib1_restrander_split.%a.txt
#SBATCH --array=1-4


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"


#mkdir for output if it doesn't exist 
mkdir -p /dcs04/hicks/data/sparthib/casey/fastqs/restrander
mkdir -p /users/sparthib/retina_lrs/code/01_fastq_processing/logs/restrander


CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
INPUT_FOLDER=/dcs04/hicks/data/sparthib/casey/fastqs
REFERENCE_FASTA=/dcs04/hicks/data/sparthib/GENCODE_FASTA.fa
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)


INPUT_FILE=/dcs04/hicks/data/sparthib/casey/fastqs/${sample}.fastq.gz
OUTPUT_FILE=/dcs04/hicks/data/sparthib/casey/fastqs/restrander/${sample}.fastq.gz
STATS_FILE=/users/sparthib/retina_lrs/code/01_fastq_processing/logs/restrander/${sample}_stats.json

cd /users/sparthib/restrander
./restrander \
    $INPUT_FILE \
    $OUTPUT_FILE \
    config/LSK114.json \
        > $STATS_FILE

echo "Memory Utilized: "
echo "**** Job ends ****"
date +"%Y-%m-%d %T"
