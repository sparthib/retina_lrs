#!/bin/bash

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --mem=200G
#SBATCH -c 10
#SBATCH --job-name=fastq_qc
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/fastq_qc.%a.txt
#SBATCH -e logs/fastq_qc.%a.txt
#SBATCH --array=1-4

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

source activate nanofilt 

CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
input_fastq=/dcs04/hicks/data/sparthib/casey/fastqs/${sample}.fastq.gz
fastq_qc_OUTPUT=/dcs04/hicks/data/sparthib/casey/fastqs_post_qc/${sample}.fastq.gz
guppy_summary_file=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $4}' $CONFIG)


mkdir -p /dcs04/hicks/data/sparthib/casey/fastqs_post_qc


gunzip -c $input_fastq | NanoFilt -q 7 --length 50 -s $guppy_summary_file | gzip > $fastq_qc_OUTPUT

echo "**** Job ends ****"
date +"%Y-%m-%d %T"