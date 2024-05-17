#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=40G
#SBATCH -c 10
#SBATCH --job-name=fastq_qc
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/02_Nanofilt/fastq_qc.%a.txt
#SBATCH -e logs/02_NanoFilt/fastq_qc.%a.txt
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

CONFIG=/users/sparthib/retina_lrs/raw_data/single_cell.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
path=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $3}' $CONFIG)
seq_sum=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $4}' $CONFIG)


input_fastq=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/01_input_fastqs/${sample}.fastq.gz
fastq_qc_output=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/02_nanofilt_processed/${sample}.fastq.gz
guppy_summary_file=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $4}' $CONFIG)

gunzip -c $input_fastq | NanoFilt -q 7 --length 50 -s ${path}${seq_sum} | gzip > $fastq_qc_output

echo "**** Job ends ****"
date +"%Y-%m-%d %T"