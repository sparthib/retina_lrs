#!/bin/bash

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --mem=250G
#SBATCH -c 10
#SBATCH --job-name=fastq_qc
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/fastq_qc.%a.txt
#SBATCH -e logs/fastq_qc.%a.txt
#SBATCH --array=1-4

CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
input_fastq=/dcs04/hicks/data/sparthib/casey/fastqs/${sample}.fastq.gz
fastq_qc_OUTPUT=/dcs04/hicks/data/sparthib/casey/fastqs_post_qc/${sample}.fastq.gz

mkdir -p /dcs04/hicks/data/sparthib/casey/fastqs_post_qc


module load conda_R/4.3

Rscript 03_fastq_qc.R -i $input_fastq -o $fastq_qc_OUTPUT 


echo "**** Job ends ****"
date +"%Y-%m-%d %T"