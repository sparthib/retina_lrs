#!/bin/bash

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH -c 10
#SBATCH --job-name=minIONQC
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/minionqc/minIONQC.%a.txt
#SBATCH -e logs/minionqc/minIONQC.%a.txt
#SBATCH --array=1-4

CONFIG=/users/sparthib/retina_lrs/raw_data/single_cell.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
path=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $3}' $CONFIG)
seq_sum=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $4}' $CONFIG)
echo $sample

input_fastq=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/01_input_fastqs/${sample}.fastq.gz
minion_qc_output=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/03_minion_qc_output/${sample}
mkdir -p $minion_qc_output
guppy_summary_file=${path}/${seq_sum}

mkdir -p $minion_qc_output

module load conda_R/4.3

Rscript minion_qc.R -i $guppy_summary_file -o $minion_qc_output -p ${SLURM_CPUS_PER_TASK}


echo "**** Job ends ****"
date +"%Y-%m-%d %T"