#!/bin/bash

#SBATCH -p shared
#SBATCH -p shared
#SBATCH --mem=25G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=counts_tpm
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/counts_tpm.%a.txt
#SBATCH -e logs/counts_tpm.%a.txt
#SBATCH --array=1-4

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"


CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
IsoQuant_dir=/dcs04/hicks/data/sparthib/casey/IsoQuant_output/${sample}/OUT
IsoQuant_tpm=$IsoQuant_dir/OUT.transcript_tpm.tsv
IsoQuant_counts=$IsoQuant_dir/OUT.transcript_counts.tsv

rm $IsoQuant_dir/tpm_counts_data.tsv
echo "id\ttpm\tcounts" > $IsoQuant_dir/tpm_counts_data.tsv
join -t $'\t' -1 1 -2 1 -o 1.1,1.2,2.2 <(sort -k1,1 $IsoQuant_tpm) <(sort -k1,1 $IsoQuant_counts) >> $IsoQuant_dir/tpm_counts_data.tsv



echo "**** Job ends ****"
date +"%Y-%m-%d %T"



