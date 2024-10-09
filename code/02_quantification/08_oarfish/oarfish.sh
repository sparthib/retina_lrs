#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH -c 10
#SBATCH --job-name=oarfish
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/oarfish.%a.txt
#SBATCH -e logs/oarfish.%a.txt
#SBATCH --array=1-12

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
echo "$sample"

bam_dir="/dcs04/hicks/data/sparthib/retina_lrs/05_bams/transcriptome/GENCODE/"
bam=$bam_dir/${sample}.bam

source activate oarfish

output_dir="/dcs04/hicks/data/sparthib/06_quantification/oarfish/$sample"
mkdir -p $output_dir
oarfish --verbose --output $output_dir --alignments $bam --model-coverage 
