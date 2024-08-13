#!/bin/bash

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --mem=40G
#SBATCH -c 10
#SBATCH --job-name=pychopper
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/pychopper/log.%a.txt
#SBATCH -e logs/pychopper/log.%a.txt
#SBATCH --array=2

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

source activate pychopper

CONFIG=/users/sparthib/retina_lrs/raw_data/single_cell.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
echo $sample


input_dir=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/02_nanofilt_processed     
output_dir=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/03_pychopper_processed/$sample
mkdir -p $output_dir

pychopper -U -k LSK114 -r $output_dir/report.pdf -u $output_dir/unclassified.fq \
-S $output_dir/stats.txt \
-w $output_dir/rescued.fq $input_dir/$sample.fastq.gz $output_dir/${sample}_full_length_output.fq

conda deactivate 
echo "**** Job ends ****"
date +"%Y-%m-%d %T"
