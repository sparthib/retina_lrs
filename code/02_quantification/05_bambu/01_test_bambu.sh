#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH -c 10
#SBATCH --job-name=test_bambu
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/bambu_quant.%a.txt
#SBATCH -e logs/bambu_quant.%a.txt
#SBATCH --time=7-00:00:00
#SBATCH --array=1-12

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echp "Array ID: ${SLURM_ARRAY_TASK_ID}"

CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
se_output=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/rc_output/$sample

for chr in {1..22} X Y M 
do
    echo "**** Processing chromosome $chr ****"
    mkdir -p $se_output/chr$chr/
done


module load conda_R/4.3.x
Rscript 01_test_bambu.R $sample

echo "**** Job ends ****"
date