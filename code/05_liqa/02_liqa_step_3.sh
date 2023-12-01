#!/bin/bash

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH -c 10
#SBATCH --job-name=liqa_diff_exp
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/liqa_diff_exp.txt
#SBATCH -e logs/liqa_diff_exp.txt

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

source activate liqa


CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
echo "${sample}"
STEP_3_OUTPUT=/dcs04/hicks/data/sparthib/casey/liqa/step3_outputs
mkdir $STEP_3_OUTPUT


##STEP 3 - Detecting differential splicing gene/isoform between conditions

liqa -task diff -condition_1 ROs_list -condition_2 RGCs_list -out $STEP_3_OUTPUT



conda deactivate

echo "**** Job ends ****"
date +"%Y-%m-%d %T"