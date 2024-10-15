#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH -c 15
#SBATCH --job-name=sq_filter
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/sq_filter.txt
#SBATCH -e logs/sq_filter.txt
#SBATCH --time=7-00:00:00


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "****"

source activate /users/sparthib/.conda/envs/SQANTI3
SQANTI_DIR=~/SQANTI3-5.2.1

INPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/high_quality/all_samples/sqanti3_qc

python $SQANTI_DIR/sqanti3_filter.py ml $INPUT_DIR/all_samples_classification.txt


echo "**** Job ends ****"
    