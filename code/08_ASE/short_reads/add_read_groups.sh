#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH -c 20
#SBATCH --job-name=bowtie
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/bowtie.txt
#SBATCH -e logs/bowtie.txt
#SBATCH -t 7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "****"


INPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/filtered_bams

samples=(SRR1091088 SRR1091091 SRR1091092)