#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH -c 2
#SBATCH --job-name=fasta_eda
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/fasta_eda.txt
#SBATCH -e logs/fasta_eda.txt
#SBATCH -t 7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "****"

ml load python/3.10.13 

# python3 extract_transcripts.py

python3 compare_genome.py

echo "**** Job ends ****"
date +"%Y-%m-%d %T"
