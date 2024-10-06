#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=150G
#SBATCH -c 20
#SBATCH --job-name=multi_exon_pcg
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/multi_exon_pcg.txt
#SBATCH -e logs/multi_exon_pcg.txt
#SBATCH --time=7-00:00:00


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

ml load python/3.10.13
python3 01_multi_exon_pcg.py 

echo "**** Job ends ****"
date +"%Y-%m-%d %T"