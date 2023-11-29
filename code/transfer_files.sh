#!/bin/bash

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --mem=2G
#SBATCH -c 1
#SBATCH --job-name=data_transfer
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/data_transfer.txt
#SBATCH -e logs/data_transfer.txt


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"


ascp cidr-zack@162.129.245.24:/193679_Zack /dcs04/hicks/data/sparthib/casey_round_2/input_fastqs


echo "****Job Ends****"
date +"%Y-%m-%d %T"