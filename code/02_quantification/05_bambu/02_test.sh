#!/bin/bash

#SBATCH -p shared
#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --job-name=your_job
#SBATCH --mail-user=email_id
#SBATCH --mail-type=ALL
#SBATCH -o logs/template_log.txt
#SBATCH -e logs/template_log.txt


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"


module load conda_R/4.3.x

echo "**** Job ends ****"
date +"%Y-%m-%d %T"