#!/bin/bash

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --mem=300G
#SBATCH -c 10
#SBATCH --job-name=scnanoseq
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/scnanoseq.txt
#SBATCH -e logs/scnanoseq.txt
#SBATCH --time=7-00:00:00


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"


ml load nextflow
ml load singularity

nextflow run nf-core/scnanoseq \
   -profile singularity \
   --input samplesheet.csv \
   --outdir <OUTDIR>

echo "**** Job ends ****"
date +"%Y-%m-%d %T"