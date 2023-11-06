#!/bin/bash

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH -c 5
#SBATCH --job-name=computeGCbias
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/computeGC_test.log
#SBATCH -e logs/computeGC_test.log


# /users/sparthib/.conda/envs/deeptools

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"


source activate deeptools 

BAM_FOLDER=/dcs04/hicks/data/sparthib/casey/bams/hRGC_chr_21.bam

computeGCBias -b $BAM_FOLDER --effectiveGenomeSize 2695000000 --genome genome.2bit -l 200 -freq freq_test.txt --region 21 --biasPlot test_gc.png


conda deactivate 

echo "**** Job ends ****"
date +"%Y-%m-%d %T"