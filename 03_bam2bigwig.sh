#!/bin/bash

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH -c 5
#SBATCH --job-name=bam2bigwig
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/bam2bigwig.log
#SBATCH -e logs/bam2bigwig.log


# /users/sparthib/.conda/envs/deeptools

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"


source activate deeptools 

BAM_FILE=/dcs04/hicks/data/sparthib/casey/bams/hRGC_chr_21_sorted.bam
TWOBIT_FILE=/dcs04/hicks/data/sparthib/GENCODE_FASTA.2bit
BIGWIGS=/dcs04/hicks/data/sparthib/casey/bigwigs/

mkdir -p $BIGWIGS

bamCoverage -b $BAM_FILE -o $BIGWIGS/hRGC_chr_21.bw

conda deactivate 

echo "**** Job ends ****"
date +"%Y-%m-%d %T"