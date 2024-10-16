#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH -c 15
#SBATCH --job-name=sq_rescue 
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/sq_rescue.txt
#SBATCH -e logs/sq_rescue.txt
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
REFERENCE_GENOME_FASTA=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa
ISOFORMS=$INPUT_DIR/all_samples_corrected.fasta
GTF=$INPUT_DIR/all_samples.filtered.gtf
REFGTF=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_assembly.gtf
REFCLASSIF=$INPUT_DIR/all_samples_classification.txt
FITLER_CLASSIFICATION=$INPUT_DIR/all_samples_MLresult_classification.txt

python $SQANTI_DIR/sqanti3_rescue.py ml --isoforms $ISOFORMS --gtf $GTF \
-g $REFGTF -f $REFERENCE_GENOME_FASTA -k $REFCLASSIF \
-d $INPUT_DIR -o post_rescue $FITLER_CLASSIFICATION

echo "**** Job ends ****"
date +"%Y-%m-%d %T"