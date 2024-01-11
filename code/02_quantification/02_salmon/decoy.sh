#!/bin/bash

#SBATCH -p shared
#SBATCH -p shared
#SBATCH --mem=150G
#SBATCH --cpus-per-task=20
#SBATCH --job-name=salmon
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/salmon.%a.txt
#SBATCH -e logs/salmon.%a.txt
#SBATCH --array=1-4

#try running for all chromosomes

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

source activate salmon


CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
echo "${sample}"
REFERENCE_FASTA=/dcs04/hicks/data/sparthib/GENCODE_FASTA.fa 

rm -r /dcs04/hicks/data/sparthib/casey/salmon_output
INPUT_FOLDER=/dcs04/hicks/data/sparthib/
mkdir -p $OUTPUT_FOLDER


# create decoy-aware transcriptome
cd ../salmon_files/
grep "^>" <(gunzip -c $INPUT_FOLDER/ENSEMBL_DNA_PRIMARY.fa.gz) | cut -d " " -f 1 > $INPUT_FOLDER/salmon_decoy.txt
sed -i.bak -e 's/>//g' $INPUT_FOLDER/salmon_decoy.txt

cat $INPUT_FOLDER/ENSEMBLE_CDNA.fa.gz $INPUT_FOLDER/ENSEMBL_DNA_PRIMARY.fa.gz > $INPUT_FOLDER/ensembl_gentrome_transcripts.fa.gz

echo "**** Job ends ****"
date +"%Y-%m-%d %T"