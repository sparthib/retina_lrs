#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=150G
#SBATCH -c 20
#SBATCH --job-name=fastq2bam
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/flair_fusion.%a.txt
#SBATCH -e logs/flair_fusion.%a.txt
#SBATCH --array=1

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

LOGS_FOLDER=/users/sparthib/retina_lrs/code/02_quantification/04_flair_fusion/logs
CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
INPUT_FOLDER=/dcs04/hicks/data/sparthib/retina_lrs/03_processed_fastqs
sample=$(awk -v Index=${SLURM_ARRAY_TASK_ID} '$1==Index {print $2}' $CONFIG)

### FLAIR-Fusion 
REFERENCE_GTF=/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf
REFERENCE_FASTA=/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa
# python ./makeShortAnno.py $REFERENCE_GTF

SHORT_GTF= /users/sparthib/FLAIR-fusion/gencode.v44.chr_patch_hapl_scaff.annotation-short.gtf
FLAIRPY_PATH=/users/sparthib/flair-2.0.0/flair.py
OUTPUT_PATH=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair_fusion/$sample
mkdir $OUTPUT_PATH
#test it on just one fastq file for now 
python3 ./19-03-2021-fasta-to-fusions-pipe.py -r ${INPUT_FOLDER}/${sample}.fastq.gz \
-f $FLAIRPY_PATH -g $REFERENCE_FASTA -t $REFERENCE_GTF \
-a $SHORT_GTF -o $OUTPUT_PATH

echo "**** Job ends ****"
date +"%Y-%m-%d %T"