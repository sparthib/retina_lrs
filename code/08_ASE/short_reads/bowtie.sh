#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH -c 20
#SBATCH --job-name=bowtie
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/bowtie.txt
#SBATCH -e logs/bowtie.txt
#SBATCH -t 7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "****"

INPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/fastq
OUTPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/sams


ml load bowtie/2.5.1
BT_INDEX_PATH=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/GRCh38_noalt_as/GRCh38_noalt_as

# Aligning paired reads
samples=(SRR1091088 SRR1091091 SRR1091092)

for sample in ${samples[@]}
do
    bowtie2 -x $BT_INDEX_PATH -p ${SLURM_CPUS_PER_TASK} \
    -1 $INPUT_DIR/${sample}_1.fastq -2 $INPUT_DIR/${sample}_2.fastq -S $OUTPUT_DIR/${sample}.sam \
    --rg-id ${sample} --rg "SM:${sample}" --rg "PL:ILLUMINA" 
done

echo "**** Job ends ****"
date +"%Y-%m-%d %T"

