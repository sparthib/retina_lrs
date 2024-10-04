#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=30G
#SBATCH -c 5
#SBATCH --job-name=filter_bam
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/filter_bam.txt
#SBATCH -e logs/filter_bam.txt
#SBATCH -t 7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "****"


INPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/sams
OUTPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/filtered_bams

ml load samtools

samples=(SRR1091088 SRR1091091 SRR1091092)

for sample in ${samples[@]}
do
    samtools view -bS -h -F 0x904 -q 20 $INPUT_DIR/${sample}.sam > $OUTPUT_DIR/${sample}.bam
    samtools sort $OUTPUT_DIR/${sample}.bam -o $OUTPUT_DIR/${sample}.sorted.bam
    samtools index $OUTPUT_DIR/${sample}.sorted.bam $OUTPUT_DIR/${sample}.sorted.bam.bai
    rm $OUTPUT_DIR/${sample}.bam
done

echo "**** Job ends ****"
date +"%Y-%m-%d %T"










samtools view -h -F 256 -q 20 file.bam