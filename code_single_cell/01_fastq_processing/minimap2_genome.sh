#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH -c 10
#SBATCH --job-name=fastq2bam
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/minimap2/genome/genome.%a.txt
#SBATCH -e logs/minimap2/genome/genome.%a.txt
#SBATCH --array=2

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

LOGS_FOLDER=/users/sparthib/retina_lrs/code_single_cell/01_fastq_processing/logs/minimap2/genome
CONFIG=/users/sparthib/retina_lrs/raw_data/single_cell.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
path=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $3}' $CONFIG)
seq_sum=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $4}' $CONFIG)
echo $sample

stranded_input_fastq=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/restrander_output/${sample}.fastq.gz
unknown_input_fastq=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/restrander_output/${sample}-unknowns.fastq.gz


LOGS_FOLDER=/users/sparthib/retina_lrs/code/01_fastq_processing/logs/bam_stats_genome_gencode_splice
REFERENCE_FASTA=/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa.gz

SAM_FOLDER=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/04_minimap2_output/genome/sams
BAM_FOLDER=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/04_minimap2_output/genome/bams

mkdir -p $SAM_FOLDER
mkdir -p $BAM_FOLDER
# 
cd ~/minimap2

#remove older sam version 
rm ${SAM_FOLDER}/${sample}_stranded.sam
./minimap2 -ax splice -y --secondary=no -t ${SLURM_CPUS_PER_TASK} $REFERENCE_FASTA $stranded_input_fastq > ${SAM_FOLDER}/${sample}_stranded.sam


# remove older bam version
rm ${BAM_FOLDER}/${sample}_stranded_sorted.bam
rm ${BAM_FOLDER}/${sample}_stranded_sorted.bam.bai

ml load samtools

samtools view -bS ${SAM_FOLDER}/${sample}_stranded.sam -o ${BAM_FOLDER}/${sample}_stranded.bam
samtools sort ${BAM_FOLDER}/${sample}_stranded.bam -o ${BAM_FOLDER}/${sample}_stranded_sorted.bam
# rm ${BAM_FOLDER}/${sample}.bam
samtools index ${BAM_FOLDER}/${sample}_stranded_sorted.bam ${BAM_FOLDER}/${sample}_stranded_sorted.bam.bai

echo "finished indexing stranded bam"
index stats ${sample}_stranded_index_stats.txt
samtools idxstats ${BAM_FOLDER}/${sample}_stranded_sorted.bam > ${LOGS_FOLDER}/${sample}_stranded_index_stats.txt


echo "finished computing stats for plotting"

echo "flagstat" > ${LOGS_FOLDER}/${sample}_stranded_bam_flagstat.txt
samtools flagstat ${BAM_FOLDER}/${sample}_stranded_sorted.bam >> ${LOGS_FOLDER}/${sample}_stranded_bam_flagstat.txt

echo "finished computing flagstats: contains percentage mapped reads"

################

#remove older sam version 
rm ${SAM_FOLDER}/${sample}_unknown.sam
./minimap2 -ax splice -y --secondary=no -t ${SLURM_CPUS_PER_TASK} $REFERENCE_FASTA $unknown_input_fastq > ${SAM_FOLDER}/${sample}_unknown.sam


# remove older bam version
rm ${BAM_FOLDER}/${sample}_unknown_sorted.bam
rm ${BAM_FOLDER}/${sample}_unknown_sorted.bam.bai

ml load samtools

samtools view -bS ${SAM_FOLDER}/${sample}_unknown.sam -o ${BAM_FOLDER}/${sample}_unknown.bam
samtools sort ${BAM_FOLDER}/${sample}_unknown.bam -o ${BAM_FOLDER}/${sample}_unknown_sorted.bam
# rm ${BAM_FOLDER}/${sample}.bam
samtools index ${BAM_FOLDER}/${sample}_unknown_sorted.bam ${BAM_FOLDER}/${sample}_unknown_sorted.bam.bai

echo "finished indexing uknown bam"
index stats ${sample}_unknown_index_stats.txt
samtools idxstats ${BAM_FOLDER}/${sample}_unknown_sorted.bam > ${LOGS_FOLDER}/${sample}_unknown_index_stats.txt


echo "finished computing stats for plotting"

echo "flagstat" > ${LOGS_FOLDER}/${sample}_unknown_bam_flagstat.txt
samtools flagstat ${BAM_FOLDER}/${sample}_unknown_sorted.bam >> ${LOGS_FOLDER}/${sample}_unknown_bam_flagstat.txt

echo "finished computing flagstats: contains percentage mapped reads"

echo "**** Job ends ****"
date +"%Y-%m-%d %T"