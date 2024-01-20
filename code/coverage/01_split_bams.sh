#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH -c 20
#SBATCH --job-name=split_bam
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/bam_stats_genome_gencode/split_bams/split.%a.txt
#SBATCH -e logs/bam_stats_genome_gencode/split.%a.txt
#SBATCH --array=1-8

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

LOGS_FOLDER=/users/sparthib/retina_lrs/code/01_fastq_processing/logs/bam_stats_genome_gencode/split_bams
mkdir -p $LOGS_FOLDER
CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
INPUT_FOLDER=/dcs04/hicks/data/sparthib/retina_lrs/03_processed_fastqs
REFERENCE_FASTA=/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa.gz
sample=$(awk -v Index=${SLURM_ARRAY_TASK_ID} '$1==Index {print $2}' $CONFIG)
echo "$sample"
BAM_FOLDER=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE
INPUT_BAM=$BAM_FOLDER/${sample}_sorted.bam
OUTPUT_BAMS=$BAM_FOLDER/${sample}
mkdir -p OUTPUT_BAMS


ml load samtools

declare -a chr_List=(
                 "chr1" 
                 "chr2" 
                 "chr3"
                 "chr4"
                 "chr5"
                 "chr6"
                 "chr7"
                 "chr8"
                 "chr9"
                 "chr10"
                 "chr11"
                 "chr12"
                 "chr13"
                 "chr14"
                 "chr15"
                 "chr16"
                 "chr17"
                 "chr18"
                 "chr19"
                 "chr20"
                 "chr21"
                 "chr22"
                 "chrX"
                 "chrY"
                 "chrM"
                )
for chr_num in ${chr_List[@]}
  do
    samtools view $INPUT_BAM $chr_num -b > $OUTPUT_BAMS/$sample/${chr_num}.bam
    
    echo "flagstat" > ${LOGS_FOLDER}/${sample}_${chr_num}_bam_flagstat.txt
    samtools flagstat $OUTPUT_BAMS/$sample/${chr_num}.bam >> ${LOGS_FOLDER}/${sample}_${chr_num}_bam_flagstat.txt

    echo "finished computing flagstats: contains percentage mapped reads"

    echo "depth of coverage" > ${LOGS_FOLDER}/${sample}_${chr_num}_depth_stats.txt
    samtools depth -a $OUTPUT_BAMS/$sample/${chr_num}.bam  | awk '{c++;s+=$3}END{print s/c}' >> ${LOGS_FOLDER}/${sample}_${chr_num}_depth_stats.txt

    echo "breadth of coverage" >> ${LOGS_FOLDER}/${sample}_${chr_num}_depth_stats.txt
    samtools depth -a $OUTPUT_BAMS/$sample/${chr_num}.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' >> ${LOGS_FOLDER}/${sample}_${chr_num}_depth_stats.txt

    echo "raw depth output" >> ${LOGS_FOLDER}/${sample}_${chr_num}_depth_stats.txt
    samtools depth -a $OUTPUT_BAMS/$sample/${chr_num}.bam  >> ${LOGS_FOLDER}/${sample}_${chr_num}_depth_stats.txt

    echo "finished computing depth stats"
  done


echo "**** Job ends ****"
date +"%Y-%m-%d %T"
    
    
  
  
