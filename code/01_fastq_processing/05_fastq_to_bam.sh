#!/bin/bash

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --mem=250G
#SBATCH -c 25
#SBATCH --job-name=fastq2bam
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/bam_stats/fastq2bam.%a.txt
#SBATCH -e logs/bam_stats/fastq2bam.%a.txt
#SBATCH --array=1-4

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

LOGS_FOLDER=/users/sparthib/retina_lrs/code/01_fastq_processing/logs/bam_stats
CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
INPUT_FOLDER=/dcs04/hicks/data/sparthib/casey/fastqs_post_qc
REFERENCE_FASTA=/dcs04/hicks/data/sparthib/GENCODE_FASTA.fa
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
SAM_FOLDER=/dcs04/hicks/data/sparthib/casey/sams
BAM_FOLDER=/dcs04/hicks/data/sparthib/casey/bams

mkdir -p $SAM_FOLDER
mkdir -p $BAM_FOLDER

cd ~/minimap2

#remove older sam version 
rm ${SAM_FOLDER}/${sample}.sam
./minimap2 -ax map-ont -N 50 --secondary=no -t ${SLURM_CPUS_PER_TASK} $REFERENCE_FASTA ${INPUT_FOLDER}/${sample}.fastq.gz > ${SAM_FOLDER}/${sample}.sam


#remove older bam version
rm ${BAM_FOLDER}/${sample}_sorted.bam
rm ${BAM_FOLDER}/${sample}_sorted.bam.bai 
ml load samtools

samtools view -bS ${SAM_FOLDER}/${sample}.sam -o ${BAM_FOLDER}/${sample}.bam
samtools sort ${BAM_FOLDER}/${sample}.bam -o ${BAM_FOLDER}/${sample}_sorted.bam
rm ${BAM_FOLDER}/${sample}.bam
samtools index ${BAM_FOLDER}/${sample}_sorted.bam ${BAM_FOLDER}/${sample}_sorted.bam.bai 

echo "finished indexing bam"
#index stats ${sample}_index_stats.txt
samtools idxstats ${BAM_FOLDER}/${sample}_sorted.bam > ${LOGS_FOLDER}/${sample}_index_stats.txt

#bam stats ${sample}_bam.stats for plotting bam stats using plot-bamstats command in samtools
samtools stats ${BAM_FOLDER}/${sample}_sorted.bam > ${LOGS_FOLDER}/${sample}_bam.stats

echo "finished computing stats for plotting"

echo "flagstat" > ${LOGS_FOLDER}/${sample}_bam_flagstat.txt
samtools flagstat ${BAM_FOLDER}/${sample}_sorted.bam >> ${LOGS_FOLDER}/${sample}_bam_flagstat.txt

echo "finished computing flagstats: contains percentage mapped reads"

echo "depth of coverage" > ${LOGS_FOLDER}/${sample}_depth_stats.txt
samtools depth -a ${BAM_FOLDER}/${sample}_sorted.bam  | awk '{c++;s+=$3}END{print s/c}' >> ${LOGS_FOLDER}/${sample}_depth_stats.txt

echo "breadth of coverage" >> ${sample}_depth_stats.txt
samtools depth -a ${BAM_FOLDER}/${sample}_sorted.bam  | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' >> ${LOGS_FOLDER}/${sample}_depth_stats.txt

echo "raw depth output" >> ${sample}_depth_stats.txt
samtools depth -a ${BAM_FOLDER}  >> ${sample}_depth_stats.txt

echo "finished computing depth stats"

echo "**** Job ends ****"
date +"%Y-%m-%d %T"