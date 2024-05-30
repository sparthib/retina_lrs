#!/bin/bash

#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=10G
#SBATCH -c 10
#SBATCH --job-name=highqualbam
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/primary_over_30/genome_splice/log.%a.txt
#SBATCH -e logs/primary_over_30/genome_splice/log.%a.txt
#SBATCH --array=1-15
#SBATCH -t 4-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

LOGS_FOLDER=/users/sparthib/retina_lrs/code/03_bam_QC/logs/primary_over_30/genome_splice
CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=${SLURM_ARRAY_TASK_ID} '$1==Index {print $2}' $CONFIG)
echo "$sample"
input_dir="/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice"
output_dir="/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only"
mkdir -p $output_dir
bam=$input_dir/${sample}_sorted.bam

ml load samtools 

samtools view -b -q 30 -F 0x800 -h $bam > $output_dir/${sample}_primary_over_30.bam
samtools sort $output_dir/${sample}_primary_over_30.bam -o $output_dir/${sample}_primary_over_30_sorted.bam
samtools index $output_dir/${sample}_primary_over_30_sorted.bam $output_dir/${sample}_primary_over_30_sorted.bam.bai

# # echo "finished indexing bam"
# samtools idxstats $primary_over_30/${sample}_sorted.bam > ${LOGS_FOLDER}/primary_over_30_${sample}_index_stats.txt
# 
# echo "finished computing stats for plotting"
# 
# echo "flagstat" > ${LOGS_FOLDER}/${sample}_bam_flagstat.txt
# samtools flagstat $primary_over_30/${sample}_sorted.bam >> ${LOGS_FOLDER}/primary_over_30_${sample}_bam_flagstat.txt

echo "**** Job ends ****"
date +"%Y-%m-%d %T"

