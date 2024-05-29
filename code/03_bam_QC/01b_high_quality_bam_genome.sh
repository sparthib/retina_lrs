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
bam_dir=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice
bam=$bam_dir/${sample}_sorted.bam
primary_over_30=$bam_dir/primary_over_30_chr_only
mkdir -p $primary_over_30

# for chr in chr{1..22} chrX chrY chrM
# do
#     echo "Processing $chr"
#     samtools view -b $bam $chr > ${output}/${sample}_${chr}.bam
#     samtools index ${output}/${sample}_${chr}.bam ${output}/${sample}_${chr}.bam.bai
# done
ml load samtools 

chr_only_bam=$bam_dir/${sample}_chromosome_level/$sample.bam
samtools view -q 30 -F 0x800 $chr_only_bam > $primary_over_30/${sample}_primary_over_30.bam
samtools sort $primary_over_30/${sample}_primary_over_30.bam -o $primary_over_30/${sample}_sorted.bam
samtools index $primary_over_30/${sample}_primary_over_30_sorted.bam $primary_over_30/${sample}_sorted.bam.bai

# echo "finished indexing bam"
samtools idxstats $primary_over_30/${sample}_sorted.bam > ${LOGS_FOLDER}/primary_over_30_${sample}_index_stats.txt

echo "finished computing stats for plotting"

echo "flagstat" > ${LOGS_FOLDER}/${sample}_bam_flagstat.txt
samtools flagstat $primary_over_30/${sample}_sorted.bam >> ${LOGS_FOLDER}/primary_over_30_${sample}_bam_flagstat.txt

echo "**** Job ends ****"
date +"%Y-%m-%d %T"

