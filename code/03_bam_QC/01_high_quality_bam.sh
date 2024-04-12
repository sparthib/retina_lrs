#!/bin/bash

#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=10G
#SBATCH -c 10
#SBATCH --job-name=highqualbam
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/primary_over_30/transcriptome/log.%a.txt
#SBATCH -e logs/primary_over_30/transcriptome/log.%a.txt
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

LOGS_FOLDER=/users/sparthib/retina_lrs/code/03_bam_QC/logs/primary_over_30/transcriptome
CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
BAM_FOLDER=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/transcriptome/GENCODE/sorted
sample=$(awk -v Index=${SLURM_ARRAY_TASK_ID} '$1==Index {print $2}' $CONFIG)
echo "$sample"



# https://wikis.utexas.edu/display/CoreNGSTools/Filtering+with+SAMTools#FilteringwithSAMTools-Filteringforhigh-qualityreads

ml load samtools

# samtools view -q 30 ${BAM_FOLDER}/${sample}.bam -o ${BAM_FOLDER}/MAPQ_filtered/${sample}.bam
mkdir -p ${BAM_FOLDER}/primary_over_30/
samtools view -q 30 -F 0x800 ${BAM_FOLDER}/${sample}_sorted.bam -o ${BAM_FOLDER}/MAPQ_FILTERED/${sample}.bam
samtools sort ${BAM_FOLDER}/MAPQ_FILTERED/${sample}.bam -o ${BAM_FOLDER}/MAPQ_FILTERED/${sample}_sorted.bam
samtools index ${BAM_FOLDER}/MAPQ_FILTERED/${sample}_sorted.bam ${BAM_FOLDER}/MAPQ_FILTERED/${sample}_sorted.bam.bai


# echo "finished indexing bam"
index stats ${sample}_primary_only_index_stats.txt
samtools idxstats ${BAM_FOLDER}/MAPQ_FILTERED/${sample}_sorted.bam > ${LOGS_FOLDER}/mapq_filtered_${sample}_index_stats.txt

#bam stats ${sample}_bam.stats for plotting bam stats using plot-bamstats command in samtools
samtools stats ${BAM_FOLDER}/MAPQ_FILTERED/${sample}_sorted.bam > ${LOGS_FOLDER}/mapq_filtered_${sample}_bam.stats

echo "finished computing stats for plotting"

echo "flagstat" > ${LOGS_FOLDER}/${sample}_bam_flagstat.txt
samtools flagstat ${BAM_FOLDER}/MAPQ_FILTERED/${sample}_sorted.bam >> ${LOGS_FOLDER}/mapq_filtered_${sample}_bam_flagstat.txt

echo "finished computing flagstats: contains percentage mapped reads"

echo "depth of coverage" > ${LOGS_FOLDER}/${sample}_depth_stats.txt
samtools depth -a ${BAM_FOLDER}/MAPQ_FILTERED/${sample}_sorted.bam  | awk '{c++;s+=$3}END{print s/c}' >> ${LOGS_FOLDER}/mapq_filtered_${sample}_depth_stats.txt

echo "breadth of coverage" >> ${sample}_depth_stats.txt
samtools depth -a ${BAM_FOLDER}/MAPQ_FILTERED/${sample}_sorted.bam  | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' >> ${LOGS_FOLDER}/mapq_filtered_${sample}_depth_stats.txt

echo "raw depth output" >> ${sample}_depth_stats.txt
samtools depth -a ${BAM_FOLDER}/MAPQ_FILTERED/${sample}_sorted.bam  >> ${LOGS_FOLDER}/mapq_filtered_${sample}_depth_stats.txt

echo "finished computing depth stats"


echo "**** Job ends ****"
date +"%Y-%m-%d %T"
