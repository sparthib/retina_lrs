#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH -c 5
#SBATCH --job-name=compute_HP
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/compute_HP.%a.txt
#SBATCH -e logs/compute_HP.%a.txt
#SBATCH --array=1-15
#SBATCH --time=7-00:00:00


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "****"


CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=${SLURM_ARRAY_TASK_ID} '$1==Index {print $2}' $CONFIG)
echo $sample

ml load samtools

longshot_output=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/01_longshot_vcfs/${sample}

total_reads=$(samtools view -c $longshot_output/${sample}.bam)
hp1_reads=$(samtools view $longshot_output/${sample}.bam | grep "HP:i:1" | wc -l)
hp2_reads=$(samtools view $longshot_output/${sample}.bam | grep "HP:i:2" | wc -l)
unphased_reads=$(samtools view $longshot_output/${sample}.bam | grep -v "HP:i:" | wc -l)

echo "Total number of reads = $total_reads"
echo "Number of HP1 reads = $hp1_reads"
echo "Number of HP2 reads = $hp2_reads"
echo "Number of unphased reads = $unphased_reads"

hp1_percent=$(echo "scale=2; ($hp1_reads / $total_reads) * 100" | bc)
hp2_percent=$(echo "scale=2; ($hp2_reads / $total_reads) * 100" | bc)
unphased_percent=$(echo "scale=2; ($unphased_reads / $total_reads) * 100" | bc)

echo "Percentage of HP1 reads = $hp1_percent"
echo "Percentage of HP2 reads = $hp2_percent"
echo "Percentage of unphased reads = $unphased_percent"


echo "**** Job ends ****"
date +"%Y-%m-%d %T"
