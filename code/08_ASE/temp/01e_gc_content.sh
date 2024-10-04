#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH -c 5
#SBATCH --job-name=gc_content
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/gc_content.%a.txt
#SBATCH -e logs/gc_content.%a.txt
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

#get chromosome of read

samtools view $longshot_output/${sample}.bam | grep "HP:i:1" | awk '{print $10}' > $longshot_output/${sample}.HP1.sequences.txt
samtools view $longshot_output/${sample}.bam | grep "HP:i:2" | awk '{print $10}' > $longshot_output/${sample}.HP2.sequences.txt
samtools view $longshot_output/${sample}.bam | grep -v "HP:i:"| awk '{print $10}' > $longshot_output/${sample}.noHP.sequences.txt

echo "**** Job ends ****"
date +"%Y-%m-%d %T"