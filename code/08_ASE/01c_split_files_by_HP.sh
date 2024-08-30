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

#get all reads that have HP:i:1 tag
samtools view -h $longshot_output/${sample}.bam | grep "HP:i:1" > /$longshot_output/${sample}.HP1.bam

#get all reads that have HP:i:2 tag
samtools view -h $longshot_output/${sample}.bam | grep "HP:i:2" > /$longshot_output/${sample}.HP1.bam

#get all reads that don't have the HP tag 
samtools view -h $longshot_output/${sample}.bam | grep -v "HP:i:" > /$longshot_output/${sample}.noHP.bam

echo "**** Job ends ****"
date +"%Y-%m-%d %T"
