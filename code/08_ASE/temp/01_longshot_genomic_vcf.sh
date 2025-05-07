#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH -c 20
#SBATCH --job-name=longshot
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/longshot.%a.txt
#SBATCH -e logs/longshot.%a.txt
#SBATCH --array=9-12
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

REFERENCE_FASTA=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa
genome_bam=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality/${sample}_primary_over_30_chr_only_sorted.bam  
source activate longshot

longshot_output=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/01_longshot_vcfs/${sample}
mkdir -p $longshot_output


longshot  --out_bam $longshot_output/${sample}.bam  --bam $genome_bam \
--ref $REFERENCE_FASTA --out $longshot_output/${sample}.vcf

echo "**** Job ends ****"
date +"%Y-%m-%d %T"
