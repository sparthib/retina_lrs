#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH -c 20
#SBATCH --job-name=split_files
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/lrRNAseqVariantCalling/split_files.%a.txt
#SBATCH -e logs/lrRNAseqVariantCalling/split_files.%a.txt
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

ml load gatk

REFERENCE_FASTA=/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa
genome_bam=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only/${sample}_primary_over_30_chr_only_sorted.bam

output_dir=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/03a_gatk_split/

gatk SplitNCigarReads \
      -R $REFERENCE_FASTA \
      -I $genome_bam \
      -O output_dir/${sample}_split.bam

echo "**** Job ends ****"
date +"%Y-%m-%d %T"
