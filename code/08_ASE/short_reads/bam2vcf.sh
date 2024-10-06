#!/bin/bash

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH -c 20
#SBATCH --job-name=bam2vcf.txt
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/bam2vcf.txt
#SBATCH -e logs/bam2vcf.txt
#SBATCH --array=1-15
#SBATCH --time=7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"


ml load bcftools/1.18

# https://samtools.github.io/bcftools/howtos/variant-calling.html

ref_fa=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa
bam_files=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/filtered_bams/bam_files.txt
output_dir=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/vcf

bcftools mpileup -Ou --threads $SLURM_CPUS_PER_TASK \
-f /dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa \
-b /dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/filtered_bams/bam_files.txt | bcftools call -mv -Ob > $output_dir/multi_sample.vcf







