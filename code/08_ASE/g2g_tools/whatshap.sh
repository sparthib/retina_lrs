#!/bin/bash

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH -c 20
#SBATCH --job-name=whatshap.txt
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/whatshap.txt
#SBATCH -e logs/whatshap.txt
#SBATCH --time=7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

ref_fa=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa
vcf_dir=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/vcf
source activate whatshap-env
genome_bam_dir=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly
transcriptome_bam_dir=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/transcriptome/ver_46

ml load bcftools 
# bcftools view -g ^miss -i 'GT[*]~"|"' $vcf_dir/SRR1091088.genotyped.vcf -o $vcf_dir/SRR1091088.phased.vcf
ml load  htslib
tabix -p vcf $vcf_dir/SRR1091088.phased.vcf

whatshap haplotag -o $vcf_dir/H9-hRGC_1_sorted_haplotagged.bam --reference $ref_fa \
$vcf_dir/SRR1091088.phased.vcf $genome_bam_dir/H9-hRGC_1_sorted.bam \
--ignore-read-groups



echo "**** Job ends ****"
date +"%Y-%m-%d %T"



