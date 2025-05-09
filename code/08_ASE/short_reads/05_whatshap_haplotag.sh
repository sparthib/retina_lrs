#!/bin/bash

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH -c 20
#SBATCH --job-name=whatshap_haplotag
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/whatshap_haplotag.txt
#SBATCH -e logs/whatshap_haplotag.txt
#SBATCH --time=7-00:00:00


### whatshap phases the variants we found using GATK with the help 
### of our long-read BAMs
### This phased VCF is further used to haplotag our BAM files. 


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

ref_fa=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa
vcf_dir=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/gvcf_ref_46
genome_bam_dir=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality
whatshap_output_dir=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output
mkdir -p $whatshap_output_dir

source activate whatshap-env

whatshap stats --gtf=$whatshap_output_dir/phased.gtf $whatshap_output_dir/phased.vcf

whatshap haplotag -o $whatshap_output_dir/H9-BRN3B_hRO_2_primary_over_30_chr_only_sorted.bam \
--reference $ref_fa $whatshap_output_dir/phased.vcf $genome_bam_dir/H9-BRN3B_hRO_2_primary_over_30_chr_only_sorted.bam

whatshap haplotag -o $whatshap_output_dir/H9-BRN3B-RO_primary_over_30_chr_only_sorted.bam \
--reference $ref_fa $whatshap_output_dir/phased.vcf $genome_bam_dir/H9-BRN3B-RO_primary_over_30_chr_only_sorted.bam

whatshap haplotag -o $whatshap_output_dir/H9-CRX_hRO_2_primary_over_30_chr_only_sorted.bam \
--reference $ref_fa $whatshap_output_dir/phased.vcf $genome_bam_dir/H9-CRX_hRO_2_primary_over_30_chr_only_sorted.bam

whatshap haplotag -o $whatshap_output_dir/H9-CRX_ROs_D45_primary_over_30_chr_only_sorted.bam \
--reference $ref_fa $whatshap_output_dir/phased.vcf $genome_bam_dir/H9-CRX_ROs_D45_primary_over_30_chr_only_sorted.bam

whatshap haplotag -o $whatshap_output_dir/H9-FT_1_primary_over_30_chr_only_sorted.bam \
--reference $ref_fa $whatshap_output_dir/phased.vcf $genome_bam_dir/H9-FT_1_primary_over_30_chr_only_sorted.bam

whatshap haplotag -o $whatshap_output_dir/H9-FT_2_primary_over_30_chr_only_sorted.bam \
--reference $ref_fa $whatshap_output_dir/phased.vcf $genome_bam_dir/H9-FT_2_primary_over_30_chr_only_sorted.bam

whatshap haplotag -o $whatshap_output_dir/H9-hRGC_1_primary_over_30_chr_only_sorted.bam \
--reference $ref_fa $whatshap_output_dir/phased.vcf $genome_bam_dir/H9-hRGC_1_primary_over_30_chr_only_sorted.bam

whatshap haplotag -o $whatshap_output_dir/H9-hRGC_2_primary_over_30_chr_only_sorted.bam \
--reference $ref_fa $whatshap_output_dir/phased.vcf $genome_bam_dir/H9-hRGC_2_primary_over_30_chr_only_sorted.bam


echo "**** Job ends ****"
date +"%Y-%m-%d %T"