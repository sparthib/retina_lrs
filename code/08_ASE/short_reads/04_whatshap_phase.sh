#!/bin/bash

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH -c 20
#SBATCH --job-name=whatshap
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/whatshap.txt
#SBATCH -e logs/whatshap.txt
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
phased_vcf_dir=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output
mkdir -p $phased_vcf_dir

source activate whatshap-env
# input vcf where vcfs were merged from individual samples= merged_variants.vcf.gz 
#corresponding output was phased.vcf instead of all_samples_phased.vcf
##--sample=SRR1091091 this tag comes right after -o if VCF has multiple samples
whatshap phase -o $phased_vcf_dir/all_samples_phased.vcf --reference=$ref_fa --ignore-read-groups $vcf_dir/all_samples_variants.vcf.gz \
$genome_bam_dir/H9-BRN3B_hRO_2_primary_over_30_chr_only_sorted.bam $genome_bam_dir/H9-BRN3B-RO_primary_over_30_chr_only_sorted.bam \
$genome_bam_dir/H9-CRX_hRO_2_primary_over_30_chr_only_sorted.bam $genome_bam_dir/H9-CRX_ROs_D45_primary_over_30_chr_only_sorted.bam \
$genome_bam_dir/H9-FT_1_primary_over_30_chr_only_sorted.bam $genome_bam_dir/H9-FT_2_primary_over_30_chr_only_sorted.bam \
$genome_bam_dir/H9-hRGC_1_primary_over_30_chr_only_sorted.bam $genome_bam_dir/H9-hRGC_2_primary_over_30_chr_only_sorted.bam 

echo "**** Job ends ****"
date +"%Y-%m-%d %T"
