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
#SBATCH --time=7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"



ref_fa=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa
INPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/filtered_bams_ref_46
OUTPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/gvcf_ref_46

mkdir -p $OUTPUT_DIR
# ml load samtools 
# samtools faidx $ref_fa

ml load gatk

# echo "🔹 Step 1: Haplotype calling to produce gVCF"
# gatk HaplotypeCaller \
#     -R "$ref_fa" \
#     -I "${INPUT_DIR}/all_samples_filtered.bam" \
#     -O "${OUTPUT_DIR}/all_samples.g.vcf.gz" \
#     -ERC GVCF
#     
# Uncomment the following lines if you want to run HaplotypeCaller for each sample individually
# echo "🔹 Step 2: Combine gVCFs from multiple samples"
# gatk CombineGVCFs \
#   -R "$ref_fa" \
#   $(for sample in "${samples[@]}"; do echo -n "-V ${OUTPUT_DIR}/${sample}.g.vcf.gz "; done) \
#   -O "${OUTPUT_DIR}/combined.g.vcf.gz"
# 
# echo "🔹 Step 3: Genotype the combined gVCF"
# gatk GenotypeGVCFs \
#   -R "$ref_fa" \
#   -V "${OUTPUT_DIR}/all_samples.g.vcf.gz" \
#   -O "${OUTPUT_DIR}/all_samples.vcf.gz"
# 
# 
# echo "🔹 Step 4a: Split VCF into SNPs and INDELs"
# 
# gatk SelectVariants \
#   -R "$ref_fa" \
#   -V "${OUTPUT_DIR}/all_samples.vcf.gz" \
#   --select-type-to-include SNP \
#   -O "${OUTPUT_DIR}/all_samples_SNPs.vcf.gz"
# 
# gatk SelectVariants \
#   -R "$ref_fa" \
#   -V "${OUTPUT_DIR}/all_samples.vcf.gz" \
#   --select-type-to-include INDEL \
#   -O "${OUTPUT_DIR}/all_samples_INDELs.vcf.gz"
#  
  echo "🔹 Step 4b: Filter SNPs and INDELs" 
 
  gatk VariantFiltration \
  -R "$ref_fa" \
  -V "${OUTPUT_DIR}/all_samples_SNPs.vcf.gz" \
   -filter "QD < 2.0" --filter-name "QD2" \
   -filter "QUAL < 30.0" --filter-name "QUAL30" \
   -filter "SOR > 3.0" --filter-name "SOR3" \
   -filter "FS > 60.0" --filter-name "FS60" \
   -filter "MQ < 40.0" --filter-name "MQ40" \
   #may not work for single samples
   # -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
   # -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  -O "${OUTPUT_DIR}/all_samples_SNPs_filtered.vcf.gz"
  
gatk VariantFiltration \
  -R "$ref_fa" \
  -V "${OUTPUT_DIR}/all_samples_INDELs.vcf.gz" \
      -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    # doens't work for single samples
    # -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O "${OUTPUT_DIR}/all_samples_INDELs_filtered.vcf.gz"
  
  
  echo "🔹 Step 4c: only keep the PASS variants from the filtered files"
gatk SelectVariants \
  -R "$ref_fa" \
  -V "${OUTPUT_DIR}/all_samples_SNPs_filtered.vcf.gz" \
  --exclude-filtered \
  -O "${OUTPUT_DIR}/all_samples_SNPs_filtered_PASS.vcf.gz"  
  
gatk SelectVariants \
  -R "$ref_fa" \
  -V "${OUTPUT_DIR}/all_samples_INDELs_filtered.vcf.gz" \
  --exclude-filtered \
  -O "${OUTPUT_DIR}/all_samples_INDELs_filtered_PASS.vcf.gz"
  
echo "🔹 Step 5: Merge SNPs and INDELs"
gatk MergeVcfs \
  -I "${OUTPUT_DIR}/all_samples_SNPs_filtered_PASS.vcf.gz" \
  -I "${OUTPUT_DIR}/all_samples_INDELs_filtered_PASS.vcf.gz" \
  -O "${OUTPUT_DIR}/all_samples_variants.vcf.gz"

echo "**** Job ends ****"
date +"%Y-%m-%d %T"




