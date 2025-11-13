#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=60G
#SBATCH -c 10
#SBATCH --job-name=filter_bam
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/filter_bam.txt
#SBATCH -e logs/filter_bam.txt
#SBATCH -t 7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "****"

# Converts SAM → BAM
# 
# Sorts and marks duplicates
# 
# Adds read groups
# 
# Performs BQSR
# 
# Filters BAMs to remove unmapped, secondary, supplementary reads, and optionally low MAPQ
# 
# Runs GATK HaplotypeCaller in GVCF mode
# 
# Combines GVCFs
# 
# Runs joint genotyping to produce a final cohort VCF

ENV_FILE="../../../.env"
if [ -f $ENV_FILE ]; then
    set -a
    source $ENV_FILE
    set +a
fi

known_sites=$retina_lrs_dir/09_ASE/H9_DNA_Seq_data/known_sites
DBSNP=${known_sites}/Homo_sapiens_assembly38.dbsnp138.vcf
MANDG_INDELS=${known_sites}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
SNPs_1000G=${known_sites}/1000G_phase1.snps.high_confidence.hg38.vcf.gz
ref_fa=$references_dirgenome/GENCODE/primary_assembly/release_46_primary_genome.fa
INPUT_DIR=$retina_lrs_dir/09_ASE/H9_DNA_Seq_data/sams_ref_46
OUTPUT_DIR=$retina_lrs_dir/09_ASE/H9_DNA_Seq_data/filtered_bams_ref_46
mkdir -p $OUTPUT_DIR

ml load gatk
ml load samtools

## first mark duplicates, then filter by MAPQ and uniquely mapped reads 

# samples=(SRR1091088 SRR1091091 SRR1091092)
# 
# mkdir -p "$OUTPUT_DIR"
# 
# for sample in "${samples[@]}"
# do
#   echo "🔄 Processing sample: $sample"
# 
# 
#   echo "🔹 Step 1a: Change SAM to BAM"
#   samtools view -@ 8 -bS ${INPUT_DIR}/${sample}.sam > ${OUTPUT_DIR}/${sample}.bam
# 
#   echo "🔹 Step 1b: Sorting BAM"
#   samtools sort -@ 8 -o ${OUTPUT_DIR}/${sample}_sorted.bam ${OUTPUT_DIR}/${sample}.bam
# 
# done 

# echo "Step 1c: Merging BAMs"
# samtools merge -@ 8 ${OUTPUT_DIR}/all_samples_merged.bam ${OUTPUT_DIR}/*_sorted.bam
#   

samtools sort -@ 8 -o "${OUTPUT_DIR}/all_samples_merged_sorted.bam" "${INPUT_DIR}"/all_samples_merged.bam

echo "🔹 Step 2: Marking duplicates"
gatk MarkDuplicates \
  I="${OUTPUT_DIR}/all_samples_merged_sorted.bam" \
  O="${OUTPUT_DIR}/all_samples_dedup.bam" \
  M="${OUTPUT_DIR}/all_samples_metrics.txt" \
  CREATE_INDEX=true 


echo "🔹 Step 3: Adding read groups to merged BAM"
gatk AddOrReplaceReadGroups \
  I="${OUTPUT_DIR}/all_samples_dedup.bam" \
  O="${OUTPUT_DIR}/all_samples_rg.bam" \
  RGID="merged_samples" \
  RGLB="lib1" \
  RGPL="illumina" \
  RGPU="unit1" \
  RGSM="merged_samples" \
  CREATE_INDEX=true


echo "🔹 Step 4a: Base recalibration (BQSR)"
gatk BaseRecalibrator \
  -R "$ref_fa" \
  -I "${OUTPUT_DIR}/all_samples_rg.bam" \
  --known-sites "$DBSNP" \
  --known-sites "$MANDG_INDELS" \
  --known-sites "$SNPs_1000G" \
  -O "${OUTPUT_DIR}/all_samples_recal_data.table" 

echo "🔹 Step 4b: Apply BQSR"
gatk ApplyBQSR \
  -R "$ref_fa" \
  -I "${OUTPUT_DIR}/all_samples_rg.bam" \
  --bqsr-recal-file "${OUTPUT_DIR}/all_samples_recal_data.table" \
  -O "${OUTPUT_DIR}/all_samples_recal.bam"


echo "🔹 Step 5: Filtering BAM (remove unmapped, secondary, supplementary; MAPQ < 20)"
samtools view -@ 8 -b -F 2308 -q 20 "${OUTPUT_DIR}/all_samples_recal.bam" > "${OUTPUT_DIR}/all_samples_filtered.bam"
samtools index "${OUTPUT_DIR}/all_samples_filtered.bam" > "${OUTPUT_DIR}/all_samples_filtered.bam.bai"


echo "**** Job ends ****"
date +"%Y-%m-%d %T"




