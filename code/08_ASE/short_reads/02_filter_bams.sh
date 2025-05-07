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
# 
known_sites=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/known_sites
DBSNP=${known_sites}/Homo_sapiens_assembly38.dbsnp138.vcf
MANDG_INDELS=${known_sites}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
SNPs_1000G=${known_sites}/1000G_phase1.snps.high_confidence.hg38.vcf.gz
ref_fa=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa
INPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/sams_ref_46
OUTPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/filtered_bams_ref_46
mkdir -p $OUTPUT_DIR

ml load gatk
ml load samtools

## first mark duplicates, then filter by MAPQ and uniquely mapped reads 

# samples=(SRR1091088 SRR1091091 SRR1091092)
samples=(SRR1091091)
mkdir -p "$OUTPUT_DIR"

for sample in "${samples[@]}"
do
  echo "🔄 Processing sample: $sample"
  
## The following steps were run already
#   echo "🔹 Step 1a: Change SAM to BAM"
#   samtools view -@ 8 -bS ${INPUT_DIR}/${sample}.sam > ${OUTPUT_DIR}/${sample}.bam
#   
#   echo "🔹 Step 1b: Sorting BAM"
#   samtools sort -@ 8 -o ${OUTPUT_DIR}/${sample}_sorted.bam ${OUTPUT_DIR}/${sample}.bam
# 
#     
#     echo "🔹 Step 2: Marking duplicates"
#   gatk MarkDuplicates \
#     I="${OUTPUT_DIR}/${sample}_sorted.bam" \
#     O="${OUTPUT_DIR}/${sample}_dedup.bam" \
#     M="${OUTPUT_DIR}/${sample}_metrics.txt" \
#     CREATE_INDEX=true
# 
#   echo "🔹 Step 3: Adding read groups"
#   gatk AddOrReplaceReadGroups \
#     I="${OUTPUT_DIR}/${sample}_dedup.bam" \
#     O="${OUTPUT_DIR}/${sample}_rg.bam" \
#     RGID="${sample}" \
#     RGLB="lib1" \
#     RGPL="illumina" \
#     RGPU="unit1" \
#     RGSM="${sample}" \
#     CREATE_INDEX=true
    
    echo "🔹 Step 4a: Base recalibration (BQSR - BaseRecalibrator)"
  gatk BaseRecalibrator \
    -R $ref_fa \
    -I ${OUTPUT_DIR}/${sample}_rg.bam \
    --known-sites "$DBSNP" \
    --known-sites "$MANDG_INDELS" \
    --known-sites "$SNPs_1000G" \
    -O "${OUTPUT_DIR}/${sample}_recal_data.table"

  echo "🔹 Step 4b: Apply BQSR"
  gatk ApplyBQSR \
    -R "$ref_fa" \
    -I "${OUTPUT_DIR}/${sample}_rg.bam" \
    --bqsr-recal-file "${OUTPUT_DIR}/${sample}_recal_data.table" \
    -O "${OUTPUT_DIR}/${sample}_recal.bam"

  echo "🔹 Step 5: Filtering BAM (remove unmapped, secondary, supplementary; MAPQ < 20)"
  samtools view -@ 8 -b -F 2308 -q 20 "${OUTPUT_DIR}/${sample}_recal.bam" > "${OUTPUT_DIR}/${sample}_filtered.bam"
  samtools -@ 8 index "${OUTPUT_DIR}/${sample}_filtered.bam"
  
done 
# samples=(SRR1091088 SRR1091091 SRR1091092)
# 
# for sample in ${samples[@]}
# do
#     samtools view -bS -h -F 0x904 -q 20 $INPUT_DIR/${sample}.sam > $OUTPUT_DIR/${sample}.bam
#     samtools sort $OUTPUT_DIR/${sample}.bam -o $OUTPUT_DIR/${sample}.sorted.bam
#     samtools index $OUTPUT_DIR/${sample}.sorted.bam $OUTPUT_DIR/${sample}.sorted.bam.bai
#     rm $OUTPUT_DIR/${sample}.bam
# done

echo "**** Job ends ****"
date +"%Y-%m-%d %T"




