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


INPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/sams_ref_46
OUTPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/filtered_bams_ref_46
mkdir -p $OUTPUT_DIR

ml load gatk
ml load samtools

## first mark duplicates, then filter by MAPQ and uniquely mapped reads 

samples=(SRR1091088 SRR1091091 SRR1091092)

mkdir -p "$OUTPUT_DIR"

for sample in "${samples[@]}"
do
  echo "🔄 Processing sample: $sample"

  echo "🔹 Step 1: Sorting SAM to BAM"
  samtools view -bS ${INPUT_DIR}/${sample}.sam > ${OUTPUT_DIR}/${sample}.bam
  gatk SortSam \
    I="${OUTPUT_DIR}/${sample}.bam" \
    O="${OUTPUT_DIR}/${sample}-sorted.bam" \
    SORT_ORDER=coordinate
    
    echo "🔹 Step 2: Marking duplicates"
  gatk MarkDuplicates \
    I="${OUTPUT_DIR}/${sample}-sorted.bam" \
    O="${OUTPUT_DIR}/${sample}-dedup.bam" \
    M="${OUTPUT_DIR}/${sample}-metrics.txt" \
    CREATE_INDEX=true

  echo "🔹 Step 3: Adding read groups"
  gatk AddOrReplaceReadGroups \
    I="${OUTPUT_DIR}/${sample}-dedup.bam" \
    O="${OUTPUT_DIR}/${sample}-rg.bam" \
    RGID="${sample}" \
    RGLB="lib1" \
    RGPL="illumina" \
    RGPU="unit1" \
    RGSM="${sample}" \
    CREATE_INDEX=true
    
    echo "🔹 Step 4a: Base recalibration (BQSR - BaseRecalibrator)"
  gatk BaseRecalibrator \
    -I "${OUTPUT_DIR}/${sample}-rg.bam" \
    -R "$REF" \
    --known-sites "$DBSNP" \
    -O "${OUTPUT_DIR}/${sample}-recal_data.table"

  echo "🔹 Step 4b: Apply BQSR"
  gatk ApplyBQSR \
    -R "$REF" \
    -I "${OUTPUT_DIR}/${sample}-rg.bam" \
    --bqsr-recal-file "${OUTPUT_DIR}/${sample}-recal_data.table" \
    -O "${OUTPUT_DIR}/${sample}-recal.bam"

  echo "🔹 Step 5: Filtering BAM (remove unmapped, secondary, supplementary; MAPQ < 20)"
  samtools view -b -F 2308 -q 20 "${OUTPUT_DIR}/${sample}-recal.bam" > "${OUTPUT_DIR}/${sample}-filtered.bam"
  samtools index "${OUTPUT_DIR}/${sample}-filtered.bam"
  
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




