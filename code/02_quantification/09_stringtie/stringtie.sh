#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=200G
#SBATCH -c 10
#SBATCH --job-name=stringtie
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/stringtie.txt
#SBATCH -e logs/stringtie.txt
#SBATCH --time=7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

BAM_FOLDER=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality
REFERENCE_GTF=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/gencode.v46.chr_patch_hapl_scaff.basic.annotation.gtf.gz
REFERENCE_FASTA=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa.gz

EP1_BRN3B_RO=$BAM_FOLDER/EP1-BRN3B-RO_primary_over_30_sorted.bam
EP1_WT_hRO_2=$BAM_FOLDER/EP1-WT_hRO_2_primary_over_30_sorted.bam
EP1_WT_ROs_D45=$BAM_FOLDER/EP1-WT_ROs_D45_primary_over_30_sorted.bam
H9_BRN3B_hRO_2=$BAM_FOLDER/H9-BRN3B_hRO_2_primary_over_30_sorted.bam
H9_BRN3B_RO=$BAM_FOLDER/H9-BRN3B-RO_primary_over_30_sorted.bam
H9_CRX_hRO_2=$BAM_FOLDER/H9-CRX_hRO_2_primary_over_30_sorted.bam
H9_CRX_ROs_D45=$BAM_FOLDER/H9-CRX_ROs_D45_primary_over_30_sorted.bam
H9_FT_1=$BAM_FOLDER/H9-FT_1_primary_over_30_sorted.bam
H9_FT_2=$BAM_FOLDER/H9-FT_2_primary_over_30_sorted.bam
H9_hRGC_1=$BAM_FOLDER/H9-hRGC_1_primary_over_30_sorted.bam
H9_hRGC_2=$BAM_FOLDER/H9-hRGC_2_primary_over_30_sorted.bam

OUTPUT_FOLDER=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/stringtie
mkdir -p $OUTPUT_FOLDER

~/stringtie/stringtie --merge -L -G $REFERENCE_GTF -o $OUTPUT_FOLDER/stringtie.gtf \
    $EP1_BRN3B_RO $EP1_WT_hRO_2 $EP1_WT_ROs_D45 $H9_BRN3B_hRO_2 $H9_BRN3B_RO $H9_CRX_hRO_2 \
    $H9_CRX_ROs_D45 $H9_FT_1 $H9_FT_2 $H9_hRGC_1 $H9_hRGC_2

echo "**** Job ends ****"
date +"%Y-%m-%d %T"