#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH -c 10
#SBATCH --job-name=stringtie
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --array=0-10 
#SBATCH --output=logs/stringtie_%a.out  # Standard output log
#SBATCH --error=logs/stringtie_%a.err   # Standard error log
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
OUTPUT_FOLDER=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/stringtie

# Define an array of sample names
SAMPLES=("EP1-BRN3B-RO" "EP1-WT_hRO_2" "EP1-WT_ROs_D45" "H9-BRN3B_hRO_2" "H9-BRN3B-RO" \
"H9-CRX_hRO_2" "H9-CRX_ROs_D45" "H9-FT_1" "H9-FT_2" "H9-hRGC_1" "H9-hRGC_2")

# Get the sample name based on the SLURM_ARRAY_TASK_ID
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Define input BAM and output GTF paths
BAM_FILE="${BAM_FOLDER}/${SAMPLE}_primary_over_30_sorted.bam"
OUTPUT_FILE="${OUTPUT_FOLDER}/${SAMPLE}.gtf"

# Run stringtie
echo "Running stringtie for sample: $SAMPLE"
~/stringtie/stringtie -L -G $REFERENCE_GTF -o $OUTPUT_FILE  $BAM_FILE