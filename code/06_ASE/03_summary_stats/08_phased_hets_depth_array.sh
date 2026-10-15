#!/bin/bash
#SBATCH --partition=shared
#SBATCH --account=jhpce
#SBATCH --job-name=ph_depth_array
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=02:00:00
#SBATCH --array=1-11
#SBATCH --output=ph_depth_array_%A_%a.log
#
# 08_phased_hets_depth_array.sh
# ---------------------------------------------------------------------------
# For each of the 11 WhatsHap input BAMs (genome-aligned long-read cDNA),
# compute per-position read depth at every phased het in the BED from 07, via
# a single streaming pass (samtools depth -b). One array task per BAM.
# `samtools depth` excludes N-skipped (spliced-out) intron spans by default,
# so this is CONTIGUOUS (allele-callable) coverage.
#
# INPUT : /users/sparthib/phased_hets.bed  (from 07)
#         BAMs in /dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality
# OUTPUT: /users/sparthib/.depth_arr/<sample>.depth  (chrom pos depth)
#
# Run after 07, before 09.
# ---------------------------------------------------------------------------
set -eo pipefail
BAMDIR=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality
BED=/users/sparthib/phased_hets.bed
OUTDIR=/users/sparthib/.depth_arr
mkdir -p $OUTDIR
module load samtools 2>/dev/null
# BAM order matches the whatshap phase command line
SAMPLES=(H9-BRN3B_hRO_2 H9-BRN3B-RO H9-CRX_hRO_2 H9-CRX_ROs_D45 H9-FT_1 H9-FT_2 H9-hRGC_1 H9-hRGC_2 EP1-BRN3B-RO EP1-WT_hRO_2 EP1-WT_ROs_D45)
S=${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}
BAM=$BAMDIR/${S}_primary_over_30_chr_only_sorted.bam
echo "task $SLURM_ARRAY_TASK_ID -> $S ($BAM)"
samtools depth -b $BED $BAM > $OUTDIR/${S}.depth
echo "lines: $(wc -l < $OUTDIR/${S}.depth)"
