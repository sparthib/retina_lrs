#!/bin/bash
# 11_salmon_quant_longreads.sh
# Salmon alignment-mode (-a) quantification of the 7 retinal-organoid long-read
# samples, against the SAME release-46 transcriptome FASTA used for the short-read
# salmon index. Purpose: quantify LR with the identical salmon EM apportionment as
# SR, removing the quantifier confound in the isoform-per-gene comparison.
#
# Input BAMs: minimap2 -ax map-ont -N 100 --secondary=no, query-grouped (unsorted),
# aligned to release_46_all_transcripts_short_header.fa (252,835 transcripts).
# Submit: sbatch --array=1-7 11_salmon_quant_longreads.sh
#SBATCH --partition=shared
#SBATCH --account=jhpce
#SBATCH --job-name=salmon_lr
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=/users/sparthib/salmon_lr_%A_%a.log
set -eo pipefail
module load salmon/1.10.1 2>/dev/null || module load Salmon/1.10.1

DD=/dcs04/hicks/data/sparthib/retina_lrs
FA=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_all_transcripts_short_header.fa
BAMDIR=$DD/05_bams/transcriptome/ver_46
OUT=$DD/10_short_reads/salmon/quants_longreads
mkdir -p $OUT

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $HOME/.lr_salmon_manifest.txt)
echo "task ${SLURM_ARRAY_TASK_ID} -> sample: $SAMPLE"

salmon quant -t $FA -l A -a $BAMDIR/${SAMPLE}.bam \
  -p 8 --noErrorModel -o $OUT/${SAMPLE}

echo "DONE $SAMPLE : $(grep -o '\"num_mapped\":[0-9]*' $OUT/${SAMPLE}/aux_info/meta_info.json 2>/dev/null || echo NA)"
