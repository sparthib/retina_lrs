#!/bin/bash
# Salmon quantification, 28 short-read samples (SRR2118xxxx / Oligo-dT re-deposit series).
# One library per sample, NO deduplication (standard for non-UMI RNA-seq; SMART-Seq v4).
# Reads a manifest: <sample>\t<R1>\t<R2>  (one line per sample).
#SBATCH --partition=shared
#SBATCH --account=jhpce
#SBATCH --job-name=salmon_quant
#SBATCH --array=1-28
#SBATCH --cpus-per-task=6
#SBATCH --mem=16G
#SBATCH --time=3:00:00
#SBATCH --output=salmon_quant_%A_%a.log
set -euo pipefail
module load Salmon/1.10.1
TGT=/dcs04/hicks/data/sparthib/retina_lrs/10_short_reads/input_fastqs
OUT=/dcs04/hicks/data/sparthib/retina_lrs/10_short_reads
IDX=$OUT/salmon/release_46_index
MAN=$OUT/salmon/salmon_quant_manifest.tsv   # sample<TAB>R1<TAB>R2
mkdir -p $OUT/salmon/quants
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MAN")
sample=$(echo "$line" | cut -f1); r1=$(echo "$line" | cut -f2); r2=$(echo "$line" | cut -f3)
echo "task $SLURM_ARRAY_TASK_ID  sample $sample"
salmon quant -i $IDX -l A \
  -1 $TGT/$r1 -2 $TGT/$r2 \
  -p 6 --validateMappings --gcBias --seqBias \
  -o $OUT/salmon/quants/$sample
