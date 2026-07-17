#!/bin/bash
# Verify gzip integrity of all short-read FASTQs (one task per file).
# NOTE: byte-size matching is NOT sufficient QC for gzipped FASTQs -- two files
# in this dataset had the exact correct size but damaged gzip streams, which
# only gzip -t detected (and which broke salmon mid-run). Always run this.
#SBATCH --partition=shared
#SBATCH --account=jhpce
#SBATCH --job-name=gzcheck
#SBATCH --array=1-56
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0:30:00
#SBATCH --output=gzcheck_%A_%a.log
set -uo pipefail
cd /dcs04/hicks/data/sparthib/retina_lrs/10_short_reads/input_fastqs
mapfile -t files < <(ls SRR211*.fastq.gz | sort)
f=${files[$((SLURM_ARRAY_TASK_ID-1))]}
gzip -t "$f" 2>/dev/null && echo "OK $f" || echo "CORRUPT $f"
