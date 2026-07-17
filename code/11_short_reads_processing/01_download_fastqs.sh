#!/bin/bash
#SBATCH --partition=transfer
#SBATCH --account=jhpce
#SBATCH --job-name=retina_ena_dl
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=%x_%j.log
#
# Download all FASTQ files for BioProject PRJNA754196 (retina organoid short-read RNA-seq)
# from ENA into the target directory, with per-file size verification.
# 56 paired-end runs -> 112 .fastq.gz files (~162 GB). 28 samples x 2 runs each.
# Re-runnable: skips files already present at the correct size.
#
# Usage:  sbatch download_fastqs.sh
# Expects fastq_file_manifest.csv in TGT (columns incl. filename,size_bytes,url).

set -uo pipefail
TGT=/dcs04/hicks/data/sparthib/retina_lrs/10_short_reads/input_fastqs
mkdir -p "$TGT"
echo "host: $(hostname)  job: ${SLURM_JOB_ID:-NA}"

MANIFEST="$TGT/fastq_file_manifest.csv"
# build url<TAB>filename<TAB>size from the manifest csv (skip header)
TSV=/tmp/dl_${SLURM_JOB_ID:-x}.tsv
awk -F, 'NR>1{print $7"\t"$5"\t"$6}' "$MANIFEST" > "$TSV"

n=$(wc -l < "$TSV"); echo "files to fetch: $n"
ok=0; skip=0; fail=0
while IFS=$'\t' read -r url fname exp; do
  dest="$TGT/$fname"
  if [ -f "$dest" ]; then
    act=$(stat -c%s "$dest" 2>/dev/null || echo 0)
    if [ "$act" = "$exp" ]; then echo "SKIP $fname ($act)"; skip=$((skip+1)); continue; fi
    echo "REDL $fname (have $act want $exp)"; rm -f "$dest"
  fi
  if wget -q --tries=5 --waitretry=15 --timeout=120 -O "$dest" "$url"; then
    act=$(stat -c%s "$dest" 2>/dev/null || echo 0)
    if [ "$act" = "$exp" ]; then echo "OK   $fname ($act)"; ok=$((ok+1));
    else echo "BADSIZE $fname (got $act want $exp)"; fail=$((fail+1)); fi
  else
    echo "FAIL $fname"; fail=$((fail+1))
  fi
done < "$TSV"

echo "=== summary ok=$ok skip=$skip fail=$fail total=$n ==="
echo "fastq present: $(ls "$TGT"/*.fastq.gz 2>/dev/null | wc -l)"
du -sh "$TGT" 2>/dev/null
rm -f "$TSV"
[ "$fail" -eq 0 ]
