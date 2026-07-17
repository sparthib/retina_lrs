#!/bin/bash
# Build salmon index for GENCODE release_46 transcriptome + tx2gene table.
# Matches the reference used by the long-read bambu quantification
# (code/03_quantification/bambu/02a_bambu_generate_rcs.R -> release_46_primary_genome.fa).
#SBATCH --partition=shared
#SBATCH --account=jhpce
#SBATCH --job-name=salmon_index
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=salmon_index_%j.log
set -euo pipefail
module load Salmon/1.10.1
PA=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly
OUT=/dcs04/hicks/data/sparthib/retina_lrs/10_short_reads
mkdir -p $OUT/salmon
# index from the clean-header transcript fasta (ENST... only -> clean quant.sf keys)
salmon index -t $PA/release_46_all_transcripts_short_header.fa \
  -i $OUT/salmon/release_46_index -k 31 -p 8
# tx2gene from the pipe-delimited full-header fasta: >ENST|ENSG|...
grep '^>' $PA/release_46_all_transcripts.fa | sed 's/^>//' \
  | awk -F'|' '{print $1"\t"$2}' > $OUT/salmon/tx2gene_release46.tsv
echo "index + tx2gene done"
