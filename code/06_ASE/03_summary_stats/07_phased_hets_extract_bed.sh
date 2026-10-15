#!/bin/bash
#SBATCH --partition=shared
#SBATCH --account=jhpce
#SBATCH --job-name=ph_extract_bed
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=00:20:00
#SBATCH --output=ph_extract_bed_%j.log
#
# 07_phased_hets_extract_bed.sh
# ---------------------------------------------------------------------------
# Extract all PHASED heterozygous sites (GT 0|1 or 1|0) from the WhatsHap
# phased VCF into a BED of 1-bp intervals. This BED feeds the per-BAM depth
# array (08) and the TSS/TTS distance analysis (11).
#
# INPUT : /dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output/all_samples_H9_and_EP1_phased.vcf.gz
# OUTPUT: /users/sparthib/phased_hets.bed   (2,041,396 positions)
#
# Run first, before 08.
# ---------------------------------------------------------------------------
set -eo pipefail
VCF=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output/all_samples_H9_and_EP1_phased.vcf.gz
BED=/users/sparthib/phased_hets.bed
module load htslib bcftools 2>/dev/null
bcftools query -f '%CHROM\t%POS\t[%GT]\n' $VCF \
  | awk -F'\t' '$3=="0|1"||$3=="1|0"{print $1"\t"($2-1)"\t"$2}' > $BED
echo "phased hets: $(wc -l < $BED)"
