
$INPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/gvcf_ref_46


bcftools view -v snps $INPUT_DIR/all_samples_variants.vcf.gz | bcftools query -f '%CHROM\n' | sort | uniq -c > $INPUT_DIR/all_samples_variants.snp.stats
bcftools view -v indels $INPUT_DIR/all_samples_variants.vcf.gz | bcftools query -f '%CHROM\n' | sort | uniq -c > $INPUT_DIR/all_samples_variants.indel.stats
bcftools stats -s - all_samples_variants.vcf.gz > all_samples_variants.stats
