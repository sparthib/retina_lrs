
#!/bin/bash

INPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output
vcf_stats_dir=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/vcf_stats
config=$vcf_stats_dir/config.toml

declare -A vcfs=(
  ["H9"]="$INPUT_DIR/all_samples_phased.vcf"
  ["H9_EP1"]="$INPUT_DIR/all_samples_H9_and_EP1_phased.vcf"
)

for sample in "${!vcfs[@]}"; do
  vcf=${vcfs[$sample]}
  outdir=$vcf_stats_dir/$sample
  mkdir -p "$outdir"

  echo "Running stats for $sample"

  ## 1. Number of variants on each chromosome DONE 

  vcfstats --vcf "$vcf" \
    --outdir "$outdir" \
    --formula 'COUNT(1) ~ CONTIG' \
    --title 'Number of variants on each chromosome' \
    --config "$config" \
    --ggs 'scale_x_discrete(name ="Chromosome", \
        limits=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10", \
        "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"]); \
        ylab("# Variants")' \
    --save

  # ## 2. Type of substitutions DONE
  vcfstats --vcf "$vcf" \
      --outdir "$outdir" \
      --formula 'COUNT(1, VARTYPE[snp]) ~ SUBST[A>T,A>G,A>C,T>A,T>G,T>C,G>A,G>T,G>C,C>A,C>T,C>G]' \
      --title 'Number of substitutions of SNPs (passed)' \
      --config "$config" \
      --save

  ## 3. Types of variants on whole genome DONE
  vcfstats --vcf "$vcf" \
      --outdir "$outdir" \
      --formula 'COUNT(1, group=VARTYPE) ~ 1' \
      --title 'Types of variants on whole genome' \
      --config "$config" \ 
      --save

  ## 4. Types of variants on each chromosome
  vcfstats --vcf "$vcf" \
      --outdir "$outdir" \
      --formula 'COUNT(1, group=VARTYPE) ~ CONTIG' \
      --title 'Types of variants on each chromosome' \
      --save \
      --config "$config" --ggs 'scale_x_discrete(name ="Chromosome", \
         limits=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10", \
         "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"]); \
         ylab("# Variants")'

  ## 5. Mean (genotype quality) and mean (depth) on each chromosome
  vcfstats --vcf "$vcf" \
      --outdir "$outdir" \
      --formula 'MEAN(GQs{0}) ~ MEAN(DEPTHs{0}, group=CHROM[chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10 \
      chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX])'\
      --title 'GQ vs depth (sample 1)' \
      --config "$config" \
      --save \
      --ggs 'scale_x_continuous(limits=(0,50)); scale_y_continuous(limits=(0,100))'

done

    

  ## 4. Alternative Allele Frequency 
  # vcfstats --vcf "$vcf" \
  #     --outdir "$outdir" \
  #     --formula 'AAF ~ CONTIG' \
  #     --title 'Allele frequency on each chromosome' \
  #     --save \
  #     --config "$config" --figtype boxplot --ggs 'scale_x_discrete(name ="Chromosome", \
  #        limits=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10", \
  #        "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"]); \
  #        ylab("# Variants")'
