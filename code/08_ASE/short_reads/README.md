# Phasing Reads

All relevant code is [here](https://github.com/sparthib/retina_lrs/tree/main/code/08_ASE/short_reads).

## Step 1:

-   3 H9 paired-end WGS samples were obtained from here: <https://www.ncbi.nlm.nih.gov/sra/?term=SRR1091092>

-   SRR1091088, SRR1091091 and SRR1091092.

## Step 2:

-   Using `bowtie` these FASTQ files were aligned to this GENCODE v 46 primary assembly reference genome. Check `01_bowtie.sh`

## Step 3:

-   QC was performed on these BAM files using the GATK pipeline here `02_filter_bams.sh`

## Step 4:

-   This bam file was then used to create VCF file in `03_bam2vcf.sh` and QC.

## Step 5:

-   The `04_whatshap_phase.sh` phases the VCF with the help of our H9 long read BAMs.

-   `05_whatshap_haplotag.sh` runs the `haplotag` command to tag the reads in our long read BAMs with their haplotype, and splits the reads into two separate BAM files, one for each haplotyp using the `split` command.

## Step 6: Number of HP reads per sample

| Sample         | HP1 reads | HP2 reads | Total aligned reads | Haplotyped Reads % |
|-----------------|---------------|---------------|---------------|---------------|
| H9-BRN3B_hRO_2 | 1099492   | 1037639   | 9272631             | 23.04              |
| H9-BRN3B-RO    | 2338537   | 2176963   | 14592201            | 30.95              |
| H9-CRX_hRO_2   | 1005074   | 963019    | 6417395             | 30.70              |
| H9-CRX_ROs_D45 | 1466690   | 1404057   | 10109987            | 28.38              |
| H9-FT_1        | 1179483   | 1121926   | 6363243             | 36.19              |
| H9-FT_2        | 1086558   | 1013173   | 6084602             | 34.58              |
| H9-hRGC_1      | 430019    | 410060    | 2585901             | 32.46              |
| H9-hRGC_2      | 1702447   | 1601349   | 12780770            | 25.82              |

## 

## Step 7: Using `featureCounts` for gene expression estimation

50 - 60% of these genome alignments were mapped to genes here `06_gene_counts.sh`

## Step 8: DE analysis on HP1 vs HP2

`07_allele_spec_expression_modeling.R` performs DGE analysis between the two haplotypes:

-    Down-regulated in H2: 614

-   NotSig 41063

-   Up-regulated in H2: 291

We may have to remove some `pseudogene biotypes`, but I've kept them all for now.

Counts matrix can be found [here](https://github.com/sparthib/retina_lrs/blob/main/processed_data/ASE/H9_DNA_Seq_data_gene_counts.tsv).

DE results [here](https://github.com/sparthib/retina_lrs/blob/main/processed_data/ASE/H9_DNA_Seq_data_DE_results.tsv).
