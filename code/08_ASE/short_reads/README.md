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

| Sample         | HP1 reads | HP2 reads | Total aligned reads |
|----------------|-----------|-----------|---------------------|
|                |           |           |                     |
| H9-BRN3B-RO    | 2338537   | 2176963   | 14592201            |
| H9-CRX_hRO_2   | 1005074   | 963019    | 6417395             |
| H9-CRX_ROs_D45 | 1466690   | 1404057   | 10109987            |
| H9-FT_1        | 1179483   | 1121926   | 6363243             |
| H9-FT_2        | 1086558   | 1013173   | 6084602             |
| H9-hRGC_1      | 430019    | 410060    | 2585901             |
| H9-hRGC_2      |           |           | 12780770            |
|                |           |           |                     |
