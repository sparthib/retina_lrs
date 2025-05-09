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

-   `05_whatshap_haplotag.sh` runs the `haplotag` command to tag the reads in our long read BAMs with their haplotype.
