# Allele-specific expression (ASE) analysis

This directory contains code and data for the analysis of allele-specific expression (ASE) in human retinal development.

# Proposed plan 

## Calling variants to get a VCF file

The goal of this section is to have as output a VCF file with the variants for each sample. If we have **matched short-read DNA-seq** for the H9 cell lines, we have a few options: 

- [GATK](): One of the most widely used tools for for calling SNVs from DNA-seq data. 
- [monopogen](https://github.com/KChen-lab/Monopogen):  Calls SNVs from single-cell sequencing, including RNA-seq, ATAC-seq, and DNA-seq technologies
- [DeepVariant](https://github.com/google/deepvariant): A deep learning-based variant caller that takes aligned reads (in BAM or CRAM format), produces pileup image tensors from them, classifies each tensor using a convolutional neural network, and finally reports the results in a standard VCF or gVCF file.

If we don't have matched short-read DNA-seq, we can use the our own **long-read RNA-seq** data to call variants. Specifially, we can use the following algorithms: 

- [Claire3RNA](https://github.com/HKU-BAL/Clair3-RNA): tool to call variants from long-read RNA (mentions PCR-cDNA and direct RNA, not direct cDNA though).

[More details here on softwares that can be used for VCF generation](https://docs.google.com/spreadsheets/d/1pfP4qprIkk7LFlwLUgiI0IDjMSSvInDqeaqIRpeLpek/edit?gid=0#gid=0)


## VC with RNA-Seq advantages
- Low cost
- Allele guaranteed to be expressed in contrast to DNA variant calling 

## Disadvantages:
- Highly non-uniform coverage (coverage is assumed to be uniform across the genome for DNA-Seq) 
    - For example, 10x coverage doesn’t really mean much for RNA-Seq because it is highly dependent on the genes that are actually expressed. 
    - I don’t understand why this is the case yet, but extremely high coverage in a region can also lead to inaccurate variant calling. 
- Allelic imbalance - can affect gene expression, lower coverage for a gene can be associated to just one allele being expressed. 
    - Can also introduce false negatives for variants (especially coupled with low sequencing depth) . 
- False positives: 
    - High sequencing error rates (especially for cDNA with ONT based on Clair3 benchmark) 
    - RNA editing events, especially A-to-I editing catalyzed by adenosine deaminases acting on RNA (ADAR), leading to A-to-G or T-to-C substitutions that can resemble genuine variants.


## Step 1: Most straightforward approach for VCF generation: 
1. Use GATK and FreeBayes for  short-read WGS data. 
2. Use monopogen for short-read single-cell data. 

## Step 2: Haplotype(phase) aligned reads using VCF

The goal of this section is to take as input (1) aligned reads e.g. a BAM file at the genome or transcriptome level and (2) the variants in the VCF and outputs variants that are phased to haplotypes (H1 or H2). However, it is important that the same reference genome/transcriptome that was used to align the reads be used again when assigning haplotypes. We can use the following tools:

### Genome-level 

- [WhatsHap](https://whatshap.readthedocs.io/en/latest/index.html): A tool to phase genomic variants using DNA-seq reads, also called read-based phasing or haplotype assembly. It is especially suitable for long reads, but works also well with short reads.

### Transcriptome-level

- If we want to phase variants to a diploid transcriptome, we can use WASP to adjust/remove reads from BAM files that are not aligned to the transcriptome. This is useful if we want to phase variants to a diploid transcriptome, e.g. for ASE analysis.
- Align reads from `oarfish` (aligns reads to a transcriptome) and then using the readIDs that were aligned in the BAM file, it will determine which haplotype the read belongs to.
