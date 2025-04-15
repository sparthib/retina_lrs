# Allele-specific expression (ASE) analysis

This directory contains code and data for the analysis of allele-specific expression (ASE) in human retinal development.

# Proposed plan 

## Calling variants to get a VCF file

The goal of this section is to have as output a VCF file with the variants for each sample. If we have **matched short-read DNA-seq** for the H9 cell lines, we have a few options: 

- [GATK](): One of the most widely used tools for for calling SNVs from DNA-seq data. 
- [monopogen](https://github.com/KChen-lab/Monopogen):  Calls SNVs from single-cell sequencing, including RNA-seq, ATAC-seq, and DNA-seq technologies
- [DeepVariant](https://github.com/google/deepvariant): A deep learning-based variant caller that takes aligned reads (in BAM or CRAM format), produces pileup image tensors from them, classifies each tensor using a convolutional neural network, and finally reports the results in a standard VCF or gVCF file.

If we don't have matched short-read DNA-seq, we can use the our own **long-read RNA-seq** data to call variants. Specifially, we can use the following algorithms: 

- [Claire3RNA](https://github.com/HKU-BAL/Clair3-RNA: too to call variants from long-read RNA (mentions PCR-cDNA and direct RNA, not direct cDNA though).
  - From the Claire3RNA paper: "RNA-seq offers advantages for variant calling and interpretation, such as that the identified allele is guaranteed to be expressed in contrast to DNA variant calling. Nevertheless, several disadvantages need to be considered such as a higher error rate than DNA-seq, which averages a 1-5% error rate. This necessitates robust variant-calling systems for distinguishing true variants from sequencing artifacts. Second, unlike the uniform read coverage typically observed in DNA-seq data [4], the coverage is uneven across genomic regions in RNA-seq, and the variable coverage poses challenges for accurate variant calling, particularly in regions characterized by inadequate read support or excessive read coverage that is not accounted for in neural network training. Additionally, RNA-seq faces RNA editing events, especially A-to-I editing catalyzed by adenosine deaminases acting on RNA (ADAR) [16], leading to A-to-G or T-to-C substitutions that can resemble genuine variants, thus introducing false positives in RNA-seq. Moreover, while phasing has been demonstrated to be advantageous for lrDNA-seq analyses [13, 14], no studies have yet investigated the impact of phasing on variant calling performance in lrRNA-seq, particularly when reads are shorter and zygosity switching can happen."
- are there others? 

## Phase variants using aligned reads

The goal of this section is to take as input (1) aligned reads e.g. a BAM file at the genome or transcriptome level and (2) the variants in the VCF and outputs variants that are phased to haplotypes (H1 or H2). However, it is important that the same reference genome/transcriptome that was used to align the reads be used again when assigning haplotypes. We can use the following tools:

### Genome-level 

- [WhatsHap](https://whatshap.readthedocs.io/en/latest/index.html): A tool to phase genomic variants using DNA-seq reads, also called read-based phasing or haplotype assembly. It is especially suitable for long reads, but works also well with short reads.

### Transcriptome-level

- If we want to phase variants to a diploid transcriptome, we can use WASP to adjust/remove reads from BAM files that are not aligned to the transcriptome. This is useful if we want to phase variants to a diploid transcriptome, e.g. for ASE analysis.
- Align reads from `oarfish` (aligns reads to a transcriptome) and then using the readIDs that were aligned in the BAM file, it will determine which haplotype the read belongs to.
