# 11_short_reads_processing

Short-read RNA-seq processing for PRJNA754196 (retinal organoid differentiation
timecourse), for the quantifier-controlled platform comparison against the
long-read data. Outputs land in `processed_data/short_reads/` (repo) and
`/dcs04/hicks/data/sparthib/retina_lrs/10_short_reads/` (full matrices, data dir).

## Dataset
- BioProject PRJNA754196: 28 samples, 2 SRA runs each (original SRR154356xx +
  re-deposit SRR2118xxxx). The two runs per sample carry the SAME read-pairs
  (verified by read-pair count + dedup-count identity); the "2x read_count"
  seen in SRA metadata is a reads-vs-pairs accounting artifact, not real depth.
- Library prep: SMART-Seq v4 Ultra Low Input (full-length cDNA, oligo-dT/polyA)
  + Nextera XT. No UMIs -> do NOT deduplicate before quantification.
- This pipeline uses the SRR2118xxxx (Oligo-dT-labelled) series, one library
  per sample.

## Reference (matches the long-read bambu quantification)
GENCODE release_46 primary assembly:
  transcriptome FASTA: references/genome/GENCODE/primary_assembly/release_46_all_transcripts_short_header.fa
  genome FASTA + GTF : references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa, release_46_primary_assembly.gtf

## Pipeline order

Preprocessing + salmon quantification (short reads)
1. 01_download_fastqs.sh          - download the FASTQs (transfer partition, array)
2. 02_make_manifests.py           - sample<->run and file manifests
3. 03_check_fastq_integrity.sh    - gzip -t on every file (run BEFORE salmon)
4. 04_salmon_index.sh             - build release_46 index + tx2gene table
5. 05_salmon_quant.sh             - 28-sample quant array (needs salmon_quant_manifest.tsv)
6. 06_tximport_gene_counts.R      - transcript -> gene + transcript raw count matrices (underscore col names)

Isoform-level comparison
7. 07_isoform_detection.R         - unique isoforms detected in the SR salmon data (detection filter)
8. 08_isoform_SR_LR_comparison.R  - SR vs LR isoform overlap (all-known and PC), LR-only isoform/gene lists
9. 09_isoforms_per_gene_SR_LR.R   - isoforms-per-gene distributions, SR vs LR (all genes and PC)
10. 10_filter_genes_isoforms_by_biotype.R - PC + expression-filtered gene/isoform matrices
                                           (filterByExpr min.count=10 gene / 2 isoform + TMM); no DE

Long-read salmon re-quantification (removes the quantifier confound)
11. 11_salmon_quant_longreads.sh  - salmon aln-mode quant of the 7 minimap2 LR BAMs vs release_46
12. 12_LRsalmon_isoform_matrix_and_comparison.R - LR-salmon tximport isoform matrix + SR-vs-LR-salmon comparison

Exon-exon junctions per read (read-length readout)
13. 13_exon_junction_per_read.py  - per-read CIGAR-N junction count over multi-exon PC genes (STAR BAMs, MAPQ=255)
14. 14_exon_junction_SR_vs_LR_plot.py - combine SR + LR per-sample junction distributions, boxplot

Related visualization/DE scripts live under `code/05_visualization/`
(08h platform DE on LR-salmon, 08f platform-DE GO, 08i correlation + length/GC,
08g LR-salmon-only isoform GO).

## Key outputs
  10_short_reads/salmon/release_46_index/          salmon index
  10_short_reads/salmon/tx2gene_release46.tsv       transcript->gene map
  10_short_reads/salmon/quants/<sample>/            per-sample quant.sf (SR)
  10_short_reads/salmon/quants_longreads/<sample>/  per-sample quant.sf (LR aln-mode)
  10_short_reads/salmon/gene_counts_matrix.csv      SR gene x 28 raw counts
  10_short_reads/salmon/transcript_counts_matrix.csv SR transcript x 28 raw counts
  10_short_reads/salmon/LRsalmon_{gene,transcript}_counts_matrix.csv  LR x 7
  10_short_reads/star_bams/<sample>/                STAR genome BAMs (junction counting)

Processed, version-controlled outputs (tables, plots, filtered matrices):
`processed_data/short_reads/` — see the README there for the section-by-section guide.

## Software versions
salmon 1.10.1; STAR 2.7.8a; minimap2 2.26-r1175; R 4.5.3 with tximport 1.37.2,
edgeR 4.8.2, limma 3.66.0.
