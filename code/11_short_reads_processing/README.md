# 11_short_reads_processing

Short-read RNA-seq processing for PRJNA754196 (retinal organoid differentiation
timecourse), for platform DE comparison against the long-read (bambu) data.

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
GENCODE release_46 primary assembly transcriptome:
  references/genome/GENCODE/primary_assembly/release_46_all_transcripts_short_header.fa

## Pipeline order
1. 01_download_fastqs.sh        - download the FASTQs (transfer partition, array)
2. 02_make_manifests.py         - sample<->run and file manifests
3. 03_check_fastq_integrity.sh  - gzip -t on every file (run BEFORE salmon)
4. 04_salmon_index.sh           - build release_46 index + tx2gene table
5. 05_salmon_quant.sh           - 28-sample quant array (needs salmon_quant_manifest.tsv)
6. 06_tximport_gene_counts.R    - transcript -> gene raw counts (underscore col names)

## Outputs
  10_short_reads/salmon/release_46_index/      salmon index
  10_short_reads/salmon/tx2gene_release46.tsv  transcript->gene map
  10_short_reads/salmon/quants/<sample>/       per-sample quant.sf
  10_short_reads/salmon/gene_counts_matrix.csv gene x 28 raw counts
