#!/usr/bin/env python3
# 13_exon_junction_per_read.py
# Count exon-exon junctions spanned per read (CIGAR N operations) across
# multi-exon protein-coding gene regions, for STAR genome-aligned short-read BAMs.
# Ported from code/02_bam_QC/01_multi_exon_pcg.py (long-read version); the
# CIGAR-N trick is aligner-agnostic. STAR MAPQ=255 = uniquely mapped -> unique-only.
#
# Regions are precomputed once (multiexon_pcg_regions.csv: seqname,start,end,gene_id)
# by an exact gene_id match on the release-46 primary-assembly GTF (multi-exon PC genes).
#
# Usage: python3 13_exon_junction_per_read.py <sample>
import pysam, sys, os
import pandas as pd
from collections import Counter

sample = sys.argv[1]
BAMOUT = "/dcs04/hicks/data/sparthib/retina_lrs/10_short_reads/star_bams"
OUTDIR = f"{BAMOUT}/junction_per_read"
bam = f"{BAMOUT}/{sample}/Aligned.sortedByCoord.out.bam"
regions = pd.read_csv(f"{BAMOUT}/multiexon_pcg_regions.csv")
max_junction = 11

bamfile = pysam.AlignmentFile(bam, "rb")
refs = set(bamfile.references)
nums = []
for contig, start, end in zip(regions['seqname'], regions['start'], regions['end']):
    if contig in refs:
        for read in bamfile.fetch(contig, int(start), int(end)):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality != 255:   # STAR: 255 = uniquely mapped
                continue
            if read.cigarstring is None:
                continue
            nums.append(read.cigarstring.count('N'))
c = Counter(nums)
vec = [c.get(i, 0) for i in range(max_junction)]   # bins 0..10; reads with >=11 junctions omitted
df = pd.DataFrame({sample: vec}).transpose()
per = df.div(df.sum(axis=1), axis=0)
os.makedirs(OUTDIR, exist_ok=True)
per.to_csv(f"{OUTDIR}/{sample}_junction_per_read.csv")
print(f"{sample}: total_reads={sum(vec)} dist={vec}")
