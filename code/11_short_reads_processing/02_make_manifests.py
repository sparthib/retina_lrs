#!/usr/bin/env python3
"""Build sample->run mapping and 112-file FASTQ manifest for BioProject PRJNA754196.

Fetches run metadata from the ENA portal API, writes:
  - sample_to_runs.csv        (28 samples x 2 runs each)
  - fastq_file_manifest.csv   (112 FASTQ files: sample, run, R1/R2, filename, size, url)

Usage:  python make_manifests.py [OUTDIR]
"""
import sys, urllib.request, csv, io

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
ACC = "PRJNA754196"
FIELDS = "run_accession,sample_accession,library_name,library_strategy,fastq_ftp,fastq_bytes,read_count"
URL = (f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={ACC}"
       f"&result=read_run&fields={FIELDS}&format=tsv")

rows = list(csv.DictReader(io.StringIO(urllib.request.urlopen(URL).read().decode()), delimiter="\t"))

# sample_to_runs.csv
samples = {}
for r in rows:
    samples.setdefault(r["library_name"], []).append(r)
with open(f"{OUT}/sample_to_runs.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["sample","sample_accession","day","replicate","n_runs",
                "run_1","run_2","run_1_reads","run_2_reads"])
    for name in sorted(samples, key=lambda s: (int(s.split("-")[0][1:]), s)):
        sub = sorted(samples[name], key=lambda r: r["run_accession"])
        runs = [r["run_accession"] for r in sub]
        day, rep = name.split("-")
        w.writerow([name, sub[0]["sample_accession"], day, rep, len(runs),
                    runs[0], runs[1] if len(runs) > 1 else "",
                    sub[0]["read_count"], sub[1]["read_count"] if len(sub) > 1 else ""])

# fastq_file_manifest.csv
with open(f"{OUT}/fastq_file_manifest.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["sample","sample_accession","run_accession","read","filename","size_bytes","url"])
    for r in rows:
        ftps = r["fastq_ftp"].split(";")
        sizes = r["fastq_bytes"].split(";")
        for i, (u, b) in enumerate(zip(ftps, sizes), start=1):
            w.writerow([r["library_name"], r["sample_accession"], r["run_accession"],
                        f"R{i}", u.split("/")[-1], b, "https://" + u])

print("wrote sample_to_runs.csv (28 samples) and fastq_file_manifest.csv (112 files)")
