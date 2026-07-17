# Short-read RNA-seq preprocessing: from raw FASTQ to counts matrix

**Data acquisition.** Short-read RNA-seq data for the retinal organoid
developmental time course were obtained from NCBI BioProject PRJNA754196
(Wahlin lab) — 28 biological samples across seven timepoints (D0, D10, D25,
D65, D100, D180, D280).

**Reference.** GENCODE release 46 primary assembly transcriptome — matching the
long-read bambu quantification so Ensembl gene IDs are directly comparable
across platforms.

**salmon index.** Decoy-free transcriptome index with salmon v1.10.1
(`salmon index -k 31`), 251,955 targets; tx2gene map from the pipe-delimited
FASTA headers.

**salmon quantification.** Each of the 28 samples quantified with
`salmon quant -l A --validateMappings --gcBias --seqBias`.

**Counts matrix (tximport).** Transcript estimates summarized to gene level
with tximport 1.37.2 (R 4.5.x) using the release-46 tx2gene map; version
suffixes stripped, columns renamed to D&lt;day&gt;_&lt;rep&gt; and
timepoint-ordered.
