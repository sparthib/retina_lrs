library(GenomicFeatures)
library(rtracklayer)
library(stringr)

# Paths
gtf_file <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_assembly.gtf"
data_dir  <- Sys.getenv("retina_lrs_dir")
ref_dir   <- Sys.getenv("references_dir")

# Output paths
igv_dir    <- file.path(data_dir, "09_ASE/H9_DNA_Seq_data/whatshap_output_phased_on_H9_and_EP1/igv")
output_bed <- file.path(igv_dir, "genes_of_interest.bed")
output_gtf <- file.path(ref_dir, "genome/GENCODE/primary_assembly/ASE_genes_of_interest.gtf")

dir.create(igv_dir,            recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_gtf), recursive = TRUE, showWarnings = FALSE)

# Genes of interest (Ensembl IDs without version numbers)
genes_of_interest <- c(
  "ENSG00000197921", "ENSG00000142937", "ENSG00000003756", "ENSG00000177879",
  "ENSG00000204356", "ENSG00000124614", "ENSG00000151893", "ENSG00000197345",
  "ENSG00000183955", "ENSG00000111775", "ENSG00000100554", "ENSG00000104129",
  "ENSG00000179918", "ENSG00000168256", "ENSG00000108953", "ENSG00000087111",
  "ENSG00000152234", "ENSG00000181588", "ENSG00000102225", "ENSG00000133169",
  "ENSG00000102144", "ENSG00000184083", "ENSG00000212747", "ENSG00000102409",
  "ENSG00000165169", "ENSG00000204272", "ENSG00000189369"
)

# Load GTF
cat("Loading GTF...\n")
gtf <- import(gtf_file)

# Strip version numbers from gene_id (e.g. ENSG00000197921.6 -> ENSG00000197921)
mcols(gtf)$gene_id_clean <- sub("\\..*", "", mcols(gtf)$gene_id)

# --- BED file: gene-level coordinates ---
gene_entries <- gtf[gtf$type == "gene" & mcols(gtf)$gene_id_clean %in% genes_of_interest]
cat("Found", length(gene_entries), "of", length(genes_of_interest), "requested genes in GTF\n")

bed_df <- data.frame(
  chrom  = as.character(seqnames(gene_entries)),
  start  = start(gene_entries) - 1,       # BED is 0-based
  end    = end(gene_entries),
  name   = mcols(gene_entries)$gene_name,
  score  = 0,
  strand = as.character(strand(gene_entries)),
  stringsAsFactors = FALSE
)

write.table(bed_df, file = output_bed, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
cat("BED file written to:", output_bed, "\n")
print(bed_df[, c("chrom", "start", "end", "name")])

# --- GTF subset: gene / transcript / exon for genes of interest ---
gtf_subset <- gtf[
  (gtf$type == "gene" | gtf$type == "transcript" | gtf$type == "exon") &
  mcols(gtf)$gene_id_clean %in% genes_of_interest
]

# Strip version numbers from all IDs
mcols(gtf_subset)$gene_id       <- sub("\\..*", "", mcols(gtf_subset)$gene_id)
mcols(gtf_subset)$transcript_id <- sub("\\..*", "", mcols(gtf_subset)$transcript_id)

export(gtf_subset, con = output_gtf, format = "gtf")
cat("GTF file written to:", output_gtf, "\n")
cat("Done.\n")
