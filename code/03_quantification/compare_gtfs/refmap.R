######

data_dir <- Sys.getenv("retina_lrs_dir")
code_dir <- Sys.getenv("retina_lrs_code")
ref_dir <- Sys.getenv("references_dir")


refmap_dir <- file.path(data_dir, "06_quantification/isoquant/high_quality/all_samples/OUT/")
refmap <- "bambu_ref_isoquant_ext.OUT.extended_annotation.gtf.refmap"

refmap_path <- paste(refmap_dir, refmap, sep="")

refmap <- read.table(refmap_path, header=TRUE, sep="\t")

#filter out ref_id where isoform doesn't start with ENST

bambu_novel_refmap <- refmap[!grepl("^ENST", refmap$ref_id), ]
nrow(bambu_novel_refmap)

bambu_novel_refmap <- refmap[!grepl("^ENST", refmap$ref_id) & refmap$class_code == "=", ]

head(bambu_novel_refmap)

#split qry_id_list based on | 

data_split <- separate(bambu_novel_refmap, qry_id_list, 
                       into = c("isoquant_gene_id", "isoquant_isoform_id"), sep = "\\|")



# nrow(bambu_novel_refmap)
#599 

#write out the refmap
# /dcs04/hicks/data/sparthib/retina_lrs/06_quantification/gffcompare
output_dir <- file.path(data_dir, "06_quantification/bambu/")
write.table(data_split, file=paste(output_dir, "bambu_isoquant_refmap.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)

